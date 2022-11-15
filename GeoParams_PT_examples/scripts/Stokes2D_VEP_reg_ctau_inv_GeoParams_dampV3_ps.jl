using ElasticArrays, Printf, GeoParams
using Plots, Plots.Measures
using ParallelStencil
using ParallelStencil.FiniteDifferences2D
@init_parallel_stencil(CUDA, Float64, 2)

plot_opts() = (aspect_ratio = 1, xlims = (0.007936507936507936, 0.9920634920634921), ylims = (0.007936507936507936, 0.9920634920634921), c = cgrad(:roma, rev=true), framestyle = :box)

@parallel function update_old!(τxx_o, Tyy_o, τxy_o, τxyv_o, Pt_old, τxx, τyy, τxy, τxyv, Pt)
    @all(τxx_o) = @all(τxx)
    @all(Tyy_o) = @all(τyy)
    @all(τxy_o) = @all(τxy)
    @all(τxyv_o) = @all(τxyv)
    @all(Pt_old) = @all(Pt)
    return nothing
end

@parallel_indices (i,j) function update_iteration_params!(η, ητ)

    @inbounds @inline function maxloc(A)
        max(
            A[i-1, j-1],
            A[i  , j-1],
            A[i+1, j-1],
            A[i-1, j  ],
            A[i  , j  ],
            A[i+1, j  ],
            A[i-1, j+1],
            A[i  , j+1],
            A[i+1, j+1],
        )
    end

    ητ[i, j] = maxloc(η)
    
    return nothing
end

function bc2!(ητ)
    nx, ny = size(ητ)  
    @parallel (1:nx) bc2_x!(ητ)
    @parallel (1:ny) bc2_y!(ητ)
end

@parallel_indices (j) function bc2_y!(A)
    nx = size(A, 1)  
    A[1 , j] = A[2   , j]
    A[nx, j] = A[nx-1, j]
    return nothing
end

@parallel_indices (i) function bc2_x!(A)
    ny = size(A, 2)  
    A[i, 1] = A[i, 2   ]
    A[i, ny] = A[i, ny-1]
    return nothing
end

@parallel_indices (i, j) function update_stress_GP_ps!(τxx, τyy, τxy, TII, τxx_o, Tyy_o, τxyv_o, εxx, εyy, εxyv, η_vep, Pt, Phasec, Phasev, MatParam, dt, lτ, r, re_mech, vdτ)

    # convinience closure (note here that i, j are captured because the closure is defined inside the loop when @parallel_indices is expanded)
    @inbounds @inline gather(A) = A[i, j], A[i + 1, j], A[i, j + 1], A[i + 1, j + 1] 

    # numerics
    θ_dτ                = lτ * (r + 2.0) / (re_mech * vdτ)
    dτ_r                = 1.0 / (θ_dτ + 1 / η_vep[i, j]) # equivalent to dτ_r = @. 1.0/(θ_dτ + η/(G*dt) + 1.0)
    # Setup up input for GeoParams.jl
    args                = (; dt=dt, P=Pt[i, j], τII_old=0.0)
    εij_p               = (εxx[i, j], εyy[i, j], gather(εxyv))
    τij_p_o             = (τxx_o[i,j], Tyy_o[i,j], gather(τxyv_o))
    phases              = (Phasec[i, j], Phasec[i, j], gather(Phasev))
    # update stress and effective viscosity
    Tij, TII[i, j],  ηᵢ = compute_τij(MatParam, εij_p, args, τij_p_o, phases)
    τxx[i,j]           += dτ_r * (-(τxx[i,j]) + Tij[1] ) /  ηᵢ # NOTE: from GP Tij = 2*η_vep * εij
    τyy[i,j]           += dτ_r * (-(τyy[i,j]) + Tij[2] ) /  ηᵢ 
    τxy[i,j]           += dτ_r * (-(τxy[i,j]) + Tij[3] ) /  ηᵢ 
    η_vep[i, j]         =  ηᵢ

    return
end

@parallel function vertex2center!(center, vertex)
    @all(center) = @av(vertex)
    return nothing
end

@parallel function center2vertex!(vertex, center)
    @inn(vertex) = @av(center)
    return nothing
end

@parallel function compute_∇V!(∇V, Vx, Vy, dx, dy) 
    @all(∇V) = @d_xa(Vx) / dx +  @d_ya(Vy) / dy
    return nothing
end

@parallel function compute_P!(P, P_old, RP, ∇V, η, K, dt, r, θ_dτ)
    @all(RP) = (-@all(∇V) - (@all(P) - @all(P_old)) / (@all(K) * dt))
    @all(P)  = @all(P) + (-@all(∇V) - (@all(P) - @all(P_old)) / (@all(K) * dt)) / (1.0 / (r / θ_dτ * @all(η)) + 1.0 /(@all(K) * dt))
    return nothing
end

@parallel function compute_strain_rate!(∇V, εxx, εyy, εxyv, Vx, Vy, dx, dy)
    @all(εxx) = @d_xa(Vx) / dx - @all(∇V) / 3.0
    @all(εyy) = @d_ya(Vy) / dy - @all(∇V) / 3.0
    @inn(εxyv) = 0.5 * (@d_yi(Vx) / dy + @d_xi(Vy) / dx)
    return nothing
end

@parallel_indices (i, j)  function update_velocities!(Vx, Vy, P, τxx, τyy, τxyv, ηdτ, ρgx, ρgy, ητ, dx, dy)

    # Again, indices i, j are captured by the closure
    @inbounds @inline d_xa(A)  = (A[i+1, j  ] - A[i  , j  ]) / dx 
    @inbounds @inline d_ya(A)  = (A[i  , j+1] - A[i  , j  ]) / dy
    @inbounds @inline d_xi(A)  = (A[i+1, j+1] - A[i  , j+1]) / dx 
    @inbounds @inline d_yi(A)  = (A[i+1, j+1] - A[i+1, j  ]) / dy
    @inbounds @inline av_xa(A) = (A[i  , j  ] + A[i+1, j  ]) * 0.5
    @inbounds @inline av_ya(A) = (A[i  , j  ] + A[i  , j+1]) * 0.5

    if i ≤ size(P, 1) - 1
        Vx[i+1, j] += (-d_xa(P) + d_xa(τxx) + d_yi(τxyv) - ρgx[i, j]) * ηdτ / av_xa(ητ)
    end
    if j ≤ size(P, 2) - 1
        Vy[i, j+1] += (-d_ya(P) + d_ya(τyy) + d_xi(τxyv) - ρgy[i, j]) * ηdτ / av_ya(ητ)
    end

    return nothing
end

@parallel function compute_residuals!(Rx, Ry, P, τxx, τyy, τxy, ρgx, ρgy, dx, dy)
    @all(Rx) = @d_xa(τxx) / dx + @d_yi(τxy) / dy - @d_xa(P) / dx + @inn_y(ρgx)
    @all(Ry) = @d_ya(τyy) / dy + @d_xi(τxy) / dx - @d_ya(P) / dy + @inn_x(ρgy)
    return nothing
end

function main(nt)
    # Physics
    do_DP = true                  # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg = 8.0e-3                # regularisation "viscosity"
    Lx, Ly = 1.0, 1.0             # domain size
    radi = 0.01                   # inclusion radius
    τ_y = 1.6                     # yield stress. If do_DP=true, τ_y stand for the cohesion: c*cos(ϕ)
    ϕ = 30.0                      # friction angle
    μ0 = 1.0                      # viscous viscosity
    G0 = 1.0                      # elastic shear modulus
    Gi = G0 / (6.0 - 4.0 * do_DP) # elastic shear modulus perturbation
    εbg = 1.0                     # background strain-rate
    Coh = τ_y / cosd(ϕ)           # cohesion; if we let Coh = Inf, we recover visco-elastic problem
    # Geoparams initialisation
    pl = DruckerPrager_regularised(; C=Coh, ϕ=ϕ, η_vp=η_reg, Ψ=0)        # non-regularized plasticity
    MatParam = (
        SetMaterialParams(;
            Name="Matrix",
            Phase=1,
            CompositeRheology=CompositeRheology(
                SetConstantElasticity(; G=G0, ν=0.5), LinearViscous(; η=μ0), pl
            ),
        ),
        SetMaterialParams(;
            Name="Inclusion",
            Phase=2,
            CompositeRheology=CompositeRheology(
                SetConstantElasticity(; G=Gi, ν=0.5), LinearViscous(; η=μ0), pl
            ),
        ),
    )
    # Numerics
    nx, ny  = 255, 255             # numerical grid resolution
    re_mech = 3π
    lτ      = min(Lx, Ly)
    dx, dy  = Lx / nx, Ly / ny
    vdτ     = min(dx, dy) / √2.1 #* 0.5 
    r       = 0.7
    θ_dτ    = lτ * (r + 2.0) / (re_mech * vdτ)
    ηdτ     = vdτ * lτ / re_mech
    # Preprocessing
    dx, dy  = Lx / nx, Ly / ny
    dt      = μ0 / G0 / 4.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = @zeros(nx, ny)
    Pt_old  = @zeros(nx, ny)
    RP      = @zeros(nx, ny)
    ∇V      = @zeros(nx, ny)
    Vx      = @zeros(nx + 1, ny)
    Vy      = @zeros(nx, ny + 1)
    εxx     = @zeros(nx, ny)
    εyy     = @zeros(nx, ny)
    Exy     = @zeros(nx, ny)
    εxyv    = @zeros(nx + 1, ny + 1)
    τxx     = @zeros(nx, ny)
    τyy     = @zeros(nx, ny)
    τxy     = @zeros(nx, ny)
    τxyv    = @zeros(nx + 1, ny + 1)
    τxx_o   = @zeros(nx, ny)
    Tyy_o   = @zeros(nx, ny)
    τxy_o   = @zeros(nx, ny)
    τxyv_o  = @zeros(nx + 1, ny + 1)
    TII     = @zeros(nx, ny)
    Rx      = @zeros(nx - 1, ny - 2)
    Ry      = @zeros(nx - 2, ny - 1)
    η_vep   = @ones(nx, ny)
    η       = @ones(nx, ny)
    ητ      = @zeros(nx, ny)
    Phasec  = Int.(@ones(nx, ny))
    Phasev  = Int.(@ones(nx + 1, ny + 1))
    ρgx     = @zeros(nx - 1, ny    )
    ρgy     = @zeros(nx    , ny - 1)
    Kb      = @fill(Inf, nx, ny)
    # Initial condition
    xc, yc                = LinRange(dx / 2, Lx - dx / 2, nx), LinRange(dy / 2, Ly - dy / 2, ny)
    xc, yc                = LinRange(dx / 2, Lx - dx / 2, nx), LinRange(dy / 2, Ly - dy / 2, ny)
    xv, yv                = LinRange(0.0, Lx, nx + 1), LinRange(0.0, Ly, ny + 1)
    Xvx                   = Data.Array([x for x in xv, y in yc])
    Yvy                   = Data.Array([y for x in xc, y in yv])
    radc                  = (xc .- Lx ./ 2) .^ 2 .+ (yc' .- Ly ./ 2) .^ 2
    radv                  = (xv .- Lx ./ 2) .^ 2 .+ (yv' .- Ly ./ 2) .^ 2
    Phasec[radc .< radi] .= 2
    Phasev[radv .< radi] .= 2
    Vx                   .=   εbg .* Xvx
    Vy                   .= .-εbg .* Yvy
    # Time loop
    t = 0.0
    maxiter = 10e3
    ncheck = 500
    ϵtol = 1e-6
    evo_t = Float64[]
    iter_evo = Int64[]
    evo_t = Float64[]
    evo_τxx = Float64[]
    errs_evo = ElasticMatrix{Float64}(undef, length(ϵtol), 0)
    for it in 1:nt
        @printf("it=%d\n", it)
        @parallel update_old!(τxx_o, Tyy_o, τxy_o, τxyv_o, Pt_old, τxx, τyy, τxy, τxyv, Pt)
        errs = 2.0 .* ϵtol
        iter = 1
        # resize!(iter_evo, 0)
        resize!(errs_evo, length(ϵtol), 0)
        while any(errs .>= ϵtol) && iter <= maxiter
            @parallel (2:nx-1, 2:ny-1) update_iteration_params!(η, ητ)
            bc2!(ητ)

            @parallel compute_∇V!(∇V, Vx, Vy, dx, dy) 
            @parallel compute_P!(Pt, Pt_old, RP, ∇V, η, Kb, dt, r, θ_dτ)
            @parallel compute_strain_rate!(∇V, εxx, εyy, εxyv, Vx, Vy, dx, dy)
            @parallel vertex2center!(Exy, εxyv)
            @parallel (1:nx, 1:ny) update_stress_GP_ps!(τxx, τyy, τxy, TII, τxx_o, Tyy_o, τxyv_o, εxx, εyy, εxyv, η_vep, Pt, Phasec, Phasev, MatParam, dt, lτ, r, re_mech, vdτ)
            @parallel center2vertex!(τxyv, τxy)

            @parallel (1:nx, 1:ny) update_velocities!(Vx, Vy, Pt, τxx, τyy, τxyv, ηdτ, ρgx, ρgy, ητ, dx, dy)

            if iter % ncheck == 0
                # update residuals
                @parallel compute_residuals!(Rx, Ry, Pt, τxx, τyy, τxyv, ρgx, ρgy, dx, dy)
                errs = maximum.((abs.(Rx), abs.(Ry), abs.(RP)))
                # push!(iter_evo, iter / max(nx, ny))
                append!(errs_evo, errs)
                @printf(
                    "  iter/nx=%.3f,errs=[ %1.3e, %1.3e, %1.3e ] \n",
                    iter / max(nx, ny),
                    errs...
                )
            end
            iter += 1
        end
        t += dt
        push!(evo_t, t)
        push!(evo_τxx, maximum(τxx))
        push!(iter_evo, iter-1)

        # visualisation
        p1 = heatmap(xc, yc, Array(Pt)'   , title="Pressure"; plot_opts()...)
        p3 = heatmap(xc, yc, Array(TII)'  , title="TII"     ; plot_opts()...)
        p2 = heatmap(xc, yc, Array(η_vep)', title="η_vep"   ; plot_opts()...)
        p4 = plot(evo_t,evo_τxx,legend=false,xlabel="time",ylabel="max(τxx)",linewidth=0,markershape=:circle,markersize=3,framestyle=:box)
        display(plot(p1,p2,p3,p4,layout=(2,2)))
        # png(plot(p1,p2,p3,p4,layout=(2,2)),@sprintf("anim/%04d.png",iframe+=1))
    end
    return 
end

nt = 12 # number of time steps
@time  main(nt);
