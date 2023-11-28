# This script implements plasticity with and w/out invariants & shows how to use GeoParams to do the local stress update
# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, GeoParams, TimerOutputs
const Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views   av2(A) = 0.25*(A[1:end-1,1:end-1].^2 .+ A[2:end,1:end-1].^2 .+ A[1:end-1,2:end].^2 .+ A[2:end,2:end].^2)
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])

# 2D Stokes routine
@views function Stokes2D_vep(doPlot, nt)
    CharDim = GEO_units(length=100km)
    do_DP   = true               # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg   = 1.5e20Pas            # regularisation "viscosity"
    # Physics
    Lx, Ly  = 100km,50km         # domain size
    radi    = 1km                # inclusion radius
    τ_y     = 35MPa              # yield stress. If do_DP=true, τ_y stand for the cohesion: c*cos(ϕ)
    ϕ       = 30.0*do_DP
    μ0      = 1e23Pas            # viscous viscosity
    G0      = 1e10Pa             # elastic shear modulus
    Gi      = G0/(6.0-4.0*do_DP) # elastic shear modulus perturbation
    εbg     = 1e-15/s            # background strain-rate
    Coh     = τ_y/cosd(ϕ)        # cohesion
    # Coh     = Inf

    # Geoparams initialisation
    pl = DruckerPrager_regularised(C=Coh, ϕ=ϕ, η_vp=η_reg, Ψ=0)        # non-regularized plasticity
    MatParam = (SetMaterialParams(Name="Matrix"   , Phase=1,
                CompositeRheology = CompositeRheology(SetConstantElasticity(G=G0, ν=0.5),LinearViscous(η=μ0), pl), CharDim=CharDim), 
                SetMaterialParams(Name="Inclusion", Phase=2,
                CompositeRheology = CompositeRheology(SetConstantElasticity(G=Gi, ν=0.5),LinearViscous(η=μ0), pl), CharDim=CharDim),
                )

    # Nondimensionalize  
    Lx, Ly  = nondimensionalize([Lx, Ly], CharDim)                
    radi    = nondimensionalize(radi, CharDim)  
    εbg     = nondimensionalize(εbg, CharDim)  
    μ0      = nondimensionalize(μ0, CharDim)   
    G0      = nondimensionalize(G0, CharDim)   
    Gi      = nondimensionalize(Gi, CharDim)   
    τ_y     = nondimensionalize(τ_y, CharDim)   
    
    # Numerics
    # nt      = 10               # number of time steps
    nx, ny  = 63, 63             # numerical grid resolution
   # nx, ny  = 128, 128           # numerical grid resolution
    
    Vdmp    = 4.0                # convergence acceleration (damping)
    Vsc     = 2.0                # iterative time step limiter
    Ptsc    = 6.0                # iterative time step limiter
    ε       = 1e-6               # nonlinear tolerence
    iterMax = 3e4                # max number of iters
    nout    = 200                # check frequency
    # Preprocessing
    dx, dy  = Lx/nx, Ly/ny
    dt      = μ0/G0/30.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = zeros(Dat, nx  ,ny  )
    P_o     = zeros(Dat, nx,  ny  )
    ∇V      = zeros(Dat, nx  ,ny  )
    Vx      = zeros(Dat, nx+1,ny  )
    Vy      = zeros(Dat, nx  ,ny+1)
    Exx     = zeros(Dat, nx  ,ny  )
    Eyy     = zeros(Dat, nx  ,ny  )
    Exy     = zeros(Dat, nx  ,ny  )
    Exyv    = zeros(Dat, nx+1,ny+1)
    Exx1    = zeros(Dat, nx  ,ny  )
    Eyy1    = zeros(Dat, nx  ,ny  )
    Exy1    = zeros(Dat, nx  ,ny  )
    Exyv1   = zeros(Dat, nx+1,ny+1)
    Txx     = zeros(Dat, nx  ,ny  )
    Tyy     = zeros(Dat, nx  ,ny  )
    Txy     = zeros(Dat, nx  ,ny  )
    Txyv    = zeros(Dat, nx+1,ny+1)
    Txx_o   = zeros(Dat, nx  ,ny  )
    Tii_o   = zeros(Dat, nx  ,ny  )
    Tyy_o   = zeros(Dat, nx  ,ny  )
    Txy_o   = zeros(Dat, nx  ,ny  )
    Txyv_o  = zeros(Dat, nx+1,ny+1)
    Tii     = zeros(Dat, nx  ,ny  )
    Eii     = zeros(Dat, nx  ,ny  )
    Eii_f   = zeros(Dat, nx  ,ny  )
    F       = zeros(Dat, nx  ,ny  )
    Fchk    = zeros(Dat, nx  ,ny  )
    Pla     = zeros(Dat, nx  ,ny  )
    λ       = zeros(Dat, nx  ,ny  )
    dQdTxx  = zeros(Dat, nx  ,ny  )
    dQdTyy  = zeros(Dat, nx  ,ny  )
    dQdTxy  = zeros(Dat, nx  ,ny  )
    Rx      = zeros(Dat, nx-1,ny  )
    Ry      = zeros(Dat, nx  ,ny-1)
    dVxdt   = zeros(Dat, nx-1,ny  )
    dVydt   = zeros(Dat, nx  ,ny-1)
    dtPt    = zeros(Dat, nx  ,ny  )
    dtVx    = zeros(Dat, nx-1,ny  )
    dtVy    = zeros(Dat, nx  ,ny-1)
    Rog     = zeros(Dat, nx  ,ny  )
    η_v     =    μ0*ones(Dat, nx, ny)
    η_vv    =    μ0*ones(Dat, nx+1, ny+1)
    η_e     = dt*G0*ones(Dat, nx, ny)
    η_ev    = dt*G0*ones(Dat, nx+1, ny+1)
    η_ve    =       ones(Dat, nx, ny)
    η_vev    =      ones(Dat, nx+1, ny+1)
    η_vep   =       ones(Dat, nx, ny)
    η_vepv  =       ones(Dat, nx+1, ny+1)
    ηc      = ones(Dat,nx,ny)
    Phasec  = ones(Int, nx  ,ny  )
    Phasev  = ones(Int, nx+1,ny+1)
    # Initial condition
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xv, yv  = LinRange(0.0, Lx, nx+1), LinRange(0.0, Ly, ny+1)
    (Xvx,Yvx) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (Xvy,Yvy) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    radc      = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv      = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    η_e[radc.<radi]    .= dt*Gi
    η_ev[radv.<radi]   .= dt*Gi
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    η_ve    .= (1.0./η_e + 1.0./η_v).^-1
    η_vev   .= (1.0./η_ev + 1.0./η_vv).^-1
    Vx      .=   εbg.*Xvx
    Vy      .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=Float64[]; evo_Txx=Float64[]
    to = TimerOutput()
    for it = 1:nt
        iter=1; err=2*ε; err_evo1=Float64[]; err_evo2=Float64[]
        Txx_o.=Txx; Tyy_o.=Tyy; Txy_o.=av(Txyv); Txyv_o.=Txyv; λ.=0.0; P_o .= Pt
        local itg
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V     .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Pt     .= Pt .- dtPt.*∇V
            # strain rates
            Exx    .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            Eyy    .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            Exyv[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            Exy    .= av(Exyv)
            # Update stresses using GeoParams
            compute_τij_stagcenter!(Txx, Tyy, Txy, Tii, η_vep, Exx, Eyy, Exyv, Pt, Txx_o, Tyy_o, Txyv_o, Phasec, Phasev, MatParam, dt) 
            Txyv[2:end-1,2:end-1].=av(Txy)      # Txyv=0 on boundaries !
            # PT timestep
            dtVx   .= min(dx,dy)^2.0./av_xa(η_vep)./4.1./Vsc
            dtVy   .= min(dx,dy)^2.0./av_ya(η_vep)./4.1./Vsc
            dtPt   .= 4.1.*η_vep./max(nx,ny)./Ptsc
            # velocities
            Rx     .= .-diff(Pt, dims=1)./dx .+ diff(Txx, dims=1)./dx .+ diff(Txyv[2:end-1,:], dims=2)./dy
            Ry     .= .-diff(Pt, dims=2)./dy .+ diff(Tyy, dims=2)./dy .+ diff(Txyv[:,2:end-1], dims=1)./dx .+ av_ya(Rog)
            dVxdt  .= dVxdt.*(1-Vdmp/nx) .+ Rx
            dVydt  .= dVydt.*(1-Vdmp/ny) .+ Ry
            Vx[2:end-1,:] .= Vx[2:end-1,:] .+ dVxdt.*dtVx
            Vy[:,2:end-1] .= Vy[:,2:end-1] .+ dVydt.*dtVy
            # convergence check
            if mod(iter, nout)==0
                norm_Rx = norm(Rx)/length(Rx); norm_Ry = norm(Ry)/length(Ry); norm_∇V = norm(∇V)/length(∇V)
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %d, iter = %d, err = %1.2e norm[Rx=%1.2e, Ry=%1.2e, ∇V=%1.2e] (Fchk=%1.2e) \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V, maximum(Fchk))
            end
            iter+=1; itg=iter
        end
        t = t + dt
        push!(evo_t, t); push!(evo_Txx, maximum(Txx))
        
        if doPlot
            # Plotting
            p1 = heatmap(xv*CharDim.length, yc*CharDim.length, Vx' , aspect_ratio=1, xlims=(0, Lx).*CharDim.length, ylims=(dy/2, Ly-dy/2).*CharDim.length, c=:inferno, title="Vx")
            # p2 = heatmap(xc, yv, Vy' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="Vy")
            p2 = heatmap(xc*CharDim.length, yc*CharDim.length, (η_vep.*CharDim.viscosity)' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2).*CharDim.length, ylims=(0, Ly).*CharDim.length, c=:inferno, title="η_vep")
            p3 = heatmap(xc*CharDim.length, yc*CharDim.length, (Tii.*CharDim.stress)' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2).*CharDim.length, ylims=(0, Ly).*CharDim.length, c=:inferno, title="τii")
            p4 = plot(evo_t.*CharDim.time, evo_Txx.*CharDim.stress , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
            plot!(evo_t.*CharDim.time, CharDim.stress*2.0.*εbg.*μ0.*(1.0.-exp.(.-evo_t.*G0./μ0)), linewidth=2.0) # analytical solution for VE loading
           # plot!(evo_t.*CharDim.time, CharDim.stress*2.0.*εbg.*μ0.*ones(size(evo_t)), linewidth=2.0)            # viscous flow stress
            if !do_DP plot!(evo_t.*CharDim.time, CharDim.stress*τ_y*ones(size(evo_t)), linewidth=2.0) end        # von Mises yield stress
            display(plot(p1, p2, p3, p4))
        end
    end
    return evo_t, evo_Txx, to
end

@time evo_t, evo_Txx_GP, to_GP = Stokes2D_vep(true, 10) # 2nd argument = timesteps  9.368483 seconds (2.73 M allocations: 13.804 GiB, 8.24% gc time)
