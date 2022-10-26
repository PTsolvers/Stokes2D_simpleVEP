# using Pkg; 
# Pkg.activate("C:\\Users\\albert\\Desktop\\GeoParams-dev\\blob\\GeoParams.jl")

# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, GeoParams
const Prec = Float64 # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end]) 
@generated function phase_viscosity(v::NTuple{N,Any}, ε̇ii, phase, args) where N
    quote
        Base.@_inline_meta
        Base.@nexprs $N i -> v[i].Phase === phase && return computeViscosity_εII(v[i].CompositeRheology[1], ε̇ii, args)
        return 0.0 
    end
end

@inline function av_c(A, i, j)
    return 0.25*(A[i,j] + A[i+1,j] + A[i,j+1] + A[i+1,j+1])
end

@inline function av_v(A, i, j)
    return 0.25*(A[i,j] + A[i-1,j] + A[i,j-1] + A[i-1,j-1])
end

# Rheology
function UpdateStressGeoParams!( ηc, ηv, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt, Phasec, Phasev )
    # Centroids
    @inbounds for j ∈ axes(ε̇xx,2), i ∈ axes(ε̇xx,1)
        τii0     =  sqrt(0.5 *(τxx0[i,j]^2 + τyy0[i,j]^2) + av_c(τxy0,i,j)^2)
        ε̇ii      =  sqrt(0.5 *( ε̇xx[i,j]^2 +  ε̇yy[i,j]^2) + av_c(ε̇xy,i,j)^2)
        args     = (; τII_old = τii0, dt=Δt)             
        η = ηc[i,j] = phase_viscosity(MatParam, ε̇ii, Phasec[i,j], args)
        τxx[i,j] = 2.0*η*ε̇xx[i,j]324
        τyy[i,j] = 2.0*η*ε̇yy[i,j]
    end
    # Vertices
    @inbounds for j ∈ 2:size(ε̇xy,2)-1, i ∈ 2:size(ε̇xy,1)-1
        τii0     = sqrt(0.5 * (av_v(τxx0,i,j)^2 + av_v(τyy0,i,j)^2) + τxy0[i,j]^2)
        ε̇ii      = sqrt(0.5 * (av_v(ε̇xx,i,j)^2 + av_v(ε̇yy,i,j)^2) + ε̇xy[i,j]^2)
        args     = (; τII_old = τii0, dt=Δt)
        η = ηv[i,j]  = phase_viscosity(MatParam, ε̇ii, Phasev[i,j], args)
        τxy[i,j] = 2.0*η*ε̇xy[i,j] 
    end
end

function viscosityGeoParams!(ηc, ηv, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt, Phasec, Phasev )
    # Centroids
    @inbounds for j ∈ axes(ε̇xx,2), i ∈ axes(ε̇xx,1)
        τii0     =  sqrt(0.5 *(τxx0[i,j]^2 + τyy0[i,j]^2) + av_c(τxy0,i,j)^2)
        ε̇ii      =  sqrt(0.5 *( ε̇xx[i,j]^2 +  ε̇yy[i,j]^2) + av_c(ε̇xy,i,j)^2)
        args     = (; τII_old = τii0, dt=Δt)             
        ηc[i,j] = phase_viscosity(MatParam, ε̇ii, Phasec[i,j], args)
    end
    # Vertices
    @inbounds for j ∈ 2:size(ε̇xy,2)-1, i ∈ 2:size(ε̇xy,1)-1
        τii0     = sqrt(0.5 * (av_v(τxx0,i,j)^2 + av_v(τyy0,i,j)^2) + τxy0[i,j]^2)
        ε̇ii      = sqrt(0.5 * (av_v(ε̇xx,i,j)^2 + av_v(ε̇yy,i,j)^2) + ε̇xy[i,j]^2)
        args     = (; τII_old = τii0, dt=Δt)
        ηv[i,j]  = phase_viscosity(MatParam, ε̇ii, Phasev[i,j], args)
    end
end
# 2D Stokes routine
@views function Stokes2D_VE_inclusion(UseGeoParams)
    # Physics
    Lx, Ly  = 1.0, 1.0  # domain size
    ξ       = 10.0      # Maxwell relaxation time
    η0      = 1.0       # viscous viscosity
    G       = 1.0       # elastic shear modulus
    εbg     = 1.0       # background strain-rate
    radi    = 0.01

    MatParam = (SetMaterialParams(Name="Matrix"   , Phase=1,
                CompositeRheology = CompositeRheology(ConstantElasticity(G=G),LinearViscous(η=η0))), 
                SetMaterialParams(Name="Inclusion", Phase=2,
                CompositeRheology = CompositeRheology(ConstantElasticity(G=G/6),LinearViscous(η=η0))),
                )

    # Numerics
    nt       = 10        # number of time steps
    ncx, ncy = 31, 31    # numerical grid resolution
    ε        = 1e-6      # nonlinear tolerence
    iterMax  = 1e4       # max number of iters
    nout     = 500       # check frequency
    # Iterative parameters -------------------------------------------
    Reopt    = 5π
    cfl      = 0.50
    ρ        = cfl*Reopt/ncx
    # Preprocessing
    Δx, Δy   = Lx/ncx, Ly/ncy
    Δt       = η0/(G*ξ + 1e-15) 
    # Array initialisation
    Pt       = zeros(Prec, ncx  ,ncy  )
    ∇V       = zeros(Prec, ncx  ,ncy  )
    Vx       = zeros(Prec, ncx+1,ncy+2)
    Vy       = zeros(Prec, ncx+2,ncy+1)
    ε̇xx      = zeros(Prec, ncx  ,ncy  )
    ε̇yy      = zeros(Prec, ncx  ,ncy  )
    ε̇xy      = zeros(Prec, ncx+1,ncy+1)    
    τxx      = zeros(Prec, ncx  ,ncy  )
    τyy      = zeros(Prec, ncx  ,ncy  )
    τxy      = zeros(Prec, ncx+1,ncy+1)
    τxx0     = zeros(Prec, ncx  ,ncy  )
    τyy0     = zeros(Prec, ncx  ,ncy  )
    τxy0     = zeros(Prec, ncx+1,ncy+1)
    Rx       = zeros(Prec, ncx+1,ncy  )
    Ry       = zeros(Prec, ncx  ,ncy+1)
    Rp       = zeros(Prec, ncx  ,ncy  )
    dVxdτ    = zeros(Prec, ncx+1,ncy  ) 
    dVydτ    = zeros(Prec, ncx  ,ncy+1)
    dPdτ     = zeros(Prec, ncx  ,ncy  )
    Δτv      = zeros(Prec, ncx+1,ncy+1)       
    Δτvx     = zeros(Prec, ncx+1,ncy  )         
    Δτvy     = zeros(Prec, ncx  ,ncy+1)         
    κΔτp     = zeros(Prec, ncx  ,ncy  )           
    Rog      = zeros(Prec, ncx  ,ncy  )
    ηc       = η0*ones(Prec, ncx, ncy)
    ηv       = η0*ones(Prec, ncx+1, ncy+1)
    Phasec   = ones(Int, ncx  ,ncy  )
    Phasev   = ones(Int, ncx+1,ncy+1)
    # For geoparams
    τii0c    = zeros(size(ε̇xx))
    τii0v    = zeros(size(ε̇xy))
    ε̇iic     = zeros(size(ε̇xx))
    ε̇iiv     = zeros(size(ε̇xy))
    # For non-geoparams version
    η, G     = 1.0, 1.0
    ηe_c     = Δt*G.*ones(Prec, ncx  ,ncy  )
    ηe_v     = Δt*G.*ones(Prec, ncx+1,ncy+1)
    ηve_c    = zeros(Prec, ncx  ,ncy  )
    ηve_v    = zeros(Prec, ncx+1,ncy+1)
    # Initialisation
    xce, yce  = LinRange(-Δx/2, Lx+Δx/2, ncx+2), LinRange(-Δy/2, Ly+Δy/2, ncy+2)
    xc, yc   = LinRange(Δx/2, Lx-Δx/2, ncx), LinRange(Δy/2, Ly-Δy/2, ncy)
    xv, yv   = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    radc     = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv     = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    # For non-geoparams version
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    ηe_c[radc.<radi]   .= Δt*G/6.0
    ηe_v[radv.<radi]   .= Δt*G/6.0
    ηve_c              .= (1.0./ηe_c .+ 1.0./η).^(-1)
    ηve_v              .= (1.0./ηe_v .+ 1.0./η).^(-1)
    ηc                 .= ηve_c
    ηv                 .= ηve_v
    # Velocity
    (Xvx,Yvx) = ([x for x=xv,y=yce], [y for x=xv,y=yce])
    (Xvy,Yvy) = ([x for x=xce,y=yv], [y for x=xce,y=yv])
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=Float64[]; evo_τxx=Float64[];
    global itg = 1
    for it = 1:15
        iter=1; err=2*ε; err_evo1=Float64[]; err_evo2=Float64[]; 
        τxx0.=τxx; τyy0.=τyy; τxy0.=τxy
        while (err>ε && iter<=iterMax)
            # BCs
            Vx[:,1]   .= Vx[:,2]     # S
            Vx[:,end] .= Vx[:,end-1] # N
            Vy[1,:]   .= Vy[2,:]     # W
            Vy[end,:] .= Vy[end-1,:] # E
            # Kinematics
            ∇V    .= diff(Vx[:,2:end-1], dims=1)./Δx .+ diff(Vy[2:end-1,:], dims=2)./Δy
            ε̇xx   .= diff(Vx[:,2:end-1], dims=1)./Δx .- 1.0/3.0*∇V
            ε̇yy   .= diff(Vy[2:end-1,:], dims=2)./Δy .- 1.0/3.0*∇V
            ε̇xy   .= 0.5.*(diff(Vx, dims=2)./Δy .+ diff(Vy, dims=1)./Δx)  
            # Stresses
            if UseGeoParams
                iter == 1 && viscosityGeoParams!(ηc, ηv, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt, Phasec, Phasev )
                τxx   .= 2 .* ηc .* ε̇xx 
                τyy   .= 2 .* ηc .* ε̇yy
                τxy   .= 2 .* ηv .* ε̇xy
                # UpdateStressGeoParams!( ηc, ηv, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt, Phasec, Phasev )
            else
                τxx   .= 2 .* ηve_c .* ( ε̇xx .+ τxx0./(2 .* ηe_c) ) 
                τyy   .= 2 .* ηve_c .* ( ε̇yy .+ τyy0./(2 .* ηe_c) )
                τxy   .= 2 .* ηve_v .* ( ε̇xy .+ τxy0./(2 .* ηe_v) )
            end
            # Residuals
            Rx[2:end-1,:] .= .-diff(Pt, dims=1)./Δx .+ diff(τxx, dims=1)./Δx .+ diff(τxy[2:end-1,:], dims=2)./Δy
            Ry[:,2:end-1] .= .-diff(Pt, dims=2)./Δy .+ diff(τyy, dims=2)./Δy .+ diff(τxy[:,2:end-1], dims=1)./Δx #.+ av_ya(Rog)
            Rp            .= .-∇V
            # PT time step -----------------------------------------------
            Δτv  .= ρ*min(Δx,Δy)^2 ./ ηv ./ 4.1 * cfl   
            @. Δτvx = (Δτv[:,1:end-1] + Δτv[:,2:end]) * 0.5
            @. Δτvy = (Δτv[1:end-1,:] + Δτv[2:end,:]) * 0.5
            @. κΔτp = cfl * ηc * Δx / Lx  
            # Calculate rate update --------------------------------------
            dVxdτ          .= (1-ρ) .* dVxdτ .+ Rx
            dVydτ          .= (1-ρ) .* dVydτ .+ Ry
            dPdτ           .= Rp
            # Update velocity and pressure -------------------------------
            Vx[:,2:end-1]  .+= Δτvx ./ ρ .* dVxdτ
            Vy[2:end-1,:]  .+= Δτvy ./ ρ .* dVydτ
            Pt             .+= κΔτp .* dPdτ
            # convergence check
            if mod(iter, nout)==0 || iter==1
                norm_Rx = norm(Rx)/sqrt(length(Rx)); norm_Ry = norm(Ry)/sqrt(length(Ry)); norm_∇V = norm(∇V)/sqrt(length(∇V))
                err = max(maximum(norm_Rx), maximum(norm_Ry), maximum(norm_∇V))
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %03d, iter = %04d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, ∇V=%1.3e] \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V)
            end
            iter+=1; global itg=iter
        end
        t = t + Δt
        # push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        # Plotting
        # p1 = heatmap(xv, yc, Vx[:,2:end-1]', aspect_ratio=1, xlims=(0, Lx), ylims=(Δy/2, Ly-Δy/2), c=:inferno, title="Vx")
        # p2 = heatmap(xc, yv, Vy[2:end-1,:]', aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        # p3 = heatmap(xc, yc, Pt' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="P")
        # # p3 = heatmap(xv, yv, τxy' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="τxy")
        # p4 = plot(evo_t, evo_τxx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
        # p4 = plot!(evo_t, 2.0.*εbg.*η0.*(1.0.-exp.(.-evo_t.*G./η0)), linewidth=2.0) # analytical solution
        # display(plot(p1, p2, p3, p4))
    end
    return
end

@time Stokes2D_VE_inclusion(false)
# 2.701631 seconds (5.39 M allocations: 3.203 GiB, 11.96% gc time) # 15 its
# 0.873773 seconds (2.76 M allocations: 1.641 GiB,  8.02% gc time) # 5 its

@time Stokes2D_VE_inclusion(true)
# 10.787914 seconds (7.92 M allocations: 5.273 GiB, 2.12% gc time) # 15 it
#  3.468200 seconds (2.44 M allocations: 1.629 GiB, 2.19% gc time) # 5 its