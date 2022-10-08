# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, GeoParams
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end]) 
# Rheology
@views function UpdateStressGeoParams!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, dt, Phasec, Phasev )
    ε̇iic  = ones(size(ε̇xx))
    ε̇iiv  = ones(size(ε̇xy))
    τii0c = sqrt.(1//2 .*(τxx0.^2 .+ τyy0.^2) .+ av(τxy0).^2)
    τii0v = zeros(size(ε̇xy))
    τii0v[2:end-1,2:end-1] = sqrt.(1//2 .*( av(τxx0).^2 .+ av(τyy0).^2) .+ τxy0[2:end-1,2:end-1].^2)
    # Centroids
    for i in eachindex(τxx)
        v      = MatParam[Phasec[i]].CompositeRheology[1]
        args   = (; τII_old = τii0c[i], dt=dt)
        η_eff  = computeViscosity_εII(v, ε̇iic[i], args)
        τxx[i] = 2*η_eff*ε̇xx[i]
        τyy[i] = 2*η_eff*ε̇yy[i]
    end
    # Vertices
    for i in eachindex(τxy)
        v      = MatParam[Phasev[i]].CompositeRheology[1]
        args   = (; τII_old = τii0v[i], dt=dt)
        η_eff  = computeViscosity_εII(v, ε̇iiv[i], args)
        τxy[i] = 2*η_eff*ε̇xy[i] 
    end
end
# 2D Stokes routine
@views function Stokes2D_ve()
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
                CompositeRheology = CompositeRheology(ConstantElasticity(G=G*0.25),LinearViscous(η=η0))),
                )

    # Numerics
    nt       = 1       # number of time steps
    ncx, ncy = 31, 31    # numerical grid resolution
    Vdmp     = 9.0       # convergence acceleration (damping)
    Ptsc     = 10.0       # iterative time step limiter
    ε        = 1e-6      # nonlinear tolerence
    iterMax  = 1e5       # max number of iters
    nout     = 500       # check frequency
    # Preprocessing
    dx, dy  = Lx/ncx, Ly/ncy
    dt      = η0/(G*ξ + 1e-15)
    # Array initialisation
    Pt      = zeros(Dat, ncx  ,ncy  )
    ∇V      = zeros(Dat, ncx  ,ncy  )
    Vx      = zeros(Dat, ncx+1,ncy  )
    Vy      = zeros(Dat, ncx  ,ncy+1)
    ε̇xx     = zeros(Dat, ncx  ,ncy  )
    ε̇yy     = zeros(Dat, ncx  ,ncy  )
    ε̇xy     = zeros(Dat, ncx+1,ncy+1)    
    τxx     = zeros(Dat, ncx  ,ncy  )
    τyy     = zeros(Dat, ncx  ,ncy  )
    τxy     = zeros(Dat, ncx+1,ncy+1)
    τxx0    = zeros(Dat, ncx  ,ncy  )
    τyy0    = zeros(Dat, ncx  ,ncy  )
    τxy0    = zeros(Dat, ncx+1,ncy+1)
    Rx      = zeros(Dat, ncx-1,ncy  )
    Ry      = zeros(Dat, ncx  ,ncy-1)
    dVxdt   = zeros(Dat, ncx-1,ncy  )
    dVydt   = zeros(Dat, ncx  ,ncy-1)
    Rog     = zeros(Dat, ncx  ,ncy  )
    Xc      =  ξ*ones(Dat, ncx, ncy)
    Xv      =  ξ*ones(Dat, ncx+1, ncy+1)
    ηc      = η0*ones(Dat, ncx, ncy)
    ηv      = η0*ones(Dat, ncx+1, ncy+1)
    Phasec    = ones(Int, ncx  ,ncy  )
    Phasev    = ones(Int, ncx+1,ncy+1)
    xc, yc  = LinRange(dx/2, Lx-dx/2, ncx), LinRange(dy/2, Ly-dy/2, ncy)
    xv, yv  = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    radc      = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv      = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    # Initialisation
    xc, yc  = LinRange(dx/2, Lx-dx/2, ncx), LinRange(dy/2, Ly-dy/2, ncy)
    xv, yv  = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    (Xvx,Yvx) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (Xvy,Yvy) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    dtVx    = min(dx,dy)^2.0./av_xa(ηc)./4.1
    dtVy    = min(dx,dy)^2.0./av_ya(ηc)./4.1
    dtPt    = 4.1*ηc/max(ncx,ncy)/Ptsc
    # Time loop
    t=0.0; evo_t=[]; evo_τxx=[]
    for it = 1:nt
        iter=1; err=2*ε; err_evo1=[]; err_evo2=[]; 
        τxx0.=τxx; τyy0.=τyy; τxy0.=τxy
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V    .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Pt    .= Pt .- dtPt.*∇V
            # strain rates
            ε̇xx   .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            ε̇yy   .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            ε̇xy[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            # stresses
            UpdateStressGeoParams!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, dt, Phasec, Phasev )
            # velocities
            Rx    .= .-diff(Pt, dims=1)./dx .+ diff(τxx, dims=1)./dx .+ diff(τxy[2:end-1,:], dims=2)./dy
            Ry    .= .-diff(Pt, dims=2)./dy .+ diff(τyy, dims=2)./dy .+ diff(τxy[:,2:end-1], dims=1)./dx .+ av_ya(Rog)
            dVxdt .= dVxdt.*(1-Vdmp/ncx) .+ Rx
            dVydt .= dVydt.*(1-Vdmp/ncy) .+ Ry
            Vx[2:end-1,:] .= Vx[2:end-1,:] .+ dVxdt.*dtVx
            Vy[:,2:end-1] .= Vy[:,2:end-1] .+ dVydt.*dtVy
            # convergence check
            if mod(iter, nout)==0
                norm_Rx = norm(Rx)/length(Rx); norm_Ry = norm(Ry)/length(Ry); norm_∇V = norm(∇V)/length(∇V)
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %d, iter = %d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, ∇V=%1.3e] \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V)
            end
            iter+=1; global itg=iter
        end
        t = t + dt
        push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        # Plotting
        p1 = heatmap(xv, yc, Vx' , aspect_ratio=1, xlims=(0, Lx), ylims=(dy/2, Ly-dy/2), c=:inferno, title="Vx")
        p2 = heatmap(xc, yv, Vy' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        p3 = plot(evo_t, evo_τxx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
            plot!(evo_t, 2.0.*εbg.*η0.*(1.0.-exp.(.-evo_t.*G./η0)), linewidth=2.0) # analytical solution
        display(plot(p1, p2, p3))
    end
    return
end

Stokes2D_ve()
