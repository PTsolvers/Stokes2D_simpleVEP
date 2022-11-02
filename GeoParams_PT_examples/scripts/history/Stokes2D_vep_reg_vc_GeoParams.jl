# Initialisation
using Plots, Printf, Statistics, LinearAlgebra
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
# Rheology
function UpdateStressVEVP!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xyc, ε̇xxv, ε̇yyv, ε̇xy, τxx0, τyy0, τxy0, τxxv0, τyyv0, τxyv0, ε̇xx1, ε̇yy1, ε̇xy1, ε̇xxv1, ε̇yyv1, ε̇xyv1, τii, τiiv, Pt, Ptv, τ_y, sinϕ, η_ve, η_vev, η_reg, dQdτxx, dQdτyy, dQdτxy, dQdτxxv, dQdτyyv, dQdτxyv, η_e, η_ev, ε̇ii, ε̇iiv )
    # visco-elastic strain rates
    ε̇xx1   .=    ε̇xx   .+ τxx0 ./2.0./η_e
    ε̇yy1   .=    ε̇yy   .+ τyy0 ./2.0./η_e  
    ε̇xy1   .=    ε̇xyc  .+ τxy0 ./2.0./η_e
    ε̇ii    .= sqrt.(0.5*(ε̇xx1.^2 .+ ε̇yy1.^2) .+ ε̇xy1.^2)
    # visco-elastic strain rates vertices
    ε̇xxv1  .=    ε̇xxv  .+ τxxv0./2.0./η_ev
    ε̇yyv1  .=    ε̇yyv  .+ τyyv0./2.0./η_ev
    ε̇xyv1  .=    ε̇xyv  .+ τxyv0./2.0./η_ev
    ε̇iiv   .= sqrt.(0.5*(ε̇xxv1.^2 .+ ε̇yyv1.^2) .+ ε̇xyv1.^2)
    # trial stress
    τxx    .= 2.0.*η_ve.*ε̇xx1
    τyy    .= 2.0.*η_ve.*ε̇yy1
    τxy    .= 2.0.*η_ve.*ε̇xy1
    τii    .= sqrt.(0.5*(τxx.^2 .+ τyy.^2) .+ τxy.^2)
    # trial stress vertices
    τxxv   .= 2.0.*η_vev.*ε̇xxv1
    τyyv   .= 2.0.*η_vev.*ε̇yyv1
    τxyv   .= 2.0.*η_vev.*ε̇xyv1
    τiiv   .= sqrt.(0.5*(τxxv.^2 .+ τyyv.^2) .+ τxyv.^2)
    # yield function
    F      .= τii .- τ_y .- Pt.*sinϕ
    Pla    .= 0.0
    Pla    .= F .> 0.0
    λ      .= Pla.*F./(η_ve .+ η_reg)
    dQdτxx .= 0.5.*τxx./τii
    dQdτyy .= 0.5.*τyy./τii
    dQdτxy .=      τxy./τii
    # yield function vertices
    Fv     .= τiiv .- τ_y .- Ptv.*sinϕ
    Plav   .= 0.0
    Plav   .= Fv .> 0.0
    λv     .= Plav.*Fv./(η_vev .+ η_reg)
    dQdτxxv.= 0.5.*τxxv./τiiv
    dQdτyyv.= 0.5.*τyyv./τiiv
    dQdτxyv.=      τxyv./τiiv
    # plastic corrections
    τxx    .= 2.0.*η_ve.*(ε̇xx1 .-      λ.*dQdτxx)
    τyy    .= 2.0.*η_ve.*(ε̇yy1 .-      λ.*dQdτyy)
    τxy    .= 2.0.*η_ve.*(ε̇xy1 .- 0.5.*λ.*dQdτxy)
    τii    .= sqrt.(0.5*(τxx.^2 .+ τyy.^2) .+ τxy.^2)
    Fchk   .= τii .- τ_y .- Pt.*sinϕ .- λ.*η_reg
    η_vep  .= τii./2.0./ε̇ii
    # plastic corrections vertices
    τxxv   .= 2.0.*η_vev.*(ε̇xxv1 .-      λv.*dQdτxxv)
    τyyv   .= 2.0.*η_vev.*(ε̇yyv1 .-      λv.*dQdτyyv)
    τxyv   .= 2.0.*η_vev.*(ε̇xyv1 .- 0.5.*λv.*dQdτxyv)
    τiiv   .= sqrt.(0.5*(τxxv.^2 .+ τyyv.^2) .+ τxyv.^2)
    Fchkv  .= τiiv .- τ_y .- Ptv.*sinϕ .- λv.*η_reg
    η_vepv .= τiiv./2.0./ε̇iiv
end
# 2D Stokes routine 
@views function Stokes2D_vep()
    do_DP   = true               # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg   = 1.2e-2             # regularisation "viscosity"
    # Physics
    Lx, Ly  = 1.0, 1.0           # domain size
    radi    = 0.01               # inclusion radius
    τ_y     = 1.6                # yield stress. If do_DP=true, τ_y stand for the cohesion: c*cos(ϕ)
    sinϕ    = sind(30)*do_DP     # sinus of the friction angle
    η0      = 1.0                # viscous viscosity
    G0      = 1.0                # elastic shear modulus
    Gi      = G0/(6.0-4.0*do_DP) # elastic shear modulus perturbation
    εbg     = 1.0                # background strain-rate
    # Numerics
    nt      = 10                 # number of time steps
    ncx, ncy  = 63, 63             # numerical grid resolution
    Vdmp    = 4.0                # convergence acceleration (damping)
    Vsc     = 2.0                # iterative time step limiter
    Ptsc    = 6.0                # iterative time step limiter
    ε       = 1e-6               # nonlinear tolerence
    iterMax = 3e4                # max number of iters
    nout    = 200                # check frequency
    # Preprocessing
    dx, dy  = Lx/ncx, Ly/ncy 
    dt      = η0/G0/4.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = zeros(Dat, ncx  ,ncy  )
    ∇V      = zeros(Dat, ncx  ,ncy  )
    Vx      = zeros(Dat, ncx+1,ncy  )
    Vy      = zeros(Dat, ncx  ,ncy+1)
    ε̇xx     = zeros(Dat, ncx  ,ncy  )
    ε̇yy     = zeros(Dat, ncx  ,ncy  )
    ε̇xyv    = zeros(Dat, ncx+1,ncy+1)
    ε̇xy     = zeros(Dat, ncx  ,ncy  )
    ε̇xx1    = zeros(Dat, ncx  ,ncy  )
    ε̇yy1    = zeros(Dat, ncx  ,ncy  )
    ε̇xy1    = zeros(Dat, ncx  ,ncy  )
    ε̇xyv1   = zeros(Dat, ncx+1,ncy+1)
    τxx     = zeros(Dat, ncx  ,ncy  )
    τyy     = zeros(Dat, ncx  ,ncy  )
    τxy     = zeros(Dat, ncx  ,ncy  )
    τxyv    = zeros(Dat, ncx+1,ncy+1)
    τxx0    = zeros(Dat, ncx  ,ncy  )
    τyy0   = zeros(Dat, ncx  ,ncy  )
    τxy0   = zeros(Dat, ncx  ,ncy  )
    τxyv0  = zeros(Dat, ncx+1,ncy+1)
    # for vertices implementation
    Ptv     = zeros(Dat, ncx+1,ncy+1)
    ε̇xxv    = zeros(Dat, ncx+1,ncy+1)
    ε̇yyv    = zeros(Dat, ncx+1,ncy+1)
    ε̇xxv1   = zeros(Dat, ncx+1,ncy+1)
    ε̇yyv1   = zeros(Dat, ncx+1,ncy+1)
    τxxv    = zeros(Dat, ncx+1,ncy+1)
    τyyv    = zeros(Dat, ncx+1,ncy+1)
    τxxv0  = zeros(Dat, ncx+1,ncy+1)
    τyyv0  = zeros(Dat, ncx+1,ncy+1)
    Fchkv   = zeros(Dat, ncx+1,ncy+1)
    Fv      = zeros(Dat, ncx+1,ncy+1)
    Plav    = zeros(Dat, ncx+1,ncy+1)
    λv      = zeros(Dat, ncx+1,ncy+1)
    dQdτxxv = zeros(Dat, ncx+1,ncy+1)
    dQdτyyv = zeros(Dat, ncx+1,ncy+1)
    dQdτxyv = zeros(Dat, ncx+1,ncy+1)
    τiiv    = zeros(Dat, ncx+1,ncy+1)
    ε̇iiv    = zeros(Dat, ncx+1,ncy+1)
    # for vertices implementation
    τii     = zeros(Dat, ncx  ,ncy  )
    ε̇ii     = zeros(Dat, ncx  ,ncy  )
    F       = zeros(Dat, ncx  ,ncy  )
    Fchk    = zeros(Dat, ncx  ,ncy  )
    Pla     = zeros(Dat, ncx  ,ncy  )
    λ       = zeros(Dat, ncx  ,ncy  )
    dQdτxx  = zeros(Dat, ncx  ,ncy  )
    dQdτyy  = zeros(Dat, ncx  ,ncy  )
    dQdτxy  = zeros(Dat, ncx  ,ncy  )
    Rx      = zeros(Dat, ncx-1,ncy  )
    Ry      = zeros(Dat, ncx  ,ncy-1)
    dVxdt   = zeros(Dat, ncx-1,ncy  )
    dVydt   = zeros(Dat, ncx  ,ncy-1)
    dtPt    = zeros(Dat, ncx  ,ncy  )
    dtVx    = zeros(Dat, ncx-1,ncy  )
    dtVy    = zeros(Dat, ncx  ,ncy-1)
    Rog     = zeros(Dat, ncx  ,ncy  )
    η_v     =    η0*ones(Dat, ncx  ,ncy  )
    η_e     = dt*G0*ones(Dat, ncx  ,ncy  )
    η_ev    = dt*G0*ones(Dat, ncx+1,ncy+1)
    η_ve    =       ones(Dat, ncx  ,ncy  )
    η_vep   =       ones(Dat, ncx  ,ncy  )
    η_vepv  =       ones(Dat, ncx+1,ncy+1)
    η_vev   =       ones(Dat, ncx+1,ncy+1) # for vertices implementation
    η_vv    =    η0*ones(Dat, ncx+1,ncy+1) # for vertices implementation
    # Initial condition
    xc, yc  = LinRange(dx/2, Lx-dx/2, ncx), LinRange(dy/2, Ly-dy/2, ncy)
    xc, yc  = LinRange(dx/2, Lx-dx/2, ncx), LinRange(dy/2, Ly-dy/2, ncy)
    xv, yv  = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    (Xvx,Yvx) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (Xvy,Yvy) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    radc      = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv      = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    Phasec    = ones(Int, ncx  ,ncy  )
    Phasev    = ones(Int, ncx+1,ncy+1)
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    η_e[radc.<radi] .= dt*Gi
    η_ev[radv.<radi].= dt*Gi
    η_ve   .= (1.0./η_e  + 1.0./η_v).^-1
    η_vev  .= (1.0./η_ev + 1.0./η_vv).^-1
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=[]; evo_τxx=[]
    for it = 1:nt
        iter=1; err=2*ε; err_evo1=[]; err_evo2=[]
        τxx0.=τxx;   τyy0.=τyy;   τxy0.=τxy;   λ.=0.0
        τxxv0.=τxxv; τyyv0.=τyyv; τxyv0.=τxyv; λv.=0.0
        local itg
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V     .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Pt     .= Pt .- dtPt.*∇V
            # strain rates
            ε̇xx    .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            ε̇yy    .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            ε̇xyv[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            ε̇xxv[2:end-1,2:end-1] .= av(ε̇xx); ε̇xxv[1,:].=ε̇xxv[2,:]; ε̇xxv[end,:].=ε̇xxv[end-1,:]; ε̇xxv[:,1].=ε̇xxv[:,2]; ε̇xxv[:,end].=ε̇xxv[:,end-1]
            ε̇yyv[2:end-1,2:end-1] .= av(ε̇yy); ε̇yyv[1,:].=ε̇yyv[2,:]; ε̇yyv[end,:].=ε̇yyv[end-1,:]; ε̇yyv[:,1].=ε̇yyv[:,2]; ε̇yyv[:,end].=ε̇yyv[:,end-1]
            Ptv[2:end-1,2:end-1]  .= av(Pt);   Ptv[1,:].= Ptv[2,:];  Ptv[end,:].= Ptv[end-1,:];  Ptv[:,1].= Ptv[:,2];  Ptv[:,end].= Ptv[:,end-1]
            ε̇xy .= av(ε̇xyv)
            # Rheology
            UpdateStressVEVP!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, ε̇xxv, ε̇yyv, ε̇xy, τxx0, τyy0, τxy0, τxxv0, τyyv0, τxyv0, ε̇xx1, ε̇yy1, ε̇xy1, ε̇xxv1, ε̇yyv1, ε̇xyv1, τii, τiiv, Pt, Ptv, τ_y, sinϕ, η_ve, η_vev, η_reg, dQdτxx, dQdτyy, dQdτxy, dQdτxxv, dQdτyyv, dQdτxyv, η_e, η_ev, ε̇ii, ε̇iiv )
            # PT timestep
            dtVx   .= min(dx,dy)^2.0./av_xa(η_vep)./4.1./Vsc
            dtVy   .= min(dx,dy)^2.0./av_ya(η_vep)./4.1./Vsc
            dtPt   .= 4.1.*η_vep./max(ncx,ncy)./Ptsc 
            # velocities
            Rx     .= .-diff(Pt, dims=1)./dx .+ diff(τxx, dims=1)./dx .+ diff(τxyv[2:end-1,:], dims=2)./dy
            Ry     .= .-diff(Pt, dims=2)./dy .+ diff(τyy, dims=2)./dy .+ diff(τxyv[:,2:end-1], dims=1)./dx .+ av_ya(Rog)
            dVxdt  .= dVxdt.*(1-Vdmp/ncx) .+ Rx
            dVydt  .= dVydt.*(1-Vdmp/ncy) .+ Ry
            Vx[2:end-1,:] .= Vx[2:end-1,:] .+ dVxdt.*dtVx
            Vy[:,2:end-1] .= Vy[:,2:end-1] .+ dVydt.*dtVy
            # convergence check
            if mod(iter, nout)==0
                norm_Rx = norm(Rx)/length(Rx); norm_Ry = norm(Ry)/length(Ry); norm_∇V = norm(∇V)/length(∇V)
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %d, iter = %d, err = %1.2e norm[Rx=%1.2e, Ry=%1.2e, ∇V=%1.2e] (Fchk=%1.2e - Fchkv=%1.2e) \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V, maximum(Fchk), maximum(Fchkv))
            end
            iter+=1; itg=iter
        end
        t = t + dt
        push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        # Plotting
        p1 = heatmap(xv, yc, Vx' , aspect_ratio=1, xlims=(0, Lx), ylims=(dy/2, Ly-dy/2), c=:inferno, title="Vx")
        # p2 = heatmap(xc, yv, Vy' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        p2 = heatmap(xc, yc, η_vep' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="η_vep")
        p3 = heatmap(xc, yc, τii' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="τii")
        p4 = plot(evo_t, evo_τxx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
            plot!(evo_t, 2.0.*εbg.*η0.*(1.0.-exp.(.-evo_t.*G0./η0)), linewidth=2.0) # analytical solution for VE loading
            plot!(evo_t, 2.0.*εbg.*η0.*ones(size(evo_t)), linewidth=2.0)            # viscous flow stress
            if !do_DP plot!(evo_t, τ_y*ones(size(evo_t)), linewidth=2.0) end        # von Mises yield stress
        display(plot(p1, p2, p3, p4))
    end
    return
end

Stokes2D_vep()
