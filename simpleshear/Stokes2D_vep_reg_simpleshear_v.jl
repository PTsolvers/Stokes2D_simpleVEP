# Initialisation
using Plots, Printf, Statistics, LinearAlgebra
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
# 2D Stokes routine
@views function Stokes2D_vep()
    do_DP   = true               # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg   = 1.2e-2             # regularisation "viscosity"
    # Physics
    Lx, Ly  = 1.0, 1.0           # domain size
    radi    = 0.01               # inclusion radius
    τ_y     = 1.6                # yield stress. If do_DP=true, τ_y stand for the cohesion: c*cos(ϕ)
    sinϕ    = sind(30)*do_DP     # sinus of the friction angle
    μ0      = 1.0                # viscous viscosity
    G0      = 1.0                # elastic shear modulus
    Gi      = G0/(6.0-4.0*do_DP) # elastic shear modulus perturbation
    εbg     = 1.0                # background strain-rate
    # Numerics
    nt      = 20                 # number of time steps
    nx, ny  = 63, 63             # numerical grid resolution
    Vdmp    = 4.0                # convergence acceleration (damping)
    Vsc     = 2.0                # iterative time step limiter
    Ptsc    = 6.0                # iterative time step limiter
    ε       = 1e-6               # nonlinear tolerence
    iterMax = 3e4                # max number of iters
    nout    = 200                # check frequency
    # Preprocessing
    dx, dy  = Lx/nx, Ly/ny
    dt      = μ0/G0/4.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = zeros(Dat, nx  ,ny  )
    ∇V      = zeros(Dat, nx  ,ny  )
    Vx      = zeros(Dat, nx+1,ny  )
    Vy      = zeros(Dat, nx  ,ny+1)
    Vxe     = zeros(Dat, nx+1,ny+2)
    Vye     = zeros(Dat, nx+2,ny+1)
    Exx     = zeros(Dat, nx  ,ny  )
    Eyy     = zeros(Dat, nx  ,ny  )
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
    Tyy_o   = zeros(Dat, nx  ,ny  )
    Txy_o   = zeros(Dat, nx  ,ny  )
    Txyv_o  = zeros(Dat, nx+1,ny+1)
    # for vertices implementation
    Ptv     = zeros(Dat, nx+1,ny+1)
    Exxv    = zeros(Dat, nx+1,ny+1)
    Eyyv    = zeros(Dat, nx+1,ny+1)
    Exxv1   = zeros(Dat, nx+1,ny+1)
    Eyyv1   = zeros(Dat, nx+1,ny+1)
    Txxv    = zeros(Dat, nx+1,ny+1)
    Tyyv    = zeros(Dat, nx+1,ny+1)
    Txxv_o  = zeros(Dat, nx+1,ny+1)
    Tyyv_o  = zeros(Dat, nx+1,ny+1)
    Fchkv   = zeros(Dat, nx+1,ny+1)
    Fv      = zeros(Dat, nx+1,ny+1)
    Plav    = zeros(Dat, nx+1,ny+1)
    λv      = zeros(Dat, nx+1,ny+1)
    dQdTxxv = zeros(Dat, nx+1,ny+1)
    dQdTyyv = zeros(Dat, nx+1,ny+1)
    dQdTxyv = zeros(Dat, nx+1,ny+1)
    Tiiv    = zeros(Dat, nx+1,ny+1)
    Eiiv    = zeros(Dat, nx+1,ny+1)
    # for vertices implementation
    Tii     = zeros(Dat, nx  ,ny  )
    Eii     = zeros(Dat, nx  ,ny  )
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
    η_v     =    μ0*ones(Dat, nx  ,ny  )
    η_e     = dt*G0*ones(Dat, nx  ,ny  )
    η_ev    = dt*G0*ones(Dat, nx+1,ny+1)
    η_ve    =       ones(Dat, nx  ,ny  )
    η_vep   =       ones(Dat, nx  ,ny  )
    η_vepv  =       ones(Dat, nx+1,ny+1)
    η_vev   =       ones(Dat, nx+1,ny+1) # for vertices implementation
    η_vv    =    μ0*ones(Dat, nx+1,ny+1) # for vertices implementation
    # Initial condition
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xv, yv  = LinRange(0.0, Lx, nx+1), LinRange(0.0, Ly, ny+1)
    (Xvx,Yvx) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (Xvy,Yvy) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    radc      = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv      = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    η_e[radc.<radi] .= dt*Gi
    η_ev[radv.<radi].= dt*Gi
    η_ve   .= (1.0./η_e  + 1.0./η_v).^-1
    η_vev  .= (1.0./η_ev + 1.0./η_vv).^-1
    Vx     .=  2.0.*εbg.*Yvx  # factor 2 such that Exy = Ebg
    Vy     .=  0.0.*Yvy
    Vx_bcN  =  2.0.*εbg.*Ly   # Vx at the top of the box
    Vx_bcS  =  2.0.*εbg.*0.0  # Vx at the bottom
    # Time loop
    t=0.0; evo_t=[]; evo_Txx=[]
    for it = 1:nt
        iter=1; err=2*ε; err_evo1=[]; err_evo2=[]
        Txx_o.=Txx;   Tyy_o.=Tyy;   Txy_o.=Txy;   λ.=0.0
        Txxv_o.=Txxv; Tyyv_o.=Tyyv; Txyv_o.=Txyv; λv.=0.0
        local itg
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V     .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Pt     .= Pt .- dtPt.*∇V
            # strain rates
            Exx    .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            Eyy    .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            Vxe[:,2:end-1] .= Vx; Vxe[:,1] .= 2.0*Vx_bcS .- Vxe[:,2]; Vxe[:,end] .= 2.0*Vx_bcN .- Vxe[:,end-1]
            Vye[2:end-1,:] .= Vy; Vye[1,:] .=               Vye[2,:]; Vye[end,:] .=               Vye[end-1,:]
            Exyv   .= 0.5.*(diff(Vxe, dims=2)./dy .+ diff(Vye, dims=1)./dx)
            Exxv[2:end-1,2:end-1] .= av(Exx); Exxv[1,:].=Exxv[2,:]; Exxv[end,:].=Exxv[end-1,:]; Exxv[:,1].=Exxv[:,2]; Exxv[:,end].=Exxv[:,end-1]
            Eyyv[2:end-1,2:end-1] .= av(Eyy); Eyyv[1,:].=Eyyv[2,:]; Eyyv[end,:].=Eyyv[end-1,:]; Eyyv[:,1].=Eyyv[:,2]; Eyyv[:,end].=Eyyv[:,end-1]
            Ptv[2:end-1,2:end-1]  .= av(Pt);   Ptv[1,:].= Ptv[2,:];  Ptv[end,:].= Ptv[end-1,:];  Ptv[:,1].= Ptv[:,2];  Ptv[:,end].= Ptv[:,end-1]
            # visco-elastic strain rates
            Exx1   .=    Exx   .+ Txx_o ./2.0./η_e
            Eyy1   .=    Eyy   .+ Tyy_o ./2.0./η_e
            Exy1   .= av(Exyv) .+ Txy_o ./2.0./η_e
            Eii    .= sqrt.(0.5*(Exx1.^2 .+ Eyy1.^2) .+ Exy1.^2)
            # visco-elastic strain rates vertices
            Exxv1  .=    Exxv  .+ Txxv_o./2.0./η_ev
            Eyyv1  .=    Eyyv  .+ Tyyv_o./2.0./η_ev
            Exyv1  .=    Exyv  .+ Txyv_o./2.0./η_ev
            Eiiv   .= sqrt.(0.5*(Exxv1.^2 .+ Eyyv1.^2) .+ Exyv1.^2)
            # trial stress vertices
            Txxv   .= 2.0.*η_vev.*Exxv1
            Tyyv   .= 2.0.*η_vev.*Eyyv1
            Txyv   .= 2.0.*η_vev.*Exyv1
            Tiiv   .= sqrt.(0.5*(Txxv.^2 .+ Tyyv.^2) .+ Txyv.^2)
            # yield function vertices
            Fv     .= Tiiv .- τ_y .- Ptv.*sinϕ
            Plav   .= 0.0
            Plav   .= Fv .> 0.0
            λv     .= Plav.*Fv./(η_vev .+ η_reg)
            dQdTxxv.= 0.5.*Txxv./Tiiv
            dQdTyyv.= 0.5.*Tyyv./Tiiv
            dQdTxyv.=      Txyv./Tiiv
            # plastic corrections vertices
            Txxv   .= 2.0.*η_vev.*(Exxv1 .-      λv.*dQdTxxv)
            Tyyv   .= 2.0.*η_vev.*(Eyyv1 .-      λv.*dQdTyyv)
            Txyv   .= 2.0.*η_vev.*(Exyv1 .- 0.5.*λv.*dQdTxyv)
            Tiiv   .= sqrt.(0.5*(Txxv.^2 .+ Tyyv.^2) .+ Txyv.^2)
            Fchkv  .= Tiiv .- τ_y .- Ptv.*sinϕ .- λv.*η_reg
            η_vepv .= Tiiv./2.0./Eiiv
            # back to centroids
            η_vep .= av(η_vepv)
            Txx   .= 2.0.*η_vep.*Exx1
            Tyy   .= 2.0.*η_vep.*Eyy1
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
                @printf("it = %d, iter = %d, err = %1.2e norm[Rx=%1.2e, Ry=%1.2e, ∇V=%1.2e] (Fchkv=%1.2e) \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V, maximum(Fchkv))
            end
            iter+=1; itg=iter
        end
        t = t + dt
        push!(evo_t, t); push!(evo_Txx, maximum(Txyv))
        # Plotting
        p1 = heatmap(xv, yc, Vx' , aspect_ratio=1, xlims=(0, Lx), ylims=(dy/2, Ly-dy/2), c=:inferno, title="Vx")
        # p2 = heatmap(xc, yv, Vy' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        p2 = heatmap(xc, yc, η_vep' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="η_vep")
        p3 = heatmap(xc, yc, av(Tiiv)' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="τii")
        p4 = plot(evo_t, evo_Txx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
            plot!(evo_t, 2.0.*εbg.*μ0.*(1.0.-exp.(.-evo_t.*G0./μ0)), linewidth=2.0) # analytical solution for VE loading
            plot!(evo_t, 2.0.*εbg.*μ0.*ones(size(evo_t)), linewidth=2.0)            # viscous flow stress
            if !do_DP plot!(evo_t, τ_y*ones(size(evo_t)), linewidth=2.0) end        # von Mises yield stress
        display(plot(p1, p2, p3, p4))
    end
    return
end

Stokes2D_vep()
