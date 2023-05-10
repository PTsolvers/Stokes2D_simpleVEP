# with hyperbolic yield
# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, SpecialFunctions
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
# 2D Stokes routine
@views function Stokes2D_vep()
    do_DP   = true               # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg   = 1.4*8.0e-3             # regularisation "viscosity"
    # Physics
    Lx, Ly  = 1.0, 1.0           # domain size
    radi    = 0.01               # inclusion radius
    yield   = :Hyperbolic        # yield function
    # yield   = :Linear        # yield function
    C       = 1.8475208614068026 # yield stress. If do_DP=true, C.*cosϕ stand for the cohesion: c*cos(ϕ)
    σt      = C/2                # tension stress
    sinϕ    = sind(30)*do_DP     # sinus of the friction angle
    cosϕ    = cosd(30)*do_DP     # cosinus of the friction angle
    sinψ    = sind(5)*do_DP      # sinus of the dilation angle
    σd      = C/2                # dilation transition function
    Pd      = C                  # dilation transition function
    μ0      = 1.0                # viscous viscosity
    G0      = 1.0                # elastic shear modulus
    Gi      = G0/(6.0-4.0*do_DP) # elastic shear modulus perturbation
    K       = 5
    εbg     = 1.0                # background strain-rate
    # Numerics
    nt      = 10                 # number of time steps
    nx, ny  = 63, 63             # numerical grid resolution
    Vdmp    = 4.0                # convergence acceleration (damping)
    Vsc     = 2.0                # iterative time step limiter
    Ptsc    = 6.0               # iterative time step limiter
    ε       = 1e-6               # nonlinear tolerence
    iterMax = 3e4                # max number of iters
    nout    = 200                # check frequency
    # Preprocessing
    dx, dy  = Lx/nx, Ly/ny
    dt      = μ0/G0/4.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = zeros(Dat, nx  ,ny  )
    Pt_o    = zeros(Dat, nx  ,ny  )
    Pt_c    = zeros(Dat, nx  ,ny  )
    Rp      = zeros(Dat, nx  ,ny  )
    ∇V      = zeros(Dat, nx  ,ny  )
    Vx      = zeros(Dat, nx+1,ny  )
    Vy      = zeros(Dat, nx  ,ny+1)
    ε̇II     = zeros(Dat, nx  ,ny  )
    Exx     = zeros(Dat, nx  ,ny  )
    Eyy     = zeros(Dat, nx  ,ny  )
    Exyv    = zeros(Dat, nx+1,ny+1)
    Exx1    = zeros(Dat, nx  ,ny  )
    Eyy1    = zeros(Dat, nx  ,ny  )
    Exy1    = zeros(Dat, nx  ,ny  )
    Txx     = zeros(Dat, nx  ,ny  )
    Tyy     = zeros(Dat, nx  ,ny  )
    Txy     = zeros(Dat, nx  ,ny  )
    Txyv    = zeros(Dat, nx+1,ny+1)
    Txx_o   = zeros(Dat, nx  ,ny  )
    Tyy_o   = zeros(Dat, nx  ,ny  )
    Txy_o   = zeros(Dat, nx  ,ny  )
    Txyv_o  = zeros(Dat, nx+1,ny+1)
    τII     = zeros(Dat, nx  ,ny  )
    τII_c   = zeros(Dat, nx  ,ny  )
    Eii     = zeros(Dat, nx  ,ny  )
    F       = zeros(Dat, nx  ,ny  )
    Fchk    = zeros(Dat, nx  ,ny  )
    Pla     = zeros(Dat, nx  ,ny  )
    λ       = zeros(Dat, nx  ,ny  )
    dλdτ    = zeros(Dat, nx  ,ny  )
    Rx      = zeros(Dat, nx-1,ny  )
    Ry      = zeros(Dat, nx  ,ny-1)
    dVxdt   = zeros(Dat, nx-1,ny  )
    dVydt   = zeros(Dat, nx  ,ny-1)
    dtPt    = zeros(Dat, nx  ,ny  )
    dtVx    = zeros(Dat, nx-1,ny  )
    dtVy    = zeros(Dat, nx  ,ny-1)
    Rog     = zeros(Dat, nx  ,ny  )
    η_v     =    μ0*ones(Dat, nx, ny)
    η_e     = dt*G0*ones(Dat, nx, ny)
    η_ev    = dt*G0*ones(Dat, nx+1, ny+1)
    η_ve    =       ones(Dat, nx, ny)
    η_vep   =       ones(Dat, nx, ny)
    c       = zeros(Dat, nx  ,ny  )
    dQdP    = zeros(Dat, nx  ,ny  )
    dQdτ    = zeros(Dat, nx  ,ny  )
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
    η_ve   .= (1.0./η_e + 1.0./η_v).^-1
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=[]; evo_Txx=[]; 

    for it = 1:nt
        iter=1; err=2*ε; err_evo1=[]; err_evo2=[]
        Txx_o.=Txx; Tyy_o.=Tyy; Txy_o.=av(Txyv); Txyv_o.=Txyv; λ.=0.0
        local itg = 0
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V     .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Rp     .= .-∇V .- (Pt.-Pt_o)./(K.*dt)
            Pt     .= Pt .+ dtPt.*Rp
            # strain rates
            Exx    .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            Eyy    .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            Exyv[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            # visco-elastic strain rates
            Exx1   .=    Exx   .+ Txx_o ./2.0./η_e
            Eyy1   .=    Eyy   .+ Tyy_o ./2.0./η_e
            Exy1   .= av(Exyv) .+ Txy_o ./2.0./η_e
            Eii    .= sqrt.(0.5*(Exx1.^2 .+ Eyy1.^2) .+ Exy1.^2)
            # trial stress
            Txx    .= 2.0.*η_ve.*Exx1
            Tyy    .= 2.0.*η_ve.*Eyy1
            Txy    .= 2.0.*η_ve.*Exy1
            τII    .= sqrt.(0.5*(Txx.^2 .+ Tyy.^2) .+ Txy.^2)
            # Plasticity
            if yield==:Hyperbolic
            # Hyperbolic yield function
                F      .= sqrt.(τII.^2 .+ (C.*cosϕ .- σt.*sinϕ).^2) .- (C.*cosϕ .+ Pt.*sinϕ)
                Pla    .= 0.0
                Pla    .= F .> 0.0
                dτ_λ    = 0.08
                c      .= 0.5 * erfc.(((Pt .- Pd) ./ σd))
                dQdP   .= -c .* sinϕ
                dQdτ   .= τII ./ sqrt.(τII^2 .+ (C.*cosϕ .- σt.*sinϕ).^2)
                τII_c  .= τII .- λ.*η_ve.*dQdτ
                Pt_c   .= Pt  .- λ*K.*dt.*dQdP
                dλdτ   .= sqrt.(τII_c.^2 .+ (C.*cosϕ .- σt.*sinϕ).^2) .- (C.*cosϕ .+ Pt_c.*sinϕ .+ λ.*η_reg)
                λ     .+= Pla.*dλdτ * dτ_λ .+ (1.0.-Pla).*0.
                η_vep  .= τII_c./2.0./Eii
                # plastic corrections
                Txx    .= 2.0.*η_vep.*Exx1
                Tyy    .= 2.0.*η_vep.*Eyy1
                Txy    .= 2.0.*η_vep.*Exy1
                Fchk   .= sqrt.(τII_c.^2 .+ (C.*cosϕ .- σt.*sinϕ).^2) .- (C.*cosϕ .+ Pt_c.*sinϕ .+ λ.*η_reg)
            elseif yield==:Linear
                # Linear yield function - PT update of λ
                F      .= τII .- C.*cosϕ .- Pt.*sinϕ# .- λ.*η_reg
                Pla    .= 0.0
                Pla    .= F .> 0.0
                dτ_λ    = 0.5
                τII_c  .= τII .- η_ve.*λ
                Pt_c   .= Pt  .+ λ.*K.*dt.*sinψ
                dλdτ   .= τII_c .- C.*cosϕ .- Pt_c.*sinϕ .- λ.*η_reg
                λ     .+= Pla.*dλdτ * dτ_λ .+ (1.0.-Pla).*0.
                η_vep  .= τII_c./2.0./Eii
                # plastic corrections
                Txx    .= 2.0.*η_vep.*Exx1
                Tyy    .= 2.0.*η_vep.*Eyy1
                Txy    .= 2.0.*η_vep.*Exy1
                τII_c  .= sqrt.(0.5*(Txx.^2 .+ Tyy.^2) .+ Txy.^2)
                Fchk   .= τII_c .- C.*cosϕ .- Pt_c.*sinϕ .- λ.*η_reg
            end
            Txyv[2:end-1,2:end-1].=av(Txy) # Txyv=0 on boundaries !
            # PT timestep
            dtVx   .= min(dx,dy)^2.0./av_xa(η_vep)./4.1./Vsc
            dtVy   .= min(dx,dy)^2.0./av_ya(η_vep)./4.1./Vsc
            dtPt   .= 4.1.*η_vep./max(nx,ny)./Ptsc
            # velocities
            Rx     .= .-diff(Pt_c, dims=1)./dx .+ diff(Txx, dims=1)./dx .+ diff(Txyv[2:end-1,:], dims=2)./dy
            Ry     .= .-diff(Pt_c, dims=2)./dy .+ diff(Tyy, dims=2)./dy .+ diff(Txyv[:,2:end-1], dims=1)./dx .+ av_ya(Rog)
            dVxdt  .= dVxdt.*(1-Vdmp/nx) .+ Rx
            dVydt  .= dVydt.*(1-Vdmp/ny) .+ Ry
            Vx[2:end-1,:] .= Vx[2:end-1,:] .+ dVxdt.*dtVx
            Vy[:,2:end-1] .= Vy[:,2:end-1] .+ dVydt.*dtVy
            # convergence check
            if mod(iter, nout)==0 || iter==1
                norm_Rx = norm(Rx)/length(Rx); norm_Ry = norm(Ry)/length(Ry); norm_Rp = norm(Rp)/length(Rp)
                err = maximum([norm_Rx, norm_Ry, norm_Rp])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %d, iter = %d, err = %1.2e norm[Rx=%1.2e, Ry=%1.2e, ∇V=%1.2e] (F=%1.2e --- Fchk=%1.2e) \n", it, itg, err, norm_Rx, norm_Ry, norm_Rp, maximum(F), maximum(Fchk))
            end
            iter+=1; itg=iter
        end
        ε̇II .= sqrt.( Exx.^2 .+ Eyy.^2 .+ av(Exyv).^2 )
        Pt  .= Pt_c # make sure corrected pressure is used for time evolution
        t    = t + dt
        push!(evo_t, t); push!(evo_Txx, maximum(Txx))
        # Plotting
        p1 = heatmap(xc, yc, τII' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(dy/2, Ly-dy/2), c=:inferno, title="τii")
        p2 = heatmap(xc, yc, Pt_c', aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(dy/2, Ly-dy/2), c=:inferno, title="η_vep")
        p3 = heatmap(xc, yc, ε̇II' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(dy/2, Ly-dy/2), c=:inferno, title="ε̇II")
        p4 = plot(evo_t, evo_Txx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
            plot!(evo_t, 2.0.*εbg.*μ0.*(1.0.-exp.(.-evo_t.*G0./μ0)), linewidth=2.0) # analytical solution for VE loading
            plot!(evo_t, 2.0.*εbg.*μ0.*ones(size(evo_t)), linewidth=2.0)            # viscous flow stress
            if !do_DP plot!(evo_t, C.*cosϕ*ones(size(evo_t)), linewidth=2.0) end        # von Mises yield stress
        display(plot(p1, p2, p3, p4))
    end
    return
end

Stokes2D_vep()
