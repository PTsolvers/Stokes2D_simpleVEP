# Initialisation
using Plots, Printf, Statistics, LinearAlgebra
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
# 2D Stokes routine
@views function Stokes2D_ve()
    # Physics
    Lx, Ly  = 10.0, 10.0  # domain size
    ξ       = 1.0         # Maxwell relaxation time
    μ0      = 1.0         # viscous viscosity
    μi      = 0.1         # viscous viscosity perturbation
    G       = 1.0         # elastic shear modulus
    ρg0     = 1.0         # density*gravity
    # Numerics
    nt      = 1           # number of time steps
    nx, ny  = 63, 63      # numerical grid resolution
    Vdmp    = 4.0         # convergence acceleration (damping)
    Ptsc    = 8.0         # iterative time step limiter
    ε       = 1e-6        # nonlinear tolerence
    iterMax = 1e5         # max number of iters
    nout    = 200         # check frequency
    # Preprocessing
    dx, dy  = Lx/nx, Ly/ny
    dt      = μ0/(G*ξ + 1e-15)
    # Array initialisation
    Pt      = zeros(Dat, nx  ,ny  )
    ∇V      = zeros(Dat, nx  ,ny  )
    Vx      = zeros(Dat, nx+1,ny  )
    Vy      = zeros(Dat, nx  ,ny+1)
    Exx     = zeros(Dat, nx  ,ny  )
    Eyy     = zeros(Dat, nx  ,ny  )
    Exy     = zeros(Dat, nx-1,ny-1)    
    Txx     = zeros(Dat, nx  ,ny  )
    Tyy     = zeros(Dat, nx  ,ny  )
    Txy     = zeros(Dat, nx+1,ny+1)
    Txx_o   = zeros(Dat, nx  ,ny  )
    Tyy_o   = zeros(Dat, nx  ,ny  )
    Txy_o   = zeros(Dat, nx+1,ny+1)
    Rx      = zeros(Dat, nx-1,ny  )
    Ry      = zeros(Dat, nx  ,ny-1)
    dVxdt   = zeros(Dat, nx-1,ny  )
    dVydt   = zeros(Dat, nx  ,ny-1)
    Rog     = zeros(Dat, nx  ,ny  )
    Xsi     =  ξ*ones(Dat, nx, ny)
    Mus     = μ0*ones(Dat, nx, ny)
    # Initialisation
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    rad          = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    Mus[rad.<1] .= μi
    Rog[rad.<1] .= ρg0
    dtVx    = min(dx,dy)^2.0./av_xa(Mus)./4.1
    dtVy    = min(dx,dy)^2.0./av_ya(Mus)./4.1
    dtPt    = 4.1*Mus/max(nx,ny)/Ptsc
    # Time loop
    err_evo1=[]; err_evo2=[]
    for it = 1:nt
        iter=1; err=2*ε;
        Txx_o.=Txx; Tyy_o.=Tyy; Txy_o.=Txy
        while (err>ε && iter<=iterMax)
            # divergence - pressure
            ∇V    .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            Pt    .= Pt .- dtPt.*∇V
            # strain rates
            Exx   .= diff(Vx, dims=1)./dx .- 1.0/3.0*∇V
            Eyy   .= diff(Vy, dims=2)./dy .- 1.0/3.0*∇V
            Exy   .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            # stresses
            Xsi   .= Mus./(G.*dt)
            Txx   .= Txx_o.*Xsi./(Xsi.+1.0) .+ 2.0.*Mus.*Exx./(Xsi.+1.0)
            Tyy   .= Tyy_o.*Xsi./(Xsi.+1.0) .+ 2.0.*Mus.*Eyy./(Xsi.+1.0)
            Txy[2:end-1,2:end-1] .= Txy_o[2:end-1,2:end-1].*av(Xsi)./(av(Xsi).+1.0) .+ 2.0.*av(Mus).*Exy./(av(Xsi).+1.0)
            # velocities
            Rx    .= .-diff(Pt, dims=1)./dx .+ diff(Txx, dims=1)./dx .+ diff(Txy[2:end-1,:], dims=2)./dy
            Ry    .= .-diff(Pt, dims=2)./dy .+ diff(Tyy, dims=2)./dy .+ diff(Txy[:,2:end-1], dims=1)./dx .+ av_ya(Rog)
            dVxdt .= dVxdt.*(1-Vdmp/nx) .+ Rx
            dVydt .= dVydt.*(1-Vdmp/ny) .+ Ry
            Vx[2:end-1,:] .= Vx[2:end-1,:] .+ dVxdt.*dtVx
            Vy[:,2:end-1] .= Vy[:,2:end-1] .+ dVydt.*dtVy
            # convergence check
            if mod(iter, nout)==0
                global max_Rx, max_Ry, max_divV
                norm_Rx = norm(Rx)/length(Rx); norm_Ry = norm(Ry)/length(Ry); norm_∇V = norm(∇V)/length(∇V)
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %d, iter = %d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, ∇V=%1.3e] \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V)
            end
            iter+=1; global itg=iter
        end
        # Plotting
        yv = LinRange(0, Ly, ny+1)
        p1 = heatmap(xc, yc, Pt' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(dy/2, Ly-dy/2), c=:inferno, title="Pressure")
        p2 = heatmap(xc, yv, Vy' , aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        p3 = heatmap(xc, yc, Tyy', aspect_ratio=1, xlims=(dx/2, Lx-dx/2), ylims=(dy/2, Ly-dy/2), c=:inferno, title="τyy")
        p4 = plot(err_evo2, log10.(err_evo1), legend=false, xlabel="# iterations", ylabel="log10(error)", linewidth=2, markershape=:circle, markersize=3, framestyle=:box, labels="max(error)")
        display(plot(p1, p2, p3, p4))
    end
end

Stokes2D_ve()
