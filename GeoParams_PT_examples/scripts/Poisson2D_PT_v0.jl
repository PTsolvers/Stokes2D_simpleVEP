# Initialisation
using Plots, Printf, Statistics, LinearAlgebra
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end]) 
@views   avh(A) = ( 0.25./A[1:end-1,1:end-1] .+ 0.25./A[1:end-1,2:end-0] .+ 0.25./A[2:end-0,1:end-1] .+ 0.25./A[2:end-0,2:end-0]).^(-1)
# 2D Poisson routine
@views function Poisson2D()
    # Physics
    Lx, Ly   = 6.0, 6.0  # domain size
    T_West   = 0.0
    T_East   = 0.0
    T_South  = 0.0
    T_North  = 0.0
    slope    = -15.
    k1       = 1
    k2       = 1e4
    # Numerics
    ncx, ncy = 51, 51    # numerical grid resolution
    ε        = 1e-6      # nonlinear tolerence
    iterMax  = 1e4       # max number of iters
    nout     = 500       # check frequency
    # Iterative parameters -------------------------------------------
    Reopt    = 0.625*π
    cfl      = 0.32
    ρ        = cfl*Reopt/ncx
    nsm      = 5
    # Reopt    = 0.625*π /3
    # cfl      = 0.3*4
    # ρ        = cfl*Reopt/ncx
    # nsm      = 10
    # Preprocessing
    Δx, Δy  = Lx/ncx, Ly/ncy
    # Array initialisation
    T         =   zeros(Dat, ncx+2, ncy+2)
    ∂T∂x      =   zeros(Dat, ncx+1, ncy  )
    ∂T∂y      =   zeros(Dat, ncx  , ncy+1)
    qx        =   zeros(Dat, ncx+1, ncy  )
    qy        =   zeros(Dat, ncx  , ncy+1)
    RT        =   zeros(Dat, ncx  , ncy  )
    dTdτ      =   zeros(Dat, ncx  , ncy  ) 
    Δτc       =   zeros(Dat, ncx  , ncy  )       
    kv        = k1*ones(Dat, ncx+1, ncy+1)
    kx        =   zeros(Dat, ncx+1, ncy+0)
    ky        =   zeros(Dat, ncx+0, ncy+1)
    H         =    ones(Dat, ncx  , ncy  )
    # Initialisation
    xc, yc    = LinRange(-Lx/2+Δx/2, Lx/2-Δx/2, ncx+0), LinRange(-Ly/2+Δy/2, Ly/2-Δy/2, ncy+0)
    xce, yce  = LinRange(-Lx/2-Δx/2, Lx/2+Δx/2, ncx+2), LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, ncy+2)
    xv, yv    = LinRange(-Lx/2, Lx/2, ncx+1), LinRange(-Ly/2, Ly/2, ncy+1)
    (xv2,yv2) = ([x for x=xv,y=yv], [y for x=xv,y=yv])
    below = yv2 .< xv2.*tand.(slope)
    kv[below] .= k2
    kx .= ( 0.5./kv[:,1:end-1] .+ 0.5./kv[:,2:end-0]).^(-1)
    ky .= ( 0.5./kv[1:end-1,:] .+ 0.5./kv[2:end-0,:]).^(-1)
    kc  = ( 0.25./kv[1:end-1,1:end-1] .+ 0.25./kv[1:end-1,2:end-0] .+ 0.25./kv[2:end-0,1:end-1] .+ 0.25./kv[2:end-0,2:end-0]).^(-1)
    for it = 1:nsm
        kv[2:end-1,2:end-1] .= av(kc)
        kc = av(kv)
    end
    # Time loop
    it=1; iter=1; err=2*ε; err_evo1=[]; err_evo2=[];
    while (err>ε && iter<=iterMax)
        # BCs
        T[:,1]   .= 2*T_South .- T[:,2]     # S
        T[:,end] .= 2*T_North .- T[:,end-1] # N
        T[1,:]   .= 2*T_West  .- T[2,:]     # W
        T[end,:] .= 2*T_East  .- T[end-1,:] # E
        # Kinematics
        ∂T∂x .= diff(T[:,2:end-1], dims=1)./Δx 
        ∂T∂y .= diff(T[2:end-1,:], dims=2)./Δy 
        # Stresses
        qx .= -kx.*∂T∂x
        qy .= -ky.*∂T∂y
        # Residuals
        RT .= .-diff(qx, dims=1)./Δx .- diff(qy, dims=2)./Δy + H
        # PT time step -----------------------------------------------
        Δτc  .= ρ*min(Δx,Δy)^2 ./ kc ./ 4.1 * cfl 
        # Calculate rate update --------------------------------------
        dTdτ          .= (1-ρ) .* dTdτ .+ RT
        # Update velocity and pressure -------------------------------
        T[2:end-1,2:end-1]  .+= Δτc ./ ρ .* dTdτ
        # convergence check
        if mod(iter, nout)==0
            norm_RT = norm(RT)/sqrt(length(RT))
            err = maximum([norm_RT])
            push!(err_evo1, err); push!(err_evo2, itg)
            @printf("it = %03d, iter = %04d, err = %1.3e norm[RT=%1.3e] \n", it, itg, err, norm_RT)
        end
        iter+=1; global itg=iter
    end
    # Plotting
    p1 = heatmap(xce, yce,         T', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="T")
    p2 = heatmap( xc, yv, log10.(ky)', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="ky")
    p3 = heatmap( xc, yc, log10.(kc)', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="kc")
    y_int = xv.*tand.(slope)
    p4 = plot(xv, y_int, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:black, label=:none)
    # lines
    for j=1:ncy+1
        y_line = yv[j]*ones(size(xv))
        p4 = plot!(xv, y_line, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:gray, label=:none)
    end
    for i=1:ncx+1
        x_line = xv[i]*ones(size(yv))
        p4 = plot!(x_line, yv, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:gray, label=:none)
    end
    display(plot(p1, p2, p3, p4))
    return
end

Poisson2D()
