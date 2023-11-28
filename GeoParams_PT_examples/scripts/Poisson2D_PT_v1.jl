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
    ncx, ncy = 11, 11    # numerical grid resolution
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
    (xc2,yc2) = ([x for x=xc,y=yc], [y for x=xc,y=yc])
    (xqx2,yqx2) = ([x for x=xv,y=yc], [y for x=xv,y=yc])
    (xqy2,yqy2) = ([x for x=xc,y=yv], [y for x=xc,y=yv])
    below = yv2 .< xv2.*tand.(slope)
    kv[below] .= k2
    kx .= ( 0.5./kv[:,1:end-1] .+ 0.5./kv[:,2:end-0]).^(-1)
    ky .= ( 0.5./kv[1:end-1,:] .+ 0.5./kv[2:end-0,:]).^(-1)
    kc  = ( 0.25./kv[1:end-1,1:end-1] .+ 0.25./kv[1:end-1,2:end-0] .+ 0.25./kv[2:end-0,1:end-1] .+ 0.25./kv[2:end-0,2:end-0]).^(-1)
    for it = 1:nsm
        kv[2:end-1,2:end-1] .= av(kc)
        kc = av(kv)
    end
    # Adapt mesh
    nc    = ncx*ncy
    XC    = xc2[:]
    YC    = yc2[:]
    num   = reshape(1:ncx*ncy, ncx, ncy)
    iS    = zeros(size(xc2)); iS[:,2:end-0] = num[:,1:end-1]
    iN    = zeros(size(xc2)); iN[:,1:end-1] = num[:,2:end-0]
    IS    = iS[:]
    IN    = iN[:]
    vertx = [ xv2[1:end-1,1:end-1][:] xv2[2:end-0,1:end-1][:] xv2[2:end-0,2:end-0][:] xv2[1:end-1,2:end-0][:] ] 
    verty = [ yv2[1:end-1,1:end-1][:] yv2[2:end-0,1:end-1][:] yv2[2:end-0,2:end-0][:] yv2[1:end-1,2:end-0][:] ] 
    for ic=1:nc
        y_int = XC[ic].*tand.(slope)
        @show y_int-YC[ic]
        @show Δy/2
        if abs(y_int-YC[ic]) < Δy/2
            println("found interface - adjust top")
            if YC[ic] < y_int
                println("adjust top")
                verty[ic,3] = vertx[ic,3][ic].*tand.(slope)
                verty[ic,4] = vertx[ic,4][ic].*tand.(slope)
            end
            if YC[ic] > y_int
                println("adjust bottom")
                verty[ic,1] = vertx[ic,1][ic].*tand.(slope)
                verty[ic,2] = vertx[ic,2][ic].*tand.(slope)
            end
        end
    end
    xqx   = zeros(nc,2)
    yqx   = zeros(nc,2)
    xqy   = zeros(nc,2)
    yqy   = zeros(nc,2)
    for ic=1:nc
        xqx[ic,1] = 0.5*( vertx[ic,1] + vertx[ic,4] )
        yqx[ic,1] = 0.5*( verty[ic,1] + verty[ic,4] )
        xqx[ic,2] = 0.5*( vertx[ic,2] + vertx[ic,3] )
        yqx[ic,2] = 0.5*( verty[ic,2] + verty[ic,3] )
        xqy[ic,1] = 0.5*( vertx[ic,1] + vertx[ic,2] )
        yqy[ic,1] = 0.5*( verty[ic,1] + verty[ic,2] )
        xqy[ic,2] = 0.5*( vertx[ic,3] + vertx[ic,4] )
        yqy[ic,2] = 0.5*( verty[ic,3] + verty[ic,4] )
    end
    # Back to 2D table
    xc3 = reshape(XC, ncx, ncy)
    yc3 = reshape(YC, ncx, ncy)
    xv3 = zeros( ncx+1, ncy+1 )
    yv3 = zeros( ncx+1, ncy+1 )
    xv3[1:end-1, 1:end-1] .= reshape(vertx[:,1], ncx, ncy)
    xv3[2:end-0, 1:end-1] .= reshape(vertx[:,2], ncx, ncy)
    xv3[2:end-0, 2:end-0] .= reshape(vertx[:,3], ncx, ncy)
    xv3[1:end-1, 2:end-0] .= reshape(vertx[:,4], ncx, ncy)
    yv3[1:end-1, 1:end-1] .= reshape(verty[:,1], ncx, ncy)
    yv3[2:end-0, 1:end-1] .= reshape(verty[:,2], ncx, ncy)
    yv3[2:end-0, 2:end-0] .= reshape(verty[:,3], ncx, ncy)
    yv3[1:end-1, 2:end-0] .= reshape(verty[:,4], ncx, ncy)
    xqx3 = zeros( ncx+1, ncy )
    yqx3 = zeros( ncx+1, ncy )
    xqy3 = zeros( ncx, ncy+1 )
    yqy3 = zeros( ncx, ncy+1 )
    xqx3[1:end-1, :] .= reshape(xqx[:,1], ncx, ncy)
    xqx3[2:end-0, :] .= reshape(xqx[:,2], ncx, ncy)
    yqx3[1:end-1, :] .= reshape(yqx[:,1], ncx, ncy)
    yqx3[2:end-0, :] .= reshape(yqx[:,2], ncx, ncy)
    xqy3[:, 1:end-1] .= reshape(xqy[:,1], ncx, ncy)
    xqy3[:, 2:end-0] .= reshape(xqy[:,2], ncx, ncy)
    yqy3[:, 1:end-1] .= reshape(yqy[:,1], ncx, ncy)
    yqy3[:, 2:end-0] .= reshape(yqy[:,2], ncx, ncy)
    @show norm(xc2.-xc3)
    @show norm(yc2.-yc3)
    @show norm(xv2.-xv3)
    @show norm(yv2.-yv3)
    @show norm(xv2.-xv3)
    @show norm(xqx3.-xqx2)
    @show norm(yqx3.-yqx2)
    @show norm(xqy3.-xqy2)
    @show norm(yqy3.-yqy2)
    Δx_∂T∂x = zeros(ncx+1,ncy)
    Δy_∂T∂y = zeros(ncx,ncy+1)
    Δx_∂T∂x[2:end-1,:] = diff(xc3, dims=1); Δx_∂T∂x[1,:] = Δx_∂T∂x[2,:]; Δx_∂T∂x[end,:] = Δx_∂T∂x[end-1,:]
    Δy_∂T∂y[:,2:end-1] = diff(yc3, dims=2); Δy_∂T∂y[:,1] = Δy_∂T∂y[:,2]; Δy_∂T∂y[:,end] = Δy_∂T∂y[:,end-1]
    Δx_∂q∂x = diff(xqx3, dims=1)
    Δy_∂q∂y = diff(yqy3, dims=2)
    # Time loop
    it=1; iter=1; err=2*ε; err_evo1=[]; err_evo2=[];
    # while (err>ε && iter<=iterMax)
    #     # BCs
    #     T[:,1]   .= 2*T_South .- T[:,2]     # S
    #     T[:,end] .= 2*T_North .- T[:,end-1] # N
    #     T[1,:]   .= 2*T_West  .- T[2,:]     # W
    #     T[end,:] .= 2*T_East  .- T[end-1,:] # E
    #     # Kinematics
    #     ∂T∂x .= diff(T[:,2:end-1], dims=1)./Δx_∂T∂x 
    #     ∂T∂y .= diff(T[2:end-1,:], dims=2)./Δy_∂T∂y 
    #     # Stresses
    #     qx .= -kx.*∂T∂x
    #     qy .= -ky.*∂T∂y
    #     # Residuals
    #     RT .= .-diff(qx, dims=1)./Δx_∂q∂x .- diff(qy, dims=2)./Δy_∂q∂y + H
    #     # PT time step -----------------------------------------------
    #     Δτc  .= ρ*min(Δx,Δy)^2 ./ kc ./ 4.1 * cfl 
    #     # Calculate rate update --------------------------------------
    #     dTdτ          .= (1-ρ) .* dTdτ .+ RT
    #     # Update velocity and pressure -------------------------------
    #     T[2:end-1,2:end-1]  .+= Δτc ./ ρ .* dTdτ
    #     # convergence check
    #     if mod(iter, nout)==0
    #         norm_RT = norm(RT)/sqrt(length(RT))
    #         err = maximum([norm_RT])
    #         push!(err_evo1, err); push!(err_evo2, itg)
    #         @printf("it = %03d, iter = %04d, err = %1.3e norm[RT=%1.3e] \n", it, itg, err, norm_RT)
    #     end
    #     iter+=1; global itg=iter
    # end
    # Plotting
    p1 = heatmap(xce, yce,         T', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="T")
    p2 = heatmap( xc, yv, log10.(ky)', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="ky")
    p3 = heatmap( xc, yc, log10.(kc)', aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), c=:turbo, title="kc")
    y_int = xv.*tand.(slope)
    p4 = plot(xv, y_int, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:black, label=:none)
    ied1  = [ 1 2 3 4] 
    ied2  = [ 2 3 4 1] 
    for ic=1:length(XC)
        for ied=1:4
            edgex = [vertx[ic, ied1[ied]]; vertx[ic, ied2[ied]]]
            edgey = [verty[ic, ied1[ied]]; verty[ic, ied2[ied]]]
            p4 = plot!( edgex, edgey, color=:red, label=:none)
        end
    end
    p4 = plot!( XC, YC , marker=:circle, linewidth=0, color=:green, label=:none)
    p4 = plot!( xqx3, yqx3, marker=:cross, linewidth=0, color=:green, label=:none)
    p4 = plot!( xqy3, yqy3, marker=:cross, linewidth=0, color=:blue, label=:none)
    # lines
    # for j=1:ncy+1
    #     y_line = yv[j]*ones(size(xv))
    #     p4 = plot!(xv, y_line, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:gray, label=:none)
    # end
    # for i=1:ncx+1
    #     x_line = xv[i]*ones(size(yv))
    #     p4 = plot!(x_line, yv, aspect_ratio=1, xlims=(-Lx/2, Lx/2), ylims=(-Ly/2, Ly/2), color=:gray, label=:none)
    # end
    # p4 = plot!( xc2, yc2, marker=:circle, linewidth=0, color=:green, label=:none)
    display(plot(p1, p2, p3, p4))
    display(plot(p4))
    return
end

Poisson2D()
