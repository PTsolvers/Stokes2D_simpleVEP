# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, GeoParams
Dat = Float64  # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end]) 
# Rheology
@views function UpdateStressGeoParams!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt )
    ε̇iic  = ones(size(ε̇xx))
    ε̇iiv  = ones(size(ε̇xy))
    τii0c = sqrt.(1//2 .*(τxx0.^2 .+ τyy0.^2) .+ av(τxy0).^2)
    τii0v = zeros(size(ε̇xy))
    τii0v[2:end-1,2:end-1] = sqrt.(1//2 .*( av(τxx0).^2 .+ av(τyy0).^2) .+ τxy0[2:end-1,2:end-1].^2)
    Phase = 1
    # Centroids
    for i in eachindex(τxx)
        v      = MatParam[Phase].CompositeRheology[1]
        args   = (; τII_old = τii0c[i], dt=Δt)
        η_eff  = computeViscosity_εII(v, ε̇iic[i], args)
        τxx[i] = 2*η_eff*ε̇xx[i]
        τyy[i] = 2*η_eff*ε̇yy[i]
    end
    # Vertices
    for i in eachindex(τxy)
        v      = MatParam[Phase].CompositeRheology[1]
        args   = (; τII_old = τii0v[i], dt=Δt)
        η_eff  = computeViscosity_εII(v, ε̇iiv[i], args)
        τxy[i] = 2*η_eff*ε̇xy[i] 
    end
end
# 2D Stokes routine
@views function Stokes2D_VE_benchmark()
    # Physics
    Lx, Ly  = 1.0, 1.0  # domain size
    ξ       = 10.0      # Maxwell relaxation time
    η0      = 1.0       # viscous viscosity
    G       = 1.0       # elastic shear modulus
    εbg     = 1.0       # background strain-rate

    MatParam = (SetMaterialParams(Name="Matrix", Phase=1,
                CompositeRheology = CompositeRheology(ConstantElasticity(G=G),LinearViscous(η=η0))),)

    # Numerics
    nt       = 1        # number of time steps
    ncx, ncy = 31, 31    # numerical grid resolution
    ε        = 1e-6      # nonlinear tolerence
    iterMax  = 1e5       # max number of iters
    nout     = 200       # check frequency
    # Iterative parameters -------------------------------------------
    Reopt    = 5π
    cfl      = 0.5
    ρ        = 1#cfl*Reopt/ncx
    # Preprocessing
    Δx, Δy  = Lx/ncx, Ly/ncy
    Δt      = η0/(G*ξ + 1e-15)
    # Array initialisation
    Pt       = zeros(Dat, ncx  ,ncy  )
    ∇V       = zeros(Dat, ncx  ,ncy  )
    Vx       = zeros(Dat, ncx+1,ncy+2)
    Vy       = zeros(Dat, ncx+2,ncy+1)
    ε̇xx      = zeros(Dat, ncx  ,ncy  )
    ε̇yy      = zeros(Dat, ncx  ,ncy  )
    ε̇xy      = zeros(Dat, ncx+1,ncy+1)    
    τxx      = zeros(Dat, ncx  ,ncy  )
    τyy      = zeros(Dat, ncx  ,ncy  )
    τxy      = zeros(Dat, ncx+1,ncy+1)
    τxx0     = zeros(Dat, ncx  ,ncy  )
    τyy0     = zeros(Dat, ncx  ,ncy  )
    τxy0     = zeros(Dat, ncx+1,ncy+1)
    Rx       = zeros(Dat, ncx+1,ncy  )
    Ry       = zeros(Dat, ncx  ,ncy+1)
    Rp       = zeros(Dat, ncx  ,ncy  )
    dVxdτ    = zeros(Dat, ncx+1,ncy  ) 
    dVydτ    = zeros(Dat, ncx  ,ncy+1)
    dPdτ     = zeros(Dat, ncx  ,ncy  )
    Δτv      = zeros(Dat, ncx+1,ncy+1)       
    Δτvx     = zeros(Dat, ncx+1,ncy  )         
    Δτvy     = zeros(Dat, ncx  ,ncy+1)         
    κΔτp     = zeros(Dat, ncx  ,ncy  )           
    Rog      = zeros(Dat, ncx  ,ncy  )
    ηc       = η0*ones(Dat, ncx, ncy)
    ηv       = η0*ones(Dat, ncx+1, ncy+1)
    # Initialisation
    xc, yc    = LinRange( Δx/2, Lx-Δx/2, ncx+0), LinRange( Δy/2, Ly-Δy/2, ncy+0)
    xce, yce  = LinRange(-Δx/2, Lx+Δx/2, ncx+2), LinRange(-Δy/2, Ly+Δy/2, ncy+2)
    xv, yv    = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    (Xvx,Yvx) = ([x for x=xv,y=yce], [y for x=xv,y=yce])
    (Xvy,Yvy) = ([x for x=xce,y=yv], [y for x=xce,y=yv])
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=[]; evo_τxx=[]
    for it = 1:nt
        iter=1; err=2*ε; err_evo1=[]; err_evo2=[]; 
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
            UpdateStressGeoParams!( τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0, MatParam, Δt )
            if iter==1
                display(τxx')
            end
            # Residuals
            Rx[2:end-1,:] .= .-diff(Pt, dims=1)./Δx .+ diff(τxx, dims=1)./Δx .+ diff(τxy[2:end-1,:], dims=2)./Δy
            Ry[:,2:end-1] .= .-diff(Pt, dims=2)./Δy .+ diff(τyy, dims=2)./Δy .+ diff(τxy[:,2:end-1], dims=1)./Δx #.+ av_ya(Rog)
            Rp            .= .-∇V
            # PT time step -----------------------------------------------
            Δτv  .= ρ*min(Δx,Δy)^2 ./ ηv ./ 4.1 * cfl 
            Δτvx .= (Δτv[:,1:end-1] .+ Δτv[:,2:end]) / 2.
            Δτvy .= (Δτv[1:end-1,:] .+ Δτv[2:end,:]) / 2.
            κΔτp .= cfl .* ηc .* Δx ./ Lx
            # Calculate rate update --------------------------------------
            dVxdτ          .= (1-ρ) .* dVxdτ .+ Rx
            dVydτ          .= (1-ρ) .* dVydτ .+ Ry
            dPdτ           .= Rp
            # Update velocity and pressure -------------------------------
            Vx[:,2:end-1]  .+= Δτvx ./ ρ .* dVxdτ
            Vy[2:end-1,:]  .+= Δτvy ./ ρ .* dVydτ
            Pt             .+= κΔτp .* dPdτ
            # convergence check
            if mod(iter, nout)==0
                norm_Rx = norm(Rx)/sqrt(length(Rx)); norm_Ry = norm(Ry)/sqrt(length(Ry)); norm_∇V = norm(∇V)/sqrt(length(∇V))
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %03d, iter = %04d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, ∇V=%1.3e] \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V)
            end
            iter+=1; global itg=iter
        end
        t = t + Δt
        push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        # Plotting
        p1 = heatmap(xv, yc, Vx[:,2:end-1]', aspect_ratio=1, xlims=(0, Lx), ylims=(Δy/2, Ly-Δy/2), c=:inferno, title="Vx")
        p2 = heatmap(xc, yv, Vy[2:end-1,:]', aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        p3 = heatmap(xc, yc, ε̇xx' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="P")
        p4 = plot(evo_t, evo_τxx , legend=false, xlabel="time", ylabel="max(τxx)", linewiΔth=0, markershape=:circle, framestyle=:box, markersize=3)
        p4 = plot!(evo_t, 2.0.*εbg.*η0.*(1.0.-exp.(.-evo_t.*G./η0)), linewiΔth=2.0) # analytical solution
        display(plot(p1, p2, p3, p4))
    end
    return
end

Stokes2D_VE_benchmark()
