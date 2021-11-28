using Plots, LinearAlgebra, Printf
# helper functions
@views     av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views  av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views  av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end])
@views maxloc(A) = max.(A[1:end-2,1:end-2],A[1:end-2,2:end-1],A[1:end-2,3:end],
                        A[2:end-1,1:end-2],A[2:end-1,2:end-1],A[2:end-1,3:end],
                        A[3:end  ,1:end-2],A[3:end  ,2:end-1],A[3:end  ,3:end])
@views   bc2!(A) = begin A[1,:] = A[2,:]; A[end,:] = A[end-1,:]; A[:,1] = A[:,2]; A[:,end] = A[:,end-1] end
# main function
@views function Stokes2D_vep()
    use_vep = true
    # phyiscs
    lx, ly  = 1.0, 1.0
    radi    = 0.1
    η0      = 1.0
    η_reg   = 1.2e-2
    G0      = 1.0
    Gi      = 0.5G0
    τ_y     = 1.6
    sinϕ    = sind(30)
    ebg     = 1.0
    dt      = η0/G0/4.0
    # numerics
    nx, ny  = 63, 63
    nt      = 10
    εnl     = 1e-6
    maxiter = 150max(nx, ny)
    nchk    = 2max(nx, ny)
    Re      = 5π
    r       = 1.0
    CFL     = 0.99/sqrt(2)
    # preprocessing
    dx, dy  = lx/nx, ly/ny
    max_lxy = max(lx, ly)
    vpdτ    = CFL*min(dx, dy)
    xc, yc  = LinRange(-(lx-dx)/2,(lx-dx)/2,nx), LinRange(-(ly-dy)/2,(ly-dy)/2,ny)
    xv, yv  = LinRange(-lx/2,lx/2,nx+1), LinRange(-ly/2,ly/2,ny+1)
    # allocate arrays
    Pr      = zeros(nx  ,ny  )
    τxx     = zeros(nx  ,ny  )
    τyy     = zeros(nx  ,ny  )
    τxy     = zeros(nx+1,ny+1)
    τxyc    = zeros(nx  ,ny  )
    τii     = zeros(nx  ,ny  )
    Eii     = zeros(nx  ,ny  )
    λ       = zeros(nx  ,ny  )
    F       = zeros(nx  ,ny  )
    Fchk    = zeros(nx  ,ny  )
    dQdTxx  = zeros(nx  ,ny  )
    dQdTyy  = zeros(nx  ,ny  )
    dQdTxy  = zeros(nx  ,ny  )
    τxx_o   = zeros(nx  ,ny  )
    τyy_o   = zeros(nx  ,ny  )
    τxyc_o  = zeros(nx  ,ny  )
    τxy_o   = zeros(nx+1,ny+1)

    τxx_r   = zeros(nx  ,ny  )
    τyy_r   = zeros(nx  ,ny  )
    τxyc_r  = zeros(nx  ,ny  )
    τxy_r   = zeros(nx+1,ny+1)

    Vx      = zeros(nx+1,ny  )
    Vy      = zeros(nx  ,ny+1)
    dVx     = zeros(nx-1,ny  )
    dVy     = zeros(nx  ,ny-1)
    Rx      = zeros(nx-1,ny  )
    Ry      = zeros(nx  ,ny-1)
    ∇V      = zeros(nx  ,ny  )
    ρg      = zeros(nx  ,ny  )
    Exx     = zeros(nx  ,ny  )
    Eyy     = zeros(nx  ,ny  )
    Exyc    = zeros(nx  ,ny  )
    Exy     = zeros(nx+1,ny+1)
    Exx_e   = zeros(nx  ,ny  )
    Eyy_e   = zeros(nx  ,ny  )
    Exyc_e  = zeros(nx  ,ny  )
    Exy_e   = zeros(nx+1,ny+1)
    Exx_τ   = zeros(nx  ,ny  )
    Eyy_τ   = zeros(nx  ,ny  )
    Exy_τ   = zeros(nx+1,ny+1)
    Exyc_τ  = zeros(nx  ,ny  )
    η_ve_τ  = zeros(nx  ,ny  )
    η_ve_τv = zeros(nx+1,ny+1)
    η_ve    = zeros(nx  ,ny  )
    η_vem   = zeros(nx  ,ny  )
    η_vev   = zeros(nx+1,ny+1)
    η_vevm  = zeros(nx+1,ny+1)
    dτ_ρ    = zeros(nx  ,ny  )
    dτ_ρv   = zeros(nx+1,ny+1)
    Gdτ     = zeros(nx  ,ny  )
    Gdτv    = zeros(nx+1,ny+1)
    η_vep   = zeros(nx  ,ny  )
    η_vepv  = zeros(nx+1,ny+1)
    η_vec   = zeros(nx  ,ny  )
    # init
    Vx = [ ebg*x for x ∈ xv, _ ∈ yc ]
    Vy = [-ebg*y for _ ∈ xc, y ∈ yv ]
    η  = fill(η0,nx,ny); ηv = fill(η0,nx+1,ny+1)
    G  = fill(G0,nx,ny); Gv = fill(G0,nx+1,ny+1)
    @. G[xc^2 + yc'^2 < radi^2]  = Gi
    @. Gv[xv^2 + yv'^2 < radi^2] = Gi
    η_e   = G.*dt; η_ev = Gv.*dt
    @. η_ve  = 1.0/(1.0/η  + 1.0/η_e)
    @. η_vec = 1.0/(1.0/η  + 1.0/η_e)
    @. η_vev = 1.0/(1.0/ηv + 1.0/η_ev)
    @. η_vep  = 1.0/(1.0/η  + 1.0/η_e)
    @. η_vepv = 1.0/(1.0/ηv + 1.0/η_ev)
    # action
    t = 0.0; evo_t = Float64[]; evo_τxx = Float64[]; niter = 0
    for it = 1:nt
        τxx_o .= τxx; τyy_o .= τyy; τxy_o .= τxy; τxyc_o .= τxyc
        err  = 2εnl; iter = 0
        while err > εnl && iter < maxiter
            if !use_vep
                η_vem[2:end-1,2:end-1]  .= maxloc(η_ve) ; bc2!(η_vem)
                η_vevm[2:end-1,2:end-1] .= maxloc(η_vev); bc2!(η_vevm)
            else
                η_vem[2:end-1,2:end-1]  .= maxloc(η_vep) ; bc2!(η_vem)
                η_vevm[2:end-1,2:end-1] .= maxloc(η_vepv); bc2!(η_vevm)
            end
            @. dτ_ρ    = vpdτ*max_lxy/Re/η_vem
            @. dτ_ρv   = vpdτ*max_lxy/Re/η_vevm
            @. Gdτ     = vpdτ^2/dτ_ρ/(r+2.0)
            @. Gdτv    = vpdτ^2/dτ_ρv/(r+2.0)
            @. η_ve_τ  = 1.0/(1.0/η + 1.0/η_e + 1.0/Gdτ)
            @. η_ve_τv = 1.0/(1.0/ηv + 1.0/η_ev + 1.0/Gdτv)
            # pressure
            ∇V    .= diff(Vx, dims=1)./dx .+ diff(Vy, dims=2)./dy
            @. Pr -= r*Gdτ*∇V
            # strain rates
            Exx   .= diff(Vx, dims=1)./dx
            Eyy   .= diff(Vy, dims=2)./dy
            Exy[2:end-1,2:end-1] .= 0.5.*(diff(Vx[2:end-1,:], dims=2)./dy .+ diff(Vy[:,2:end-1], dims=1)./dx)
            Exyc  .= av(Exy)
            # viscoelastic strain rates
            @. Exx_e  = Exx  + τxx_o /2.0/η_e
            @. Eyy_e  = Eyy  + τyy_o /2.0/η_e
            @. Exy_e  = Exy  + τxy_o /2.0/η_ev
            @. Exyc_e = Exyc + τxyc_o/2.0/η_e
            # viscoelastic pseudo-transient strain rates
            @. Exx_τ  = Exx_e  + τxx /2.0/Gdτ
            @. Eyy_τ  = Eyy_e  + τyy /2.0/Gdτ
            @. Exy_τ  = Exy_e  + τxy /2.0/Gdτv
            @. Exyc_τ = Exyc_e + τxyc/2.0/Gdτ
            # stress update
            @. τxx    = 2.0*η_ve_τ *Exx_τ
            @. τyy    = 2.0*η_ve_τ *Eyy_τ
            @. τxy    = 2.0*η_ve_τv*Exy_τ
            @. τxyc   = 2.0*η_ve_τ *Exyc_τ
            # stress and strain rate invariants
            @. τii    = sqrt(0.5*(τxx^2 + τyy^2) + τxyc*τxyc)
            @. Eii    = sqrt(0.5*(Exx_τ^2 + Eyy_τ^2) + Exyc_τ*Exyc_τ)
            # yield function
            @. F      = τii - τ_y - Pr.*sinϕ
            @. λ      = max(F,0.0)/(η_ve_τ + η_reg)
            @. dQdTxx = 0.5*τxx /τii
            @. dQdTyy = 0.5*τyy /τii
            @. dQdTxy =     τxyc/τii
            # plastic correction
            @. τxx    = 2.0*η_ve_τ *(Exx_τ  - λ*dQdTxx)
            @. τyy    = 2.0*η_ve_τ *(Eyy_τ  - λ*dQdTyy)
            @. τxyc   = 2.0*η_ve_τ *(Exyc_τ - 0.5*λ*dQdTxy)
            τxy[2:end-1,2:end-1] .= 2.0 .* η_ve_τv[2:end-1,2:end-1].*(Exy_τ[2:end-1,2:end-1] .- 0.5 .* av(λ.*dQdTxy))
            @. τii    = sqrt(0.5*(τxx^2 + τyy^2) + τxyc*τxyc)
            @. Fchk   = τii - τ_y - Pr*sinϕ - λ*η_reg
            @. η_vep  = τii / 2.0 / Eii * 19.3 # nx, ny = 63, 63 
            # @. η_vep  = τii / 2.0 / Eii * 19.3 * 1.99 # nx, ny = 127, 127
            η_vepv[2:end-1,2:end-1] .= av(η_vep); bc2!(η_vep)
            # velocity update
            dVx .= av_xa(dτ_ρ) .* (.-diff(Pr, dims=1)./dx .+ diff(τxx, dims=1)./dx .+ diff(τxy[2:end-1,:], dims=2)./dy)
            dVy .= av_ya(dτ_ρ) .* (.-diff(Pr, dims=2)./dy .+ diff(τyy, dims=2)./dy .+ diff(τxy[:,2:end-1], dims=1)./dx .+ av_ya(ρg))
            @. Vx[2:end-1,:] = Vx[2:end-1,:] + dVx
            @. Vy[:,2:end-1] = Vy[:,2:end-1] + dVy
            if iter % nchk == 0
                @. τxx_r   = 2.0*η_vec*Exx_e
                @. τyy_r   = 2.0*η_vec*Eyy_e
                @. τxy_r   = 2.0*η_vev*Exy_e
                @. τxyc_r  = 2.0*η_vec*Exyc_e
                @. τii     = sqrt(0.5*(τxx_r^2 + τyy_r^2) + τxyc_r*τxyc_r)
                # yield function
                @. F      = τii - τ_y - Pr.*sinϕ
                @. λ      = max(F,0.0)/(η_vec + η_reg)
                @. dQdTxx = 0.5*τxx_r /τii
                @. dQdTyy = 0.5*τyy_r /τii
                @. dQdTxy =     τxyc_r/τii
                # plastic correction
                @. τxx_r  = 2.0*η_vec*(Exx_e  - λ*dQdTxx)
                @. τyy_r  = 2.0*η_vec*(Eyy_e  - λ*dQdTyy)
                @. τxyc_r = 2.0*η_vec*(Exyc_e - 0.5*λ*dQdTxy)
                τxy_r[2:end-1,2:end-1] .= 2.0 .* η_vev[2:end-1,2:end-1].*(Exy_e[2:end-1,2:end-1] .- 0.5 .* av(λ.*dQdTxy))
                Rx .= .-diff(Pr, dims=1)./dx .+ diff(τxx_r, dims=1)./dx .+ diff(τxy_r[2:end-1,:], dims=2)./dy
                Ry .= .-diff(Pr, dims=2)./dy .+ diff(τyy_r, dims=2)./dy .+ diff(τxy_r[:,2:end-1], dims=1)./dx .+ av_ya(ρg)
                norm_Rx = norm(Rx)/sqrt(length(Rx)); norm_Ry = norm(Ry)/sqrt(length(Ry)); norm_∇V = norm(∇V)/sqrt(length(∇V))
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                @printf("it = %d, iter = %d, err = %1.2e norm[Rx=%1.2e, Ry=%1.2e, ∇V=%1.2e] (Fchk=%1.2e) \n", it, iter, err, norm_Rx, norm_Ry, norm_∇V, maximum(Fchk))
            end
            iter += 1; niter += 1
        end
        println(norm(Exyc_τ))
        t += dt; push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        p1 = heatmap(xc,yc,τii',aspect_ratio=1,xlims=(-lx/2,lx/2),ylims=(-ly/2,ly/2),title="τii")
        # p3 = heatmap(xc,yc,η_vep',aspect_ratio=1,xlims=(-lx/2,lx/2),ylims=(-ly/2,ly/2),title="τii")
        p2 = plot(evo_t, evo_τxx , legend=false, xlabel="time", ylabel="max(τxx)", linewidth=0, markershape=:circle, framestyle=:box, markersize=3)
             plot!(evo_t, 2.0.*ebg.*η0.*(1.0.-exp.(.-evo_t.*G0./η0)), linewidth=2.0) # analytical solution for VE loading
             plot!(evo_t, 2.0.*ebg.*η0.*ones(size(evo_t)), linewidth=2.0)            # viscous flow stress
        display(plot(p1,p2))
    end
    println(niter)
    return
end
# action
Stokes2D_vep()
