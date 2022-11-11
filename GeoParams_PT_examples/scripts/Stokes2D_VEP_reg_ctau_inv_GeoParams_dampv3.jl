using ElasticArrays, Printf, GeoParams
using Plots,Plots.Measures
default(size=(800,500),framestyle=:box,label=false,grid=false,lw=6,labelfontsize=11,tickfontsize=11,titlefontsize=11)

@inline amean(a,b) = 0.5*(a + b)
@inline hmean(a,b) = 2.0/(1.0/a + 1.0/b)
@inline amean4(a,b,c,d) = 0.25*(a+b+c+d)
@inline hmean4(a,b,c,d) = 4.0/(1.0/a+1.0/b+1.0/c+1.0/d)
const av  = amean
const av4 = amean4
@views amean1(A)  = 0.5.*(A[1:end-1] .+ A[2:end])
@views avx(A)     = av.(A[1:end-1,:], A[2:end,:])
@views avy(A)     = av.(A[:,1:end-1], A[:,2:end])
@views avxy(A)    = av4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])
@views ameanx(A)  = amean.(A[1:end-1,:], A[2:end,:])
@views ameany(A)  = amean.(A[:,1:end-1], A[:,2:end])
@views ameanxy(A) = amean4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])
@views hmeanx(A)  = hmean.(A[1:end-1,:], A[2:end,:])
@views hmeany(A)  = hmean.(A[:,1:end-1], A[:,2:end])
@views hmeanxy(A) = hmean4.(A[1:end-1,1:end-1],A[2:end,1:end-1],A[1:end-1,2:end],A[2:end,2:end])
@views maxloc(A)  = max.(A[1:end-2,1:end-2],A[1:end-2,2:end-1],A[1:end-2,3:end],
                         A[2:end-1,1:end-2],A[2:end-1,2:end-1],A[2:end-1,3:end],
                         A[3:end  ,1:end-2],A[3:end  ,2:end-1],A[3:end  ,3:end])
@views bc2!(A)    = begin A[[1,end],:]=A[[2,end-1],:]; A[:,[1,end]]=A[:,[2,end-1]]; end

@views function update_old!(Txx_o,Tyy_o,Txy_o,Txyv_o,Pr_old,Txx,Tyy,Txy,Txyv,Pr)
    copyto!(Txx_o , Txx)
    copyto!(Tyy_o , Tyy)
    copyto!(Txy_o , Txy)
    copyto!(Txyv_o, Txyv)
    copyto!(Pr_old, Pr)
    # Pr      .= Pr_c
    # λ       .= 0.0
    return
end


@views function update_iteration_params!(η,ητ)
    ητ[2:end-1,2:end-1] .= maxloc(η); bc2!(ητ)
    return
end

# @views function update_stresses!((;Pr,Pr_old,Pr_c,dPr,K,Txx,Tyy,Txy,dτxx,dτyy,dTxyv,Txyv,Txx_o,Tyy_o,Txy_o,dQdτxx,dQdτyy,dQdTxyv,Fchk,TII,F,λ,Exx,Eyy,Exy,Exyv,Exx_ve,Eyy_ve,Exy_ve,EII_ve,η_vep,Vx,Vy,∇V,η,G,dτ_r),τ_y,sinϕ,sinψ,η_reg,relλ,dt,re_mech,vdτ,lτ,r,dx,dy,iter)
#     θ_dτ    = lτ*(r+2.0)/(re_mech*vdτ)
#     dτ_r   .= 1.0./(θ_dτ .+ η./(G.*dt) .+ 1.0)
#     ∇V     .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
#     # dPr    .= .-∇V
#     dPr    .= .-∇V .- (Pr .- Pr_old)./K./dt
#     # Pr    .+= (r/θ_dτ).*η.*dPr
#     Pr    .+= dPr./(1.0./(r/θ_dτ.*η) .+ 1.0./K./dt)
#     Exx    .= diff(Vx,dims=1)./dx .- ∇V./3.0
#     Eyy    .= diff(Vy,dims=2)./dy .- ∇V./3.0
#     Exyv[2:end-1,2:end-1] .= 0.5*(diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)
#     Exy    .= ameanxy(Exyv)
#     # visco-elastic strain rates
#     Exx_ve .= Exx .+ 0.5.*Txx_o./(G.*dt)
#     Eyy_ve .= Eyy .+ 0.5.*Tyy_o./(G.*dt)
#     Exy_ve .= Exy .+ 0.5.*Txy_o./(G.*dt)
#     EII_ve .= sqrt.(0.5.*(Exx_ve.^2 .+ Eyy_ve.^2) .+ Exy_ve.^2)
#     # stress increments
#     dτxx   .= (.-(Txx .- Txx_o).*η./(G.*dt) .- Txx .+ 2.0.*η.*Exx).*dτ_r
#     dτyy   .= (.-(Tyy .- Tyy_o).*η./(G.*dt) .- Tyy .+ 2.0.*η.*Eyy).*dτ_r
#     dTxyv   .= (.-(Txy .- Txy_o).*η./(G.*dt) .- Txy .+ 2.0.*η.*Exy).*dτ_r
#     TII    .= sqrt.(0.5.*((Txx.+dτxx).^2 .+ (Tyy.+dτyy).^2) .+ (Txy.+dTxyv).^2)
#     # # yield function
#     F      .= TII .- τ_y .- Pr.*sinϕ
#     if iter>100
#     λ      .= (1.0 .- relλ).*λ .+ relλ.*(max.(F,0.0)./(dτ_r.*η .+ η_reg .+ K.*dt.*sinϕ.*sinψ))
#     dQdτxx .= 0.5.*(Txx.+dτxx)./TII
#     dQdτyy .= 0.5.*(Tyy.+dτyy)./TII
#     dQdTxyv .=      (Txy.+dTxyv)./TII
#     end
#     Pr_c   .= Pr .+ K.*dt.*λ.*sinψ
#     Txx   .+= (.-(Txx .- Txx_o).*η./(G.*dt) .- Txx .+ 2.0.*η.*(Exx .-      λ.*dQdτxx)).*dτ_r
#     Tyy   .+= (.-(Tyy .- Tyy_o).*η./(G.*dt) .- Tyy .+ 2.0.*η.*(Eyy .-      λ.*dQdτyy)).*dτ_r
#     Txy   .+= (.-(Txy .- Txy_o).*η./(G.*dt) .- Txy .+ 2.0.*η.*(Exy .- 0.5.*λ.*dQdTxyv)).*dτ_r
#     Txyv[2:end-1,2:end-1] .= ameanxy(Txy)
#     TII    .= sqrt.(0.5.*(Txx.^2 .+ Tyy.^2) .+ Txy.^2)
#     Fchk   .= TII .- τ_y .- Pr_c.*sinϕ .- λ.*η_reg
#     η_vep  .= TII ./ 2.0 ./ EII_ve
#     return
# end


function update_stress_GP2!(Txx, Tyy, Txy, Txyv, Tii, Txx_o, Tyy_o, Txy_o, Txyv_o, Exx, Eyy, Exyv, η, η_vep, Pt, Phasec, Phasev, MatParam, dt, G, lτ, r, re_mech, vdτ)
    θ_dτ    = lτ*(r+2.0)/(re_mech*vdτ)
    for j in axes(Exx,2), i in axes(Exx,1)
        # dτ_r = 1.0/(θ_dτ + η[i,j]/(G[i,j]*dt) + 1.0)
        dτ_r = 1.0/(θ_dτ + 1/η_vep[i,j]) # equivalent to dτ_r = @. 1.0/(θ_dτ + η/(G*dt) + 1.0)
        args = (; dt=dt, P=Pt[i,j], τII_old=0.0)
        # # gather strain rate
        εij_v = (
            Exyv[i  ,j  ], 
            Exyv[i+1,j  ], 
            Exyv[i  ,j+1], 
            Exyv[i+1,j+1]
        ) # gather vertices around ij center
        
        εij_p = (Exx[i,j], Eyy[i,j], εij_v)
        # # gather deviatoric stress
        # τij_v = (
        #     Txyv_o[i  ,j  ], 
        #     Txyv_o[i+1,j  ], 
        #     Txyv_o[i  ,j+1], 
        #     Txyv_o[i+1,j+1]
        # ) # gather vertices around ij center
        # τij_p_o = (Txx_o[i,j], Tyy_o[i,j], τij_v)
        τij_v = (
            Txy_o[i,j], 
            Txy_o[i,j], 
            Txy_o[i,j], 
            Txy_o[i,j]
        )
        τij_p_o = (Txx_o[i,j], Tyy_o[i,j], τij_v)
        # gathermaterial phases
        phases_v = (
            Phasev[i  ,j  ], 
            Phasev[i+1,j  ], 
            Phasev[i  ,j+1], 
            Phasev[i+1,j+1]
        ) # gather vertices around ij center
        phases = (Phasec[i,j], Phasec[i,j], phases_v)
        # update stress and effective viscosity
        Tij, Tii[i,j], η_vep[i,j] = compute_τij(MatParam, εij_p, args, τij_p_o, phases)
        Exy = sum(εij_v)*0.25
        xyv_o = sum(τij_v)*0.25

        # need to compute τij/η_ve(p?) in
        # Txx[i,j] += dτ_r * (-(Txx[i,j] - Txx_o[i,j])/(G[i,j]*dt) - Txx[i,j]/η[i,j] + 2*Exx[i,j])
        # Tyy[i,j] += dτ_r * (-(Tyy[i,j] - Tyy_o[i,j])/(G[i,j]*dt) - Tyy[i,j]/η[i,j] + 2*Eyy[i,j])
        # Txy[i,j] += dτ_r * (-(Txy[i,j] -      xyv_o)/(G[i,j]*dt) - Txy[i,j]/η[i,j] + 2*sum(εij_v)/4)

        Txx[i,j] += dτ_r * (-Txx[i,j] /(G[i,j]*dt) - Txx[i,j]/η[i,j] + Tij[1] / η_vep[i,j])
        Tyy[i,j] += dτ_r * (-Tyy[i,j] /(G[i,j]*dt) - Tyy[i,j]/η[i,j] + Tij[2] / η_vep[i,j])
        Txy[i,j] += dτ_r * (-Txy[i,j] /(G[i,j]*dt) - Txy[i,j]/η[i,j] + Tij[3] / η_vep[i,j])

        # Txx[i,j] += dτ_r * (-(Txx[i,j])/(G[i,j]*dt) - Txx[i,j]/η[i,j] + 2*ε_eff[1])
        # Tyy[i,j] += dτ_r * (-(Tyy[i,j])/(G[i,j]*dt) - Tyy[i,j]/η[i,j] + 2*ε_eff[2])
        # Txy[i,j] += dτ_r * (-(Txy[i,j])/(G[i,j]*dt) - Txy[i,j]/η[i,j] + 2*ε_eff[3])

        # Txx[i,j] += dτ_r * (-(Txx[i,j])/(G[i,j]*dt) - Txx[i,j]/η[i,j] + Tij[1] / η_vep[i,j])
        # Tyy[i,j] += dτ_r * (-(Tyy[i,j])/(G[i,j]*dt) - Tyy[i,j]/η[i,j] + Tij[2] / η_vep[i,j])
        # Txy[i,j] += dτ_r * (-(Txy[i,j])/(G[i,j]*dt) - Txy[i,j]/η[i,j] + Tij[3] / η_vep[i,j])

        # Txx[i,j] += dτ_r * (-(Txx[i,j]) + Tij[1] ) / η_vep[i,j] # NOTE: from GP Tij = 2*η_vep * εij
        # Tyy[i,j] += dτ_r * (-(Tyy[i,j]) + Tij[2] ) / η_vep[i,j] 
        # Txy[i,j] += dτ_r * (-(Txy[i,j]) + Tij[3] ) / η_vep[i,j] 
    end
    Txyv[2:end-1,2:end-1] .= ameanxy(Txy)

end

@views function update_stresses!(Pr,Pr_old,dPr,K,Txx,Tyy,Txy,Txyv,Txx_o,Tyy_o,Txy_o,TII,Exx,Eyy,Exy,Exyv,Exx_ve,Eyy_ve,Exy_ve,EII_ve,η_vep,Vx,Vy,∇V,η,G,dτ_r,dt,re_mech,vdτ,lτ,r,dx,dy)
    θ_dτ    = lτ*(r+2.0)/(re_mech*vdτ)
    dτ_r   .= 1.0./(θ_dτ .+ η./(G.*dt) .+ 1.0)
    ∇V     .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    # dPr    .= .-∇V
    dPr    .= .-∇V .- (Pr .- Pr_old)./K./dt
    # Pr    .+= (r/θ_dτ).*η.*dPr
    Pr    .+= dPr./(1.0./(r/θ_dτ.*η) .+ 1.0./K./dt)
    Exx    .= diff(Vx,dims=1)./dx .- ∇V./3.0
    Eyy    .= diff(Vy,dims=2)./dy .- ∇V./3.0
    Exyv[2:end-1,2:end-1] .= 0.5*(diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)
    Exy    .= ameanxy(Exyv)
    # visco-elastic strain rates
    Exx_ve .= Exx .+ 0.5.*Txx_o./(G.*dt)
    Eyy_ve .= Eyy .+ 0.5.*Tyy_o./(G.*dt)
    Exy_ve .= Exy .+ 0.5.*Txy_o./(G.*dt)
    EII_ve .= sqrt.(0.5.*(Exx_ve.^2 .+ Eyy_ve.^2) .+ Exy_ve.^2)
    # stress increments
    Txx   .+= (.-(Txx .- Txx_o)./(G.*dt) .- Txx./η .+ 2.0.*Exx).*dτ_r
    Tyy   .+= (.-(Tyy .- Tyy_o)./(G.*dt) .- Tyy./η .+ 2.0.*Eyy).*dτ_r
    Txy   .+= (.-(Txy .- Txy_o)./(G.*dt) .- Txy./η .+ 2.0.*Exy).*dτ_r
    Txyv[2:end-1,2:end-1] .= ameanxy(Txy)
    TII   .= sqrt.(0.5.*((Txx).^2 .+ (Tyy).^2) .+ (Txy).^2)
    η_vep .= TII ./ 2.0 ./ EII_ve
    return
end

@views function update_P_strain!(Pr,Pr_old,dPr,K,Exx,Eyy,Exy,Exyv,Vx,Vy,∇V,η,G,dτ_r,dt,re_mech,vdτ,lτ,r,dx,dy)
    θ_dτ    = lτ*(r+2.0)/(re_mech*vdτ)
    dτ_r   .= 1.0./(θ_dτ .+ η./(G.*dt) .+ 1.0)
    ∇V     .= diff(Vx,dims=1)./dx .+ diff(Vy,dims=2)./dy
    # dPr    .= .-∇V
    dPr    .= .-∇V .- (Pr .- Pr_old)./K./dt
    # Pr    .+= (r/θ_dτ).*η.*dPr
    Pr    .+= dPr./(1.0./(r/θ_dτ.*η) .+ 1.0./K./dt)
    Exx    .= diff(Vx,dims=1)./dx .- ∇V./3.0
    Eyy    .= diff(Vy,dims=2)./dy .- ∇V./3.0
    Exyv[2:end-1,2:end-1] .= 0.5*(diff(Vx[2:end-1,:],dims=2)./dy .+ diff(Vy[:,2:end-1],dims=1)./dx)
    Exy    .= ameanxy(Exyv)
    return
end

@views function update_velocities!(Vx,Vy,Pr_c,Txx,Tyy,Txyv,ητ,ρgx,ρgy,vdτ,lτ,re_mech,dx,dy)
    ηdτ = vdτ*lτ/re_mech
    Vx[2:end-1,:] .+= (diff(.-Pr_c.+Txx,dims=1)./dx .+ diff(Txyv[2:end-1,:],dims=2)./dy .- ρgx).*ηdτ./avx(ητ)
    Vy[:,2:end-1] .+= (diff(.-Pr_c.+Tyy,dims=2)./dy .+ diff(Txyv[:,2:end-1],dims=1)./dx .- ρgy).*ηdτ./avy(ητ)
    return
end

@views function compute_residuals!(r_Vx,r_Vy,Pr_c,Txx,Tyy,Txyv,ρgx,ρgy,dx,dy)
    r_Vx .= diff(.-Pr_c[:,2:end-1].+Txx[:,2:end-1],dims=1)./dx .+ diff(Txyv[2:end-1,2:end-1],dims=2)./dy .- ρgx[:,2:end-1]
    r_Vy .= diff(.-Pr_c[2:end-1,:].+Tyy[2:end-1,:],dims=2)./dy .+ diff(Txyv[2:end-1,2:end-1],dims=1)./dx .- ρgy[2:end-1,:]
    return
end

function main(UseGP, nt)
    pl_correction = :native_naive
    pl_correction = :native_inv1
    pl_correction = :native_inv2
    pl_correction = :native_inv3

    #gp_correction = :loop
    gp_correction = :native_gp
    # gp_correction = :native_gp_dilation
    
    do_DP   = true               # do_DP=false: Von Mises, do_DP=true: Drucker-Prager (friction angle)
    η_reg   = 8.0e-3             # regularisation "viscosity"
    # Physics
    Lx, Ly  = 1.0, 1.0           # domain size
    radi    = 0.01               # inclusion radius
    τ_y     = 1.6                # yield stress. If do_DP=true, τ_y stand for the cohesion: c*cos(ϕ)
    ϕ       = 30.0
    sinϕ    = sind(ϕ)*do_DP      # sinus of the friction angle
    μ0      = 1.0                # viscous viscosity
    G0      = 1.0                # elastic shear modulus
    Gi      = G0/(6.0-4.0*do_DP) # elastic shear modulus perturbation
    εbg     = 1.0                # background strain-rate
    Coh     = τ_y/cosd(ϕ)        # cohesion
    # Coh     = Inf

    # Geoparams initialisation
    pl = DruckerPrager_regularised(C=Coh, ϕ=ϕ, η_vp=η_reg, Ψ=0)        # non-regularized plasticity
    MatParam = (SetMaterialParams(Name="Matrix"   , Phase=1,
                CompositeRheology = CompositeRheology(SetConstantElasticity(G=G0, ν=0.5),LinearViscous(η=μ0), pl)), 
                SetMaterialParams(Name="Inclusion", Phase=2,
                CompositeRheology = CompositeRheology(SetConstantElasticity(G=Gi, ν=0.5),LinearViscous(η=μ0), pl)),
                )
    # Numerics
    # nt      = 10               # number of time steps
    nx, ny  = 63, 63             # numerical grid resolution
    Vdmp    = 4.0                # convergence acceleration (damping)
    Vsc     = 2.0                # iterative time step limiter
    Ptsc    = 6.0                # iterative time step limiter
    ε       = 1e-6               # nonlinear tolerence
    iterMax = 3e4                # max number of iters
    nout    = 200                # check frequency
    re_mech = 3π
    lτ      = min(Lx,Ly)
    dx, dy  = Lx/nx, Ly/ny
    vdτ     = 0.5*min(dx,dy)/√2.1
    r       = 0.7
    # Preprocessing
    dx, dy  = Lx/nx, Ly/ny
    dt      = μ0/G0/4.0 # assumes Maxwell time of 4
    # Array initialisation
    Pt      = zeros(nx  ,ny  )
    dPt     = zeros(nx  ,ny  )
    Pt_c    = zeros(nx  ,ny  )
    P_o     = zeros(nx,  ny  )
    ∇V      = zeros(nx  ,ny  )
    Vx      = zeros(nx+1,ny  )
    Vy      = zeros(nx  ,ny+1)
    Exx     = zeros(nx  ,ny  )
    Eyy     = zeros(nx  ,ny  )
    Exy     = zeros(nx  ,ny  )
    Exx_ve  = zeros(nx  ,ny  )
    Eyy_ve  = zeros(nx  ,ny  )
    Exy_ve  = zeros(nx  ,ny  )
    EII_ve  = zeros(nx  ,ny  )
    Exyv    = zeros(nx+1,ny+1)
    Txx     = zeros(nx  ,ny  )
    Tyy     = zeros(nx  ,ny  )
    Txy     = zeros(nx  ,ny  )
    Txyv    = zeros(nx+1,ny+1)
    Txx_o   = zeros(nx  ,ny  )
    TII_o   = zeros(nx  ,ny  )
    Tyy_o   = zeros(nx  ,ny  )
    Txy_o   = zeros(nx  ,ny  )
    Txyv_o  = zeros(nx+1,ny+1)
    TII     = zeros(nx  ,ny  )
    EII     = zeros(nx  ,ny  )
    EII_f   = zeros(nx  ,ny  )
    F       = zeros(nx  ,ny  )
    Fchk    = zeros(nx  ,ny  )
    Pla     = zeros(nx  ,ny  )
    λ       = zeros(nx  ,ny  )
    dQdTxx  = zeros(nx  ,ny  )
    dQdTyy  = zeros(nx  ,ny  )
    dQdTxy  = zeros(nx  ,ny  )    
    r_Vx    = zeros(nx-1,ny-2)
    r_Vy    = zeros(nx-2,ny-1)
    η_v     =    μ0*ones(nx, ny)
    η_e     = dt*G0*ones(nx, ny)
    η_ev    = dt*G0*ones(nx+1, ny+1)
    η_vep   =       ones(nx, ny)
    ητ      = zeros(nx  ,ny  )
    dτ_r    = zeros(nx  ,ny  )

    η       = ones(nx,ny)
    Phasec  = ones(Int, nx  ,ny  )
    Phasev  = ones(Int, nx+1,ny+1)
    # Initial condition
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xc, yc  = LinRange(dx/2, Lx-dx/2, nx), LinRange(dy/2, Ly-dy/2, ny)
    xv, yv  = LinRange(0.0, Lx, nx+1), LinRange(0.0, Ly, ny+1)
    Xvx     = [x for x=xv,y=yc]
    Yvy     = [y for x=xc,y=yv]
    radc    = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv    = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    η_e[radc.<radi]    .= dt*Gi
    η_ev[radv.<radi]   .= dt*Gi
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    η_vep    .= (1.0./η_e + 1.0./η_v).^-1
    Vx      .=   εbg.*Xvx
    Vy      .= .-εbg.*Yvy

    Kb = fill(Inf, nx, ny)
    G  = fill(G0,nx,ny)
    G[radc.<radi] .= Gi
    ρgx = zeros(nx-1,ny)
    ρgy = zeros(nx,ny-1)
    Pt_old = deepcopy(Pt)
    ϵtol = 1e-6

    # Time loop
    t=0.0; evo_t=Float64[]; evo_Txx=Float64[]

    iter_evo = Float64[]
    errs_evo = ElasticMatrix{Float64}(undef,length(ϵtol),0)
    opts = (aspect_ratio=1, xlims=extrema(xc), ylims=extrema(yc), c=:turbo, framestyle=:box)
    # mask = copy(fields.phase); mask[mask.<0.7].=NaN
    t = 0.0; evo_t=[]; evo_τxx=[]
    # ispath("anim")&&rm("anim",recursive=true);mkdir("anim");iframe = -1
    # time loop
    # nt = 10
    maxiter = 10e3
    ncheck = 100
    for it = 1:nt
        @printf("it=%d\n",it)
        update_old!(Txx_o,Tyy_o,Txy_o,Txyv_o,Pt_old,Txx,Tyy,Txy,Txy,Pt)
        errs = 2.0.*ϵtol; iter = 1
        resize!(iter_evo,0); resize!(errs_evo,length(ϵtol),0)
        while any(errs .>= ϵtol) && iter <= maxiter
            update_iteration_params!(η,ητ)

            if UseGP
                update_P_strain!(Pt,Pt_old,dPt,Kb,Exx,Eyy,Exy,Exyv,Vx,Vy,∇V,η,G,dτ_r,dt,re_mech,vdτ,lτ,r,dx,dy)
                update_stress_GP2!(Txx, Tyy, Txy, Txyv, TII, Txx_o, Tyy_o, Txy_o, Txyv_o, Exx, Eyy, Exyv, η, η_vep, Pt, Phasec, Phasev, MatParam, dt,G, lτ, r, re_mech, vdτ)
            else
                update_stresses!(Pt,Pt_old,dPt,Kb,Txx,Tyy,Txy,Txyv,Txx_o,Tyy_o,Txy_o,TII,Exx,Eyy,Exy,Exyv,Exx_ve,Eyy_ve,Exy_ve,EII_ve,η_vep,Vx,Vy,∇V,η,G,dτ_r,dt,re_mech,vdτ,lτ,r,dx,dy)
            end                   
            update_velocities!(Vx,Vy,Pt,Txx,Tyy,Txyv,ητ,ρgx,ρgy,vdτ,lτ,re_mech,dx,dy)
            
            # update_velocities!(fields,vdτ,lτ,re_mech,dx,dy)
            if iter % ncheck == 0
                # update residuals
                compute_residuals!(r_Vx,r_Vy,Pt,Txx,Tyy,Txyv,ρgx,ρgy,dx,dy)

                errs = maximum.((abs.(r_Vx),abs.(r_Vy),abs.(dPt)))
                push!(iter_evo,iter/max(nx,ny));append!(errs_evo,errs)
                @printf("  iter/nx=%.3f,errs=[ %1.3e, %1.3e, %1.3e ] \n",iter/max(nx,ny),errs...)
            end
            iter += 1
        end
        t += dt
        push!(evo_t,t); push!(evo_τxx,maximum(Txx))
        # visualisation
        # mask .= fields.phase; mask[mask.<0.7].=NaN
        # fields.Vmag .= sqrt.(ameanx(fields.Vx).^2 + ameany(fields.Vy).^2)
        # fields.TII  .= sqrt.(0.5.*(fields.Txx.^2 .+ fields.Tyy.^2) .+ fields.Txy.^2)
        # # p1=heatmap(xc,yc,log10.(fields.η)',title="log10(η)";opts...)
        # p1=heatmap(xc,yc,mask' .* fields.Pr',title="Pressure";opts...)
        # p3=heatmap(xc,yc,mask' .* fields.TII',title="TII";opts...)
        # p2=heatmap(xc,yc,mask' .* fields.η_vep',title="η_vep";opts...)
        # # p2=heatmap(xc,yc,#=mask' .*=# fields.F',title="F";opts...)
        # # p3=heatmap(xc[2:end-1],yc[2:end-1],mask[2:end-1,2:end-1]' .* ameany(fields.r_Vy)',title="Vmag";opts...)
        # p4=heatmap(xc,yc,mask' .* fields.Vmag',title="Vmag";opts...)
        # # p4=plot(evo_t,evo_τxx,legend=false,xlabel="time",ylabel="max(Txx)",linewidth=0,markershape=:circle,markersize=3,framestyle=:box)
        # display(plot(p1,p2,p3,p4,layout=(2,2)))
        # png(plot(p1,p2,p3,p4,layout=(2,2)),@sprintf("anim/%04d.png",iframe+=1))

        TII  .= sqrt.(0.5.*(Txx.^2 .+ Tyy.^2) .+ Txy.^2)
        p2=heatmap(xc,yc, Txy',title="η_vep")
        display(p2)
    end
    return evo_t, evo_τxx, iter_evo, TII, η_vep
end

@time evo_t, evo_τxx_GP, iter_evo_GP, TII, η_vep = main(true, 5)
@time evo_t, evo_τxx, iter_evo, TII, η_vep = main(false, 5)

sol = @. 2.0*1*1*(1.0-exp(-evo_t*1/1))
scatter(evo_t, sol)
plot!(evo_t, evo_τxx_GP, color=:red)
plot!(evo_t, evo_τxx, color=:black)

plot(evo_t, @.((evo_τxx_GP-evo_τxx)/evo_τxx*100), color=:black)
