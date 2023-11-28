# Initialisation
using Plots, Printf, Statistics, LinearAlgebra, GeoParams
const Dat = Float64 # Precision (double=Float64 or single=Float32)
# Macros
@views    av(A) = 0.25*(A[1:end-1,1:end-1].+A[2:end,1:end-1].+A[1:end-1,2:end].+A[2:end,2:end])
@views av_xa(A) =  0.5*(A[1:end-1,:].+A[2:end,:])
@views av_ya(A) =  0.5*(A[:,1:end-1].+A[:,2:end]) 

@inline function av_c(A, i, j)
    return 0.25*(A[i,j] + A[i+1,j] + A[i,j+1] + A[i+1,j+1])
end

# Average square of the 4 vertex points
@inline function av2_c(A, i, j)
    return 0.25*(A[i,j]^2 + A[i+1,j]^2 + A[i,j+1]^2 + A[i+1,j+1]^2)
end

@inline function av_v(A, i, j)
    return 0.25*(A[i,j] + A[i-1,j] + A[i,j-1] + A[i-1,j-1])
end

@inline function av2_v(A, i, j)
    return 0.25*(A[i,j]^2 + A[i-1,j]^2 + A[i,j-1]^2 + A[i-1,j-1]^2)
end


@views function average_c!(Ac::AbstractArray{N,T}, Av::AbstractArray{N,T}) where {N,T} # extrapolate vertex -> center
 
    for j ∈ axes(Ac,2), i ∈ axes(Ac,1)
        Ac[i,j] = av_c(Av,i,j)
    end

    return nothing
end


@views function cen2ver!(Av::Array{N,T}, Ac::Array{N,T})  where {N,T} # extrapolate center -> vertex
    average_c!( Av[2:end-1,2:end-1], Ac)
    Av[1,:].=Av[2,:]; 
    Av[end,:].=Av[end-1,:]; 
    Av[:,1].=Av[:,2]; 
    Av[:,end].=Av[:,end-1]
    return nothing
end

@generated function phase_viscosity(v::NTuple{N,Any}, ε̇ii, phase, args) where N
    quote
        Base.@_inline_meta
        Base.@nexprs $N i -> v[i].Phase === phase && return computeViscosity_εII(v[i].CompositeRheology[1], ε̇ii, args)        
    end
end
# Rheology
@views function UpdateStressGeoParams!( ηc, ηv,  τii, τxx, τyy, τxy, τxy_c, ε̇xx, ε̇yy, ε̇xy, ε̇xy_c, ε̇iic, ε̇iiv, τxx0, τyy0, τxy0, τii0c, τii0v, Pt, MatParam, Δt, Phasec, Phasev )
    # ε̇iic                   .= sqrt.(1//2 .*(ε̇xx.^2 .+ ε̇yy.^2) .+ av(ε̇xy).^2)
    # ε̇iiv[2:end-1,2:end-1]  .= sqrt.(1//2 .*( av(ε̇xx).^2 .+ av(ε̇yy).^2) .+ ε̇xy[2:end-1,2:end-1].^2)
    # τii0c                  .= sqrt.(1//2 .*(τxx0.^2 .+ τyy0.^2) .+ av(τxy0).^2)
    # τii0v[2:end-1,2:end-1] .= sqrt.(1//2 .*( av(τxx0).^2 .+ av(τyy0).^2) .+ τxy0[2:end-1,2:end-1].^2)

    average_c!(ε̇xy_c, ε̇xy)  # average vertices -> centers

    # Centroids
    @inbounds for j ∈ axes(ε̇xx,2), i ∈ axes(ε̇xx,1)
        # compute second invariants from surrounding points
        τii0     =  sqrt(0.5 *(τxx0[i,j]^2 + τyy0[i,j]^2) + av2_c(τxy0,i,j))
        ε̇ii      =  sqrt(0.5 *( ε̇xx[i,j]^2 +  ε̇yy[i,j]^2) + av2_c(ε̇xy,i,j))

        args       = (; τII_old = τii0, dt=Δt, P=Pt[i,j])             
        
        η = ηc[i,j] = phase_viscosity(MatParam, ε̇ii, Phasec[i,j], args)
        τxx[i,j]   = 2*η*ε̇xx[i,j]
        τyy[i,j]   = 2*η*ε̇yy[i,j]
        τxy_c[i,j] = 2*η*ε̇xy_c[i,j]

        τii[i,j]   = 2*η*ε̇ii                # mostly for debugging
    end
    cen2ver!(ηv, ηc)    # extrapolate from centers -> vertices

    # Vertices
   # @inbounds for j ∈ 2:size(ε̇xy,2)-1, i ∈ 2:size(ε̇xy,1)-1
   #     τxy[i,j] = 2*ηv[i,j]*ε̇xy[i,j] 
   # end
    cen2ver!(τxy, τxy_c)    # extrapolate from centers -> vertices

    
    return nothing
end

@views function UpdateStressNative_VE!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0)
    
    # purely viscoelastic correction
    τxx   .= 2 .* ηve_c .* ( ε̇xx .+ τxx0./(2 .* ηe_c) ) 
    τyy   .= 2 .* ηve_c .* ( ε̇yy .+ τyy0./(2 .* ηe_c) )
    τxy   .= 2 .* ηve_v .* ( ε̇xy .+ τxy0./(2 .* ηe_v) )

    return nothing
end

@views function UpdateStressNative_VEP_invariants!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0,
    ε̇ii, τii, τii0, ε̇xy_c, τxy_c,ε̇xy2_c, τxy2_c, F, τ_y, Pt, sinϕ, η_reg, Pla, λ, ε̇ii_pl, Fchk, ηc, ηv)
    
    # visco-elastic strain rates
    average_c!(τxy2_c, τxy0.^2)     # average squares!
    average_c!(ε̇xy2_c, ε̇xy.^2 )
    
    τii0   .=  sqrt.( 0.5.*(τxx0.^2 .+ τyy0.^2) .+ τxy2_c)  # old stress invariant @ cente
    ε̇ii    .=  sqrt.( 0.5.*(ε̇xx.^2 .+ ε̇yy.^2) .+ ε̇xy2_c)    # strain rate invariant @ center
    τii    .=  2.0.*ηve_c.*(ε̇ii .+ τii0./(2.0.*ηe_c) )      # trial stress
   
    average_c!(ε̇xy_c, ε̇xy)

    # Compute plasticity @ center
    F      .= τii .- τ_y .- Pt.*sinϕ
    Pla    .= 0.0
    Pla    .= F .> 0.0
    λ      .= Pla.*F./(ηve_c .+ η_reg)
    ε̇ii_pl .= λ*0.5                     # 2nd invariant of plastic strainrate (check!)
    τii    .= 2.0.*ηve_c.*(ε̇ii .+ τii0./(2.0*ηe_c)  .- ε̇ii_pl) # updated stress @ center
    Fchk   .= τii .- τ_y .- Pt.*sinϕ .- λ.*η_reg

    # Update stress components
    ηc      .=  τii./(2.0*ε̇ii)
    τxx     .=  2.0.*ηc.*ε̇xx 
    τyy     .=  2.0.*ηc.*ε̇yy 
    τxy_c   .=  2.0.*ηc.*ε̇xy_c 

    cen2ver!(ηv, ηc)                # extrapolate from centers -> vertices
    cen2ver!(τxy, τxy_c)                # extrapolate from centers -> vertices

#    τxy     .=  2.0.*ηv.*ε̇xy        # update stress @ vertexes

    return nothing
end


@views function UpdateStressNative_VEP!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0,
    ε̇xx1, ε̇yy1, ε̇xy1, ε̇xy1_c, τxy0_c, ε̇ii, τxy_c, τii, F, τ_y, Pt, sinϕ, η_reg, Pla, λ, dQdTxx, dQdTyy, dQdTxy, Fchk, ηvep_c, ηvep_v)
    
    # visco-elastic strain rates
    ε̇xx1   .=    ε̇xx  .+ τxx0  ./2.0./ηe_c
    ε̇yy1   .=    ε̇yy  .+ τyy0  ./2.0./ηe_c
    ε̇xy1   .=    ε̇xy  .+ τxy0  ./2.0./ηe_v
    ε̇xy1_c .= av(ε̇xy) .+ τxy0_c./2.0./ηe_c
    ε̇ii    .= sqrt.(0.5*(ε̇xx1.^2 .+ ε̇yy1.^2) .+ ε̇xy1_c.^2)
    
    # trial stress
    τxx    .= 2.0.*ηve_c.*ε̇xx1
    τyy    .= 2.0.*ηve_c.*ε̇yy1
    τxy_c  .= 2.0.*ηve_c.*ε̇xy1_c
    τii    .= sqrt.(0.5*(τxx.^2 .+ τyy.^2) .+ τxy_c.^2)
    
    
    # yield function
    F      .= τii .- τ_y .- Pt.*sinϕ
    Pla    .= 0.0
    Pla    .= F .> 0.0
    λ      .= Pla.*F./(ηve_c .+ η_reg)
    dQdTxx .= 0.5.*τxx./τii
    dQdTyy .= 0.5.*τyy./τii
    dQdTxy .=    τxy_c./τii 

    # plastic corrections
    τxx    .= 2.0.*ηve_c.*(ε̇xx1   .-      λ.*dQdTxx)
    τyy    .= 2.0.*ηve_c.*(ε̇yy1   .-      λ.*dQdTyy)
    τxy_c  .= 2.0.*ηve_c.*(ε̇xy1_c .- 0.5.*λ.*dQdTxy)
    τii    .= sqrt.(0.5*(τxx.^2 .+ τyy.^2) .+ τxy_c.^2)
    Fchk   .= τii .- τ_y .- Pt.*sinϕ .- λ.*η_reg
    
    ηvep_c .= τii./2.0./ε̇ii
    cen2ver!(ηvep_v, ηvep_c)

    τxy   .= 2.0.*ηvep_v.*ε̇xy1

    # purely viscoelastic correction
   # τxx   .= 2 .* ηve_c .* ( ε̇xx .+ τxx0./(2 .* ηe_c) ) 
   # τyy   .= 2 .* ηve_c .* ( ε̇yy .+ τyy0./(2 .* ηe_c) )
   # τxy   .= 2 .* ηve_v .* ( ε̇xy1 )

    return nothing
end


# 2D Stokes routine
@views function Stokes2D_VE_inclusion(UseGeoParams, doPlot = false)
    # Physics
    do_DP   = true
    Lx, Ly  = 1.0, 1.0  # domain size
    ξ       = 10.0      # Maxwell relaxation time
    η0      = 1.0       # viscous viscosity
    G0      = 1.0       # elastic shear modulus
    εbg     = 1.0       # background strain-rate
    radi    = 0.01
    τ_y     = 1.6
    Gi      = G0/(6.0-4.0*do_DP)      # inclusion shear modulus
    η_reg   = 2*1.2e-2            # regularisation "viscosity"
    ϕ       = 30*do_DP    
    sinϕ    = sind(ϕ)*do_DP      
    Coh     = τ_y/cosd(ϕ)      # cohesion
    
    #pl = DruckerPrager(C=Coh, ϕ=ϕ)        # non-regularized plasticity
    #pl = Parallel(DruckerPrager(C=Coh, ϕ=ϕ), LinearViscous(η=η_reg))
    pl = DruckerPrager_regularised(C=Coh, ϕ=ϕ, η_vp=η_reg)        # non-regularized plasticity
    
    MatParam = (SetMaterialParams(Name="Matrix"   , Phase=1,
                CompositeRheology = CompositeRheology(ConstantElasticity(G=G0),LinearViscous(η=η0),pl)), 
                SetMaterialParams(Name="Inclusion", Phase=2,
                CompositeRheology = CompositeRheology(ConstantElasticity(G=Gi),LinearViscous(η=η0),pl)),
                )

    # Numerics
    nt       = 20        # number of time steps
    # ncx, ncy = 31, 31    # numerical grid resolution
    ncx, ncy = 63, 63    # numerical grid resolution
    
    ε        = 1e-6      # nonlinear tolerence
    iterMax  = 1e5       # max number of iters
    nout     = 500       # check frequency
    # Iterative parameters -------------------------------------------
    Reopt    = 4π
    cfl      = 0.50
    ρ        = cfl*Reopt/ncx
    # Preprocessing
    Δx, Δy   = Lx/ncx, Ly/ncy
    Δt       = η0/(G0*ξ + 1e-15) 
    # Array initialisation
    Pt       = zeros(Dat, ncx  ,ncy  )
    ∇V       = zeros(Dat, ncx  ,ncy  )
    Vx       = zeros(Dat, ncx+1,ncy+2)
    Vy       = zeros(Dat, ncx+2,ncy+1)
    ε̇xx      = zeros(Dat, ncx  ,ncy  )
    ε̇yy      = zeros(Dat, ncx  ,ncy  )
    ε̇xy      = zeros(Dat, ncx+1,ncy+1)    
    
    ε̇xx1     = zeros(Dat, ncx  ,ncy  )
    ε̇yy1     = zeros(Dat, ncx  ,ncy  )
    ε̇xy1     = zeros(Dat, ncx+1,ncy+1)    
    ε̇xy1_c   = zeros(Dat, ncx  ,ncy  )
    
    τxx      = zeros(Dat, ncx  ,ncy  )
    τyy      = zeros(Dat, ncx  ,ncy  )
    τxy      = zeros(Dat, ncx+1,ncy+1)
    τxy_c    = zeros(Dat, ncx  ,ncy  )
    τxy2_c   = zeros(Dat, ncx  ,ncy  )
    
    τxx0     = zeros(Dat, ncx  ,ncy  )
    τyy0     = zeros(Dat, ncx  ,ncy  )
    τxy0     = zeros(Dat, ncx+1,ncy+1)
    τxy0_c   = zeros(Dat, ncx  ,ncy  )
    
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
    Phasec   = ones(Int, ncx  ,ncy  )
    Phasev   = ones(Int, ncx+1,ncy+1)

    # For native plastic version 
    F        = zeros(Dat, ncx  ,ncy  )
    Fchk     = zeros(Dat, ncx  ,ncy  )
    τxy2_c   = zeros(Dat, ncx  ,ncy  )
    τxy_c    = zeros(Dat, ncx  ,ncy  )
    ε̇xy1_c   = zeros(Dat, ncx  ,ncy  )
    ε̇xy2_c   = zeros(Dat, ncx  ,ncy  )
    ε̇xy_c    = zeros(Dat, ncx  ,ncy  )
    τii      = zeros(size(ε̇xx))
    ε̇ii      = zeros(size(ε̇xx))
    λ        = zeros(size(ε̇xx))
    Pla      = zeros(size(ε̇xx))
    dQdTxx   = zeros(size(ε̇xx))
    dQdTyy   = zeros(size(ε̇xx))
    dQdTxy   = zeros(size(ε̇xx))
    ηvep_v   = ones(Dat, ncx+1, ncy+1)
    ηvep_c   = ones(Dat, ncx  , ncy  )
    
    # For geoparams
    τii0     = zeros(size(ε̇xx))
    ε̇ii_pl   = zeros(size(ε̇xx))
    
    τii0c    = zeros(size(ε̇xx))
    τii0v    = zeros(size(ε̇xy))
    ε̇iic     = zeros(size(ε̇xx))
    ε̇iiv     = zeros(size(ε̇xy))
    
    # For non-geoparams version
    η, G     = 1.0, 1.0
    ηe_c     = Δt*G.*ones(Dat, ncx  ,ncy  )
    ηe_v     = Δt*G.*ones(Dat, ncx+1,ncy+1)
    ηve_c    = zeros(Dat, ncx  ,ncy  )
    ηve_v    = zeros(Dat, ncx+1,ncy+1)
    η_vep    = ones(Dat, ncx  ,ncy  )
    η_vep_v  = ones(Dat, ncx+1,ncy+1)

    # Initialisation
    xce, yce  = LinRange(-Δx/2, Lx+Δx/2, ncx+2), LinRange(-Δy/2, Ly+Δy/2, ncy+2)
    xc, yc   = LinRange(Δx/2, Lx-Δx/2, ncx), LinRange(Δy/2, Ly-Δy/2, ncy)
    xv, yv   = LinRange(0.0, Lx, ncx+1), LinRange(0.0, Ly, ncy+1)
    radc     = (xc.-Lx./2).^2 .+ (yc'.-Ly./2).^2
    radv     = (xv.-Lx./2).^2 .+ (yv'.-Ly./2).^2
    # For non-geoparams version
    Phasec[radc.<radi] .= 2
    Phasev[radv.<radi] .= 2
    ηe_c[radc.<radi]   .= Δt*(G0/(6.0-4.0*do_DP))
    ηe_v[radv.<radi]   .= Δt*(G0/(6.0-4.0*do_DP))
    ηve_c              .= (1.0./ηe_c .+ 1.0./η).^(-1)
    ηve_v              .= (1.0./ηe_v .+ 1.0./η).^(-1)
    ηc                 .= ηve_c
    ηv                 .= ηve_v

    ηc1 = copy(ηc)
    ηv1 = copy(ηv)


    # for checking 
    ηc_check = zeros(size(ηc))
    ηv_check = zeros(size(ηv))
    τii_check      = zeros(size(τii))
    τxx_check      = zeros(size(τxx))
    τyy_check      = zeros(size(τyy))
    τxy_check      = zeros(size(τxy))
    
    # Velocity
    (Xvx,Yvx) = ([x for x=xv,y=yce], [y for x=xv,y=yce])
    (Xvy,Yvy) = ([x for x=xce,y=yv], [y for x=xce,y=yv])
    Vx     .=   εbg.*Xvx
    Vy     .= .-εbg.*Yvy
    # Time loop
    t=0.0; evo_t=[]; evo_τxx=[];
    global itg = 1
    global it_total = 0
    
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
            if UseGeoParams
                UpdateStressGeoParams!( ηc, ηv, τii, τxx, τyy, τxy, τxy_c, ε̇xx, ε̇yy, ε̇xy, ε̇xy_c, ε̇iic, ε̇iiv, τxx0, τyy0, τxy0, τii0c, τii0v, Pt, MatParam, Δt, Phasec, Phasev )

                if 1==0
                    # for debugging: do the exact same with the "native" routine and check that the values are identical

                    # Do the same calculation with native routine
                    UpdateStressNative_VEP_invariants!(ηe_c, ηe_v, ηve_c, ηve_v, τxx_check, τyy_check, τxy_check, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0,
                        ε̇ii, τii_check, τii0, ε̇xy_c, τxy_c, ε̇xy2_c, τxy2_c, F, τ_y, Pt, sinϕ, η_reg, Pla, λ, ε̇ii_pl, Fchk, ηc_check, ηv_check)

                    #UpdateStressNative_VE!(ηe_c, ηe_v, ηve_c, ηve_v, τxx_check, τyy_check, τxy_check, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0)
                
                    error_τii = norm(τii .- τii_check)
                    error_τxx = norm(τxx .- τxx_check)
                    error_τxy = norm(τxy .- τxy_check)
                    error_ηc  = norm(ηc .- ηc_check)
                    error_ηv  = norm(ηv .- ηv_check)
                    
                    if norm(error_τii)>1e-10
                        @show error_τii, error_τxx, error_τxy, error_ηc, error_ηv τii[100]-τii_check[100], F[100]
                        error("stop - difference too large")
                    end

                    if mod(iter, nout)==0
                        @show error_τii, error_τxx, error_τxy, error_ηc, error_ηv τii[100]-τii_check[100]
                    end
                end
            else
                #UpdateStressNative_VE!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0)
                #UpdateStressNative_VEP!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0,
                #                        ε̇xx1, ε̇yy1, ε̇xy1, ε̇xy1_c, ε̇xy2_c, τxy0_c, ε̇ii, τxy_c, τxy2_c, τii, F, τ_y, Pt, sinϕ, η_reg, Pla, λ, dQdTxx, dQdTyy, dQdTxy, Fchk, ηvep_c, ηvep_v,  ηc, ηv)
                
                UpdateStressNative_VEP_invariants!(ηe_c, ηe_v, ηve_c, ηve_v, τxx, τyy, τxy, ε̇xx, ε̇yy, ε̇xy, τxx0, τyy0, τxy0,
                                            ε̇ii, τii, τii0, ε̇xy_c, τxy_c, ε̇xy2_c, τxy2_c, F, τ_y, Pt, sinϕ, η_reg, Pla, λ, ε̇ii_pl, Fchk, ηc, ηv)

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
            if mod(iter, nout)==0 || iter==1
                norm_Rx = norm(Rx)/sqrt(length(Rx)); norm_Ry = norm(Ry)/sqrt(length(Ry)); norm_∇V = norm(∇V)/sqrt(length(∇V))
                err = maximum([norm_Rx, norm_Ry, norm_∇V])
                push!(err_evo1, err); push!(err_evo2, itg)
                @printf("it = %03d, iter = %04d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, ∇V=%1.3e] Fchk = %1.3e \n", it, itg, err, norm_Rx, norm_Ry, norm_∇V, maximum(Fchk))
            end
            iter+=1; global itg = iter; 
        end
       # error("stop here")

        global it_total += itg
        t = t + Δt
        push!(evo_t, t); push!(evo_τxx, maximum(τxx))
        if doPlot
            # Plotting
            #p1 = heatmap(xv, yc, Vx[:,2:end-1]', aspect_ratio=1, xlims=(0, Lx), ylims=(Δy/2, Ly-Δy/2), c=:inferno, title="Vx")
            p1 = heatmap(xv, yv, τxy', aspect_ratio=1, xlims=(0, Lx), ylims=(Δy/2, Ly-Δy/2), c=:inferno, title="τxy")
            
            #p2 = heatmap(xc, yv, Vy[2:end-1,:]', aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="Vy")
        #    p2 = heatmap(xv, yv, η_vep_v', aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="Vy")
            #p2 = heatmap(xv, yv, ηv' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="η_vep")
            p2 = heatmap(xv, yv, ηv' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="η_vep")
        
            p3 = heatmap(xc, yc, Pt' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="P")
            #p3 = heatmap(xv, yv, τxy' , aspect_ratio=1, xlims=(Δx/2, Lx-Δx/2), ylims=(0, Ly), c=:inferno, title="τxy")
            p4 = plot(evo_t, evo_τxx , xlabel="time", ylabel="max(τxx)", legend=false, linewidth=0, markershape=:circle, markersize=3)
            p4 = plot!(evo_t, 2.0.*εbg.*η0.*(1.0.-exp.(.-evo_t.*G./η0)), linewidth=2.0) # analytical solution
            display(plot(p1, p2, p3, p4))
        end
    end
    @show it_total
    return
end

for i=1:1
    println("step $i")
    doPlots = true
    @time Stokes2D_VE_inclusion(false, doPlots)
    # @time Stokes2D_VE_inclusion(true, doPlots)
end