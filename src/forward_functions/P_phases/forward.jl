# To simplify these forward functions, consider defining a set of base velocity functions for different parameterisations
# 1. Isotropic
# 2. Harmonic
# 3. Thomsen
# 4. Invariant Thomsen
# Then, within each forward function, you can call the above methods with the appropriate inputs
# This should reduce some redundancy and probability of errors

# NEW ANISOTROPIC FUNCTIONS
function tt_P_dlnVp_AniMag_SymAzm_SymRad_Thomsen(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    # Extract Views to Perturbational Fields
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    AniMag = get_field("AniMag",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    SymAzm = get_field("SymAzm",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    SymRad = get_field("SymRad",raytmp.fields,rays2nodes,rays_outdom,nray,model)

    # Populate raytmp with appropriate slownesses
    p_slowness_thomsen_dlnVp_AniMag_SymAzm_SymElv!(raytmp,nray,vnox,rays2nodes,dlnVp,AniMag,SymAzm,SymRad; tf_sin_sym_elv = true)

    # Compute ray segment travel-times
    average_velocity(raytmp,rays2nodes,nray) # Ray segment averaged velocity
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]

    return nothing
end
function tt_P_dlnVp_AniMag_SymAzm_SymElv_Thomsen(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    # Extract Views to Perturbational Fields
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    AniMag = get_field("AniMag",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    SymAzm = get_field("SymAzm",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    SymElv = get_field("SymElv",raytmp.fields,rays2nodes,rays_outdom,nray,model)

    # Populate raytmp with appropriate slownesses
    p_slowness_thomsen_dlnVp_AniMag_SymAzm_SymElv!(raytmp,nray,vnox,rays2nodes,dlnVp,AniMag,SymAzm,SymElv; tf_sin_sym_elv = false)

    # Compute ray segment travel-times
    average_velocity(raytmp,rays2nodes,nray) # Ray segment averaged velocity
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]

    return nothing
end
function p_slowness_thomsen_dlnVp_AniMag_SymAzm_SymElv!(raytmp,nray,vnox,rays2nodes,dlnVp,AniMag,SymAzm,SymElv; tf_sin_sym_elv = false)
    # Hard-coded Thomsen Ratios
    k_ϵ, k_δ = 1.0, 1.0 # Elliptical (±, slow/fast symmetry)
    # k_ϵ, k_δ = -1.0, -1.3101 # Oceanic Lithosphere from Russell et al. 2019

    # Compute phase velocities
    vp_0 = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]] # Reference model Vp
    RayAzm = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]] # Ray azimuth
    RayElv = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]] # Ray elevation
    q16_15, q4_15 = (16.0/15.0), (4.0/15.0) # Fractions for computing invariant isotropic velocity
    @inbounds for i in eachindex(dlnVp)
        # Directional components
        sym_elv_i = tf_sin_sym_elv ? asin(SymElv[i]) : SymElv[i]
        cosx = symmetry_axis_cosine(SymAzm[i], sym_elv_i, RayAzm[i], RayElv[i])
        cosx2 = cosx^2
        sinx2 = 1.0 - cosx2
        sinx4 = sinx2^2
        # Interpret dlnVp as a perturbation to the invariant isotropic velocity
        vp_i = vp_0[i]*(1.0 + dlnVp[i]) 
        # Thomsen Parameters
        ϵ, δ = k_ϵ*AniMag[i], k_δ*AniMag[i] # Define Thomsen parameters as scalar multiples of anisotropic magnitude
        α = vp_i/sqrt(1.0 + q16_15*ϵ + q4_15*δ) # The Thomsen P-velocity (i.e. symmetry axis velocity)
        # P-slowness: Hexagonal Weak Elastic Anisotropy
        raytmp.u[i] = 1.0/(α*sqrt(1.0 + 2.0*δ*sinx2*cosx2 + 2.0*ϵ*sinx4))
    end
    return nothing
end
function symmetry_axis_cosine(symmetry_azimuth, symmetry_elevation, propagation_azimuth, propagation_elevation)
    cosΔλ = cos(propagation_azimuth - symmetry_azimuth)
    sinϕp, cosϕp = sincos(propagation_elevation)
    sinϕs, cosϕs = sincos(symmetry_elevation)
    # Cosine of angle between propagation direction and symmetry axis
    cosθ = cosΔλ*cosϕp*cosϕs + sinϕp*sinϕs
    return cosθ
end

# Start Gianmarco Forward Functions in v1.5
# -- P-wave anisotropic travel time with relative perturbation dlnVp and radial anisotropy (spherical parametrization)
function tt_P_dlnVp_fp(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_dlnVp_fp(raytmp,nray,vnox,rays2nodes,dlnVp,fp)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_fp(raytmp,nray,vnox,rays2nodes,dlnVp,fp)
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(sin(θ[i]))^2-1.0))*(1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# -- P-wave anisotropic travel time with relative perturbation dlnVp and hexagonal anisotropy (MAV parametrization)
function tt_P_dlnVp_fp_psi_v3(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    v3 = get_field("v3",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_P_dlnVp_fp_psi_v3(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,v3)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_fp_psi_v3(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,v3)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        gamma = asin(v3[i])
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(θ[i])*cos(gamma)*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma))^2-1.0))*(1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# Note! The tt_aniso_P_th amd tt_aniso_P_th_ell functions are deprecated by tt_P_dlnVp_fp_psi_gamma_thomsen
# Did not include tt_aniso_P_th_crack because parameters are too specific to a particular DEM model
# End Gianmarco Forward Functions in v1.5


# -- P-wave isotropic travel time with absolute Vp
function tt_P_Vp(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    Vp = get_field("Vp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_Vp(raytmp,nray,vnox,rays2nodes,Vp)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_Vp(raytmp,nray,vnox,rays2nodes,Vp)
    @inbounds for i in eachindex(Vp)
        raytmp.u[i] = 1.0/(Vp[i])
    end
end

# -- P-wave isotropic travel time with relative perturbation dlnVp to reference velocity field
function tt_P_dlnVp(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_dlnVp(raytmp,nray,vnox,rays2nodes,dlnVp)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp(raytmp,nray,vnox,rays2nodes,dlnVp)
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# -- P-wave isotropic travel time with relative perturbation dlnVp to free 1D Vp field
function tt_P_dlnVp_Vp(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp = get_field("Vp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_Vp_dlnVp(raytmp,nray,vnox,rays2nodes,Vp,dlnVp)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_Vp_dlnVp(raytmp,nray,vnox,rays2nodes,Vp,dlnVp)
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + dlnVp[i])*Vp[i])
    end
end

# -- P-wave isotropic travel time with relative perturbation dlnVs to reference velocity field and ratio Vp/Vs
function tt_P_dlnVs_Vp2Vs(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp2Vs = get_field("Vp2Vs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_dlnVs_Vp2Vs(raytmp,nray,vnox,rays2nodes,dlnVs,Vp2Vs)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVs_Vp2Vs(raytmp,nray,vnox,rays2nodes,dlnVs,Vp2Vs)
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVs)
        raytmp.u[i] = 1.0/((1.0 + dlnVs[i])*v_1D_S[i]*Vp2Vs[i])
    end
end

# -- P-wave anisotropic travel time with absolute Vp and azimuthal anisotropy (spherical parametrization)
function tt_P_Vp_fp_psi(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    Vp = get_field("Vp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_Vp_fp_psi(raytmp,nray,vnox,rays2nodes,Vp,fp,psi)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_Vp_fp_psi(raytmp,nray,vnox,rays2nodes,Vp,fp,psi)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(Vp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(ϕ[i]-psi[i]))^2-1.0))*(Vp[i]))
    end
end

# -- P-wave anisotropic travel time with relative perturbation dlnVp and azimuthal anisotropy (spherical parametrization)
function tt_P_dlnVp_fp_psi(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_P_dlnVp_fp_psi(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_fp_psi(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(ϕ[i]-psi[i]))^2-1.0))*(1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# -- P-wave anisotropic travel time with relative perturbation dlnVp and hexagonal anisotropy (spherical parametrization)
function tt_P_dlnVp_fp_psi_gamma(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_P_dlnVp_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,gamma)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(θ[i])*cos(gamma[i])*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma[i]))^2-1.0))*(1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# -- P-wave anisotropic travel time with dlnVp and hexagonal anisotropy (spherical -> Thomsen re-parametrization)
function tt_P_dlnVp_fp_psi_gamma_thomsen(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_P_dlnVp_fp_psi_gamma_thomsen(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_fp_psi_gamma_thomsen(raytmp,nray,vnox,rays2nodes,dlnVp,fp,psi,gamma;
    k_ϵ = -1.0, k_δ = -1.0, tf_elv_is_prj = true)
    # Note! tf_elv_prj should be set to true when sampling the projection of the fast axis instead of the elevation angle
    # Note! k_v referes to the ratio v/fp where v is one of the Thomsen anisotropic parameters δ or ϵ
    # For elliptical anisotropy, k_δ = k_ϵ
    # Assuming that fp >= 0, reasonable values for the mantle (taken from Russell et al., 2019) are:
    #   k_ϵ = -1.0; k_δ = -1.3101
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    q16_15, q4_15 = (16.0/15.0), (4.0/15.0) # Fractions for computing invariant isotropic velocity
    @inbounds for i in eachindex(dlnVp)
        selv = tf_elv_is_prj ? asin(gamma[i]) : gamma[i]
        sinγ, cosγ = sincos(selv)
        sinθ, cosθ = sincos(θ[i])
        cosx2 = (cos(ϕ[i] - psi[i])*cosθ*cosγ + sinθ*sinγ)^2
        sinx2 = 1.0 - cosx2
        sinx4 = sinx2^2

        ϵ, δ = k_ϵ*fp[i], k_δ*fp[i] # Define Thomsen parameters as scalar multiples of fp
        αiso = v_1D_P[i]*(1.0 + dlnVp[i]) # Invariant isotropic velocity; interpret dlnVp as perturbation to invariant isotropic velocity
        α = αiso/sqrt(1.0 + q16_15*ϵ + q4_15*δ) # The Thomsen P-velocity (i.e. symmetry axis velocity)

        raytmp.u[i] = 1.0/(α*sqrt(1.0 + 2.0*δ*sinx2*cosx2 + 2.0*ϵ*sinx4)) # Approximate Thomsen (Hexagonal) Anisotropy (more accurate)
        # raytmp.u[i] = 1.0/(α*(1.0 + δ*sinx2*cosx2 + ϵ*sinx4)) # Approximate Weak Thomsen (Hexagonal) Anisotropy
        # raytmp.u[i] = 1.0/((1.0 + fp[i]*(2.0*cosx2 - 1.0))*(1.0 + dlnVp[i])*v_1D_P[i]) # Elliptical Anisotropy
    end
end

# -- P-wave anisotropic travel time with relative perturbation dlnVp to free 1D Vp and hexagonal anisotropy (spherical parametrization)
function tt_P_dlnVp_Vp_fp_psi_gamma(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp = get_field("Vp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_P_dlnVp_Vp_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,Vp,dlnVp,fp,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVp_Vp_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,Vp,dlnVp,fp,psi,gamma)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(θ[i])*cos(gamma[i])*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma[i]))^2-1.0))*(1.0 + dlnVp[i])*Vp[i])
    end
end

# -- P-wave anisotropic travel time with relative perturbation dlnVs, Vp2Vs and hexagonal anisotropy (spherical parametrization)
function tt_P_dlnVs_Vp2Vs_fp_psi_gamma(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp2Vs = get_field("Vp2Vs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fp = get_field("fp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_P_dlnVs_Vp2Vs_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,Vp2Vs,fp,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_P_dlnVs_Vp2Vs_fp_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,Vp2Vs,fp,psi,gamma)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + fp[i]*(2*(cos(θ[i])*cos(gamma[i])*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma[i]))^2-1.0))*(1.0 + dlnVs[i])*v_1D_S[i]*Vp2Vs[i])
    end
end
