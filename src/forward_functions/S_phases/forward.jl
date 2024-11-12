# -- S-wave isotropic travel time with absolute Vs
function tt_S_Vs(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    Vs = get_field("Vs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_S_Vs(raytmp,nray,vnox,rays2nodes,Vs)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_Vs(raytmp,nray,vnox,rays2nodes,Vs)
    @inbounds for i in eachindex(Vs)
        raytmp.u[i] = 1.0/(Vs[i])
    end
end

# -- S-wave isotropic travel time with relative perturbation dlnVs to reference velocity field
function tt_S_dlnVs(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_S_dlnVs(raytmp,nray,vnox,rays2nodes,dlnVs)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_dlnVs(raytmp,nray,vnox,rays2nodes,dlnVs)
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVs)
        raytmp.u[i] = 1.0/((1.0 + dlnVs[i])*v_1D_S[i])
    end
end

# -- S-wave isotropic travel time with absolute Vp and ratio Vp/Vs
function tt_S_Vp_Vp2Vs(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    Vp = get_field("Vp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp2Vs = get_field("Vp2Vs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_S_Vp_Vp2Vs(raytmp,nray,vnox,rays2nodes,Vp,Vp2Vs)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_Vp_Vp2Vs(raytmp,nray,vnox,rays2nodes,Vp,Vp2Vs)
    @inbounds for i in eachindex(Vp)
        raytmp.u[i] = Vp2Vs[i]/Vp[i]
    end
end

# -- S-wave isotropic travel time with relative perturbation dlnVp and ratio Vp/Vs
function tt_S_dlnVp_Vp2Vs(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    Vp2Vs = get_field("Vp2Vs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_S_dlnVp_Vp2Vs(raytmp,nray,vnox,rays2nodes,dlnVp,Vp2Vs)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt)  - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_dlnVp_Vp2Vs(raytmp,nray,vnox,rays2nodes,dlnVp,Vp2Vs)
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = Vp2Vs[i]*1.0/((1.0 + dlnVp[i])*v_1D_P[i])
    end
end

# -- S-wave isotropic travel time with relative perturbation dlnVp and ratio dlnVs/dlnVp
function tt_S_dlnVp_r(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs) 
    dlnVp = get_field("dlnVp",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    r = get_field("r",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    c_S_dlnVp_r(raytmp,nray,vnox,rays2nodes,dlnVp,r)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt)  - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_dlnVp_r(raytmp,nray,vnox,rays2nodes,dlnVp,r)
    v_1D_P = @view vnox[6,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVp)
        raytmp.u[i] = 1.0/((1.0 + dlnVp[i]*r[i])*v_1D_S[i])
    end
end

# -- S-wave isotropic travel time with relative perturbation dlnVs and hexagonal anisotropy (spherical)
function tt_S_dlnVs_fsh_psi_gamma(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fsh = get_field("fsh",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_S_dlnVs_fsh_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_dlnVs_fsh_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma; k_fsv = 0.0)
    # Note! k_fsv = fsv/fsh. Set k_fsv = 0 for elliptical anisotropy. A reasonable choice for
    # mantle hexagonal anisotropy is k_fsv = -1.0/4.75
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    ζ = @view vnox[10,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVs)
        cos2α = 2*(cos(θ[i])*cos(gamma[i])*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma[i]))^2-1.0
        β = atan(-sin(ϕ[i]-psi[i])*cos(gamma[i])/(cos(ϕ[i]-psi[i])*sin(θ[i])*cos(gamma[i])-cos(θ[i])*sin(gamma[i]))) - ζ[i]
        fsv = k_fsv*fsh
        raytmp.u[i] = ((sin(β)^2)/(1.0+fsh[i]*cos2α) + (cos(β)^2)*(1.0+(fsv))/((1.0+fsh[i])*(1.0+(fsv)*(2(cos2α)^2-1.0))))/(v_1D_S[i]*(1.0+dlnVs[i]))
    end
end

# -- S-wave isotropic travel time with dlnVs and hexagonal anisotropy (spherical -> Thomsen re-parametrization)
function tt_S_dlnVs_fsh_psi_gamma_thomsen(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fsh = get_field("fsh",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_S_dlnVs_fsh_psi_gamma_thomsen(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt) - obs.ref_t[obs.ray2obs[nray]]
end
function c_S_dlnVs_fsh_psi_gamma_thomsen(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma;
    k_γ = -1.0, k_η = 0.0, k_ϵ = 0.0, vp2vs = 1.8201, tf_elv_is_prj = false)
    # Note! tf_elv_prj should be set to true when sampling the projection of the fast axis instead of the elevation angle
    # Note! vp2vs refers to the Vp/Vs ratio of the invariant isotropic velocitites
    # Note! k_v referes to the ratio v/fsh where v is one of the Thomsen anisotropic parameters γ, η = ϵ - δ, or ϵ
    # For elliptical anisotropy, η = 0.0 and the values of k_ϵ and vp2vs are inconsequential to the shear velocity calculation
    # Assuming that fsh >= 0, reasonable values for the mantle (taken from Russell et al., 2019) are:
    #   k_γ = -1.0; k_η = 0.5018; k_ϵ = -1.6185
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    ζ = @view vnox[10,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    q2_3, q2_15, q16_15, q4_15 = (2.0/3.0), (2.0/15.0), (16.0/15.0), (4.0/15.0) # Fractions for computing invariant isotropic velocities
    @inbounds for i in eachindex(dlnVs)
        selv = tf_elv_is_prj ? asin(gamma[i]) : gamma[i]
        cosx, Δζ = symmetry_axis_cosine(psi[i], selv, ϕ[i], θ[i], ζ[i])
        cosx2 = cosx^2
        sinx2 = 1.0 - cosx2

        # Retrieve Thomsen's anisotropic parameters assuming linear scaling with fsh
        γ, η = k_γ*fsh[i], k_η*fsh[i] # η = ϵ - δ
        ϵ = k_ϵ*γ
        δ = ϵ - η
        # Compute invariant isotropic velocities; interpret dlnVs as perturbation to invariant isotropic velocity
        βiso = v_1D_S[i]*(1.0 + dlnVs[i])
        αiso = βiso*vp2vs # Assume linear scaling with βiso
        # Thomsen velocities
        α = αiso/sqrt(1.0 + q16_15*ϵ + q4_15*(ϵ - η))
        β = sqrt( ((βiso^2) - q2_15*η*(α^2))/(1.0 + q2_3*γ) )

        # Quasi-shear slownesses
        us1 = 1.0/(β*sqrt(1.0 + 2.0*qa2_b2*η*cosx2*sinx2))
        us2 = 1.0/(β*sqrt(1.0 + 2.0*γ*sinx2))

        # Effective anisotropic shear slowness
        raytmp.u[i] = us2 + (us1 - us2)*(cos(Δζ)^2)
    end
end
function symmetry_axis_cosine(symmetry_azimuth, symmetry_elevation, propagation_azimuth, propagation_elevation, qt_polarization)
    sinΔλ, cosΔλ = sincos(propagation_azimuth - symmetry_azimuth) 
    sinϕp, cosϕp = sincos(propagation_elevation)
    sinϕs, cosϕs = sincos(symmetry_elevation)
    # Cosine of angle between propagation direction and symmetry axis
    cosθ = cosΔλ*cosϕp*cosϕs + sinϕp*sinϕs
    # Angle between polarization vector and projection of symmetry axis in QT-plane (i.e. ray-normal plane)
    # Do not return cos(ζ). The sign of this angle is important for splitting intensity.
    ζ = atan(-sinΔλ*cosϕs, cosΔλ*sinϕp*cosϕs - cosϕp*sinϕs) - qt_polarization
    return cosθ, ζ
end

# -- S-wave splitting intensity with dlnVs and hexagonal anisotropy (spherical parametrization)
function si_S_dlnVs_fsh_psi_gamma(raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred,obs)
    dlnVs = get_field("dlnVs",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    fsh = get_field("fsh",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    psi = get_field("psi",raytmp.fields,rays2nodes,rays_outdom,nray,model)
    gamma = get_field("gamma",raytmp.fields,rays2nodes,rays_outdom,nray,model) 
    c_si_dlnVs_fsh_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma)
    average_velocity(raytmp,rays2nodes,nray)
    for (i,j) in enumerate(rays2nodes[1,nray]:rays2nodes[2,nray]-1)
        raytmp.dt[i] = vnox[11,j] * raytmp.u_path[i]
    end
    pred[obs.ray2obs[nray]] = sum(raytmp.dt)
end
function c_si_dlnVs_fsh_psi_gamma(raytmp,nray,vnox,rays2nodes,dlnVs,fsh,psi,gamma; k_fsv = 0.0)
    ϕ = @view vnox[8,rays2nodes[1,nray]:rays2nodes[2,nray]]
    θ = @view vnox[9,rays2nodes[1,nray]:rays2nodes[2,nray]]
    ζ = @view vnox[10,rays2nodes[1,nray]:rays2nodes[2,nray]]
    v_1D_S = @view vnox[7,rays2nodes[1,nray]:rays2nodes[2,nray]]
    @inbounds for i in eachindex(dlnVs)
        cos2α = 2*(cos(θ[i])*cos(gamma[i])*cos(ϕ[i]-psi[i])+sin(θ[i])*sin(gamma[i]))^2-1.0
        β = atan(-sin(ϕ[i]-psi[i])*cos(gamma[i])/(cos(ϕ[i]-psi[i])*sin(θ[i])*cos(gamma[i])-cos(θ[i])*sin(gamma[i]))) - ζ[i]
        fsv = k_fsv*fsh[i]
        raytmp.u[i] = 0.5*sin(2β)*((1.0)/(1.0+fsh[i]*cos2α) - (1.0+(fsv))/((1.0+fsh[i])*(1.0+(fsv)*(2(cos2α)^2-1.0))))/(v_1D_S[i]*(1.0+dlnVs[i]))
    end
end


