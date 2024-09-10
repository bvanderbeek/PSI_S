


"""
    vqp_thomsen(pazm, pelv, sazm, selv, α_iso, ϵ, δ)
    Approximate (but not weak) Thomsen qP phase-velocity expression

    pazm: Propagation azimuth (radians)
    pelv: Propagation elevation (radians)
    sazm: Symmetry axis azimuth (radians)
    selv: Symmetry axis elevation (radians)
    vp1D: reference 1D P-velocity
    dlnVp: Invariant isotropic P-velocity perturbation!
    ϵ: Thomsen parameter, 0.5*(c11 - c33)/c33
    δ: Thomsen parameter, (c13 - c33 + 2*c44)/c33
"""
function vqp_thomsen(pazm, pelv, sazm, selv, vp1D, dlnVp, ϵ, δ)
    sin_selv, cos_selv = sincos(selv)
    sin_pelv, cos_pelv = sincos(pelv)
    cosx2 = (cos(pazm - sazm)*cos_pelv*cos_selv + sin_pelv*sin_selv)^2
    sinx2 = 1.0 - cosx2
    sinx4 = sinx2^2

    # Thomsen α-parameter from invariant isotropic velocity
    α_iso = vp1D*(1.0 + dlnVp)
    α = α_iso/sqrt(1.0 + (16/15)*ϵ + (4/15)*δ)

    return α*sqrt(1.0 + 2.0*δ*sinx2*cosx2 + 2.0*ϵ*sinx4)
end

"""
    parametric_thomsen_delta(ϵ; p1 = 0.6742, p2 = -0.8169, q1 = 0.04419)
    Return Thomsen δ parameter as a function of ϵ. Fit parameters taken from
    Manuele's DEM model for penny-shaped fractures with aspect ratio 100.
"""
function parametric_thomsen_delta(ϵ; p1 = 0.6742, p2 = -0.8169, q1 = 0.04419)
    return ϵ > 0.0 ? ϵ*(p1*ϵ + p2)/(ϵ + q1) : ϵ # Defaults to elliptical
end

"""
    parameteric_thomsen_melt(f, dlnVp_0)
    Return Thomsen parameters (ϵ, δ) and *total* velocity perturbation
    (dlnVp) as a function of melt fraction (f)

    f: Melt fraction
    dlnVp_0: Initial isotropic velocity perturbation
"""
function parameteric_thomsen_melt(f, dlnVp_0)
    # dlnVp
    p_dvp = (-0.4363, -0.1672, -2.759e-5)
    q_dvp = (2.693e-2)
    dlnVp = dlnVp_0 + ((p_dvp[1]*f^2 + p_dvp[2]*f + p_dvp[3])./(f + q_dvp[1]))
    # ϵ
    p_e = (-2.1911, 2.1881, -1.4518e-3)
    ϵ = p_e[1]*(f^2) + p_e[2]*f + p_e[3]
    # δ
    p1 = -0.3302
    q1 = 0.01573
    δ = ϵ*p1/(f + q1)

    return ϵ, δ, dlnVp
end

