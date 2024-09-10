 
function build_obslist()
    IP_obs = Vector{obsinfo}()
    # Positional Inputs:
    # 1. Initialized observation vector
    # 2. Observation type (not used for teleseismic delay inversions)
    # 3. Relative path to data file
    # 4. Demean data?
    # 5. Solve for event statics?
    # 6. Solve for station statics?
    # 7. Initial noise guess for data
    # 8. Forward function for making predictions
    # 8.1. "tt_P_dlnVp", "tt_P_dlnVp_fp_psi_gamma", "tt_P_dlnVp_fp_psi_gamma_thomsen"
    add_obs(IP_obs,"teleseismic_P_delays",string((@__DIR__),"/input/compressional_delays_melt.dat"),true,true,false,0.06,"tt_P_dlnVp")

    return IP_obs
end

function build_fieldslist()
    IP_fields = Vector{fieldinfo}()
    # Positional Inputs:
    # 1. Initialized parameter fields vector
    # 2. Field name
    # 2.1. dlnVp
    # 3. "Prior"
    # 3.1. "uniform"
    # 4. [min., max.] Latitude limits (deg.)
    # 5. [min., max.] Longitude limits (deg.)
    # 6. [min., max.] Elevation limits (km)
    # 7. [min., max.] Time limits for 4D
    # 8. Extrapolation value when squeezing is enables
    add_field(IP_fields,"dlnVp","uniform",[-0.05,0.05],[-17.0,-17.0],[-121.4,-102.6],[-600.0,0.0],[0.0,0.0],0.0)
    #add_field(IP_fields,"fp","uniform",[0.0,0.06],[36.0,54.0],[-134.0,-112.00],[-500.0,3.136],[0.0,0.0],0.0)
    #add_field(IP_fields,"psi","uniform",[-2.0*pi,2.0*pi],[36.0,54.0],[-134.0,-112.00],[-500.0,3.136],[0.0,0.0],0.0)
    #add_field(IP_fields,"gamma","uniform",[0.0,0.5*pi],[36.0,54.0],[-134.0,-112.00],[-500.0,3.136],[0.0,0.0],0.0)

    return IP_fields 
end


