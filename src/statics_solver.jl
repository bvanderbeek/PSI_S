function static_generalized_inverse(IP,observables)
    number_of_measurements = 0
    number_of_statics = 0
    for obs in observables.Obs
        e_flag = obs.solve_evt
        s_flag = obs.solve_sta
        (e_flag || s_flag) && (number_of_measurements += length(obs.prdval))
        e_flag && (number_of_statics += length(obs.evtids))
        s_flag && (number_of_statics += length(obs.staids))
    end

    obs_progression = 0
    statics_progression = 0
    I, J, V = Int64[], Int64[], Float64[]
    Is, Js, Vs = Int64[], Int64[], Float64[]

    for obs in observables.Obs
        events_progression = length(obs.evtids)
        e_flag = obs.solve_evt
        s_flag = obs.solve_sta
        for i in eachindex(obs.obsval)
            if e_flag 
                event_static_id = obs.obs2evt[i] + statics_progression
                push!(I,i+obs_progression)
                push!(J,event_static_id)
                push!(V,1.0)
            end
            if s_flag
                station_static_id = obs.obs2sta[i] + statics_progression
                e_flag && (station_static_id += events_progression)
                push!(I,i+obs_progression)
                push!(J,station_static_id)
                push!(V,1.0)
                push!(Is,i+obs_progression)
                push!(Js,station_static_id)
                push!(Vs,1.0)
            end
        end
        e_flag && (statics_progression += length(obs.evtids))
        s_flag && (statics_progression += length(obs.staids))
        (e_flag || s_flag) && (obs_progression += length(obs.obsval))
    end
    G = sparse(I,J,V)
    if length(Vs) != 0
        station_statics_diag = sparse(Is,Js,Vs)
        station_damp = station_statics_diag'*station_statics_diag
    else
        station_damp = sparse(zeros(Float64,number_of_statics,number_of_statics))
    end
    if length(V) != 0
        damping_factor = (1.0/IP.SC.ss_uncertainty)^2
        Gg = qr(G'G + damping_factor*station_damp)
    else
        Gg = qr(sparse(zeros(Float64,1,1)))
        G = sparse(zeros(Float64,1,1))
    end

    return Gg, G'
end

function solve_statics(model,observables,Gg,Gt,residuals,estimated_statics)
    residuals .= 0.0
    estimated_statics .= 0.0
    global_idx = 0
    for obs in observables.Obs
        e_flag = obs.solve_evt
        s_flag = obs.solve_sta
        if e_flag || s_flag
            for i in eachindex(obs.obsval)
                global_idx += 1
                residuals[global_idx] = obs.obsval[i]-obs.prdval[i]
            end
        end
    end
    
    estimated_statics .= (Gg)\(Gt*residuals)

    statics_progression = 0
    for obsid in eachindex(observables.Obs)
        obs = observables.Obs[obsid]
        events_progression = length(obs.evtids)
        e_flag = obs.solve_evt
        s_flag = obs.solve_sta
        for i in eachindex(obs.obsval)
            if e_flag 
                event_static_id = obs.obs2evt[i] + statics_progression
                model.dataspace.Obs[obsid].estatics[obs.obs2evt[i]] = estimated_statics[event_static_id]
            end
            if s_flag
                station_static_id = obs.obs2sta[i] + statics_progression
                e_flag && (station_static_id += events_progression)
                model.dataspace.Obs[obsid].sstatics[obs.obs2sta[i]] = estimated_statics[station_static_id]
            end
        end
        e_flag && (statics_progression += length(obs.evtids))
        s_flag && (statics_progression += length(obs.staids))
    end
end
