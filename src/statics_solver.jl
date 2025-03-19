function static_generalized_inverse(IP,observables;up_σ=false,σss=Float64[])

    evtid2st = Vector{Dict{Int64,Int64}}()
    staid2st = Vector{Dict{Int64,Int64}}()

    number_of_measurements = 0
    number_of_statics = 0
    for obs in observables.Obs
        e_flag = (obs.solve_evt != 0)
        s_flag = (obs.solve_sta != 0)
        (e_flag || s_flag) && (number_of_measurements += length(obs.prdval))
        e_flag && (number_of_statics += length(obs.evtids))
        s_flag && (number_of_statics += length(obs.staids))
    end

    obs_progression = 0
    statics_progression = 0
    I, J, V = Int64[], Int64[], Float64[]
    
    Is = sparse(zeros(Float64,number_of_statics,number_of_statics))
    Id, Jd, Vd = Int64[], Int64[], Float64[]

    for iobs in eachindex(observables.Obs)
        obs = observables.Obs[iobs]
        if up_σ
            noise_estimate = σss[iobs]
        else
            noise_estimate = obs.noise_guess
        end
        events_progression = length(obs.evtids)
        e_flag = (obs.solve_evt != 0)
        s_flag = (obs.solve_sta != 0)
        push!(evtid2st,Dict{Int64,Int64}())
        push!(staid2st,Dict{Int64,Int64}())
        for i in eachindex(obs.obsval)
            if e_flag 
                event_static_id = obs.obs2evt[i] + statics_progression
                push!(I,i+obs_progression)
                push!(J,event_static_id)
                push!(V,1.0)
                evtid2st[end][obs.evtids[obs.obs2evt[i]]] = event_static_id
            end
            if s_flag
                station_static_id = obs.obs2sta[i] + statics_progression
                e_flag && (station_static_id += events_progression)
                push!(I,i+obs_progression)
                push!(J,station_static_id)
                push!(V,1.0)
                Is[station_static_id,station_static_id] = (1.0/IP.SC.ss_uncertainty*noise_estimate)^2
                staid2st[end][obs.staids[obs.obs2sta[i]]] = station_static_id
            end
            push!(Id,i+obs_progression)
            push!(Jd,i+obs_progression)
            push!(Vd,1/(noise_estimate)^2)
        end
        e_flag && (statics_progression += length(obs.evtids))
        s_flag && (statics_progression += length(obs.staids))
        (e_flag || s_flag) && (obs_progression += length(obs.obsval))
    end
    N = [obs_progression]

    for iobs in eachindex(observables.Obs)
        obs = observables.Obs[iobs]
        e_flag = (obs.solve_evt != 0)
        s_flag = (obs.solve_sta != 0)

        if e_flag
            for jobs in eachindex(observables.Obs)
                (iobs == jobs) && continue
                if observables.Obs[iobs].solve_evt == observables.Obs[jobs].solve_evt
                    constrain_mobs_statics(evtid2st,iobs,jobs,I,J,V,Id,Jd,Vd,N)
                end
            end
        end
        if s_flag
            for jobs in eachindex(observables.Obs)
                (iobs == jobs) && continue
                if observables.Obs[iobs].solve_sta == observables.Obs[jobs].solve_sta
                    constrain_mobs_statics(staid2st,iobs,jobs,I,J,V,Id,Jd,Vd,N)
                end
            end
        end
    end

    G = sparse(I,J,V)
    Cd = sparse(Id,Jd,Vd)

    if length(V) != 0
        Gg = qr(G'*Cd*G + Is)
        Gt = G'*Cd
        statics_status = true
    else
        Gg = qr(sparse(zeros(Float64,1,1)))
        G = sparse(zeros(Float64,1,1))
        Gt = sparse(zeros(Float64,1,1))
        statics_status = false
    end

    return Gg, Gt, statics_status
end

function constrain_mobs_statics(id2st,iobs,jobs,I,J,V,Id,Jd,Vd,N)
    for (i,j) in id2st[iobs]
        if haskey(id2st[jobs],i)
            N[1] += 1
            push!(I,N[1])
            push!(J,j)
            push!(V,1.0)
            push!(I,N[1])
            push!(J,id2st[jobs][i])
            push!(V,-1.0)
            push!(Id,N[1])
            push!(Jd,N[1])
            push!(Vd,1/(0.001)^2)   # -- 2 ooms lower than residuals' uncertainty...
        end
    end
end

function solve_statics(model,observables,Gg,Gt,residuals,estimated_statics)
    residuals .= 0.0
    estimated_statics .= 0.0
    global_idx = 0
    for iobs in eachindex(observables.Obs)
        obs = observables.Obs[iobs]
        e_flag = (obs.solve_evt != 0)
        s_flag = (obs.solve_sta != 0)
        if e_flag || s_flag
            for i in eachindex(obs.obsval)
                global_idx += 1
                residuals[global_idx] = obs.obsval[i]-obs.prdval[i]
                # σd = model.dataspace.Obs[iobs].noise[1]
                # Cd[global_idx,global_idx] = 1/σd^2
            end
        end
    end

    estimated_statics .= (Gg)\(Gt*residuals)

    statics_progression = 0
    for obsid in eachindex(observables.Obs)
        obs = observables.Obs[obsid]
        events_progression = length(obs.evtids)
        e_flag = (obs.solve_evt != 0)
        s_flag = (obs.solve_sta != 0)
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

function reinitialize_statics(MarkovChains,observables,IP)
    # -- reinitializes noise for statics evaluation
    σss = zeros(Float64,length(observables.Obs))

    nchains = length(MarkovChains)
    chain_models = Int64(round(IP.RayTracingInit.sub_its / IP.MCS.saveint))
    nsamples = 0
    for chain in eachindex(MarkovChains)
        MarkovChain = MarkovChains[chain]
        for model in MarkovChain[end-chain_models+1:end]
        nsamples += 1
            for iobs in eachindex(observables.Obs)
                obs = observables.Obs[iobs]
                σss[iobs] += model.dataspace.Obs[iobs].noise[1]
            end
        end
    end
    @. σss = σss / nsamples
    Gg, Gt, statics_status = static_generalized_inverse(IP,observables;up_σ=true,σss=σss)
    print("\nreinitializing noise for statics evauation... ",nsamples," ",σss)
    return Gg, Gt
end