
function vereornox(rays,rnodes)
    vnox = zeros(Float64,14,length(rnodes.n2r))
    vnox[1:4,:] .= rnodes.c
    vnox[5,:] .= rnodes.r[1,:]
    for nray in eachindex(rays)
        nodes = rnodes.r2n[nray]
        for (i,node) in enumerate(nodes)
            vnox[6,node] = rays[nray].v_1D_P[i]
            vnox[7,node] = rays[nray].v_1D_S[i]
            vnox[8,node] = rays[nray].ϕ[i]
            vnox[9,node] = rays[nray].θ[i]
            vnox[10,node] = rays[nray].ζ[i]
            if node != nodes[end]
                vnox[11,node] = rays[nray].L[i]
            end
        end
    end
    vnox = SharedArray(vnox)
    nodes2rays = SharedArray(rnodes.n2r)
    rays2nodes = zeros(Int64,2,length(rnodes.r2n))
    for nray in eachindex(rays)
        rays2nodes[1,nray] = rnodes.r2n[nray][begin]
        rays2nodes[2,nray] = rnodes.r2n[nray][end]
    end
    rays2nodes = SharedArray(rays2nodes)
    rays_outdom = zeros(Int64,2*length(rays[begin].outdomain),length(rnodes.r2n))
    for nray in eachindex(rays)
        for fieldid in eachindex(rays[nray].outdomain)
            rays_outdom[(fieldid-1)*2 + 1,nray] = rays[nray].outdomain[fieldid][1]
            rays_outdom[(fieldid-1)*2 + 2,nray] = rays[nray].outdomain[fieldid][2]
        end
    end
    rays_outdom = SharedArray(rays_outdom)
    return vnox, nodes2rays, rays2nodes, rays_outdom
end


# This function is intended to deprecate run_RJMCMC
function run_psi_s(parameter_file::String, ichunk)
    # Assumes that chunks are run serially and chunk 1 saves the IPConst file
    tf_save = ichunk > 1 ? false : true
    println("Progress: Building input structures...")
    P, _, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager = build_inputs(parameter_file; tf_save = tf_save, tf_serial = false)
    run_psi_s(ichunk, P, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager)
    return nothing
end
function run_psi_s(ichunk, P::Dict, IP_fields::Vector{fieldinfo}, IP::IPConst, rnodes, rays, evtsta, observables, LocalRaysManager)

    println("Progress: Making output directory...")
    run_directory = P["out_directory"]*"/"*P["run_id"]
    mkpath(run_directory) # Make output path if it doesn't exist
    nchains = IP.MCS.nchains # -- number of gianMarkov chains

    println("Progress: Building voronoi diagrams...")
    @time vnox,nodes2rays,rays2nodes,rays_outdom = vereornox(rays,rnodes)

    # -- divides iterations in sub-samples; accounts for parallel tempering and ray-tracing
    println("Progress: Initializing sub-samples...")
    sub_samples, PT_samples = 1, 1
    iterations = IP.MCS.niterations
    if LocalRaysManager.status 
        sub_samples = floor(Int64, IP.MCS.niterations / LocalRaysManager.sub_its)
        if sub_samples == 0
            sub_samples = 1
        else        
            iterations = LocalRaysManager.sub_its
        end
    end
    if IP.PT.status
        PT_samples = floor(Int64, iterations / IP.PT.intsw_its)
        iterations = IP.PT.intsw_its
    end

    # -- generalized inverse matrices for static inversions and PT utilities
    println("Progress: Generalized inverse for statics...")
    Gg, Gt, statics_status = static_generalized_inverse(IP,observables)
    if !statics_status
        print("\nno static corrections\n")
    end
    temperatures = Float64[]
    temp_ind = Int64[]
    obs_lengths = Int64[]
    obs_noise = zeros(Float64,nchains,length(observables.Obs))

    # -- distributed arrays for parallel chains
    println("Progress: Initializing parallel chains...")
    @time ObjsInChains, MarkovChains = initialize_parallel_chains(IP,IP_fields,rnodes,evtsta,observables,rays,Gg,Gt,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
    # -- tempered utilities
    [push!(temperatures,ObjsInChains[n].model.T[1]) for n in eachindex(ObjsInChains)]
    [push!(temp_ind,i) for i in eachindex(temperatures)]
    [push!(obs_lengths,length(obs.obsval)) for obs in ObjsInChains[begin].observables.Obs]
    for i in eachindex(ObjsInChains)
        [obs_noise[i,j] = ObjsInChains[i].model.dataspace.Obs[j].noise[1] for j in eachindex(ObjsInChains[i].model.dataspace.Obs)]
    end
    temperatures = SharedArray(temperatures)
    temp_ind = SharedArray(temp_ind)
    T_misfits = SharedArray(zeros(Float64,nchains))
    obs_lengths = SharedArray(obs_lengths)
    obs_noise = SharedArray(obs_noise)

    IP.PT.status && print("\ntemperatures: ",temperatures)
    TMatrix = zeros(Float64,nchains,nchains)
    DMatrix = zeros(Float64,nchains,nchains) 

    println("Progress: Defining actions...")
    actions = [1,2,3,4]
    if (IP.MCS.hierarchical == true)
		push!(actions,5)
	end
    fields = Int64[]
    for i in eachindex(MarkovChains[1][end].fields)
        if MarkovChains[1][end].fields[i].vlims[1] != MarkovChains[1][end].fields[i].vlims[2]
            push!(fields,i)
        end
    end
    actions, fields = SharedArray(actions), SharedArray(fields)

    println("Progress: Starting rj-mcmc sampling...")
    for ii in 1:sub_samples
        for jj in 1:PT_samples
            it0 = (ii-1)*LocalRaysManager.sub_its + (jj-1)*IP.PT.intsw_its + 1
            @sync @distributed for chain in 1:nchains
                ObjsInChain, MarkovChain = localpart(ObjsInChains)[1], localpart(MarkovChains)[1]
                ObjsInChain.model.T[1] = temperatures[chain]
                (nchains == 1) && (chain = ichunk)
                RJMCMC(ObjsInChain, MarkovChain, vnox, nodes2rays, rays2nodes, rays_outdom, IP, chain, iterations; it0=it0, actions=actions, fields=fields, statics_status=statics_status)
                T_misfits[chain] = ObjsInChain.model.misfit[1]
                [obs_noise[chain,j] = ObjsInChain.model.dataspace.Obs[j].noise[1] for j in eachindex(ObjsInChain.model.dataspace.Obs)]
            end
            if IP.PT.status && ((it0 + iterations) >= IP.PT.pt_pause)
               @time parallel_tempering(T_misfits,obs_lengths,obs_noise,temperatures,temp_ind,nchains,TMatrix,DMatrix)
            end 
	    end
        if LocalRaysManager.status && (ii < sub_samples)
            rnodes = raytracing(rays,LocalRaysManager,MarkovChains,observables,evtsta,IP;it=ii)
            vnox,nodes2rays,rays2nodes,rays_outdom = vereornox(rays,rnodes)
            Gg, Gt = reinitialize_statics(MarkovChains,observables,IP)
            reinitialize_models(ObjsInChains,rnodes,rays,observables,Gg,Gt,IP,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
        end
    end

    println("Progress: Saving Markov chains...")
    for chain in 1:nchains
        MarkovChain = MarkovChains[chain]
        chain = nchains * (ichunk-1) + chain
        (nchains == 1) && (chain = ichunk)
        save(string(run_directory, "/chain.", chain, ".jld"), "MarkovChain", MarkovChain)
    end

    @. TMatrix = TMatrix / DMatrix

    if ichunk == 1
        save(string(run_directory, "/psi_s_utilities.jld"), "rnodes", rnodes, "observables", observables, "TransitionMatrix", TMatrix, "evtsta", evtsta)
    end

    return nothing
end

function run_RJMCMC(wrk_dir, name, chain_id)
    # -- loading input
    @time input = load(string(wrk_dir, "/output/", name, "/", name,".jld"),"IP","IP_fields","rnodes","rays","evtsta","observables","LocalRaysManager")
    IP, IP_fields, rnodes, rays, evtsta, observables, LocalRaysManager = input[1], input[2], input[3], input[4], input[5], input[6], input[7]
    nchains = IP.MCS.nchains # -- number of gianMarkov chains

    @time vnox,nodes2rays,rays2nodes,rays_outdom = vereornox(rays,rnodes)

    # -- divides iterations in sub-samples; accounts for parallel tempering and ray-tracing
    sub_samples, PT_samples = 1, 1
    iterations = IP.MCS.niterations
    if LocalRaysManager.status 
        sub_samples = floor(Int64, IP.MCS.niterations / LocalRaysManager.sub_its)
        if sub_samples == 0
            sub_samples = 1
        else        
            iterations = LocalRaysManager.sub_its
        end
    end
    if IP.PT.status
        PT_samples = floor(Int64, iterations / IP.PT.intsw_its)
        iterations = IP.PT.intsw_its
    end

    # -- generalized inverse matrices for static inversions and PT utilities
    Gg, Gt, statics_status = static_generalized_inverse(IP,observables)
    if !statics_status
        print("\nno static corrections\n")
    end
    temperatures = Float64[]
    temp_ind = Int64[]
    obs_lengths = Int64[]
    obs_noise = zeros(Float64,nchains,length(observables.Obs))

    # -- distributed arrays for parallel chains
    @time ObjsInChains, MarkovChains = initialize_parallel_chains(IP,IP_fields,rnodes,evtsta,observables,rays,Gg,Gt,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
    # -- tempered utilities
    [push!(temperatures,ObjsInChains[n].model.T[1]) for n in eachindex(ObjsInChains)]
    [push!(temp_ind,i) for i in eachindex(temperatures)]
    [push!(obs_lengths,length(obs.obsval)) for obs in ObjsInChains[begin].observables.Obs]
    for i in eachindex(ObjsInChains)
        [obs_noise[i,j] = ObjsInChains[i].model.dataspace.Obs[j].noise[1] for j in eachindex(ObjsInChains[i].model.dataspace.Obs)]
    end
    temperatures = SharedArray(temperatures)
    temp_ind = SharedArray(temp_ind)
    T_misfits = SharedArray(zeros(Float64,nchains))
    obs_lengths = SharedArray(obs_lengths)
    obs_noise = SharedArray(obs_noise)

    IP.PT.status && print("\ntemperatures: ",temperatures)
    TMatrix = zeros(Float64,nchains,nchains)
    DMatrix = zeros(Float64,nchains,nchains) 

    actions = [1,2,3,4]
    if (IP.MCS.hierarchical == true)
		push!(actions,5)
	end
    fields = Int64[]
    for i in eachindex(MarkovChains[1][end].fields)
        if MarkovChains[1][end].fields[i].vlims[1] != MarkovChains[1][end].fields[i].vlims[2]
            push!(fields,i)
        end
    end
    actions, fields = SharedArray(actions), SharedArray(fields)

    for ii in 1:sub_samples
        for jj in 1:PT_samples
            it0 = (ii-1)*LocalRaysManager.sub_its + (jj-1)*IP.PT.intsw_its + 1
            @sync @distributed for chain in 1:nchains
                ObjsInChain, MarkovChain = localpart(ObjsInChains)[1], localpart(MarkovChains)[1]
                ObjsInChain.model.T[1] = temperatures[chain]
                (nchains == 1) && (chain = chain_id)
                RJMCMC(ObjsInChain, MarkovChain, vnox, nodes2rays, rays2nodes, rays_outdom, IP, chain, iterations; it0=it0, actions=actions, fields=fields, statics_status=statics_status)
                T_misfits[chain] = ObjsInChain.model.misfit[1]
                [obs_noise[chain,j] = ObjsInChain.model.dataspace.Obs[j].noise[1] for j in eachindex(ObjsInChain.model.dataspace.Obs)]
            end
            if IP.PT.status && ((it0 + iterations) >= IP.PT.pt_pause)
               @time parallel_tempering(T_misfits,obs_lengths,obs_noise,temperatures,temp_ind,nchains,TMatrix,DMatrix)
            end 
	    end
        if LocalRaysManager.status && (ii < sub_samples)
            rnodes = raytracing(rays,LocalRaysManager,MarkovChains,observables,evtsta,IP;it=ii)
            vnox,nodes2rays,rays2nodes,rays_outdom = vereornox(rays,rnodes)
            Gg, Gt = reinitialize_statics(MarkovChains,observables,IP)
            reinitialize_models(ObjsInChains,rnodes,rays,observables,Gg,Gt,IP,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
        end
    end

    for chain in 1:nchains
        MarkovChain = MarkovChains[chain]
        chain = nchains * (chain_id-1) + chain
        (nchains == 1) && (chain = chain_id)
        save(string(wrk_dir, "/output/", name, "/", name, ".", chain, ".jld"), "MarkovChain", MarkovChain)
    end

    @. TMatrix = TMatrix / DMatrix

    save(string(wrk_dir, "/output/", name, "/", name, "_utilities.jld"), "rays",rays,"rnodes", rnodes,"observables",observables,"TransitionMatrix",TMatrix,"evtsta",evtsta)

end

function RJMCMC(ObjsInChain, MarkovChain, vnox, nodes2rays, rays2nodes, rays_outdom, IP, chain, iterations; it0=1, actions=[1], fields=[1], statics_status=false)

    (; model, observables, update) = ObjsInChain
    (; Gg, Gt) = update
    σ = IP.MCS.pert_size
    if IP.PT.temp_pert
        σ = IP.MCS.pert_size * sqrt(model.T[1])
    end

    print("\nbegin: ",model.rms," - ",model.misfit,"\n")
    for i = it0:(it0+iterations-1)
        print_progress(IP,model,i,chain)

        field = rand(fields)
        if i <= IP.MCS.pertv_initrange
            action = 1
        else
            action = rand(actions)
        end
        if statics_status && ((i == IP.SC.statics_pause) || (((i % (IP.SC.statics_its)) .== 0) && (i >= IP.SC.statics_pause)))
            action = 6
        end

        update.field .= field
        update.action .= action
        clear_updtmp(update, model, observables, field)

        # -- perturb model --
        if action == 1
            pertv(model.fields[field],σ,vnox,nodes2rays,rays2nodes,update,IP)
        elseif action == 2
            pertp(model.fields[field],σ,vnox,nodes2rays,rays2nodes,update,IP)
        elseif action == 3
            birth(model.fields[field],σ,vnox,nodes2rays,rays2nodes,update,IP)
        elseif action == 4
            death(model.fields[field],σ,vnox,nodes2rays,rays2nodes,update,IP)
        elseif action == 5
            obs_id = rand(eachindex(observables.Obs))
            perturb_noise(update.noise[obs_id],σ,observables.Obs[obs_id],obs_id,update,IP)
            update.term[1] = update.term[1]/model.T[1]
        elseif action == 6
            print("\nstatics inversion running...\n")
            print(model.rms,"\n")
            @time solve_statics(model,observables,Gg,Gt,update.statics_residuals,update.statics_values)
	    end
        
        # -- validate model --
        update_predictions(model,observables,vnox,rays2nodes,rays_outdom,update)
        evaluate_model(model,observables,update,IP,vnox,nodes2rays,rays2nodes,rays_outdom)
        if action == 6
            print(model.rms,"\n")
        end
        if (((i % IP.MCS.saveint) .== 0))
            chainmodel = output_model(model)
            push!(MarkovChain, chainmodel)
        end
    end
    print("\nend: ",model.rms," - ",model.misfit,"\n")
end

function pertv(voronoi,σ,vnox,nodes2rays,rays2nodes,update,IP;ind=0) 
    if ind == 0
        ind = rand(1:voronoi.n[1])    
    end    
    update.ind .= ind
    old_v = voronoi.v[ind]
    dv = randn()*σ*(voronoi.vlims[2] - voronoi.vlims[1]) 
    voronoi.v[ind] = old_v + dv
    if IP.DR.status
        if update.trial_sample.n[1] == 1
            dvi = update.trial_sample.dv1
        elseif update.trial_sample.n[1] == 2
            dvi = update.trial_sample.dv2
        end
        dvi .= 0.0
        dvi .= dv
    end

    for nray in voronoi.nuclei2rays[ind]
        update.influence.rays[nray] = true            # -- the influence has not been changed, this may be not necessary
    end

    if voronoi.prior == "uniform"
        update.term .= 0.0
    elseif voronoi.prior == "normal"
        update.term .= -0.5*((voronoi.v[ind]-voronoi.vlims[3])^2-(old_v-voronoi.vlims[3])^2)/(voronoi.vlims[4])^2
        return
    elseif occursin("half-normal",voronoi.prior)
        # BPV ONLY VALID WHEN MEAN = 0.0!
        update.term .= -0.5*((voronoi.v[ind]-voronoi.vlims[3])^2-(old_v-voronoi.vlims[3])^2)/(voronoi.vlims[4])^2
        if voronoi.prior == "right-half-normal"
            if voronoi.v[ind] < voronoi.vlims[3]
		        update.term .= -Inf
            end
            return
        elseif voronoi.prior == "left-half-normal"
            if voronoi.v[ind] > voronoi.vlims[3]
		        update.term .= -Inf
            end
            return
        end
    end
    if !check_boundaries(voronoi.v[ind],voronoi.vlims)
        update.term .= -Inf
    end
end

function pertp(voronoi,σ,vnox,nodes2rays,rays2nodes,update,IP;ind=0) 
    if ind == 0
        ind = rand(1:voronoi.n[1])   
    end         
    update.ind .= ind
    # -- isolates the nuclei/rays/nodes influenced by the perturbation in the precedent diagram
    old_influence(voronoi,nodes2rays,rays2nodes,update)  
    # -- perturbs nucleus' position       
    r = sqrt(voronoi.c[1,ind]^2+voronoi.c[2,ind]^2+voronoi.c[3,ind]^2)
    θ = asin(voronoi.c[3,ind]/r)
    φ = atan(voronoi.c[2,ind],voronoi.c[1,ind])
    dr = randn() * σ * (voronoi.slims[3][2]-voronoi.slims[3][1])
    dθ = randn() * σ * (voronoi.slims[1][2]-voronoi.slims[1][1])
    dφ = randn() * σ * (voronoi.slims[2][2]-voronoi.slims[2][1])
    dt = 0.0
    r += dr
    θ += dθ
    φ += dφ

    # -- update diagram
    voronoi.c[3,ind] = r * sin(θ)
    voronoi.c[2,ind] = r * cos(θ) * sin(φ)
    voronoi.c[1,ind] = r * cos(θ) * cos(φ)
    voronoi.r[1,ind] = r
    if IP.B4DI.status 
        dt = randn() * σ * (voronoi.slims[4][2]-voronoi.slims[4][1]) 
        voronoi.c[4,ind] = dt + voronoi.c[4,ind]
    end
    if IP.DR.status
        if update.trial_sample.n[1] == 1
            ds = update.trial_sample.ds1
        elseif update.trial_sample.n[1] == 2
            ds = update.trial_sample.ds2
        end
        ds .= 0.0
        ds .= [dθ,dφ,dr,dt]
    end

    # -- isolates the nuclei/rays/nodes influenced by the perturbation in the current diagram
    new_influence(voronoi,nodes2rays,rays2nodes,update)      
    # -- backup map voronoi2rays (modified in the interpolation step)
    save_n2r(voronoi,update)
    # -- interpolation step, running on the isolated set of nuclei/nodes
    interpolate_diagram(voronoi,vnox,nodes2rays,update) 
    # -- computes priors ratio term
    t = 0.0
    IP.B4DI.status && (t = voronoi.c[4,ind])
    if !(check_boundaries(θ,voronoi.slims[1]) &&
        check_boundaries(φ,voronoi.slims[2]) &&
        check_boundaries(r,voronoi.slims[3]) &&
        check_boundaries(t,voronoi.slims[4])
       )
       update.term .= -Inf
    end
end

function birth(voronoi,σ,vnox,nodes2rays,rays2nodes,update,IP)
    if voronoi.n[1] == IP.MCS.max_nuclei 
        update.term[1] = -Inf
        return
    end

    newθ = (rand() * (voronoi.slims[1][2]-voronoi.slims[1][1])) + voronoi.slims[1][1]      # -- generates a new nucleus
    newφ = (rand() * (voronoi.slims[2][2]-voronoi.slims[2][1])) + voronoi.slims[2][1]
    newr = (rand() * (voronoi.slims[3][2]-voronoi.slims[3][1])) + voronoi.slims[3][1]
    newx = newr * cos(newθ) * cos(newφ)
    newy = newr * cos(newθ) * sin(newφ)
    newz = newr * sin(newθ)

    if check_1dim(voronoi)
        ind = v_dist(newr,voronoi)
    else
        if IP.B4DI.status
            newt = (rand() * (voronoi.slims[4][2]-voronoi.slims[4][1])) + voronoi.slims[4][1]
            ind = v_dist(newx,newy,newz,newt,voronoi)
        else
            ind = v_dist(newx,newy,newz,voronoi)    
        end
    end
    update.ind .= ind
    old_influence(voronoi,nodes2rays,rays2nodes,update)                             # -- recovers influence of the old diagram
    oldv = voronoi.v[ind]

    newv = randn()*σ*(voronoi.vlims[2] - voronoi.vlims[1]) + oldv 
    birth_term = log.(sqrt(2*pi)*σ) + (((newv - oldv)^2)/(2*(σ*(voronoi.vlims[2] - voronoi.vlims[1]))^2))
    
    voronoi.n[1] += 1
    voronoi.c[1,voronoi.n[1]] = newx
    voronoi.c[2,voronoi.n[1]] = newy
    voronoi.c[3,voronoi.n[1]] = newz
    voronoi.r[1,voronoi.n[1]] = newr
    IP.B4DI.status && (voronoi.c[4,voronoi.n[1]] = newt)
    voronoi.v[voronoi.n[1]] = newv

    update.ind[1] = voronoi.n[1]           
    new_influence(voronoi,nodes2rays,rays2nodes,update)          # -- saves all the rays,nodes and nuclei in the NEW diagram that may suffer this change
    save_n2r(voronoi,update)
    interpolate_diagram(voronoi,vnox,nodes2rays,update)

    if IP.MCS.nuclei_prior == "uniform"
        ratio_nuclei_prior = 0.0
    elseif IP.MCS.nuclei_prior == "log-uniform"
        ratio_nuclei_prior = log(voronoi.n[1]-1) - log(voronoi.n[1])
    end
    if voronoi.prior == "uniform"
        update.term .=  ratio_nuclei_prior + birth_term
    elseif voronoi.prior == "normal"
        σ0 = voronoi.vlims[4]
        birth_term += log(voronoi.vlims[2]-voronoi.vlims[1])
        update.term .= ratio_nuclei_prior - log(sqrt(2*pi)*σ0) - (((newv - voronoi.vlims[3])^2)/(2*(σ0)^2)) + birth_term
    elseif occursin("half-normal",voronoi.prior)
        σ0 = voronoi.vlims[4]
        birth_term += log(voronoi.vlims[2]-voronoi.vlims[1])
        update.term .= ratio_nuclei_prior - log(sqrt(2*pi)*σ0) - (((newv - voronoi.vlims[3])^2)/(2*(σ0)^2)) + log(2) + birth_term
    end
    
    if voronoi.prior == "normal"
        return
    elseif voronoi.prior == "right-half-normal"
        if newv < voronoi.vlims[3]
            update.term .= -Inf
        end
        return
    elseif voronoi.prior == "left-half-normal"
        if newv > voronoi.vlims[3]
            update.term .= -Inf
        end
        return
    end
    if !check_boundaries(newv,voronoi.vlims)
        update.term .= -Inf
    end
end

function death(voronoi,σ,vnox,nodes2rays,rays2nodes,update,IP)
    min_n = 0
    if check_1dim(voronoi)
        min_n = 1
    else
        if IP.B4DI.status
            min_n = 6
        else
            min_n = 4
        end
    end
    if voronoi.n[1] == min_n
        update.term .= -Inf 
        return 
    end

    ind = rand(1:voronoi.n[1])           
    oldv = voronoi.v[ind]
    oldx = voronoi.c[1,ind]
    oldy = voronoi.c[2,ind]
    oldz = voronoi.c[3,ind]
    oldr = voronoi.r[1,ind]
    IP.B4DI.status && (oldt = voronoi.c[4,ind])
    last = voronoi.n[1]
    if ind != last
        swap(voronoi,ind,last)  
        swap(update.old_voronoi,ind,last)                                          
    end     
    update.ind .= last                                        
    old_influence(voronoi,nodes2rays,rays2nodes,update)     
    voronoi.c[1,last] = 0.0
    voronoi.c[2,last] = 0.0
    voronoi.c[3,last] = 0.0
    voronoi.r[1,last] = 0.0
    IP.B4DI.status && (voronoi.c[4,last] = 0.0)
    voronoi.v[last] = 0.0
    voronoi.n[1] -= 1
    save_n2r(voronoi,update)
    interpolate_diagram(voronoi,vnox,nodes2rays,update)

    if check_1dim(voronoi)
        ind = v_dist(oldr,voronoi)
    else
        if IP.B4DI.status
            ind = v_dist(oldx,oldy,oldz,oldt,voronoi)
        else
            ind = v_dist(oldx,oldy,oldz,voronoi)
        end
    end
    newv = voronoi.v[ind]

    death_term = - log((sqrt(2*pi)*σ)) - (((oldv - newv)^2)/(2*(σ*(voronoi.vlims[2] - voronoi.vlims[1]))^2))

    if IP.MCS.nuclei_prior == "uniform"
        ratio_nuclei_prior = 0.0
    elseif IP.MCS.nuclei_prior == "log-uniform"
        ratio_nuclei_prior = log(voronoi.n[1]+1) - log(voronoi.n[1])
    end
    if voronoi.prior == "uniform"
        update.term .=  ratio_nuclei_prior + death_term
    elseif voronoi.prior == "normal"
        σ0 = voronoi.vlims[4]
        death_term -= log(voronoi.vlims[2]-voronoi.vlims[1])
        update.term .= ratio_nuclei_prior + log(sqrt(2*pi)*σ0) + (((oldv - voronoi.vlims[3])^2)/(2*(σ0)^2)) + death_term
    elseif occursin("half-normal",voronoi.prior)
        σ0 = voronoi.vlims[4]
        death_term -= log(voronoi.vlims[2]-voronoi.vlims[1])
        update.term .= ratio_nuclei_prior + log(sqrt(2*pi)*σ0) + (((oldv - voronoi.vlims[3])^2)/(2*(σ0)^2)) - log(2) + death_term
    end
end

function perturb_noise(noise,σ,observable,obsid,update,IP)
    update.obs_id .= obsid
    old_noise = noise[1]
    pert_noise = noise .+ randn()*σ*observable.noise_guess
    if pert_noise[1] <= 0.0
        update.term .= -Inf
        pert_noise[1] = noise[1]
    else
        noise .= pert_noise
        update.term .= length(observable.obsval)*(log(old_noise) - log(noise[1]))
    end
end

function update_predictions(model,observables,vnox,rays2nodes,rays_outdom,update)
    if update.term[1] == -Inf
        return
    end
    for i in eachindex(observables.Obs)
        predict(update, vnox, rays2nodes, rays_outdom, model, update.pred[i], observables.Obs[i], getfield(Main, observables.Obs[i].forward_name))
    end
end

function predict(update, vnox, rays2nodes, rays_outdom, model, pred_tmp, obs, forward::Function)
    for nray in eachindex(update.influence.rays)
        if update.influence.rays[nray]
            if haskey(obs.ray2obs,nray)
                update.raytmp.fields .= 0.0
                update.raytmp.u .= 0.0
                update.raytmp.u_path .= 0.0
                update.raytmp.dt .= 0.0
                forward(update.raytmp,vnox,rays2nodes,rays_outdom,nray,model,pred_tmp,obs)
            end
        end
    end
end

function evaluate_model(model,observables,update,IP,vnox,nodes2rays,rays2nodes,rays_outdom)
    p = IP.MCS.Lnorm
    update.misfit .= 0.0
    for i in eachindex(observables.Obs)
        events_corrections = @view model.dataspace.Obs[i].estatics[observables.Obs[i].obs2evt]
        stations_corrections =  @view model.dataspace.Obs[i].sstatics[observables.Obs[i].obs2sta]
        @. update.buffer[i] = abs(observables.Obs[i].obsval - 
                                    update.pred[i] - 
                                    events_corrections - 
                                    stations_corrections
                                    )^p
        update.rms[i] = (sum(update.buffer[i])/length(observables.Obs[i].obsval))^(1/p)
        update.misfit[1] += sum(update.buffer[i])/(update.noise[i][1]^p)
    end
    r = log(rand())
    acceptance_ratio = ((1/p)*(model.misfit[1] - update.misfit[1])/model.T[1] + update.term[1])
    if IP.DR.status && (update.trial_sample.n[1] == 2)
        acceptance_ratio = acceptance_ratio_delayed_rejection(model,update,IP)
    end
    if r <= acceptance_ratio
        obs_id = update.obs_id[1]
        model.rms .= update.rms
        model.misfit .= update.misfit
        model.accepted[1] += 1
        if update.action[1] < 5
            [observables.Obs[i].prdval .= update.pred[i] for i in eachindex(observables.Obs)] 
        elseif update.action[1] == 5
            model.dataspace.Obs[obs_id].noise .= update.noise[obs_id]
        end
    else
        if update.action[1] < 5
            copy_model(model,update)
        end
        if IP.DR.status
            if (update.action[1] != 1 && update.action[1] != 2) || (update.trial_sample.n[1] == 2)
                return
            end
            ind = update.ind[1]
            action = update.action[1]
            update.trial_sample.misfit .= update.misfit
            update.trial_sample.term .= update.term
            clear_updtmp(update,model,observables,update.field[1])
            delayed_rejection(model,update,observables,vnox,nodes2rays,rays2nodes,rays_outdom,IP,ind,action)
        end
    end
end

function temp_obj(model,rays,rnodes,observables,Gg,Gt,IP)
    buffer = Vector{Vector{Float64}}()
    [push!(buffer,zeros(Float64,length(obs.obsval))) for obs in observables.Obs]
    pred_tmp = Vector{Vector{Float64}}()
    [push!(pred_tmp,copy(observables.Obs[i].prdval)) for i in eachindex(observables.Obs)]
    noise_tmp = Vector{Vector{Float64}}()
    [push!(noise_tmp, [0.0]) for i in eachindex(observables.Obs)]
    statics_residuals = zeros(Float64,length(axes(Gt,2)))   # -- buffers for statics solver
    statics_values = zeros(Float64,length(axes(Gt,1)))
    influence_rays = Array{Bool}(undef,length(rnodes.r2n))
    influence_nodes = Array{Bool}(undef,length(rnodes.n2r))
    influence_nuclei = Array{Bool}(undef,IP.MCS.max_nuclei)
    neighbours = Set{Int64}()
    new_neighborhood = Set{Int64}()
    influence = InfluenceConst(influence_rays,influence_nodes,influence_nuclei,neighbours,new_neighborhood)
    n = zeros(Int64,1)
    c = zeros(Float64,4,IP.MCS.max_nuclei)
    r = zeros(Float64,1,IP.MCS.max_nuclei)
    v = zeros(Float64,IP.MCS.max_nuclei)
    nodes2nuclei = zeros(Int64,length(rnodes.n2r))
    misfit = zeros(Int64,1)
    rms = zeros(Float64,length(model.rms))
    nuclei2rays = Vector{Set{Int64}}()
    old_voronoi = Voronoi("","",n,c,r,v,nuclei2rays,nodes2nuclei,[0.0],[[0.0]],0.0)
    max_L = 0
    for nray in eachindex(rays)
        max_L = max(max_L,length(rays[nray].L)+1)
    end
    fields = zeros(Float64,length(model.fields),max_L)
    u = zeros(Float64,max_L)
    u_path = zeros(Float64,max_L-1)
    dt = zeros(Float64,max_L-1)
    raytmp = RayTmpConst(fields,u,u_path,dt)

    trial_sample = TrialSampleConst([1],[0.0],[0.0],zeros(Float64,4),zeros(Float64,4),[0.0],[0.0])

    update = UpdateConst([0], [0], [0], [0], [0.0], buffer,pred_tmp,noise_tmp,statics_residuals,statics_values,Gg,Gt,influence,old_voronoi,misfit,rms,raytmp,trial_sample)

    return update
end

function initialize_parallel_chains(IP,IP_fields,rnodes,evtsta,observables,rays,Gg,Gt,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
    parallel_chains = Vector{ObjsInChainConst}()
    MarkovChains = Vector{Vector{ModelConst}}()
    for i in 1:nchains
        model = initialize_model(IP,IP_fields,rnodes,evtsta,observables,rays,Gg,Gt,vnox,nodes2rays,rays2nodes,rays_outdom)
        update = temp_obj(model,rays,rnodes,observables,Gg,Gt,IP)
        if IP.PT.status && (i > IP.PT.ncold)
            n = i - IP.PT.ncold + 1
            x = log10(IP.PT.maxT)
            model.T .= 10^(x*(n-1)/((nchains-IP.PT.ncold)))
        end
        push!(MarkovChains,Vector{ModelConst}())
        chainmodel = output_model(model)
        push!(MarkovChains[end],chainmodel)
        ObjsInChain = ObjsInChainConst(
            model,
            deepcopy(observables),
            update
        )
        push!(parallel_chains,ObjsInChain)
    end
    distributed_parallel_chains = distribute(parallel_chains)
    distributed_markov_chains = distribute(MarkovChains)
    parallel_chains = nothing
    MarkovChains = nothing
    return distributed_parallel_chains, distributed_markov_chains        
end

function reinitialize_models(ObjsInChains,rnodes,rays,observables,Gg,Gt,IP,nchains,vnox,nodes2rays,rays2nodes,rays_outdom)
    @sync @distributed for chain in 1:nchains
        ObjsInChain = localpart(ObjsInChains)[1]
        model_chain = ObjsInChain.model
        update = temp_obj(model_chain,rays,rnodes,observables,Gg,Gt,IP)
        observables_chain = ObjsInChain.observables
        for voronoi in model_chain.fields
            deleteat!(voronoi.nodes2nuclei, 1:length(voronoi.nodes2nuclei))
            for i in eachindex(voronoi.nuclei2rays)
                [delete!(voronoi.nuclei2rays[i],j) for j in voronoi.nuclei2rays[i]]
            end
            map_voro2space2(voronoi,rnodes)
        end
        for i in eachindex(observables_chain.Obs)
            observables_chain.Obs[i].obsval .= observables.Obs[i].obsval 
        end
        initialize_fitstats(model_chain,observables_chain,rnodes,rays,IP,Gg,Gt,vnox,nodes2rays,rays2nodes,rays_outdom)
        ObjsInChain.update = update
    end
end

function clear_updtmp(update, model, observables, field)
    # -- clear temporary variables
    for i in eachindex(observables.Obs)
        update.pred[i] .= observables.Obs[i].prdval
        update.noise[i] .= model.dataspace.Obs[i].noise
    end
    update.ind .= 0
    update.obs_id .= 0
    update.term .= 0.0
    [update.buffer[i] .= 0.0 for i in eachindex(update.buffer)]
    # -- clear influence
    update.influence.rays .= false
    update.influence.nodes .= false
    update.influence.nuclei .= false
    [delete!(update.influence.neighbours,i) for i in update.influence.neighbours]
    [delete!(update.influence.new_neighborhood,i) for i in update.influence.new_neighborhood]
    # -- clear old voronoi
    voronoi = model.fields[field]
    old_voronoi = update.old_voronoi
    old_voronoi.n[1] = voronoi.n[1]
    old_voronoi.v .= 0.0
    old_voronoi.r .= 0.0
    old_voronoi.c .= 0.0
    for i in axes(voronoi.c,2)[1:voronoi.n[1]]
        old_voronoi.v[i] = voronoi.v[i]
        old_voronoi.r[1,i] = voronoi.r[1,i]
        for j in axes(voronoi.c,1)
            old_voronoi.c[j,i] = voronoi.c[j,i]
        end
    end
    deleteat!(old_voronoi.nuclei2rays,1:length(old_voronoi.nuclei2rays))
    old_voronoi.nodes2nuclei .= voronoi.nodes2nuclei
    update.misfit .= 0.0
    update.rms .= 0.0
    update.trial_sample.n[1] = 1
end

function print_progress(IP,model,i,chain_id)
    if i % IP.MCS.printon .== 0
            print(
                  string("chain #", string(chain_id), " at ",
                  string(round(100*(i/IP.MCS.niterations); digits=4)), "; rms of ", 
                  string(round.(1000*model.rms; digits=4)),"[ms], accepted: ",round(100*(model.accepted[1]/(i))),"%\n")
                  )
    end
end
