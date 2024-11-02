# This project has been developed with the idea that every user should be able to implement 
# observations and models of interest (compatible with the seismological core of this code), without modifying
# the inner parts of this package. Here is the factory that builds all the objects needed keeping 
# some generality in the approach and leaving space for future developments.

function build_obslist(B)
    # Add observations
    IP_obs = Vector{obsinfo}()
    for D in B # Loop over datasets to fit
        add_obs(IP_obs, D["type"], D["file"], D["tf_event_demean"], D["tf_event_statics"],
            D["tf_station_statics"], D["data_uncertainty"], D["forward_function"])
    end
    return IP_obs
end

function build_fieldslist(F)
    # Add parameter fields
    IP_fields = Vector{fieldinfo}()
    for D in F # Loop over parameter fields to estimate
        add_field(IP_fields, D["field"], D["prior"], D["pdf_parameters"], D["lat_lims"],
            D["lon_lims"], D["elv_lims"], D["time_lims"], D["extrap_value"])
    end
    return IP_fields 
end

function build_inputs(parameter_file::String; tf_save = false, tf_serial = false)
    # Load parameter file
    P = TOML.parsefile(parameter_file)
    if tf_serial
        # For serial calculations we run 1 chain at a time
        num_chains = P["MonteCarloSolver"]["num_chains"]*P["MonteCarloSolver"]["num_chunks"]
        P["MonteCarloSolver"]["num_chains"] = 1
        P["MonteCarloSolver"]["num_chunks"] = num_chains
    end
    # Build inputs
    IP_obs = build_obslist(P["DataFields"])
    IP_fields = build_fieldslist(P["ParameterFields"])
    IP = IPConst(P)
    rnodes, rays, evtsta, observables, LocalRaysManager = initialize_IP(IP, IP_obs, IP_fields)
    # Save inputs
    if tf_save
        # Make output directory
        run_directory = P["out_directory"]*"/"*P["run_id"]
        mkpath(run_directory)
        # Save copy of parameter file
        cp(parameter_file, run_directory*"/parameters.toml")
        save(
            string(run_directory, "/IPConst.jld"), "IP", IP,
            "IP_fields", IP_fields,
            "rnodes", rnodes,
            "rays", rays,
            "evtsta", evtsta,
            "observables", observables,
            "LocalRaysManager", LocalRaysManager
        )
    end

    return P, IP_obs, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager
end

# -- this script creates the setup for the inverse problem (IP).
function initialize_IP(IP, IP_obs, IP_fields)

    # -- conversion deg to rad for domain's boundaries
    @. IP.lims.lat = deg2rad(IP.lims.lat)
    @. IP.lims.lon = deg2rad(IP.lims.lon)

    # -- rotation reference frame for azimuth and elevation
    φ, θ = 0.5*(IP.lims.lon[begin]+IP.lims.lon[end]), 0.5*(IP.lims.lat[begin]+IP.lims.lat[end])
    Rx, Ry = rotation_matrices(φ,θ)
    lonm, latm = mean(IP.lims.lon), mean(IP.lims.lat)
    Px, Py, Pz = [1 * cos(latm) * cos(lonm)], [1 * cos(latm) * sin(lonm)], [1 * sin(latm)]
    rotate_basis(Py,Pz,Px,Rx,Ry)
    print("\n",Px," ",Py," ",Pz)

    # -- extraction Voronoi diagrams boundaries
    voronoi_slims = []
    for fieldid in eachindex(IP_fields)
        push!(voronoi_slims,[])
        [push!(voronoi_slims[fieldid],zeros(Float64,2)) for i in 1:3]
        field = IP_fields[fieldid]
        R1, R2 = field.depthlims[1], field.depthlims[2]
        θ1, θ2 = (field.θlims[1]), (field.θlims[2])
        φ1, φ2 = (field.φlims[1]), (field.φlims[2])
        voronoi_slims[fieldid][1][1], voronoi_slims[fieldid][1][2] = θ1, θ2
        voronoi_slims[fieldid][2][1], voronoi_slims[fieldid][2][2] = φ1, φ2
        voronoi_slims[fieldid][3][1], voronoi_slims[fieldid][3][2] = R1, R2
    end

    
    # -- populates event and station objects
	ievt = readdlm(IP.InputEvt)   # -- reads the events file
	ista = readdlm(IP.InputSta)   # -- reads the stations file
    nevt = length(ievt[:,1])    # -- number of events
    nsta = length(ista[:,1])    # -- number of stations
    evts = Vector{evt}()
    stas = Vector{sta}()
    for i in eachindex(ievt[:,1])
        T0 = 0.0
        IP.B4DI.status && (T0 = ievt[i,5])
        evt_i = evt(ievt[i,1],deg2rad(ievt[i,2]),deg2rad(ievt[i,3]),-ievt[i,4],T0)
        push!(evts,evt_i)
    end
    for i in eachindex(ista[:,1])
        sta_i = sta(ista[i,1],deg2rad(ista[i,2]),deg2rad(ista[i,3]),ista[i,4])
        push!(stas,sta_i)
    end
    evtsta = evtstaConst(evts,stas)         # -- events&stations structure is created

    # -- populates observables objects
    phases = collect_phases(IP_obs)         # -- this array collects all the seismic phases that appear in the observables provided 
    nphases = length(phases)

    evtsta2rays = zeros(Int64,nevt,nsta,nphases)    # -- [evt_id,sta_id,phase matrix] to keep trace of rays already registered
    paths_list = Vector{Vector}()

    obsvec = Vector{ObsConst}()         # -- vector of observables used in the misfit function (predictive)
    descriptive = Vector{ObsConst}()    # -- vector of observables NOT used in the misfit function (descriptive, e.g. polarization)

    obslist = Dict{String, Int64}()     # -- dictionary observables' names -> indeces in the vectors created
    obs_id = [0]
    newray = 0

    # -- loop over input observables
    for obs in IP_obs
        if obs.forward_name != ""       # -- check if obs is predictive
            obs_id[1] += 1
            obslist[obs.name] = obs_id[1]
        end
        iobs = readdlm(obs.file)        # -- reads observable's data file

        evts_obs = Set{Int64}()         # -- collects events associated with the observable
        stas_obs = Set{Int64}()         # -- collects stations associated with the observable

        obs2evt = Vector{Int64}()       # -- associate every observable's measurement to an event
        obs2sta = Vector{Int64}()       # -- associate every observable's measurement to a station

        obsvalues = Vector{Float64}()   # -- measurements array
        prdvalues = Vector{Float64}()   # -- predictions array

        ray2obs = Dict{Int64, Int64}()  # -- rays to observed values

        for i in axes(iobs,1)           # -- loops over all the observed values for observable obs

            evt_id = Int(iobs[i,1])
            sta_id = Int(iobs[i,2])
            phase_name = iobs[i,4]
            phase_id = phases[phase_name]   # -- convert phase to integer ID using the matrix phases
            push!(obsvalues,iobs[i,3])      # -- push measurement's value
            push!(prdvalues,0.0)            # -- initializing prediction value to 0

            # -- check if ray already exists and associate ray-ID to measured value
            if (evtsta2rays[evt_id,sta_id,phase_id] == 0)
                newray += 1 
                evtsta2rays[evt_id,sta_id,phase_id] = newray    # new ray is registered in the ray-list
                push!(paths_list,[evt_id,sta_id,phase_name])
                ray2obs[newray] = i                             # associates the ray to the observable measurement
            else
                oldray = evtsta2rays[evt_id,sta_id,phase_id]    # for this combination a ray has already been registered in the ray-list
                ray2obs[oldray] = i                             # this existing ray is associated also with the new measurement
            end

            push!(evts_obs,evt_id)
            push!(stas_obs,sta_id)

            push!(obs2evt,evt_id)
            push!(obs2sta,sta_id)
        end

        # -- evts_obs and stas_obs are Set (unique) of events and stations' IDs registered for the observable obs. These sets typically do not
        #    collect ALL the possible evt and sta IDS, they are vectorized (so now they are ordered), and in the following steps two maps 
        #    are defined to point from obs -> vector(evt_id), vector(sta_id)

        evts_obs = collect(evts_obs)
        sort!(evts_obs)
        stas_obs = collect(stas_obs)
        sort!(stas_obs)

        # -- Absolute evt and sta IDS are not needed, so in the maps from obs to events and stations the indeces of the local vectors are stored.
        # -- creates maps from measur. to the array of evts and stas (only storing evts and stas of this observable!!)
        for i in eachindex(evts_obs)
            for j in eachindex(obs2evt)
                (obs2evt[j] == evts_obs[i]) && (obs2evt[j] = i)
            end
        end
        for i in eachindex(stas_obs)
            for j in eachindex(obs2sta)
                (obs2sta[j] == stas_obs[i]) && (obs2sta[j] = i)
            end
        end
        

        # -- initialize static structures
        # obs_evts = zeros(Float64,length(evts_obs))
        # obs_stas = zeros(Float64,length(stas_obs))

        observable = ObsConst(
            obs.name,
            Symbol(obs.forward_name),
            obsvalues,
            prdvalues,
            zeros(Float64,length(prdvalues)),
            evts_obs,
            stas_obs,
            obs2evt,
            obs2sta,
            obs.noise_guess,
            ray2obs,
            obs.demean,
            obs.eventstatics,
            obs.stationstatics
        )

        if obs.forward_name != ""           # -- obs is pushed in predictive/descriptive vector
            push!(obsvec,observable)
        else
            push!(descriptive,observable)
        end
    end

    observables = ObservablesConst(         # -- House of the Observables
        length(obsvec),
        obslist,
        obsvec,
        )

    
    # -- the code compares the number of local earthquakes with the number of stations, and decides which ones to use as sources and receivers.
    # -- is relocation status is ON, stations are always the sources
    refm = readdlm(IP.velocitymodel,' ',Float64,'\n')   # -- reads the reference 1D velocity model
    if IP.RayTracingInit.allTauP
        local_evts, local_stats = Set{Int64}(), Set{Int64}()
    else
        local_evts, local_stats = indomain(evtsta,paths_list,IP.lims)  # -- is local ray-tracing required? (at least one event inside the inversion domain) 
    end
    ray2source_receiver = Dict{Int64,Vector{Int64}}()   # associates every ray-ID to [source_id,station_id,phase]
    pair2ray = Dict{Vector{Int64}, Int64}()         # associates every pair source-receiver to a ray-ID
    source2receivers = Dict{Int64,Vector{Int64}}()  # associates every source to the set of receivers
    nnodes, fw_level, pert, carving, sub_its, relocations = IP.RayTracingInit.nnodes, IP.RayTracingInit.fw_level, IP.RayTracingInit.perturb, IP.RayTracingInit.carving, IP.RayTracingInit.sub_its, IP.EQRLOC.relocate

    if  ((length(local_evts) <= length(local_stats)) && !relocations) || IP.B4DI.status
        for evt_id in local_evts
            receivers = Int64[]
            for i in eachindex(paths_list)
                if paths_list[i][1] == evt_id
                    sta_id = paths_list[i][2]
                    (paths_list[i][3] == "P") && (phase = 1)
                    (paths_list[i][3] == "S") && (phase = 2)
                    ray2source_receiver[i] = [evt_id,sta_id,phase]
                    push!(receivers,sta_id)
                    pair2ray[[evt_id,sta_id,phase]] = i
                end
            end
            source2receivers[evt_id] = receivers
        end
    else
        for sta_id in local_stats
            receivers = Int64[]
            for i in eachindex(paths_list)
                if paths_list[i][2] == sta_id
                    evt_id = paths_list[i][1]
                    (paths_list[i][3] == "P") && (phase = 1)
                    (paths_list[i][3] == "S") && (phase = 2)
                    ray2source_receiver[i] = [sta_id,evt_id,phase]
                    push!(receivers,evt_id)
                    pair2ray[[sta_id,evt_id,phase]] = i
                end
            end
            source2receivers[sta_id] = receivers
        end
    end

    # -- RayTracing status
    RT_status = false
    (length(local_evts) != 0) && (RT_status = true)

    source_nodes, receiv_nodes = Dict{Int64, Int64}(), Dict{Int64, Int64}()
    LocalRaysManager = LocalRaysManagerConst(RT_status,nnodes,fw_level,pert,carving,ray2source_receiver,source2receivers,pair2ray,local_evts,local_stats,source_nodes,receiv_nodes,sub_its,relocations)    
    
    paths = Array{pathConst,1}(undef,length(paths_list))
    rays = Array{RayConst,1}(undef,length(paths_list))
    initialize_paths(paths,paths_list,LocalRaysManager,IP,evtsta,observables,ievt,ista,refm)

    rnodes = build_rays(paths,rays,observables,evtsta,LocalRaysManager,IP,refm)
    voronoi_domain_pierce(rays,rnodes,voronoi_slims,IP)

    for obs in descriptive
        if obs.obsname == "polarization"
            for nray in eachindex(rays)
                if haskey(obs.ray2obs,nray)
                    @. rays[nray].ζ = obs.obsval[obs.ray2obs[nray]]
                    # correct_polarization(rnodes,rays[nray],nray,evtsta)
                end
            end
        end
    end

    return rnodes, rays, evtsta, observables, LocalRaysManager

end

function build_rays(paths,rays,observables,evtsta,LocalRaysManager,IP,refm)
    rayn_x = Vector{Float64}()
    rayn_y = Vector{Float64}()
    rayn_z = Vector{Float64}()
    rayn_r = Vector{Float64}()
    rayn_t = Vector{Float64}()
    vps = Vector{Float64}()
    vss = Vector{Float64}()
    rays2nodes = Vector{Vector{Int64}}()   # -- associates each ray to the corresponding nodes
    nodes2rays = Vector{Int64}()           # -- associates each node to the corresponding ray

    for nray in eachindex(paths)
        evt = paths[nray].evt
        sta = paths[nray].sta
        phase = paths[nray].phase
        lat = paths[nray].lat 
        lon = paths[nray].lon 
        rad = paths[nray].rad 
        times = paths[nray].t

        # -- evaluates where the rays are piercing the inversion domain
        lat, lon, rad, timepath = domain_pierce(lat,lon,rad,times,IP)

        x = @. rad * cos(lat) * cos(lon)
        y = @. rad * cos(lat) * sin(lon)
        z = @. rad * sin(lat)

        x, y, z = raylin_interpol(IP,x,y,z)
        dx, dy, dz = diff(x), diff(y), diff(z)
        L = @. sqrt(dx^2+dy^2+dz^2)

        rads = @. sqrt(x^2+y^2+z^2)

        pre_ind = length(vps) + 1
        [push!(vps,ref_V1D(rads[l],refm[:,[1,2]])) for l in eachindex(rads)]
        [push!(vss,ref_V1D(rads[l],refm[:,[1,3]])) for l in eachindex(rads)]
        lst_ind = length(vps)

        ray = RayConst(                     # -- see rayConst structure in structures.jl for details
            evt,
            sta,
            phase,
            L, 
            zeros(Float64,length(x)),
            zeros(Float64,length(x)), 
            zeros(Float64,length(dx)),
            zeros(Float64,length(dx)),
            zeros(Float64,length(x)), 
            [timepath],
            []
        )

        for obs in observables.Obs
            if haskey(obs.ray2obs,nray)
                obs.ref_t[obs.ray2obs[nray]] = ray.t_1D[1]
            end
        end

        for i in eachindex(x)[begin:end-1]
            φ, θ = atan(y[i],x[i]), asin(z[i]/sqrt(x[i]^2+y[i]^2+z[i]^2))
            Rx, Ry = rotation_matrices(φ,θ)
            R = Ry*Rx
            b = inv(R)*[dy[i],dz[i],dx[i]]
            dy[i], dz[i], dx[i] = b[1], b[2], b[3]
        end

        @. ray.ϕ = atan(dz,dy)          # ray-path azimuth
        @. ray.θ = asin(dx/ray.L)       # ray-path elevation

        push!(ray.ϕ,ray.ϕ[end])
        push!(ray.θ,ray.θ[end])

        @. ray.v_1D_P = vps[pre_ind:lst_ind]    # -- 1D velocity along the ray-path
        @. ray.v_1D_S = vss[pre_ind:lst_ind]

        rays[nray] = ray
        
        [push!(rayn_x, x[l]) for l in eachindex(x)]
        [push!(rayn_y, y[l]) for l in eachindex(y)]
        [push!(rayn_z, z[l]) for l in eachindex(z)]
        [push!(rayn_t, evtsta.evts[evt].T0) for l in eachindex(x)]
        [push!(rayn_r, rads[l]) for l in eachindex(rads)]

        push!(rays2nodes,range(length(rayn_x)-length(x)+1,length(rayn_x),step=1))
        [push!(nodes2rays,nray) for j in eachindex(x)]

    end

    coords = permutedims(hcat(rayn_x,rayn_y,rayn_z,rayn_t))
    rad_coords = permutedims(hcat(rayn_r))

    rnodes = NodesConst(
        coords, 
        rad_coords, 
        rays2nodes,
        nodes2rays
        )

    return rnodes

end


function initialize_model(IP,IP_fields,rnodes,evtsta,observables,rays,Gg,Gt,vnox,nodes2rays,rays2nodes,rays_outdom)

    refm = readdlm(IP.velocitymodel,' ',Float64,'\n')   # -- reads the reference 1D velocity model
    
    vorovec = Vector{Voronoi}()      # -- vector of Voronoi diagram to build
    nfields = length(IP_fields)
    nrays = length(rnodes.r2n)
    ini_nuclei = IP.MCS.ini_nuclei
    diagram_id = 0
    fieldslist = Dict{String, Int64}()

    for field in IP_fields
        diagram_id += 1
        if !IP.B4DI.status 
            field.timelims[1] = 0.0 
            field.timelims[2] = 0.0
        end
        R1, R2 = field.depthlims[1], field.depthlims[2]
        θ1, θ2 = (field.θlims[1]), (field.θlims[2])
        φ1, φ2 = (field.φlims[1]), (field.φlims[2])
        t1, t2 = field.timelims[1], field.timelims[2]

        nodes2nuclei = Int64[]
        fieldslist[field.name] = diagram_id

        if IP.B4DI.status
            coords = zeros(Float64,4,IP.MCS.max_nuclei)
        else
            coords = zeros(Float64,3,IP.MCS.max_nuclei)
        end
        rad_coords = zeros(Float64,1,IP.MCS.max_nuclei)
        v = zeros(Float64,IP.MCS.max_nuclei)

        for i in 1:ini_nuclei
            θs = (θ2 - θ1)*rand() .+ θ1
            φs = (φ2 - φ1)*rand() .+ φ1
            rs = (R2 - R1)*rand() .+ R1
            t = (t2 - t1)*rand() .+ t1
            # x, y, z = @cartesian(θs,φs,rs)
            x, y, z = geo_to_cartesian(θs,φs,rs)
            vs = field.init_val
            if IP.MCS.rand_init_values
                if field.prior == "uniform"
                    vs = (field.vlims[2] - field.vlims[1])*rand() .+ field.vlims[1]
                elseif field.prior == "normal"
                    vs = field.vlims[4]*randn() .+ field.vlims[3]
                elseif field.prior == "right-half-normal"
                    vs = field.vlims[4]*randn() .+ field.vlims[3]
                    while vs < field.vlims[3]
                        vs = field.vlims[4]*randn() .+ field.vlims[3]
                    end
                elseif field.prior == "left-half-normal"
                    vs = field.vlims[4]*randn() .+ field.vlims[3]
                    while vs > field.vlims[3]
                        vs = field.vlims[4]*randn() .+ field.vlims[3]
                    end
                else
                    print("\nprior not recognized!")
                    return
                end
            end
            coords[1,i] = x
            coords[2,i] = y
            coords[3,i] = z
            IP.B4DI.status && (coords[4,i] = t)
            rad_coords[1,i] = rs
            v[i] = vs
        end

        nuclei2rays = Vector{Set{Int64}}()
        [push!(nuclei2rays,Set{Int64}()) for i in 1:IP.MCS.max_nuclei]

        voronoi = Voronoi(      # -- builds Voronoi diagram using all the information provided
            field.name,
            field.prior,
            [ini_nuclei],
            coords,
            rad_coords,
            v,
            nuclei2rays,
            nodes2nuclei,
            field.vlims,
            [[θ1, θ2],[φ1, φ2],[R1, R2], [t1, t2]],
            field.init_val
    )

        # if check_1dim(voronoi)
        #     for ind in axes(voronoi.r,2)[1:voronoi.n[1]]
        #         if voronoi.fieldname == "Vp"
        #             voronoi.v[ind] = ref_V1D(voronoi.r[1,ind],refm[:,[1,2]])
        #         elseif voronoi.fieldname == "Vs"
        #             voronoi.v[ind] = ref_V1D(voronoi.r[1,ind],refm[:,[1,3]])
        #         end
        #     end
        # end
        map_voro2space2(voronoi,rnodes)  # -- for each v-cell identifies the ray-paths and ray-nodes inside the influence area.
        push!(vorovec,voronoi)
    end

    dataspace = DataSpaceConst(
            Vector{DataPConst}()
        )
    for obs in observables.Obs
        estatics = zeros(Float64,length(obs.evtids))
        sstatics = zeros(Float64,length(obs.staids))
        noise = [obs.noise_guess]
        datap = DataPConst(
            estatics,
            sstatics,
            noise
        )
        push!(dataspace.Obs,datap)
    end

    model = ModelConst(
        nfields, 
        fieldslist, 
        vorovec, 
        zeros(Float64,1), 
        zeros(Float64,length(observables.obslist)),
        [1.0],
        [1],
        dataspace
        )

    # -- initialize rms and misfit (expand initial predictions computation (to do))
    initialize_fitstats(model,observables,rnodes,rays,IP,Gg,Gt,vnox,nodes2rays,rays2nodes,rays_outdom; firstcall=true)

    return model


end

function initialize_paths(paths,paths_list,LocalRaysManager,IP,evtsta,observables,ievt,ista,refm)
    nrays = length(axes(paths_list,1))
    s_r = LocalRaysManager.ray2source_receiver
    Rmax = IP.RayTracingInit.Rmax
    if (length(s_r) != 0)
        grid = instance_grid(evtsta, LocalRaysManager, IP)
        raytracing_dijkstra(grid, observables, LocalRaysManager, paths, evtsta, IP)
    end
    for i in axes(paths_list,1)
        print("progress: ",round(i/nrays*100; digits=2),"%\r")
        if haskey(s_r,i)
            continue
        else
            PathObj = buildPathObj(IP.RayTracingInit.TauP_Model)
            event, station, phase = paths_list[i][1], paths_list[i][2], paths_list[i][3]
            evt_depth = (Rmax-R) + ievt[event,4]
            sta_depth = (Rmax-R) - ista[station,4]
            d,r,t,ϕ,λ,e = taup_path!(PathObj,String(phase),ievt[event,2],ievt[event,3],evt_depth,ista[station,2],ista[station,3],sta_depth)
            lon = deg2rad.(λ)
            lat = deg2rad.(ϕ)
            rad = Rmax .- (r)
            paths[i] = pathConst(event,station,phase,lat,lon,rad,t)
        end
    end
end

function raylin_interpol(IP,x,y,z)

    dx = diff(x)
    dy = diff(y)
    dz = diff(z)

    dLi = zeros(Float64,length(dx))
    dLq = IP.RayTracingInit.dray
    @. dLi = sqrt((dx^2) + (dy^2) + (dz^2))
    Li  = vcat(0,cumsum(dLi))
    n_fine = convert(Int64,trunc((Li[end]-Li[1])/dLq, digits=0))
    if n_fine == 0
        return x, y, z
    end
    Lq   = range(Li[1],Li[end],length=1 + n_fine)
    if IP.RayTracingInit.itp == "linear"
        itpx  = LinearInterpolation(Li,x)
        itpy  = LinearInterpolation(Li,y)
        itpz  = LinearInterpolation(Li,z)
    elseif IP.RayTracingInit.itp == "spline"
        itpx  = Spline1D(Li,x)
        itpy  = Spline1D(Li,y)
        itpz  = Spline1D(Li,z)
    else
        print("\ninterpolation scheme not recognized")
        return
    end
    x_int = [itpx(p) for p in Lq]
    y_int = [itpy(p) for p in Lq]
    z_int = [itpz(p) for p in Lq] 
    
    return x_int, y_int, z_int
    
end


function ref_V1D(rad,refm)

    rr = refm[:,1]
    for i in eachindex(rr)
        if i == 1
            continue
        end
        if (rr[i-1] > rad >= rr[i]) 
            return refm[i-1,2] + (rad-rr[i-1])*(refm[i,2]-refm[i-1,2])/(rr[i]-rr[i-1])
        end
    end
    if rad >= rr[1]
        return refm[1,2]
    end
    if rad < rr[end]
        return refm[end,2]
    end
    return refm[end,2]

end


function domain_pierce(lat, lon, rad, times, IP)
    lat1,lat2,lon1,lon2,rad1,rad2 = 0.0,0.0,0.0,0.0,0.0,0.0
    rad1,rad2 = R + IP.lims.depth[1],R + IP.lims.depth[2]
    lat1,lat2 = IP.lims.lat[1],IP.lims.lat[2]
    lon1,lon2 = IP.lims.lon[1],IP.lims.lon[2]
    global_idx = 0
    for i in eachindex(lat)
        global_idx += 1
        if ((lat1 <= lat[i] <= lat2) && (lon1 <= lon[i] <= lon2) && (rad1 <= rad[i] <= rad2))
            break
        end
    end
    timepath = times[end] - times[global_idx]
    
    return lat[global_idx:end], lon[global_idx:end], rad[global_idx:end], timepath
end

function voronoi_domain_pierce(rays,rnodes,voronoi_slims,IP)
    for nray in eachindex(rays)
        ray = rays[nray]
        nodes = rnodes.r2n[nray]
        r = rnodes.r[nodes]
        r = round.(r,digits=1)
        for field in voronoi_slims
            push!(ray.outdomain,zeros(Int64,2))
            if !IP.MCS.squeezing
                continue
            end
            if field[3][1] <= r[begin] <= field[3][2]
                continue
            else
                ray.outdomain[end][1] = 1
                for i in eachindex(r)[2:end]
                    if field[3][1] < r[i] < field[3][2]
                        ray.outdomain[end][2] = i
                        break
                    end
                end
            end
        end
    end
end

function check_coordinates(field,θ,φ,r)
    return  check_boundaries(θ,field[1]) &&
            check_boundaries(φ,field[2]) &&
            check_boundaries(r,field[3])
end

function check_boundaries(v,vlims)
    if vlims[1] == vlims[2]
        return true
    end
    if ((v < vlims[1]) || (v > vlims[2]))
        return false     # -- out of boundaries
    else
        return true
    end
end

function check_1dim(voronoi)
    if (voronoi.slims[1][1] == voronoi.slims[1][2]) && (voronoi.slims[2][1] == voronoi.slims[2][2])
        return true
    else
        return false
    end
end

function collect_phases(IP_obs)
    phases_set = Set{String}()
    for obs in IP_obs
        iobs = readdlm(obs.file)
        for j in eachindex(iobs[:,1])
            phase = iobs[j,end]
            push!(phases_set,phase)
        end
    end
    phases_vec = collect(phases_set)
    phases_dict = Dict{String, Int64}()
    for i in eachindex(phases_vec)
        phases_dict[phases_vec[i]] = i
    end

    return phases_dict
end

function add_field(list,name,prior,_vlims,_θlims,_φlims,_depthlims,_timelims,_init_val)
    added_field = fieldinfo(
        name,
        prior,
        _vlims,
        deg2rad.(_θlims),
        deg2rad.(_φlims),
        R .+ (_depthlims),
        _timelims,
        _init_val
    )
    push!(list,added_field)
end

function add_obs(list,name,file,demean,evt_stats,sta_stats,noise,forward_name)
    added_field = obsinfo(
        name,
        file,
        demean,
        evt_stats,
        sta_stats,
        noise,
        forward_name
    )
    push!(list,added_field)
end

function initialize_fitstats(model,observables,rnodes,rays,IP,Gg,Gt,vnox,nodes2rays,rays2nodes,rays_outdom; firstcall=false)
    n = IP.MCS.Lnorm
    misfit = 0.0
    obs_id = 0

    update = temp_obj(model,rays,rnodes,observables,Gg,Gt,IP)
    update.influence.rays .= true
    update_predictions(model,observables,vnox,rays2nodes,rays_outdom,update)
    [observables.Obs[i].prdval .= update.pred[i] for i in eachindex(observables.Obs)]

    # -- initialize statics
    for ii in eachindex(observables.Obs)
        obs = observables.Obs[ii]
        if obs.demean && firstcall
            for i in eachindex(model.dataspace.Obs[ii].estatics)
                residuals = Float64[]
                for j in eachindex(obs.obsval)
                    (obs.obs2evt[j] == i) && (push!(residuals,obs.obsval[j]-obs.prdval[j]))
                end
                if length(residuals) > 0
                     for j in eachindex(obs.obsval)
                        if obs.obs2evt[j] == i
                            obs.obsval[j] = obs.obsval[j] - mean(residuals)
                        end
                     end
                end
            end
        end
    end

    for obsid in eachindex(observables.Obs)
        obs = observables.Obs[obsid]
        obs_id += 1
        sum_sq = sum(@. abs(obs.obsval - 
        obs.prdval - 
        model.dataspace.Obs[obsid].estatics[obs.obs2evt] - 
        model.dataspace.Obs[obsid].sstatics[obs.obs2sta]
        )^n)
        rms = sqrt(sum_sq/length(obs.obsval))
        misfit += (sum_sq / (model.dataspace.Obs[obsid].noise[1])^n)
        model.rms[obs_id] = rms
    end
    model.misfit .= misfit

end

function raytime_1D(ray)
    if ray.phase == "P"
        @. ray.u = 1.0/ray.v_1D_P
        average_velocity(ray.u_path,ray.u)
    elseif ray.phase == "S"
        @. ray.u = 1.0/ray.v_1D_S
        average_velocity(ray.u_path,ray.u)
    else
        return 0.0
    end
    return sum(@. ray.L * ray.u_path)
end

function indomain(evtsta,paths_list,lims)
    evt_set = Set{Int64}()
    sta_set = Set{Int64}()
    for i in eachindex(paths_list)
        evt = evtsta.evts[paths_list[i][1]]
        lat, lon, depth = evt.lat, evt.lon, evt.depth
        if (check_boundaries(lat, lims.lat) && (check_boundaries(lon, lims.lon)) && (check_boundaries(depth, lims.depth)))
            push!(evt_set,evt.id)
            for j in eachindex(paths_list)
                (paths_list[j][1] == evt.id) && push!(sta_set,paths_list[j][2])
            end
        end
    end
    return evt_set, sta_set
end

function correct_polarization(rnodes,ray,nray,evtsta)
    error("Does polarization need to be corrected?") # it seems it does, at least for SPECFEM data
    # Polarization angle is defined in ray-aligned QTL coordinates and should be constant along the ray!
    # To compute S-wave velocities and splitting intensities, we need the angle between the S-wave
    # polarization azimuth (in the QT-plane) and the projection of the hexagonal symmetry axis into the QT-plane.
    # The polarization azimuth only needs to be corrected if we are computing this angle in some other coordinate
    # system.
    # nodes = rnodes.r2n[nray][1]:rnodes.r2n[nray][end]
	# sas = Vector{Float64}()
	# a, f = 6371.0, 0.0
    # n1,n2 = nodes[1],nodes[end]
    # lat2, lon2 =  asin(rnodes.c[3,n1]/rnodes.r[1,n1]), atan(rnodes.c[2,n1],rnodes.c[1,n1])
	# for k in nodes
    #     lat1,lon1 = asin(rnodes.c[3,k]/rnodes.r[1,k]),atan(rnodes.c[2,k],rnodes.c[1,k])
	# 	dist, az, baz = Geodesics.inverse(((lon1,lat1,lon2,lat2))...,a,f)
	# 	if ((3pi/2 - az) > pi)
	# 		push!(sas,(3pi/2 - az) - pi)
	# 	else
	# 		push!(sas,(3pi/2 - az))
	# 	end
	# end
    # # BPV lazily removed polarization 'correction'. No adjustment for ray orientation.
	# for k = 1:length(nodes)
	# 	ray.ζ[k] = ray.ζ[k]  + sas[k] - ray.ϕ[k]
	# end
end

function rotation_matrices(φ,θ)
    latrot = -θ
    lonrot = φ
    Rx = zeros(Float64,3,3)
    Rx[1,1] = 1
    Rx[2,2] = cos(latrot)
    Rx[2,3] = -sin(latrot)
    Rx[3,2] = sin(latrot)
    Rx[3,3] = cos(latrot)
    Ry = zeros(Float64,3,3) 
    Ry[1,1] = cos(lonrot)
    Ry[1,3] = sin(lonrot)
    Ry[2,2] = 1
    Ry[3,1] = -sin(lonrot)
    Ry[3,3] = cos(lonrot)
    return Rx,Ry
end

function rotate_basis(v1,v2,v3,Rx,Ry)
    R = Ry*Rx
    for i in eachindex(v1)
        b = inv(R)*[v1[i],v2[i],v3[i]]
        v1[i], v2[i], v3[i] = b[1], b[2], b[3]
    end
end


