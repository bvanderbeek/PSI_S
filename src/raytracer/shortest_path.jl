struct ShortestPathConst
    previous::Vector{Int64}
    distance::Vector{Float64}
end

function tracing_vfield(Vp1D,Vs1D,gVp,gVs,vfields,active_fields)
    Vp = copy(Vp1D)
    Vs = copy(Vs1D)
    if haskey(active_fields,"Vp") 
        Vp .= vfields[active_fields["Vp"]]
    end 
    if haskey(active_fields,"dlnVp") 
        dlnVp = vfields[active_fields["dlnVp"]]
        @. Vp = Vp * (1.0 + dlnVp)
    end

    if haskey(active_fields,"Vs") 
        Vs .= vfields[active_fields["Vs"]]
    end 
    if haskey(active_fields,"dlnVs") 
        dlnVs = vfields[active_fields["dlnVs"]]
        @. Vs = Vs * (1.0 + dlnVs)    
    end
    if haskey(active_fields,"Vp2Vs") 
        Vp2Vs = vfields[active_fields["Vp2Vs"]]
        @. Vs = Vp / Vp2Vs
    end
    if haskey(active_fields,"η") 
        dlnVp = vfields[active_fields["dlnVp"]]
        η = vfields[active_fields["η"]]
        @. Vs = Vs * (1.0 + dlnVp*η)
    end
    @. gVp = gVp + Vp
    @. gVs = gVs + Vs
end

function raytracing_dijkstra(gr, observables, LocalRaysManager, paths, evtsta, IP; firstcall=true, it=0, chains=Vector{Vector{ModelConst}}, aniso_status=false)
    print("\nray-tracing is running...\n")
    print("\n",maximum(gr.Vp)," ",minimum(gr.Vp),"\n")
    print("\n",maximum(gr.Vs)," ",minimum(gr.Vs),"\n")
    relocation_status = false
    if IP.EQRLOC.relocate
        if (IP.EQRLOC.relocations_init) && (it==0)
            relocation_status = true
        elseif ((it % IP.EQRLOC.relocations_its) == 0) && (it >= IP.EQRLOC.relocations_pause)
            relocation_status = true
        end
    end

    length(LocalRaysManager.source2receivers) == length(LocalRaysManager.local_evts) ? rev = true : rev = false
    if IP.B4DI.status && !firstcall
        velocity_4D_field = instance_velocity_grid(observables,evtsta,LocalRaysManager,IP,chains)
        node2node = grids_map(gr,velocity_4D_field)
    end
    DShPa = Vector{Vector{ShortestPathConst}}()
    if relocation_status
        Ns = maximum(LocalRaysManager.source2receivers)[1]
        for i in 1:Ns
            push!(DShPa,Vector{ShortestPathConst}(undef,2))
        end
    end

   # lowmen = [zeros(Int64,length(gr.x)) for i in 1:Threads.nthreads()]      # -- lower interval nodes collector

    sr = collect((LocalRaysManager.source2receivers))
    Threads.@threads for npair in eachindex(sr)
        t_id = Threads.threadid()
        (source,receivers) = sr[npair]
        if relocation_status
            D = DShPa[source]
        else
            D = Vector{ShortestPathConst}(undef,2)
        end
        if IP.B4DI.status && !firstcall
            grid = copy_grid(gr)
            interpolate_4Dfield(grid, velocity_4D_field, evtsta, source, node2node)
        else
            grid = gr
        end
        source_node = LocalRaysManager.source_nodes[source]
        phases = Set{Int64}()
        for receiver in receivers
            (haskey(LocalRaysManager.pair2ray,[source,receiver,1])) && push!(phases,1)
            (haskey(LocalRaysManager.pair2ray,[source,receiver,2])) && push!(phases,2)
        end
        visited = zeros(Bool,length(grid.x)) # -- visited nodes
        mask_topography(grid,visited)
        phase_visited = zeros(Bool,length(grid.x))
        if LocalRaysManager.carving
            rec_nodes = [LocalRaysManager.receiv_nodes[receiver] for receiver in receivers]
            carve_grid(visited,grid,source_node,rec_nodes)
        end
        for phase in phases 
            phase_visited .= visited
            # D[phase] = dijkstra_interval(grid,source_node,phase,phase_visited,lowmen[t_id];aniso_status=aniso_status)
            D[phase] = dijkstra(grid,source_node,phase,phase_visited;aniso_status=aniso_status)
        end
        if !relocation_status
            get_path(D,grid,source,receivers,LocalRaysManager,paths,rev,evtsta)
        end
    end 

    if relocation_status
        if IP.B4DI.status
            print("\nrelocations not allowed in 4D imaging...\n")
        else
            EQrelocation(gr,DShPa,IP,LocalRaysManager,observables,evtsta;it=it)
            for npair in eachindex(collect((LocalRaysManager.source2receivers)))
                (source,receivers) = collect((LocalRaysManager.source2receivers))[npair]
                get_path(DShPa[source],gr,source,receivers,LocalRaysManager,paths,rev,evtsta)
            end
            # test_reloc_lonlat(gr,evtsta,LocalRaysManager,paths,it)
            # test_reloc_londepth(gr,evtsta,LocalRaysManager,paths,it)
        end
    end
    
end

function get_path(D,gr,source,receivers,LocalRaysManager,paths,rev,evtsta)
    source_node = LocalRaysManager.source_nodes[source]
    phases = Set{Int64}()
    for receiver in receivers
        (haskey(LocalRaysManager.pair2ray,[source,receiver,1])) && push!(phases,1)
        (haskey(LocalRaysManager.pair2ray,[source,receiver,2])) && push!(phases,2)
    end
    for receiver in receivers
        for phase in phases
            if haskey(LocalRaysManager.pair2ray,[source,receiver,phase])
                rayid = LocalRaysManager.pair2ray[[source,receiver,phase]]
                receiv_node = LocalRaysManager.receiv_nodes[receiver]
                # print("\n",source," ",receiver," ",)
                p = shortest_path(D[phase],source_node,receiv_node; rev = rev)
                if rev
                    ievt, ista = source, receiver
                else
                    ievt, ista = receiver, source
                end
                (phase == 1) ? ph = "P" : ph = "S"
                paths[rayid] = pathConst(ievt, ista, ph, gr.θ[p], gr.φ[p], gr.r[p], [0.0], [D[phase].distance[receiv_node]])
            end
        end
    end
end

function shortest_path(D,source,receiver; rev = true)
    prev = D.previous
    path = Int[receiver]
    ipath = prev[receiver]
    # print(D.distance[receiver],"\n")
    while ipath ∉ path
        (ipath == 0) && break
        push!(path, ipath)
        ipath = prev[ipath]
    end
    rev ? (return reverse(path)) : (return path)
end

function raytracing(rays,LocalRaysManager,MarkovChains,observables,evtsta,IP;it=0)
    # -- here at the moment, but dreaming of IP parameters
    aniso_status = false

    grid = instance_grid(observables, evtsta, LocalRaysManager, IP; aniso_status=aniso_status)
    !IP.B4DI.status && fill_grid(grid,MarkovChains,evtsta,observables,IP;aniso_status=aniso_status) # -- if 4D is not active -> one velocity field evaluation for all the events

    refm = readdlm(IP.velocitymodel,' ',Float64,'\n')   # -- reads the reference 1D velocity model

    paths = Array{pathConst,1}(undef,length(rays))
    if IP.B4DI.status
        raytracing_dijkstra(grid, observables, LocalRaysManager, paths, evtsta, IP; firstcall=false, chains = MarkovChains)
    else
        raytracing_dijkstra(grid, observables, LocalRaysManager, paths, evtsta, IP; firstcall=false, it=it, aniso_status=aniso_status)
    end

    rnodes = build_rays(paths,rays,observables,evtsta,LocalRaysManager,IP,refm)

    voronoi_slims = []
    voronoi_fields = MarkovChains[begin][end].fields
    for fieldid in eachindex(voronoi_fields)
        push!(voronoi_slims,[])
        [push!(voronoi_slims[fieldid],zeros(Float64,2)) for i in 1:3]
        field = voronoi_fields[fieldid]
        voronoi_slims[fieldid][1][1], voronoi_slims[fieldid][1][2] = field.slims[1][1], field.slims[1][2]
        voronoi_slims[fieldid][2][1], voronoi_slims[fieldid][2][2] = field.slims[2][1], field.slims[2][2]
        voronoi_slims[fieldid][3][1], voronoi_slims[fieldid][3][2] = field.slims[3][1], field.slims[3][2]
    end
    voronoi_domain_pierce(rays,rnodes,voronoi_slims,IP)
    return rnodes
end

function fill_grid(grid, MarkovChains, evtsta, observables, IP; aniso_status = false)

    grid.σ .= 0.0
    [grid.estatics[iobs] .= 0.0 for iobs in eachindex(grid.estatics)]
    [grid.sstatics[iobs] .= 0.0 for iobs in eachindex(grid.sstatics)]

    nsamples = 0
    points = permutedims(hcat(grid.x, grid.y, grid.z))
    points_radial = permutedims(hcat(grid.r))
    fieldslist = MarkovChains[begin][end].fieldslist
    nchains = length(MarkovChains)
    tracing_fields = ["Vp","dlnVp","Vs","dlnVs","Vp2Vs","η"]
    active_fields = Dict{String,Int64}()
    vfields = Vector{Vector{Float64}}()
    for fieldname in fieldslist
        (fieldname[1] ∉ tracing_fields) && continue    
        push!(vfields,zeros(Float64,length(grid.x)))
    end
    Vp1D, Vs1D = copy(grid.Vp), copy(grid.Vs)
    grid.Vp .= 0.0
    grid.Vs .= 0.0
    chain_models = Int64(round(IP.RayTracingInit.sub_its / IP.MCS.saveint))
    print("\n",chain_models,"\n")
    for chain in eachindex(MarkovChains)
        MarkovChain = MarkovChains[chain]
        for model in MarkovChain[end:end]#[end-chain_models+1:end]
            [(vfields[i] .= 0.0) for i in eachindex(vfields)]
            fieldid = 0
            if model.T[1] != 1
                continue
            end
            nsamples += 1
            for iobs in eachindex(observables.Obs)
                obs = observables.Obs[iobs]
                grid.σ[iobs] += model.dataspace.Obs[iobs].noise[1]
                for i in eachindex(model.dataspace.Obs[iobs].estatics)
                    ievt = obs.evtids[i]
                    grid.estatics[iobs][ievt] += model.dataspace.Obs[iobs].estatics[i]
                end
                for i in eachindex(model.dataspace.Obs[iobs].sstatics)
                    ista = obs.staids[i]
                    grid.sstatics[iobs][ista] += model.dataspace.Obs[iobs].sstatics[i]
                end
            end
            for i in eachindex(model.fields)
                voronoi = model.fields[i]
                fieldname = voronoi.fieldname
                (fieldname ∉ tracing_fields) && continue
                fieldid += 1
                if check_1dim(voronoi)
                    inds = NN_interpolation(points_radial, voronoi.r[:,begin:voronoi.n[1]])
                else
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                end
                for j in eachindex(grid.x)
                    r = grid.r[j]
                    if IP.MCS.squeezing
                        if voronoi.slims[3][1] <= r <= voronoi.slims[3][2]
                            vfields[fieldid][j] = voronoi.v[inds[j]]
                        else
                            vfields[fieldid][j] = voronoi.ref_value
                        end
                    else
                        vfields[fieldid][j] = voronoi.v[inds[j]]
                    end
                end
                active_fields[fieldname] = fieldid
            end
            tracing_vfield(Vp1D, Vs1D, grid.Vp, grid.Vs, vfields, active_fields)
        end
    end
    @. grid.Vp = grid.Vp / nsamples
    @. grid.Vs = grid.Vs / nsamples
    # -- anisotropic ray-tracing --
    if aniso_status
        aniso_fields(grid,MarkovChains,IP,IP.MCS)
    end
    # -----------------------------
    @. grid.σ = grid.σ / nsamples
    for iobs in eachindex(observables.Obs)
        @. grid.estatics[iobs] = grid.estatics[iobs] / nsamples
        @. grid.sstatics[iobs] = grid.sstatics[iobs] / nsamples
    end
    print("\nnoise: ",grid.σ)
end

function fill_4Dgrid(grid, MarkovChains, IP)
    nsamples = 0
    points = permutedims(hcat(grid.x, grid.y, grid.z, ones(length(grid.x))*0.0))
    points_radial = permutedims(hcat(grid.r, ones(length(grid.x))*0.0))
    fieldslist = MarkovChains[begin][end].fieldslist
    nchains = length(MarkovChains)
    tracing_fields = ["Vp","dlnVp","Vs","dlnVs","Vp2Vs","η"]
    
    for frame in eachindex(grid.tp)
        timestep = grid.tp[frame]
        points[4,:] .= timestep 
        points_radial[2,:] .= timestep 
        active_fields = Dict{String,Int64}()
        vfields = Vector{Vector{Float64}}()
        for fieldname in fieldslist
            (fieldname[1] ∉ tracing_fields) && continue    
            push!(vfields,zeros(Float64,length(grid.x)))
        end
        Vp1D, Vs1D = copy(grid.Vp[frame]), copy(grid.Vs[frame])
        grid.Vp[frame] .= 0.0
        grid.Vs[frame] .= 0.0
        for chain in eachindex(MarkovChains)
            [(vfields[i] .= 0.0) for i in eachindex(vfields)]
            MarkovChain = MarkovChains[chain]
            fieldid = 0
            if MarkovChain[end].T[1] != 1
                continue
            end
            nsamples += 1
            for i in eachindex(MarkovChain[end].fields)
                voronoi = MarkovChain[end].fields[i]
                fieldname = voronoi.fieldname
                (fieldname ∉ tracing_fields) && continue
                fieldid += 1
                if check_1dim(voronoi)
                    inds = NN_interpolation(points_radial, voronoi.r[:,begin:voronoi.n[1]])
                else
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                end
                for j in eachindex(grid.x)
                    x, y, z = grid.x[j], grid.y[j], grid.z[j]
                    r = sqrt(x^2+y^2+z^2)
                    θ = asin(z/r)
                    φ = atan(y,x)
                    ((r < grid.rp[begin]) && (r = grid.rp[begin]))
                    ((r > grid.rp[end]) && (r = grid.rp[end]))
                    if IP.MCS.squeezing
                        if voronoi.slims[3][1] <= r <= voronoi.slims[3][2]
                            vfields[fieldid][j] = voronoi.v[inds[j]]
                        else
                            vfields[fieldid][j] = voronoi.ref_value
                        end
                    else
                        vfields[fieldid][j] = voronoi.v[inds[j]]
                    end
                end
                active_fields[fieldname] = fieldid
            end
            tracing_vfield(Vp1D, Vs1D, grid.Vp[frame], grid.Vs[frame], vfields, active_fields)
        end
        @. grid.Vp[frame] = grid.Vp[frame] / nsamples
        @. grid.Vs[frame] = grid.Vs[frame] / nsamples
    end
end

function interpolate_4Dfield(grid, velocity_4D_field, evtsta, source, node2node)
    T0 = evtsta.evts[source].T0
    frame = v_dist_ind(T0,velocity_4D_field.tp)
    for ind in eachindex(grid.Vp)
        vel_ind = node2node[ind]
        grid.Vp[ind] = velocity_4D_field.Vp[frame][vel_ind]
        grid.Vs[ind] = velocity_4D_field.Vs[frame][vel_ind]
    end
end

function aniso_fields(grid,MarkovChains,IP,MCS)
    print("\ncalculating average anisotropic fields...\n")
    nfields = MarkovChains[begin][end].nfields
    fieldsname = MarkovChains[begin][end].fieldslist
    nchains = length(MarkovChains)
    points = permutedims(hcat(grid.x, grid.y, grid.z))
    chain_models = Int64(round(IP.RayTracingInit.sub_its / MCS.saveint))
    nsamples = 0
    ε_map = zeros(Float64,length(grid.x))
    δ_map = zeros(Float64,length(grid.x))
    𝛙_map = zeros(Float64,length(grid.x))
    v3_map = zeros(Float64,length(grid.x))
    n1_map = zeros(Float64,length(grid.x))
    n2_map = zeros(Float64,length(grid.x))
    n3_map = zeros(Float64,length(grid.x))
    db_map = zeros(Float64,length(grid.x))
    v_ecef = [zeros(Float64,3) for i in eachindex(grid.x)]
    T = Vector{Matrix{Float64}}()
    [push!(T,zeros(Float64,3,3)) for i in eachindex(grid.x)]
    for chain in eachindex(MarkovChains)
        MarkovChain = MarkovChains[chain]
        for model in MarkovChain[end:end]#[end:-5:end-chain_models+1]
            nsamples += 1
            for i in eachindex(model.fields)
                voronoi = model.fields[i]
                fieldname = voronoi.fieldname
                if fieldname == "ε"
                    ε_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [ε_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                    δ_map .= ε_map
                elseif fieldname == "δ"
                    δ_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [δ_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                elseif fieldname == "𝛙"
                    𝛙_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [𝛙_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                elseif fieldname == "v3"
                    v3_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [v3_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                end
            end
            [updateT(T[i],6.0,0.0,ε_map[i],δ_map[i],𝛙_map[i],v3_map[i]) for i in eachindex(grid.x)]
        end
    end
    slims = MarkovChains[begin][begin].fields[1].slims
    midlat, midlon = mean([slims[1][1],slims[1][2]]), mean([slims[2][1],slims[2][2]])
    Threads.@threads for i in eachindex(T)
        v, b = zeros(Float64,3), zeros(Float64,3)
        @. T[i] = T[i] / nsamples
        Ti = T[i]
        vecs = eigvecs(Ti)
        vals = abs.(eigvals(Ti))
        inds = sortperm(vals)
        ind_a, ind_b, ind_c = reverse(inds)
        λ1, λ2, λ3 = vals[ind_a], vals[ind_b], vals[ind_c]
        crit = sqrt(0.5)*(sqrt((λ1-λ2)^2+(λ2-λ3)^2+(λ3-λ1)^2))/(sqrt((λ1)^2+(λ2)^2+(λ3)^2))
        db_map[i] = crit
        v_ecef[i] .= [(vals[ind_a])*vecs[1,ind_a],(vals[ind_a])*vecs[2,ind_a],-(vals[ind_a])*vecs[3,ind_a]]
        b[1], b[2], b[3] = ecef_vector(v_ecef[i], rad2deg(midlat), rad2deg(midlon))
        #b[1], b[2], b[3] = v[1], v[2], v[3]
	    nb = sqrt(b[1]^2+b[2]^2+b[3]^2)
        grid.n1[i] = b[1]/nb
        grid.n2[i] = b[2]/nb
        grid.n3[i] = b[3]/nb

        grid.ε[i] = sqrt(b[1]^2+b[2]^2+b[3]^2)
        grid.δ[i] = grid.ε[i]
        # p1 = 0.6742
        # p2 = -0.8169
        # q1 = 0.04419
        # grid.δ[i] = grid.ε[i]*(p1*grid.ε[i] + p2)/(grid.ε[i] + q1)
        if crit < 0.75
            grid.ε[i] = 0.0
            grid.δ[i] = 0.0
        end
        v_ecef[i] .= v_ecef[i]/norm(v_ecef[i])
        grid.vecef1[i] = v_ecef[i][1]
        grid.vecef2[i] = v_ecef[i][2]
        grid.vecef3[i] = v_ecef[i][3]
    end

    # -- recompute Thomsen parameters
    grid.ε .= 0.0
    grid.δ .= 0.0
    nsamples = 0
    for chain in eachindex(MarkovChains)
        MarkovChain = MarkovChains[chain]
        for model in MarkovChain[end:end]#[end:-5:end-chain_models+1]
            nsamples += 1
            for i in eachindex(model.fields)
                voronoi = model.fields[i]
                fieldname = voronoi.fieldname
                if fieldname == "ε"
                    ε_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [ε_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                elseif fieldname == "δ"
                    δ_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [δ_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                elseif fieldname == "𝛙"
                    𝛙_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [𝛙_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                elseif fieldname == "v3"
                    v3_map .= 0.0
                    inds = NN_interpolation(points, voronoi.c[:,begin:voronoi.n[1]])
                    [v3_map[i] = voronoi.v[inds[i]] for i in eachindex(inds)]
                end
            end
            for i in eachindex(grid.ε)
                if db_map[i] > 0.75
                    γ = asin(v3_map[i])
                    v1, v2, v3 = cos(𝛙_map[i])*cos(γ), sin(𝛙_map[i])*cos(γ), sin(γ)
                    grid.ε[i] += ε_map[i] * abs(v1*v_ecef[i][1]+v2*v_ecef[i][2]+v3*v_ecef[i][3])
                end
            end
        end
    end

    @. grid.ε = grid.ε / nsamples
    #@. grid.δ = grid.δ / nsamples
    # -- elliptical anisotropy
    #@. grid.δ = grid.ε
    # -- low-aspect ratio (100) cracks anisotropy
    # p1 = 0.6742
    # p2 = -0.8169
    # q1 = 0.04419
    @. grid.δ = grid.ε#*(p1*grid.ε + p2)/(grid.ε + q1)

    return true

end

function updateT(T,v_1D_P,dlnVp,ε,δ,azi,v3)

    if (δ != ε) && (δ/(δ-ε) > 0) && (sqrt(0.5*(δ/(δ-ε))) <= 1)
        θt = asin(sqrt(0.5*(δ/(δ-ε))))
    else
        θt = 0.0
    end
    sinx2, cosx2 = sin(θt)^2, cos(θt)^2
    sinx4 = sinx2^2
    q16_15, q4_15 = (16.0/15.0), (4.0/15.0) 
    αiso = v_1D_P*(1.0 + dlnVp) 
    α = αiso/sqrt(1.0 + q16_15*ε + q4_15*δ) 
    α1 = α
    α2 = α*sqrt(1+2.0*ε)
    α3 = α*sqrt(1.0 + 2.0*δ*sinx2*cosx2 + 2.0*ε*sinx4)
    αmin, αmax = min(α1,α2,α3), max(α1,α2,α3)

    f = (αmax-αmin)/(αmax+αmin)

    γ = asin(v3)
    𝛙 = azi

    i = pi/2 - γ
    ϕ = 3*pi/2 - 𝛙

    T[1,1] += f*sin(i)^2*sin(ϕ)^2
    T[1,2] += f*sin(i)^2*sin(ϕ)*cos(ϕ)
    T[1,3] += f*sin(i)*cos(i)*sin(ϕ)

    T[2,1] += f*sin(i)^2*sin(ϕ)*cos(ϕ)
    T[2,2] += f*sin(i)^2*cos(ϕ)^2
    T[2,3] += f*sin(i)*cos(i)*cos(ϕ)

    T[3,1] += f*sin(i)*cos(i)*sin(ϕ)
    T[3,2] += f*sin(i)*cos(i)*cos(ϕ)
    T[3,3] += f*cos(i)^2
end

function ecef_vector(east_north_elv, latitude, longitude)
    # ECEF components for vector at (0°N, 0°E)
    w = (east_north_elv[3], east_north_elv[1], east_north_elv[2])
    # Rotate to geographic position
    c = π/180.0
    R = rotation_matrix((-c*latitude, c*longitude), (2, 3))
    sx = R[1,1]*w[1] + R[1,2]*w[2] + R[1,3]*w[3]
    sy = R[2,1]*w[1] + R[2,2]*w[2] + R[2,3]*w[3]
    sz = R[3,1]*w[1] + R[3,2]*w[2] + R[3,3]*w[3]

    return sx, sy, sz
end

function rotation_matrix(α, n)
    sinα, cosα = sincos(α)
    if n == 1
        return [1.0 0.0 0.0; 0.0 cosα -sinα; 0.0 sinα cosα]
    elseif n == 2
        return [cosα 0.0 sinα; 0.0 1.0 0.0; -sinα 0.0 cosα]
    elseif n == 3
        return [cosα -sinα 0.0; sinα cosα 0.0; 0.0 0.0 1.0]
    else
        error("Requested rotation axis out-of-bounds!")
    end
end

function rotation_matrix(α::Tuple, n::Tuple)
    R = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(n)
        R = rotation_matrix(α[i],n[i])*R
    end
    return R
end

# using Plots
# function buildV(gr)
#     Vp = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     Vs = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     for n in eachindex(gr.Vp)
#         i, j, k = CartesianIndex(gr,n)
#         Vp[i,j,k] = gr.Vp[n]
#         Vs[i,j,k] = gr.Vs[n]
#     end
#     latg = rad2deg.(collect(range(gr.θ[1],gr.θ[end],length(Vp[:,1,1]))))
#     long = rad2deg.(collect(range(gr.φ[1],gr.φ[end],length(Vp[1,:,1]))))
#     dg = collect(range(gr.r[1],6371,length(Vp[1,1,:]))) .- R
#     return latg,long,dg,Vp,Vs
# end

# function test_reloc_lonlat(gr,evtsta,LocalRaysManager,paths,it)
#     latg,long,dg,Vp,Vs = buildV(gr)

#     evtfile = readdlm("input/evt.dat")
#     stafile = readdlm("input/sta.dat")
#     truelat, truelon, truedepth = Float64[], Float64[], Float64[]
#     loclat, loclon, locdepth = Float64[], Float64[], Float64[]
#     for i in axes(evtfile,1)
#         push!(truelat,evtfile[i,2])
#         push!(truelon,evtfile[i,3])
#         push!(truedepth,-evtfile[i,4])
#         push!(loclon,rad2deg(evtsta.evts[i].lon))
#         push!(locdepth,(evtsta.evts[i].depth))
#         push!(loclat,rad2deg(evtsta.evts[i].lat))
#     end

#     p1 = heatmap(latg,long,Vp[:,:,1],clim=(6,8),c=cgrad(:rainbow, rev=true),aspect_ratio=:equal,dpi=600)
#     # scatter!(stafile[:,3],stafile[:,2],color=:black,markershape=:utriangle,markersize=2,legend=false)
#     for npair in eachindex(collect((LocalRaysManager.source2receivers)))
#         (source,receivers) = collect((LocalRaysManager.source2receivers))[npair]
#         for receiver in receivers
#             if haskey(LocalRaysManager.pair2ray,[source,receiver,1])
#                 rayid = LocalRaysManager.pair2ray[[source,receiver,1]]
#                 θ, φ, r = paths[rayid].lat, paths[rayid].lon, paths[rayid].rad
#                 if source == 1
#                     # plot!(rad2deg.(φ),r .- R,color=:black,alpha=0.3)
#                     plot!(rad2deg.(φ),rad2deg.(θ),color=:black,alpha=0.3)
#                 end
#             end
#         end
#     end
#     scatter!(truelon,truelat,color=:red,markershape=:star8,legend=false)
#     scatter!(loclon,loclat,color=:green,markershape=:star8,legend=false)
#     loc_dists = Float64[]
#     for i in eachindex(loclon)
#         plot!([truelon[i],loclon[i]],[truelat[i],loclat[i]],color=:white,alpha=0.4)
#         lat_eq,lon_eq,r_eq = deg2rad(truelat[i]), deg2rad(truelon[i]), R+truedepth[i]
#         x1,y1,z1 = @cartesian(lat_eq,lon_eq,r_eq)
#         lat_eq,lon_eq,r_eq = deg2rad(loclat[i]), deg2rad(loclon[i]), R+locdepth[i]
#         x2,y2,z2 = @cartesian(lat_eq,lon_eq,r_eq)
#         push!(loc_dists,sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2))
#     end
#     Vp2Vs = @. Vp / Vs
#     p2 = heatmap(latg,long,Vp2Vs[:,:,1],clim=(1.6,2.2),c=cgrad(:roma, rev=true),aspect_ratio=:equal,dpi=600)
#     plot(p1,p2)
#     name = "output/synth_PS/figs/reloc/relocations_$it.png"
#     Plots.savefig(name)

#     histogram(loc_dists,bins=30,color=:grey)
#     Plots.savefig("output/synth_PS/figs/reloc/discr_$it.png")
# end

# function test_reloc_londepth(gr,evtsta,LocalRaysManager,paths,it)
#     latg,long,dg,Vp,Vs = buildV(gr)

#     evtfile = readdlm("input/evt.dat")
#     truelat, truelon, truedepth = Float64[], Float64[], Float64[]
#     loclat, loclon, locdepth = Float64[], Float64[], Float64[]
#     for i in axes(evtfile,1)
#         push!(truelat,evtfile[i,2])
#         push!(truelon,evtfile[i,3])
#         push!(truedepth,-evtfile[i,4])
#         push!(loclon,rad2deg(evtsta.evts[i].lon))
#         push!(locdepth,(evtsta.evts[i].depth))
#         push!(loclat,rad2deg(evtsta.evts[i].lat))
#     end
    
#     # -- convert degrees to distances
#     truelon_dist = @. R * deg2rad(truelon)
#     loclon_dist = @. R * deg2rad(loclon)

#     dist, depth = R * deg2rad.(long), dg .- R
#     p1 = heatmap(dist,dg,Vp[1,:,:]',clim=(4,6),c=cgrad(:rainbow, rev=true),aspect_ratio=:equal,dpi=1200)
#     for npair in eachindex(collect((LocalRaysManager.source2receivers)))
#         (source,receivers) = collect((LocalRaysManager.source2receivers))[npair]
#         for receiver in receivers
#             if haskey(LocalRaysManager.pair2ray,[source,receiver,1])
#                 rayid = LocalRaysManager.pair2ray[[source,receiver,1]]
#                 θ, φ, r = paths[rayid].lat, paths[rayid].lon, paths[rayid].rad
#                 if source == 1
#                     # plot!(rad2deg.(φ),r .- R,color=:black,alpha=0.3)
#                     lonray, depthray = R * φ, r .- R
#                     plot!(lonray,depthray,color=:black,alpha=0.3)
#                 end
#             end
#         end
#     end
#     scatter!(truelon_dist,truedepth,color=:red,markershape=:star8,legend=false)
#     scatter!(loclon_dist,locdepth,color=:green,markershape=:star8,legend=false)
#     loc_dists = Float64[]
#     for i in eachindex(loclon)
#         plot!([truelon_dist[i],loclon_dist[i]],[truedepth[i],locdepth[i]],color=:white,alpha=0.4)
#         lat_eq,lon_eq,r_eq = deg2rad(truelat[i]), deg2rad(truelon[i]), R+truedepth[i]
#         x1,y1,z1 = @cartesian(lat_eq,lon_eq,r_eq)
#         lat_eq,lon_eq,r_eq = deg2rad(loclat[i]), deg2rad(loclon[i]), R+locdepth[i]
#         x2,y2,z2 = @cartesian(lat_eq,lon_eq,r_eq)
#         push!(loc_dists,sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2))
#     end
#     # p2 = heatmap(dist,dg,Vs[1,:,:]',clim=(2,3.25),c=cgrad(:rainbow, rev=true),aspect_ratio=:equal,dpi=1200)
#     # plot(p1,p2,layout=(2,1))
#     name = "output/vert_PS/figs/reloc/relocations_$it.png"
#     Plots.savefig(name)

#     histogram(loc_dists,bins=30,color=:grey)
#     Plots.savefig("output/vert_PS/figs/reloc/discr_$it.png")

# end

