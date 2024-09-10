# using PyPlot
struct ShortestPathConst
    previous::Vector{Int64}
    distance::Vector{Float64}
end

# function get_field_spm(fields,fieldsname,name)
#     return [fields[i][fieldsname[name]] for i in eachindex(fields)]
# end

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

function buildV(gr)
    V = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
    for n in eachindex(gr.Vp)
        i, j, k = CartesianIndex(gr,n)
        V[i,j,k] = gr.Vp[n]
    end
    V = V[1,:,:]
    lg = rad2deg.(collect(range(gr.φ[1],gr.φ[end],length(V[:,1]))))
    dg = collect(range(gr.r[1],6371,length(V[1,:]))) .- R
    return lg,dg,V
end

function raytracing_dijkstra(gr, observables, LocalRaysManager, paths, evtsta, IP; firstcall=true, it=0, chains=Vector{Vector{ModelConst}})
    print("\nray-tracing is running...\n")
    print("\n",maximum(gr.Vp)," ",minimum(gr.Vp),"\n")
    # lg,dg,V = buildV(gr)
    relocation_status = IP.EQRLOC.relocate && ((it % IP.EQRLOC.relocations_its) == 0)

    length(LocalRaysManager.source2receivers) == length(LocalRaysManager.local_evts) ? rev = true : rev = false
    if IP.B4DI.status && !firstcall
        velocity_4D_field = instance_velocity_grid(evtsta,LocalRaysManager,IP,chains)
        node2node = grids_map(gr,velocity_4D_field)
    end

    DShPa = Vector{Vector{ShortestPathConst}}()
    if relocation_status
        Ns = maximum(LocalRaysManager.source2receivers)[1]
        for i in 1:Ns
            push!(DShPa,Vector{ShortestPathConst}(undef,2))
        end
    end

    sr = collect((LocalRaysManager.source2receivers))
    Threads.@threads for npair in eachindex(sr)
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
        visited = zeros(Bool,length(gr.x)) # -- visited nodes
        phase_visited = zeros(Bool,length(gr.x))
        if LocalRaysManager.carving
            rec_nodes = [LocalRaysManager.receiv_nodes[receiver] for receiver in receivers]
            carve_grid(visited,grid,source_node,rec_nodes)
        end
        for phase in phases 
            phase_visited .= visited
            D[phase] = dijkstra_interval(grid,source_node,phase,phase_visited)
            #@time D[phase] = dijkstra(grid,source_node,phase,phase_visited)
            # @time D[phase] = bfm(grid,source_node,phase)
        end
        if !relocation_status
            get_path(D,grid,source,receivers,LocalRaysManager,paths,rev)
        end
    end 

    if relocation_status
        if IP.B4DI.status
            print("\nrelocations not allowed in 4D imaging...\n")
        else
            EQrelocation(gr,DShPa,IP,LocalRaysManager,observables,evtsta)
            for npair in eachindex(collect((LocalRaysManager.source2receivers)))
                (source,receivers) = collect((LocalRaysManager.source2receivers))[npair]
                get_path(DShPa[source],gr,source,receivers,LocalRaysManager,paths,rev)
            end
        end
    end

    # evtfile = readdlm("../../InputData/joint_1D_3D/evt.dat")
    # truelon, truedepth = Float64[], Float64[]
    # loclon, locdepth = Float64[], Float64[]
    # for i in axes(evtfile,1)
    #     push!(truelon,evtfile[i,3])
    #     push!(truedepth,-evtfile[i,4])
    #     push!(loclon,rad2deg(evtsta.evts[i].lon))
    #     push!(locdepth,(evtsta.evts[i].depth))
    # end

    # fig,ax = PyPlot.subplots()
    # levels = range(4,7,length=100)
    # ax.contourf(lg,dg,V',vmin=4,vmax=7,cmap="gist_rainbow",levels=levels,extend="both")

    # for npair in eachindex(collect((LocalRaysManager.source2receivers)))
    #     (source,receivers) = collect((LocalRaysManager.source2receivers))[npair]
    #     for receiver in receivers
    #         if haskey(LocalRaysManager.pair2ray,[source,receiver,1])
    #             rayid = LocalRaysManager.pair2ray[[source,receiver,1]]
    #             θ, φ, r = paths[rayid].lat, paths[rayid].lon, paths[rayid].rad
    #             if source == 1
    #                 ax.plot(rad2deg.(φ),r .- R,color="black",alpha=0.3)
    #             end
    #         end
    #     end
    # end
    # ax.scatter(truelon,truedepth,color="red",marker="*",s=50)
    # ax.scatter(loclon,locdepth,color="green",marker="*",s=50)
    # for i in eachindex(loclon)
    #     ax.plot([truelon[i],loclon[i]],[truedepth[i],locdepth[i]],color="white",alpha=0.4)
    # end
    # name = "relocations_$it.png"
    # PyPlot.savefig(name)
end

function get_path(D,gr,source,receivers,LocalRaysManager,paths,rev)
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
                p = shortest_path(D[phase],source_node,receiv_node; rev = rev)
                if rev
                    ievt, ista = source, receiver
                else
                    ievt, ista = receiver, source
                end
                (phase == 1) ? ph = "P" : ph = "S"
                paths[rayid] = pathConst(ievt, ista, ph, gr.θ[p], gr.φ[p], gr.r[p], [0.0])
            end
        end
    end
end

function shortest_path(D,source,receiver; rev = true)
    prev = D.previous
    path = Int[receiver]
    ipath = prev[receiver]
    while ipath ∉ path
        (ipath == 0) && break
        push!(path, ipath)
        ipath = prev[ipath]
    end
    rev ? (return reverse(path)) : (return path)
end

function raytracing(rays,LocalRaysManager,MarkovChains,observables,evtsta,IP;it=0)

    grid = instance_grid(evtsta, LocalRaysManager, IP)

    !IP.B4DI.status && fill_grid(grid,MarkovChains,evtsta,IP) # -- if 4D is not active -> one velocity field evaluation for all the events

    refm = readdlm(IP.velocitymodel,' ',Float64,'\n')   # -- reads the reference 1D velocity model

    paths = Array{pathConst,1}(undef,length(rays))
    if IP.B4DI.status
        raytracing_dijkstra(grid, observables, LocalRaysManager, paths, evtsta, IP; firstcall=false, chains = MarkovChains)
    else
        raytracing_dijkstra(grid, observables, LocalRaysManager, paths, evtsta, IP; firstcall=false, it=it)
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

function fill_grid(grid, MarkovChains, evtsta, IP)
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
        for model in MarkovChain[end-chain_models+1:end]
            [(vfields[i] .= 0.0) for i in eachindex(vfields)]
            fieldid = 0
            if model.T[1] != 1
                continue
            end
            nsamples += 1
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
end

function fill_4Dgrid(grid, MarkovChains)
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


