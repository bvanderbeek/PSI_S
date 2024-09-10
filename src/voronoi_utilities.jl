function neighbours(voronoi,ind,update) # -- compute Voronoi neighbours via Delaunay triangulation
    if check_1dim(voronoi)
        if voronoi.n[1] == 1
            return
        end
        coords = @view voronoi.r[:,1:voronoi.n[1]]
        inds = 1:voronoi.n[1]
        sorter = vcat(coords,inds')
        sorter = sorter[:,sortperm(sorter[1, :])]
        inds = Int64.(sorter[2,:])
        i = 0
        for j in eachindex(inds)
            if inds[j] == ind
                i = j
            end
        end
        if i == 1
            push!(update.influence.neighbours,inds[i+1])
        elseif i == voronoi.n[1]
            push!(update.influence.neighbours,inds[i-1])
        else
            push!(update.influence.neighbours,inds[i+1])
            push!(update.influence.neighbours,inds[i-1])
        end
    else
        del = MiniQhull.delaunay(@view voronoi.c[:,1:voronoi.n[1]])
        for tri in axes(del,2)
            for node1 in axes(del,1)
                if del[node1,tri] == ind
                    for node2 in axes(del,1)
                        if del[node2,tri] != ind
                            push!(update.influence.neighbours,del[node2,tri])
                        end
                    end
                end
            end
        end
    end
end

function map_voro2space2(voronoi,rnodes)
    [push!(voronoi.nodes2nuclei,0) for i in eachindex(rnodes.n2r)]
    v_dist_grid2(rnodes,voronoi)
end

function map_voro2space(voronoi,vnox)
    [push!(voronoi.nodes2nuclei,0) for i in axes(vnox,2)]
    v_dist_grid(vnox,voronoi)
end

function NN_interpolation(points, nuclei)
    kdtree = KDTree(nuclei)
    inds = nn(kdtree, points)[1]
    return inds
end


function old_influence(voronoi,nodes2rays,rays2nodes,update)
    neighbours(voronoi,update.ind[1],update)
    nray = 0
    for nray in voronoi.nuclei2rays[update.ind[1]]
        update.influence.rays[nray] = true
        for i in rays2nodes[1,nray]:rays2nodes[2,nray]
            if voronoi.nodes2nuclei[i] == update.ind[1]
                update.influence.nodes[i] = true
            end
        end
    end
    for i in (update.influence.neighbours)
        update.influence.nuclei[i] = true
    end
    update.influence.nuclei[update.ind[1]] = true
    push!(update.influence.new_neighborhood,update.ind[1])
end

function new_influence(voronoi,nodes2rays,rays2nodes,update)
    neighbours(voronoi,update.ind[1],update)
    nray = 0
    for nray in voronoi.nuclei2rays[update.ind[1]]
        update.influence.rays[nray] = true
        for i in rays2nodes[1,nray]:rays2nodes[2,nray]
            if voronoi.nodes2nuclei[i] == update.ind[1]
                update.influence.nodes[i] = true
            end
        end
    end
    for j in update.influence.neighbours
        for nray in voronoi.nuclei2rays[j]
            update.influence.rays[nray] = true
            for i in rays2nodes[1,nray]:rays2nodes[2,nray]
                if voronoi.nodes2nuclei[i] == j
                    update.influence.nodes[i] = true
                end
            end
        end
    end

    for i in (update.influence.neighbours)
        update.influence.nuclei[i] = true
        push!(update.influence.new_neighborhood,i)
    end
    update.influence.nuclei[update.ind[1]] = true
    push!(update.influence.new_neighborhood,update.ind[1])
end

function save_n2r(voronoi,update)
    for i in eachindex(update.influence.nuclei)
        if update.influence.nuclei[i]
            push!(update.old_voronoi.nuclei2rays,Set{Int64}())
            [push!(update.old_voronoi.nuclei2rays[end],j) for j in voronoi.nuclei2rays[i]]
        end
    end
end

function swap(voronoi,ind,last)
    for i in eachindex(voronoi.nodes2nuclei)
        if voronoi.nodes2nuclei[i] == last
            voronoi.nodes2nuclei[i] = ind
        elseif voronoi.nodes2nuclei[i] == ind
            voronoi.nodes2nuclei[i] = last
        end
    end
    for i in axes(voronoi.c,1)
        x1, x2 = voronoi.c[i,ind], voronoi.c[i,last]
        voronoi.c[i,last], voronoi.c[i,ind] = x1, x2
    end
    voronoi.v[ind], voronoi.v[last] = voronoi.v[last], voronoi.v[ind]
    voronoi.r[1,ind], voronoi.r[1,last] = voronoi.r[1,last], voronoi.r[1,ind]
    if !isempty(voronoi.nuclei2rays)
        voronoi.nuclei2rays[ind], voronoi.nuclei2rays[last] = voronoi.nuclei2rays[last], voronoi.nuclei2rays[ind]
    end
end

function interpolate_diagram(voronoi,vnox,nodes2rays,update)
    for i in update.influence.new_neighborhood
        [delete!(voronoi.nuclei2rays[i],j) for j in voronoi.nuclei2rays[i]]
    end
    v_dist_grid(vnox,voronoi,nodes2rays;infl=true,inodes=update.influence.nodes,inuclei=update.influence.nuclei)
end

function v_dist(r, voronoi)
    min_ind = zero(Int64)
    min_dist = Inf
	for i in eachindex(voronoi.v)[1:voronoi.n[1]]
		distance = (r - voronoi.r[1,i])^2 
        if distance < min_dist
           min_dist = distance
           min_ind = i
        end
    end
    return min_ind
end

function v_dist(x, y, z, voronoi)
    min_ind = zero(Int64)
    min_dist = Inf
	for i in eachindex(voronoi.v)[1:voronoi.n[1]]
		distance = (x - voronoi.c[1,i])^2 + (y - voronoi.c[2,i])^2 + (z - voronoi.c[3,i])^2
        if distance < min_dist
           min_dist = distance
           min_ind = i
        end
    end
    return min_ind
end

function v_dist(x, y, z, t, voronoi)
    min_ind = zero(Int64)
    min_dist = Inf
	for i in eachindex(voronoi.v)[1:voronoi.n[1]]
		distance = (x - voronoi.c[1,i])^2 + (y - voronoi.c[2,i])^2 + (z - voronoi.c[3,i])^2 + (t - voronoi.c[4,i])^2
        if distance < min_dist
           min_dist = distance
           min_ind = i
        end
    end
    return min_ind
end

function v_dist_grid(vnox, voronoi, nodes2rays; infl=false, inodes=[0], inuclei=[0])
    if check_1dim(voronoi)
        coords1 = @view vnox[5:5,:]
        coords2 = voronoi.r
    else
        if length(axes(voronoi.c,1)) == 4
            coords1 = @view vnox[1:4,:]
        else
            coords1 = @view vnox[1:3,:]
        end
        coords2 = voronoi.c
    end
    nray = 0
    nvoro = 0
    for i in axes(coords1,2)
        if infl && !inodes[i]
            continue
        end
        min_dist = Inf
        for j in axes(coords2,2)[begin:voronoi.n[1]]
            if infl && !inuclei[j]
                continue
            end
            distance = 0.0
            for k in axes(coords2,1)
                distance += (coords1[k,i] - coords2[k,j])^2
            end
            if distance < min_dist
                min_dist = distance
                voronoi.nodes2nuclei[i] = j
            end
        end
        if (nray != nodes2rays[i]) || (nvoro != voronoi.nodes2nuclei[i])
            push!(voronoi.nuclei2rays[voronoi.nodes2nuclei[i]],nodes2rays[i])
            nray = nodes2rays[i]
            nvoro = voronoi.nodes2nuclei[i]
        end
    end
end

function v_dist_grid2(rnodes, voronoi; infl=false, inodes=[0], inuclei=[0])
    if check_1dim(voronoi)
        coords1, coords2 = rnodes.r, voronoi.r
    else
        if length(axes(voronoi.c,1)) == 4
            coords1 = @view rnodes.c[1:4,:]
        else
            coords1 = @view rnodes.c[1:3,:]
        end
        coords2 = voronoi.c
    end
    for i in axes(coords1,2)
        if infl && !inodes[i]
            continue
        end
        min_dist = Inf
        for j in axes(coords2,2)[begin:voronoi.n[1]]
            if infl && !inuclei[j]
                continue
            end
            distance = 0.0
            for k in axes(coords2,1)
                distance += (coords1[k,i] - coords2[k,j])^2
            end
            if distance < min_dist
                min_dist = distance
                voronoi.nodes2nuclei[i] = j
            end
        end
        push!(voronoi.nuclei2rays[voronoi.nodes2nuclei[i]],rnodes.n2r[i])
    end
end
