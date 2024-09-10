function bfm(gr,source,phase)
    if phase == 1
        V = gr.Vp
    elseif phase == 2
        V = gr.Vs
    end

    fw_influence = -gr.fw_level:gr.fw_level
    n = length(gr.x)
    dist = zeros(Float64,length(gr.x))  # -- distances from source
    for i in eachindex(gr.x)
        dist[i] = Inf
    end
    dist[source] = 0.0
    dist0 = deepcopy(dist)
    p = zeros(Int64,length(gr.x))
    Q = falses(n)
    init_Q!(Q,source,gr,fw_influence)

    it = 1

    @inbounds while sum(Q) != 0
        # relax edges (parallel process)
        relax!(dist, p, dist0, Q, gr, V, fw_influence)

        # pop queue (serial-but-fast process)
        # fillfalse!(Q)
        fill!(Q, false)

        # update nodal queue (parallel process)
        update_Q!(Q, dist, dist0, gr, fw_influence)

        # update old distance vector (TODO parallel version)
        copyto!(dist0, dist)

        # update iteration counter
        it += 1
    end

    println("Converged in $it iterations")

    return ShortestPathConst(p, dist)
end

function relax!(dist, p, dist0, Q, gr, V, fw_influence)

    # iterate over queue. Unfortunately @threads can't iterate 
    # over a Set, so we need to collect() it. This yields an 
    # allocation, but it's worth it in this case as it saves 
    # a decent number of empty iterations and removes a layer
    # of branching
    Threads.@threads for i in findall(Q)
        _relax!(p, dist, dist0, i, gr, V, fw_influence)
    end
end

function _relax!(p, dist, dist0, nw, gr, V, fw_influence)
    redundancyQ = Set{Int32}()
    i, j, k = CartesianIndex(gr,nw)
    vmain = V[nw]
    xw, yw, zw = gr.x[nw], gr.y[nw], gr.z[nw]
    di = dist0[nw]
    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    ik = 0
    Threads.@threads for k3 in fw_influence
        ik += 1
        k3 += k
        (k3 < 1 || k3 > nz) && continue 
        ij = 0
        @inbounds for k2 in fw_influence
            ij += 1
            k2 += j
            (k2 < 1 || k2 > ny) && continue 
            ii = 0
            @inbounds for k1 in fw_influence
                ii += 1
                k1 += i
                (k1 < 1 || k1 > nx) && continue 
                nn = LinearIndex(gr, k1, k2, k3)
                # !(nn ∉ redundancyQ) && continue
                # push!(redundancyQ, nn)
                dGi = dist0[nn]
                tempDist = dGi + 2*sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)/(vmain+V[nn])
                if tempDist < di
                    di = tempDist
                    p[nw] = nn 
                end
            end
        end
    end
    dist[nw] = di
end

function init_Q!(Q,nw,gr,fw_influence)
    Q[nw] = true
    i, j, k = CartesianIndex(gr,nw)
    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    ik = 0
    Threads.@threads for k3 in fw_influence
        ik += 1
        k3 += k
        (k3 < 1 || k3 > nz) && continue 
        ij = 0
        @inbounds for k2 in fw_influence
            ij += 1
            k2 += j
            (k2 < 1 || k2 > ny) && continue 
            ii = 0
            @inbounds for k1 in fw_influence
                ii += 1
                k1 += i
                (k1 < 1 || k1 > nx) && continue 
                nn = LinearIndex(gr, k1, k2, k3)
                Q[nn] = true
            end
        end
    end
end

function update_Q!(Q, dist, dist0, gr, fw_influence)
    # update queue: if new distance is smaller 
    # than the previous one, add to adjecent 
    # to the queue nodes
    Threads.@threads for i in eachindex(Q)
        @inbounds if (dist[i] < Inf) && (dist[i] < dist0[i])
            init_Q!(Q,i,gr,fw_influence)
        end
    end
end