function iso_rt_dist(gr,V,distance,previous,Q,main,nn)
    u_mean = 0.5*(1/V[main]+1/V[nn])
    tempDist = distance[main] + sqrt((gr.x[main]-gr.x[nn])^2+(gr.y[main]-gr.y[nn])^2+(gr.z[main]-gr.z[nn])^2)*u_mean
    update_dist(tempDist,distance,previous,Q,main,nn)
end

function aniso_rt_dist(gr,V,n1,n2,n3,ε,δ,distance,previous,Q,main,nn)
    dist = sqrt((gr.x[main]-gr.x[nn])^2+(gr.y[main]-gr.y[nn])^2+(gr.z[main]-gr.z[nn])^2)
    cosx2 = (((gr.x[main]-gr.x[nn])*n1[main] + (gr.y[main]-gr.y[nn])*n2[main] + (gr.z[main]-gr.z[nn])*n3[main])/dist)^2
    sinx2 = 1 - cosx2
    sinx4 = sinx2^2
    α_iso = V[main]
    α = α_iso/sqrt(1.0 + (16/15)*ε[main] + (4/15)*δ[main])
    v1 = α*sqrt(1.0 + 2.0*δ[main]*sinx2*cosx2 + 2.0*ε[main]*sinx4)

    cosx2 = (((gr.x[main]-gr.x[nn])*n1[nn] + (gr.y[main]-gr.y[nn])*n2[nn] + (gr.z[main]-gr.z[nn])*n3[nn])/dist)^2
    sinx2 = 1 - cosx2
    sinx4 = sinx2^2
    α_iso = V[nn]
    α = α_iso/sqrt(1.0 + (16/15)*ε[nn] + (4/15)*δ[nn])
    v2 = α*sqrt(1.0 + 2.0*δ[nn]*sinx2*cosx2 + 2.0*ε[nn]*sinx4)
    u_mean = 0.5*(1/v1+1/v2)
    tempDist = distance[main] + dist*u_mean
    update_dist(tempDist,distance,previous,Q,main,nn)
end

function update_dist(tempDist,distance,previous,Q,main,nn)
    if tempDist < distance[nn]
        distance[nn] = tempDist
        previous[nn] = main 
        Q[nn] = tempDist
    end
end

function dijkstra(gr,source,phase,visited;aniso_status=false)

    if phase == 1
        V = gr.Vp
    elseif phase == 2
        V = gr.Vs
    end
    if aniso_status
        ε = gr.ε
        δ = gr.δ
        n1 = gr.n1
        n2 = gr.n2
        n3 = gr.n3
    end

    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    fw_influence = -gr.fw_level:gr.fw_level

    Q = PriorityQueue{Int64,Float64}() # -- queue
    distance = zeros(Float64,length(gr.x))
    previous = zeros(Int64,length(gr.x))
    for i in eachindex(gr.x)
        distance[i] = Inf
    end
    Q[source] = 0.0
    distance[source] = 0.0

    reg_nodes = gr.nnodes[1]*gr.nnodes[2]*gr.nnodes[3]
    while !isempty(Q)
        main = dequeue!(Q)
        visited[main] = true 
        if main > reg_nodes
            i, j, k = CartesianIndex(gr,gr.sr2n[main])
        else      
            i, j, k = CartesianIndex(gr,main)
        end
        for nn in gr.n2sr[main]
            if aniso_status
                aniso_rt_dist(gr,V,n1,n2,n3,ε,δ,distance,previous,Q,main,nn)
            else
                iso_rt_dist(gr,V,distance,previous,Q,main,nn)
            end
        end
        ik = 0
        @inbounds for k3 in fw_influence
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
                    (visited[nn]) && continue
                    if aniso_status
                        aniso_rt_dist(gr,V,n1,n2,n3,ε,δ,distance,previous,Q,main,nn)
                    else
                        iso_rt_dist(gr,V,distance,previous,Q,main,nn)
                    end                    
                end
            end
        end
    end
    return ShortestPathConst(previous, distance)
end
