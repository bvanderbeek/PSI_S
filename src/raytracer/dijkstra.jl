function dijkstra_interval(gr,source,phase,visited)

    if phase == 1
        V = gr.Vp
    elseif phase == 2
        V = gr.Vs
    end

    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    fw_influence = -gr.fw_level:gr.fw_level

    lowmen = zeros(Int64,length(gr.x))      # -- lower interval nodes collector

    # visited = zeros(Bool,length(gr.x))      # -- visited nodes
    distance = zeros(Float64,length(gr.x))  # -- distances from source
    previous = zeros(Int64,length(gr.x))
    for id in eachindex(gr.x)
        distance[id] = Inf
    end
    distance[source] = 0.0
    dtmin = Inf
    biarray = [Inf,Inf]
    for i in eachindex(gr.x)
        biarray[1] = 1/V[i]
        biarray[2] = dtmin
        dtmin = minimum(biarray)  
    end
    dtmin = gr.dmin*dtmin
    maxI = 0.0

    minQ, maxQ = 1, length(gr.x)

    while minQ <= maxQ
        maxI = maxI + dtmin
        while ((minQ <= maxQ) && visited[minQ])
            minQ = minQ + 1
        end
        while ((maxQ >= minQ) && visited[maxQ])
            maxQ = maxQ - 1
        end

        numQ = 0
        for indQ in minQ:maxQ
          if (!visited[indQ]) && (distance[indQ] < maxI)
             numQ = numQ + 1
             lowmen[numQ] = indQ
             visited[indQ] = true
          end
        end

        for indx in 1:numQ
            nw = lowmen[indx]
            xw, yw, zw = gr.x[nw], gr.y[nw], gr.z[nw]
            dmain = distance[nw]
            vmain = V[nw]
            i, j, k  = CartesianIndex(gr,nw)
            ik = 0
            @inbounds for k3 in fw_influence    # -- starts loop over the forward-star
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
                        tempDist = dmain + 2*sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)/(vmain+V[nn])
                        if tempDist < distance[nn]
                           distance[nn] = tempDist
                           previous[nn] = nw 
                        end
                    end
                end
            end
        end
    end

    return ShortestPathConst(previous, distance)

end


function dijkstra(gr,source,phase,visited)

    if phase == 1
        V = gr.Vp
    elseif phase == 2
        V = gr.Vs
    end

    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    fw_influence = -gr.fw_level:gr.fw_level

    Q = PriorityQueue{Int64,Float64}() # -- queue
    distance = zeros(Float64,length(gr.x))
    previous = zeros(Int64,length(gr.x))
    for id in eachindex(gr.x)
        distance[id] = Inf
    end
    Q[source] = 0.0
    distance[source] = 0.0

    
    while !isempty(Q)
        main = dequeue!(Q)
        xw, yw, zw = gr.x[main], gr.y[main], gr.z[main]
        visited[main] = true        
        i, j, k = CartesianIndex(gr,main)
        dmain = distance[main]
        vmain = V[main]
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
                    tempDist = dmain + 2*sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)/(vmain+V[nn])
                    if tempDist < distance[nn]
                       distance[nn] = tempDist
                       previous[nn] = main 
                       Q[nn] = tempDist
                    end
                end
            end
        end
    end
    return ShortestPathConst(previous, distance)
end