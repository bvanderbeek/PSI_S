# function dijkstra_interval(gr,source,phase,visited,lowmen;aniso_status=false)

#     if phase == 1
#         V = gr.Vp
#     elseif phase == 2
#         V = gr.Vs
#     end
#     if aniso_status
#         ε = gr.ε
#         δ = gr.δ
#         n1 = gr.n1
#         n2 = gr.n2
#         n3 = gr.n3
#     end
#     nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
#     fw_influence = -gr.fw_level:gr.fw_level

#     lowmen .= 0      # -- lower interval nodes collector
#     distance = zeros(Float64,length(gr.x))  # -- distances from source
#     previous = zeros(Int64,length(gr.x)) 
#     for id in eachindex(gr.x)
#         distance[id] = Inf
#     end
#     distance[source] = 0.0
#     dtmin = Inf
#     biarray = [Inf,Inf]
#     for i in eachindex(gr.x)
#         biarray[1] = 1/V[i]
#         biarray[2] = dtmin
#         dtmin = minimum(biarray)  
#     end
#     dtmin = gr.dmin*dtmin
#     maxI = 0.0

#     minQ, maxQ = 1, length(gr.x)

#     while minQ <= maxQ
#         maxI = maxI + dtmin
#         while ((minQ <= maxQ) && visited[minQ])
#             minQ = minQ + 1
#         end
#         while ((maxQ >= minQ) && visited[maxQ])
#             maxQ = maxQ - 1
#         end

#         numQ = 0
#         for indQ in minQ:maxQ
#           if (!visited[indQ]) && (distance[indQ] < maxI)
#              numQ = numQ + 1
#              lowmen[numQ] = indQ
#              visited[indQ] = true
#           end
#         end

#         for indx in 1:numQ
#             nw = lowmen[indx]
#             xw, yw, zw = gr.x[nw], gr.y[nw], gr.z[nw]
#             dmain = distance[nw]
#             vmain = V[nw]
#             if aniso_status
#                 εmain = ε[nw]
#                 δmain = δ[nw]
#                 n1main = n1[nw]
#                 n2main = n2[nw]
#                 n3main = n3[nw]
#             end
#             i, j, k  = CartesianIndex(gr,nw)
#             ik = 0
#             @inbounds for k3 in fw_influence    # -- starts loop over the forward-star
#                 ik += 1
#                 k3 += k
#                 (k3 < 1 || k3 > nz) && continue 
#                 ij = 0
#                 @inbounds for k2 in fw_influence
#                     ij += 1
#                     k2 += j
#                     (k2 < 1 || k2 > ny) && continue 
#                     ii = 0
#                     @inbounds for k1 in fw_influence
#                         ii += 1
#                         k1 += i
#                         (k1 < 1 || k1 > nx) && continue 
#                         nn = LinearIndex(gr, k1, k2, k3)
#                         (visited[nn]) && continue
#                         if aniso_status
#                             dist = sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)
#                             cosx2 = (((xw-gr.x[nn])*n1main + (yw-gr.y[nn])*n2main + (zw-gr.z[nn])*n3main)/dist)^2
#                             sinx2 = 1 - cosx2
#                             sinx4 = sinx2^2
#                             α_iso = vmain
#                             α = α_iso/sqrt(1.0 + (16/15)*εmain + (4/15)*δmain)
#                             v1 = α*sqrt(1.0 + 2.0*δmain*sinx2*cosx2 + 2.0*εmain*sinx4)

#                             cosx2 = (((xw-gr.x[nn])*n1[nn] + (yw-gr.y[nn])*n2[nn] + (zw-gr.z[nn])*n3[nn])/dist)^2
#                             sinx2 = 1 - cosx2
#                             sinx4 = sinx2^2
#                             α_iso = V[nn]
#                             α = α_iso/sqrt(1.0 + (16/15)*ε[nn] + (4/15)*δ[nn])
#                             v2 = α*sqrt(1.0 + 2.0*δ[nn]*sinx2*cosx2 + 2.0*ε[nn]*sinx4)
#                             tempDist = dmain + 2*dist/(v1+v2)
#                         else
#                             tempDist = dmain + 2*sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)/(vmain+V[nn])
#                         end
#                         if tempDist < distance[nn]
#                            distance[nn] = tempDist
#                            previous[nn] = nw 
#                         end
#                     end
#                 end
#             end
#         end
#     end

#     return ShortestPathConst((previous), (distance))

# end


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
        if aniso_status
            εmain = ε[main]
            δmain = δ[main]
            n1main = n1[main]
            n2main = n2[main]
            n3main = n3[main]
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
                        dist = sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)
                        cosx2 = (((xw-gr.x[nn])*n1main + (yw-gr.y[nn])*n2main + (zw-gr.z[nn])*n3main)/dist)^2
                        sinx2 = 1 - cosx2
                        sinx4 = sinx2^2
                        α_iso = vmain
                        α = α_iso/sqrt(1.0 + (16/15)*εmain + (4/15)*δmain)
                        v1 = α*sqrt(1.0 + 2.0*δmain*sinx2*cosx2 + 2.0*εmain*sinx4)

                        cosx2 = (((xw-gr.x[nn])*n1[nn] + (yw-gr.y[nn])*n2[nn] + (zw-gr.z[nn])*n3[nn])/dist)^2
                        sinx2 = 1 - cosx2
                        sinx4 = sinx2^2
                        α_iso = V[nn]
                        α = α_iso/sqrt(1.0 + (16/15)*ε[nn] + (4/15)*δ[nn])
                        v2 = α*sqrt(1.0 + 2.0*δ[nn]*sinx2*cosx2 + 2.0*ε[nn]*sinx4)
                        u_mean = 0.5*(1/v1+1/v2)
                        tempDist = dmain + dist*u_mean
                    else
                        u_mean = 0.5*(1/vmain+1/V[nn])
                        tempDist = dmain + sqrt((xw-gr.x[nn])^2+(yw-gr.y[nn])^2+(zw-gr.z[nn])^2)*u_mean
                    end                    
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