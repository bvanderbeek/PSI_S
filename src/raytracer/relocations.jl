# -- NON-LINEAR MAXIMUM LIKELIHOOD EARTHQUAKE RELOCATOR
# using PyPlot
function EQrelocation(gr,D,IP,LocalRaysManager,observables,evtsta)

    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    dn = Int64(round(5*gr.nnodes[3]/(gr.r[end]-gr.r[begin]))) + 1
    fw_influence = -dn:dn

    evt_rms = Dict{Int64,Float64}()
    errormap = zeros(Float64,length(gr.x))
    misfitmap = zeros(Float64,length(gr.x))
    staticmap = zeros(Float64,length(gr.x))
    total_residuals = Float64[]
    for inevt in eachindex(collect(LocalRaysManager.local_evts))
        nevt = collect(LocalRaysManager.local_evts)[inevt]
        ind = LocalRaysManager.receiv_nodes[nevt] 
        print("progress: ",round(inevt/length(collect(LocalRaysManager.local_evts))*100; digits=2),"%\r")
        i,j,k = CartesianIndex(gr,ind)
        errormap .= Inf
        staticmap .= Inf
        misfitmap .= Inf
        rays = Set{Int64}()
        ik = 0
        @inbounds for k3 in fw_influence
            ik += 1
            k3 += k
            (k3 < 1 || k3 > gr.nnodes[3]) && continue 
            ij = 0
            @inbounds for k2 in fw_influence
                ij += 1
                k2 += j
                (k2 < 1 || k2 > gr.nnodes[2]) && continue 
                ii = 0
                @inbounds for k1 in fw_influence
                    ii += 1
                    k1 += i
                    (k1 < 1 || k1 > gr.nnodes[1]) && continue 
                    nn = LinearIndex(gr, k1, k2, k3)
                    if gr.r[nn] > R
                        continue
                    end
                    residuals = Float64[]
                    noise_res = Float64[]
                    for receiver in LocalRaysManager.local_stats
                        for ph in [1,2]
                            if haskey(LocalRaysManager.pair2ray,[receiver,nevt,ph])
                                nray = LocalRaysManager.pair2ray[[receiver,nevt,ph]]
                                for Obs in observables.Obs
                                    if haskey(Obs.ray2obs,nray)
                                        if Obs.obsname == "local_traveltimes_P" || Obs.obsname == "local_traveltimes_S"
                                            push!(residuals,Obs.obsval[Obs.ray2obs[nray]]-D[receiver][ph].distance[nn])
                                            push!(noise_res,Obs.noise_guess)
                                            push!(rays,nray)
                                        end
                                    end
                                end
                            else 
                                continue
                            end
                        end
                    end
                    mr = mean(residuals)
                    dm_residuals = residuals .- mr
                    errormap[nn] = mean(@. dm_residuals^2)
                    misfitmap[nn] = sum(@. (dm_residuals/noise_res)^2)
                    staticmap[nn] = mr
                end
            end
        end
        min,ind = findmin(misfitmap)
        evt_rms[nevt] = errormap[ind]
        T0 = evtsta.evts[nevt].T0
        evtsta.evts[nevt] = evt(nevt,gr.θ[ind],gr.φ[ind],gr.r[ind]-R,T0)
        LocalRaysManager.receiv_nodes[nevt] = ind
        for nray in rays
            for Obs in observables.Obs
                if (Obs.obsname == "local_traveltimes_P" || Obs.obsname == "local_traveltimes_S") && haskey(Obs.ray2obs,nray)
                    Obs.obsval[Obs.ray2obs[nray]] = Obs.obsval[Obs.ray2obs[nray]] - staticmap[ind]
                end
            end
        end   
        for receiver in LocalRaysManager.local_stats
            for ph in [1,2]
                if haskey(LocalRaysManager.pair2ray,[receiver,nevt,ph])
                    nray = LocalRaysManager.pair2ray[[receiver,nevt,ph]]
                    for Obs in observables.Obs
                        if haskey(Obs.ray2obs,nray)
                            if Obs.obsname == "local_traveltimes_P" || Obs.obsname == "local_traveltimes_S"
                                push!(total_residuals,Obs.obsval[Obs.ray2obs[nray]]-D[receiver][ph].distance[ind])
                            end
                        end
                    end
                else 
                    continue
                end
            end
        end    
    end
    
    print("\n\nlocations_RMS: ",sqrt(mean(@. total_residuals^2)),"\n")    
    # save("evt_rms.jld","evt2rms",evt_rms)
end