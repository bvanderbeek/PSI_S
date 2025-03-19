# -- NON-LINEAR MAXIMUM LIKELIHOOD EARTHQUAKE RELOCATOR
# using PyPlot
function EQrelocation(gr,D,IP,LocalRaysManager,observables,evtsta;it=0)

    debug = false
    time_shift = Vector{Float64}()
    eq_residuals = Vector{Vector{Float64}}()
    for i in eachindex(evtsta.evts)
        push!(time_shift,0.0)
        push!(eq_residuals,Vector{Float64}()) 
    end

    nx, ny, nz = gr.nnodes[1], gr.nnodes[2], gr.nnodes[3]
    dn = IP.EQRLOC.fw_reloc + 1
    fw_influence = -dn:dn

    evt_rms = Dict{Int64,Float64}()
    evt_azmgap = Dict{Int64,Float64}()
    errormap = zeros(Float64,length(gr.x))
    misfitmap = zeros(Float64,length(gr.x))
    staticmap = zeros(Float64,length(gr.x))
    for inevt in eachindex(collect(LocalRaysManager.local_evts))
        active_stations = Set{Int64}()
        if debug 
            all_st = Vector{Float64}()
            all_res = Vector{Vector{Float64}}()
            all_double = Vector{Dict{Int64,Vector{Float64}}}()
            all_p = Vector{Vector{Float64}}()
        end
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
                    # if gr.r[nn] > R
                    #     continue
                    # end
                    residuals = Float64[]
                    noise_res = Float64[]
                    double_phase_tts = Dict{Int64,Vector{Float64}}()
                    double_phase_noise = Dict{Int64,Vector{Float64}}()
                    for receiver in LocalRaysManager.local_stats
                        for ph in [1,2]
                            if haskey(LocalRaysManager.pair2ray,[receiver,nevt,ph])
                                push!(active_stations,receiver)
                                nray = LocalRaysManager.pair2ray[[receiver,nevt,ph]]
                                for iobs in eachindex(observables.Obs)
                                    obs = observables.Obs[iobs]
                                    if haskey(obs.ray2obs,nray)
                                        if (obs.obsname == "local_traveltimes_P") || (obs.obsname == "local_traveltimes_S")
                                            σ = gr.σ[iobs]
                                        #     σ = gr.σP[1]
                                        # elseif Obs.obsname == "local_traveltimes_S"
                                        #     σ = gr.σS[1]
                                        else
                                            continue
                                        end
                                        push!(residuals,obs.obsval[obs.ray2obs[nray]]-D[receiver][ph].distance[nn])
                                        # -- add statics
                                        residuals[end] += (-gr.estatics[iobs][nevt]-gr.sstatics[iobs][receiver])
                                        # --------------
                                        # print("\n",ph," ",Obs.obsval[Obs.ray2obs[nray]]," ",D[receiver][ph].distance[nn])
                                        push!(noise_res,σ)
                                        push!(rays,nray)
                                        if (ph == 1) && (obs.obsname == "local_traveltimes_P")
                                            double_phase_tts[receiver] = [obs.obsval[obs.ray2obs[nray]]]
                                            double_phase_noise[receiver] = [σ]
                                        elseif (ph == 2) && (obs.obsname == "local_traveltimes_S")
                                            if haskey(double_phase_tts,receiver)
                                                push!(double_phase_tts[receiver],obs.obsval[obs.ray2obs[nray]])
                                                push!(double_phase_noise[receiver],σ)
                                            end
                                        end
                                    end
                                end
                            else 
                                continue
                            end
                        end
                    end

                    if it == 0
                        if IP.EQRLOC.wadati_init
                            mr, p_wd = wadati_OT(nevt,double_phase_tts)
                        else
                            mr, p_wd = weighted_average(residuals,noise_res), [0,0]
                        end
                    else
                        mr, p_wd = weighted_average(residuals,noise_res), [0,0]
                    end
                    dm_residuals = residuals .- mr
                    errormap[nn] = sqrt(mean(@. dm_residuals^2))
                    misfitmap[nn] = sum(@. (dm_residuals/noise_res)^2)
                    # errormap[nn] = sqrt(mean(@. abs(dm_residuals)))
                    # misfitmap[nn] = sum(@. (abs(dm_residuals)/noise_res)^1)
                    staticmap[nn] = mr
                    if debug 
                        push!(all_st,mr)
                        push!(all_res,residuals)
                        push!(all_double,double_phase_tts)
                        push!(all_p,p_wd)
                    end
                end
            end
        end
        min,ind = findmin(misfitmap)
        evt_rms[nevt] = (errormap[ind])
        T0 = evtsta.evts[nevt].T0
        evtsta.evts[nevt] = evt(nevt,gr.θ[ind],gr.φ[ind],gr.r[ind]-R,T0+staticmap[ind])
        LocalRaysManager.receiv_nodes[nevt] = ind
        for nray in rays
            for Obs in observables.Obs
                if (Obs.obsname == "local_traveltimes_P" || Obs.obsname == "local_traveltimes_S") && haskey(Obs.ray2obs,nray)
                    Obs.obsval[Obs.ray2obs[nray]] = Obs.obsval[Obs.ray2obs[nray]] - staticmap[ind]
                end
            end
        end   

        # # -- azimuthal gap
        # azm_gap = calculate_azimuthal_gap(evtsta,nevt,active_stations)
        # evt_azmgap[nevt] = azm_gap

        # # -- calculates hypocenter's uncertainty
        # Ntts = length(all_res[ind])
        # p = Vector{Vector{Float64}}()
        # max_dT = 0
        # for i in eachindex(misfitmap)
        #     if misfitmap[i] <= (misfitmap[ind] + Ntts)
        #         # push!(p,[gr.x[i]-gr.x[ind],gr.y[i]-gr.y[ind],gr.z[i]-gr.z[ind]])
        #         # push!(p,[gr.y[i]-gr.y[ind],gr.z[i]-gr.z[ind]]) # -- lon-lat test
        #         push!(p,[gr.x[i]-gr.x[ind],gr.y[i]-gr.y[ind]]) # -- lon-depth test
        #         max_dT = maximum([max_dT,abs(staticmap[i]-staticmap[ind])])
        #     end
        # end
        # max_dT = round(max_dT;digits=2)
        # hull = convex_hull(p)
        # # G = zeros(Float64,length(hull),6)
        # G = zeros(Float64,length(hull),3) # -- 2D tests
        # for (i,c) in enumerate(hull)
        #     # x, y, z = c[1], c[2], c[3]
        #     # G[i,1], G[i,2], G[i,3], G[i,4], G[i,5], G[i,6] = x^2, y^2, z^2, 2*x*y, 2*x*z, 2*y*z
        #     # -- 2D test (x-y <-> y-z)
        #     x, y = c[1], c[2]
        #     G[i,1], G[i,2], G[i,3] = x^2, y^2, 2*x*y
        # end        

        # if det(G'G) != 0
        #     p = inv(G'G)*G'*ones(length(axes(G,1)))
        #     # -- ellipsoid matrix
        #     # A = [
        #     #         p[1] p[4] p[5]
        #     #         p[4] p[2] p[6]
        #     #         p[5] p[6] p[3]
        #     # ]

        #     # -- 2D test (ellipse matrix)
        #     A = [
        #             p[1] p[3]
        #             p[3] p[2]
        #     ]

        #     λ = abs.(eigvals(A))
        #     axs = @. 1/sqrt(λ)
        #     vecs = eigvecs(A)
        # else
        #     axs = ones(3)*0.2
        #     vecs = [
        #         1 0 0
        #         0 1 0
        #         0 0 1
        #     ]
        # end

        # if debug
        #     time_shift[nevt] = all_st[ind]
        #     eq_residuals[nevt] = all_res[ind]
        #     # -- plot error metrics
        #     # plot_eqerr_lonlat(gr,errormap,misfitmap,ind,nevt,Ntts,axs,vecs,max_dT,azm_gap,time_shift,eq_residuals,it,active_stations)
        #     # plot_eqerr_londepth(gr,errormap,misfitmap,ind,nevt,Ntts,axs,vecs,max_dT,azm_gap,time_shift,eq_residuals,it,active_stations,all_double[ind],all_p[ind])
        # end

    end

end

function weighted_average(residuals,noise_res)
    w,µ = 0.0, 0.0
    for i in eachindex(residuals)
        w += 1/(noise_res[i])^2
        µ += residuals[i]/(noise_res[i])^2
    end
    return µ/w
end

function wadati_OT(nevt,double_phase_tts)
    double_phase_tts = [PS_tts for (i,PS_tts) in double_phase_tts]
    G = zeros(Float64,length(double_phase_tts),2)
    d = zeros(Float64,length(double_phase_tts))
    for i in axes(G,1)
        G[i,1], G[i,2] = double_phase_tts[i][1], 1
        d[i] = double_phase_tts[i][2]-double_phase_tts[i][1]
    end
    p = inv(G'G)*G'd 
    t0 = -p[2]/p[1]

    return t0, p
end

function calculate_azimuthal_gap(evtsta,nevt,active_stations)
    elat, elon = rad2deg(evtsta.evts[nevt].lat), rad2deg(evtsta.evts[nevt].lon)
    azm = Float64[]
    for ista in active_stations
        slat, slon = rad2deg(evtsta.stas[ista].lat), rad2deg(evtsta.stas[ista].lon)
        Δ, α = inverse_geodesic(elat, elon, slat, slon)
        if α < 0
            α += 360
        end
        push!(azm,α)
    end
    sort!(azm)
    azm_gap = diff(azm)
    push!(azm_gap,360-(azm[end]-azm[begin]))
    return maximum(azm_gap)
end

# using Plots
# function plot_eqerr_lonlat(gr,errormap,misfitmap,ind,nevt,Ntts,λ,vecs,max_dT,azm_gap,time_shift,eq_residuals,it,active_stations)
#     # (nevt != 6) && return
#     latg = rad2deg.(collect(range(gr.θ[1],gr.θ[end],gr.nnodes[1])))
#     long = rad2deg.(collect(range(gr.φ[1],gr.φ[end],gr.nnodes[2])))
#     dg = collect(range(gr.r[1],6371,gr.nnodes[3])) .- R

#     rms = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     misfit = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     for n in eachindex(gr.Vp)
#         i, j, k = CartesianIndex(gr,n)
#         rms[i,j,k] = errormap[n]*1000
#         misfit[i,j,k] = misfitmap[n]
#     end

#     min_rms = round(errormap[ind];digits=3)*1000
#     azm_gap = round(azm_gap;digits=0)
#     p1 = heatmap(long,latg,rms[:,:,1],c=cgrad(:magma, rev=false),clim=(0,1000),title="RMS: $(min_rms) [ms], azm_gap: $(azm_gap)°",titlefont=font(8,"Computer Modern"),aspect_ratio=:equal,dpi=600)
#     evtfile = readdlm("input/evt.dat")
#     stafile = readdlm("input/sta.dat")
#     evtlat_true, evtlon_true, evtdepth_true = evtfile[nevt,2], evtfile[nevt,3], -evtfile[nevt,4]
#     evtlat, evtlon, evtdepth = rad2deg(gr.θ[ind]), rad2deg(gr.φ[ind]), gr.r[ind]-R
#     scatter!([evtlon_true],[evtlat_true],color=:white,markersize=3,legend=false)
#     scatter!([evtlon],[evtlat],color=:red,markersize=3,legend=false)
#     stat_ids = collect(active_stations)
#     scatter!(stafile[stat_ids,3],stafile[stat_ids,2],color=:black,markershape=:utriangle,markersize=2,legend=false)
#     max_dT = max_dT*1000
#     p2 = heatmap(long,latg,log10.(misfit[:,:,1]),c=cgrad(:magma, rev=false),clim=(1,3),aspect_ratio=:equal,title="RMS: σT: $(max_dT) [ms]",titlefont=font(8,"Computer Modern"),dpi=600)
#     var = log10(misfitmap[ind] + Ntts)
#     p2 = contour!(long,latg,log10.(misfit[:,:,1]),color=:white,levels=[var])
#     scatter!([evtlon_true],[evtlat_true],color=:white,markersize=3,legend=false)
#     scatter!([evtlon],[evtlat],color=:red,markersize=3,legend=false)

#     v1, v2 = vecs[:,1], vecs[:,2] 
#     λ1, λ2 = λ[1], λ[2]
#     t = range(0,2π,length=100)

#     ellipse_points_y = @. λ1 * cos(t)*v1[1] + λ2 * sin(t)*v2[1]
#     ellipse_points_z = @. λ1 * cos(t)*v1[2] + λ2 * sin(t)*v2[2]

#     y = @. ellipse_points_y + gr.y[ind]
#     z = @. ellipse_points_z + gr.z[ind]

#     ϕ = @. rad2deg(atan(y,6371*cos(y/6371)))
#     θ = @. rad2deg(asin(z/6371.0))

#     plot!(ϕ,θ,color=:red,legend=false)
#     # plot!([gr.y[ind],gr.y[ind]+λ[1]*vecs[1,1]],[gr.z[ind],gr.z[ind]+λ[1]*vecs[2,1]],color=:red,legend=false)
#     # plot!([gr.y[ind],gr.y[ind]+λ[2]*vecs[1,2]],[gr.z[ind],gr.z[ind]+λ[2]*vecs[2,2]],color=:red,legend=false)

#     plot(p1,p2,layout = (1, 2))

#     Plots.savefig("output/synth_PS/figs/reloc/hypocenters/hypo_$(nevt)_error.png")

#     histogram(time_shift,bins=30,color=:grey)
#     Plots.savefig("output/synth_PS/figs/reloc/time_shift_$(it).png")
#     for i in eachindex(eq_residuals)
#         histogram(eq_residuals[i],bins=30,color=:grey)
#         Plots.savefig("output/synth_PS/figs/reloc/eq_residuals/residuals_$(i).png")
#     end
# end

# function plot_eqerr_londepth(gr,errormap,misfitmap,ind,nevt,Ntts,λ,vecs,max_dT,azm_gap,time_shift,eq_residuals,it,active_stations,double_phase_tts,p_wd)
#     # (nevt != 6) && return
#     latg = rad2deg.(collect(range(gr.θ[1],gr.θ[end],gr.nnodes[1])))
#     long = rad2deg.(collect(range(gr.φ[1],gr.φ[end],gr.nnodes[2])))
#     dg = collect(range(gr.r[1],6371,gr.nnodes[3])) .- R

#     rms = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     misfit = zeros(Float64,gr.nnodes[1],gr.nnodes[2],gr.nnodes[3])
#     for n in eachindex(gr.Vp)
#         i, j, k = CartesianIndex(gr,n)
#         rms[i,j,k] = errormap[n]*1000
#         misfit[i,j,k] = misfitmap[n]
#     end

#     dist, depth = R * deg2rad.(long), dg .- R
#     min_rms = round(errormap[ind];digits=3)*1000
#     azm_gap = round(azm_gap;digits=0)
#     p1 = heatmap(dist,dg,rms[1,:,:]',c=cgrad(:magma, rev=false),clim=(0,1000),title="RMS: $(min_rms) [ms], azm_gap: $(azm_gap)°",titlefont=font(8,"Computer Modern"),aspect_ratio=:equal,dpi=600)
#     evtfile = readdlm("input/evt.dat")
#     evtlat, evtlon, evtdepth_true = evtfile[nevt,2], evtfile[nevt,3], -evtfile[nevt,4]
#     evtlon_dist_true = R * deg2rad(evtlon)
#     evtlat, evtlon, evtdepth = gr.θ[ind],gr.φ[ind],gr.r[ind]-R
#     evtlon_dist = R * evtlon
#     scatter!([evtlon_dist_true],[evtdepth_true],color=:white,markersize=3,legend=false)
#     scatter!([evtlon_dist],[gr.r[ind]-R],color=:red,markersize=3,legend=false)
#     max_dT = max_dT*1000
#     p2 = heatmap(dist,dg,log10.(misfit[1,:,:])',c=cgrad(:magma, rev=false),clim=(1,3),aspect_ratio=:equal,title="RMS: σT: $(max_dT) [ms]",titlefont=font(8,"Computer Modern"),dpi=600)
#     var = log10(misfitmap[ind] + Ntts)
#     p2 = contour!(dist,dg,log10.(misfit[1,:,:])',color=:white,levels=[var])
#     scatter!([evtlon_dist_true],[evtdepth_true],color=:white,markersize=3,legend=false)
#     scatter!([evtlon_dist],[gr.r[ind]-R],color=:red,markersize=3,legend=false)

#     v1, v2 = vecs[:,1], vecs[:,2] 
#     λ1, λ2 = λ[1], λ[2]
#     t = range(0,2π,length=100)

#     ellipse_points_x = @. λ1 * cos(t)*v1[1] + λ2 * sin(t)*v2[1]
#     ellipse_points_y = @. λ1 * cos(t)*v1[2] + λ2 * sin(t)*v2[2]

#     x = @. ellipse_points_x + (gr.x[ind]-R)
#     y = @. ellipse_points_y + gr.y[ind]

#     plot!(y,x,color=:red,legend=false)
#     plot!([gr.y[ind],gr.y[ind]+λ[1]*vecs[2,1]],[gr.x[ind]-R,gr.x[ind]+λ[1]*vecs[1,1]-R],color=:red,legend=false)
#     plot!([gr.y[ind],gr.y[ind]+λ[2]*vecs[2,2]],[gr.x[ind]-R,gr.x[ind]+λ[2]*vecs[1,2]-R],color=:red,legend=false)

#     plot(p1,p2,layout = (2, 1))

#     Plots.savefig("output/vert_PS/figs/reloc/hypocenters/hypo_$(nevt)_error.png")

#     histogram(time_shift,bins=30,color=:grey)
#     Plots.savefig("output/vert_PS/figs/reloc/time_shift_$(it).png")
#     for i in eachindex(eq_residuals)
#         histogram(eq_residuals[i],bins=30,color=:grey)
#         Plots.savefig("output/vert_PS/figs/reloc/eq_residuals/residuals_$(i).png")
#     end

#     P_tts = [double_phase_tts[i][1] for i in eachindex(double_phase_tts)]
#     S_tts = [double_phase_tts[i][2] for i in eachindex(double_phase_tts)]
#     t0 = round(-p_wd[2]/p_wd[1];digits=1)
#     scatter(P_tts,S_tts-P_tts,color=:black,legend=false,title="Vp/Vs: $(round(p_wd[1]+1;digits=2)), t0: $(t0) [s]",titlefont=font(8,"Computer Modern"))
#     x = [minimum(P_tts),maximum(P_tts)]
#     y = @. p_wd[1]*x + p_wd[2]
#     plot!(x,y,color=:red,legend=false)
#     Plots.savefig("output/vert_PS/figs/reloc/wadati/eq_$(nevt).png")
# end
