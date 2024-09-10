struct GridConst
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    θ::Vector{Float64}
    φ::Vector{Float64}
    r::Vector{Float64}
    Vp::Vector{Float64}
    Vs::Vector{Float64}
    fw_level::Int64
    nnodes::Vector{Int64}
    nxny::Int64
    dmin::Float64
    perturbed::Bool
end

struct Velocity4DGridConst
    θp::Vector{Float64}
    φp::Vector{Float64}
    rp::Vector{Float64}
    tp::Vector{Float64}
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    r::Vector{Float64}
    Vp::Vector{Vector{Float64}}
    Vs::Vector{Vector{Float64}}
    nnodes::Vector{Int64}
    nxny::Int64
    upscale::Float64
end

function copy_grid(grid)
    return GridConst(
        grid.x,
        grid.y,
        grid.z,
        grid.θ,
        grid.φ,
        grid.r,
        copy(grid.Vp),
        copy(grid.Vs),
        grid.fw_level,
        grid.nnodes,
        grid.nxny,
        grid.dmin,
        grid.perturbed,
    )
end

function instance_grid(evtsta, raytracer, IP)
    lims = IP.lims
    refm = readdlm(IP.velocitymodel,' ',Float64,'\n')   # -- reads the reference 1D velocity model

    source = ""
    length(raytracer.source2receivers) == length(raytracer.local_evts) ? source = "evt" : source = "sta"
    nnodes, fw_level, perturb = raytracer.nnodes, raytracer.fw_level, raytracer.perturb
    nn1, nn2, nn3 = nnodes[1], nnodes[2], nnodes[3]
    θmin, θmax = (lims.lat[1]), (lims.lat[2])
    φmin, φmax = (lims.lon[1]), (lims.lon[2])
    rmin, rmax = R + lims.depth[1], R + lims.depth[2]

    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()
    θ = Vector{Float64}()
    φ = Vector{Float64}()
    r = Vector{Float64}()
    vp = Vector{Float64}()
    vs = Vector{Float64}()

    θp = ifelse(
        θmin == θmax,
        [θmin],
        range(θmin,θmax,length=nn1)
    )

    φp = ifelse(
        φmin == φmax,
        [φmin],
        range(φmin,φmax,length=nn2)
    )

    rp = ifelse(
        rmin == rmax,
        [rmin],
        range(rmin,rmax,length=nn3)
    )

    if length(θp) > 1 
        θp1, θp2 = θp[1], θp[2] 
    else
        θp1, θp2 = θp[1], θp[1]
    end
    if length(φp) > 1 
        φp1, φp2 = φp[1], φp[2]
    else
        φp1, φp2 = φp[1], φp[2]
    end
    if length(rp) > 1 
        rp1, rp2 = rp[1], rp[2]
    else
        rp1, rp2 = rp[1], rp[1]
    end

    dθ, dφ, dr = θp2 - θp1, φp2 - φp1, rp2 - rp1
    # x2 = @cartesian(θp2,φp2,rp2)
    # x1 = @cartesian(θp1,φp1,rp1)
    x2 = geo_to_cartesian(θp2,φp2,rp2)
    x1 = geo_to_cartesian(θp1,φp1,rp1)
    dx, dy, dz = abs(x2[1] - x1[1]), abs(x2[2] - x1[2]), abs(x2[3] - x1[3])

    dmin = minimum([dx,dy,dz])
    if length(θp) == 1
        dmin = minimum([dx,dy])
    elseif length(φp) == 1
        dmin = minimum([dx,dz])
    elseif length(rp) == 1
        dmin = minimum([dy,dz])
    end

    gr = GridConst(x, y, z, θ, φ, r, vp, vs, fw_level, [nn1,nn2,nn3], nn1*nn2, dmin, perturb)

    for  k in eachindex(rp), j in eachindex(φp), i in eachindex(θp)
        push!(gr.θ,θp[i])
        push!(gr.φ,φp[j])
        push!(gr.r,rp[k])
        # xt, yt, zt = @cartesian(θp[i],φp[j],rp[k])
        xt, yt, zt = geo_to_cartesian(θp[i],φp[j],rp[k])
        push!(gr.x,xt)
        push!(gr.y,yt)
        push!(gr.z,zt)
        push!(gr.Vp,0.0)
        push!(gr.Vs,0.0)
    end

    # -- perturbs primary grid's nodal positions
    if perturb
        σs = [1,2,3,4]
        for i in eachindex(gr.x)
            σ = rand(σs,3)
            # -- if uniform_noise
                pert_θ = rand()*dθ/σ[1] - dθ/σ[1]/2
                pert_φ = rand()*dφ/σ[2] - dφ/σ[2]/2
                pert_r = rand()*dr/σ[3] - dr/σ[3]/2
                if θp[begin] < (gr.θ[i] + pert_θ) < θp[end]
                    gr.θ[i] += pert_θ
                end
                if φp[begin] < (gr.φ[i] + pert_φ) < φp[end]
                    gr.φ[i] += pert_φ
                end
                if rp[begin] < (gr.r[i] + pert_r) < rp[end]
                    gr.r[i] += pert_r
                end
                # gr.x[i], gr.y[i], gr.z[i] = @cartesian(gr.θ[i],gr.φ[i],gr.r[i])
                gr.x[i], gr.y[i], gr.z[i] = geo_to_cartesian(gr.θ[i],gr.φ[i],gr.r[i])
        end
    end
    
    # -- moves grid points to concide with events and stations' positions
    for id in raytracer.local_evts
        evt = evtsta.evts[id]
        # xe, ye, ze = @cartesian(evt.lat, evt.lon, R + evt.depth)
        xe, ye, ze = geo_to_cartesian(evt.lat, evt.lon, R + evt.depth)
        ind = closest_point(gr,xe,ye,ze)
        gr.x[ind], gr.y[ind], gr.z[ind] = xe, ye, ze
        gr.θ[ind], gr.φ[ind], gr.r[ind] = evt.lat, evt.lon, R + evt.depth
        source == "evt" ? raytracer.source_nodes[id] = ind : raytracer.receiv_nodes[id] = ind
    end
    for id in raytracer.local_stats
        sta = evtsta.stas[id]
        # xs, ys, zs = @cartesian(sta.lat, sta.lon, R + sta.elevation)
        xs, ys, zs = geo_to_cartesian(sta.lat, sta.lon, R + sta.elevation)
        ind = closest_point(gr,xs,ys,zs)
        gr.x[ind], gr.y[ind], gr.z[ind] = xs, ys, zs
        gr.θ[ind], gr.φ[ind], gr.r[ind] = sta.lat, sta.lon, R + sta.elevation
        source == "evt" ? raytracer.receiv_nodes[id] = ind : raytracer.source_nodes[id] = ind
    end

    # minimum_dist(grid)

    for i in eachindex(gr.r)
        gr.Vp[i] = ref_V1D(gr.r[i],refm[:,[1,2]])
        gr.Vs[i] = ref_V1D(gr.r[i],refm[:,[1,3]])
    end

    return gr
end

function LinearIndex(gr, i, j, k)
    nx = gr.nnodes[1]
    return (k-1)*gr.nxny + (j-1)*nx + i
end

function CartesianIndex(gr, I)
    nx = gr.nnodes[1]
    k = div(I,gr.nxny,RoundUp)
    j = div(I-gr.nxny*(k-1), nx, RoundUp)
    i = mod(I-1, nx)+1
    return (i, j, k)
end

function closest_point(gr, px, py, pz; system = :cartesian)
    n = length(gr.x)
    dist = Inf
    di = 0.0
    index = -1

    v_x, v_y, v_z = ifelse(
        system == :cartesian,
        (px, py, pz),
        geo_to_cartesian(px, py, pz)
    )

    @inbounds for i in 1:n

        di = compute_distance(gr.x[i], gr.y[i], gr.z[i], v_x, v_y, v_z)

        if di < dist
            index = i
            dist = di
        end

    end
    return index
end

function compute_distance(x1,y1,z1,x2,y2,z2)
    return sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)
end 

function instance_velocity_grid(evtsta,raytracer, IP, chains)
    lims = IP.lims
    refm = readdlm(string("../../VelocityModels/",IP.velocitymodel),' ',Float64,'\n')   # -- reads the reference 1D velocity model

    vel_nnodes, time_steps = 64000, 50

    nnodes, upscale = upscaling(raytracer,vel_nnodes)
    nn1, nn2, nn3 = nnodes[1], nnodes[2], nnodes[3]
    θmin, θmax = (lims.lat[1]), (lims.lat[2])
    φmin, φmax = (lims.lon[1]), (lims.lon[2])
    rmin, rmax = R + lims.depth[1], R + lims.depth[2]
    tmin, tmax = evtt_extremes(evtsta)

    x = Vector{Float64}()
    y = Vector{Float64}()
    z = Vector{Float64}()
    r = Vector{Float64}()
    Vp = Vector{Vector{Float64}}()
    Vs = Vector{Vector{Float64}}()

    θp = ifelse(
        θmin == θmax,
        [θmin],
        range(θmin,θmax,length=nn1)
    )

    φp = ifelse(
        φmin == φmax,
        [φmin],
        range(φmin,φmax,length=nn2)
    )

    rp = ifelse(
        rmin == rmax,
        [rmin],
        range(rmin,rmax,length=nn3)
    )

    tp = range(tmin,tmax,length=time_steps)
    for frame in eachindex(tp)
        push!(Vp,Float64[])
        push!(Vs,Float64[])
    end

    gr = Velocity4DGridConst(θp, φp, rp, tp, x, y, z, r, Vp, Vs, nnodes, nn1*nn2, upscale)

    for  k in eachindex(rp), j in eachindex(φp), i in eachindex(θp)
        # xt, yt, zt = @cartesian(θp[i],φp[j],rp[k])
        xt, yt, zt = geo_to_cartesian(θp[i],φp[j],rp[k])
        push!(gr.x,xt)
        push!(gr.y,yt)
        push!(gr.z,zt)
        push!(gr.r,rp[k])
        Vp_ref = ref_V1D(rp[k],refm[:,[1,2]])
        Vs_ref = ref_V1D(rp[k],refm[:,[1,3]])
        for frame in eachindex(tp)
            push!(gr.Vp[frame],Vp_ref)
            push!(gr.Vs[frame],Vs_ref)
        end
    end

    fill_4Dgrid(gr,chains)

    return gr

end

function grids_map(grid1,grid2)
    node2node = Int64[]
    for ind in eachindex(grid1.Vp)
        i = v_dist_ind(grid1.θ[ind],grid2.θp)
        j = v_dist_ind(grid1.φ[ind],grid2.φp)
        k = v_dist_ind(grid1.r[ind],grid2.rp)
        vel_ind = LinearIndex(grid2, i, j, k)
        push!(node2node,vel_ind)
    end
    return node2node
end

function v_dist_ind(value,points)
    min_ind = zero(Int64)
    min_dist = Inf
	@inbounds for i in eachindex(points)
		distance = abs(value - points[i])
        if distance < min_dist
           min_dist = distance
           min_ind = i
        end
    end
	return min_ind
end

function upscaling(raytracer,nnodes)
    nn1, nn2, nn3 = raytracer.nnodes[1], raytracer.nnodes[2], raytracer.nnodes[3]
    upscale = ((nn1*nn2*nn3)/nnodes)^(1/3)
    nn1 = Int(floor(nn1/upscale))
    nn2 = Int(floor(nn2/upscale))
    nn3 = Int(floor(nn3/upscale))
    return [nn1, nn2, nn3], upscale
end

function evtt_extremes(evtsta)
    T0s = Float64[]
    for evt in evtsta.evts
        push!(T0s,evt.T0)
    end
    return minimum(T0s), maximum(T0s)
end