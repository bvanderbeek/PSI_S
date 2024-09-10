function gnomonic_proj(proj,θ,φ,θmid,φmid)
    cos_Δ = sin(φmid)*sin(φ) + cos(φmid)*cos(φ)*cos(θ-θmid)
    proj[1] = R*cos(φ)*sin(θ-θmid)/cos_Δ
    proj[2] = R*(cos(φmid)*sin(φ) - sin(φmid)*cos(φ)*cos(θ-θmid))/cos_Δ
end


function carve_grid(visited,grid,source,receivers)
    θmid = 0.5*(grid.θ[begin]+grid.θ[end])
    φmid = 0.5*(grid.φ[begin]+grid.φ[end])
    p = Vector{Vector{Float64}}()

    proj = zeros(Float64,2)
    gnomonic_proj(proj,grid.θ[source],grid.φ[source],θmid,φmid)
    push!(p,proj)
    for receiver in receivers
        rproj = zeros(Float64,2)
        gnomonic_proj(rproj,grid.θ[receiver],grid.φ[receiver],θmid,φmid)
        push!(p,rproj)
    end

    npoints = 8
    radius = 3
    radius_φ = deg2rad(1.0 / (111.319 * cos(grid.φ[source])) * radius)
    radius_θ = deg2rad(1.0 / 110.574 * radius)
    dΩ = 2*pi/npoints
    Ω = 0.0
    for i in 1:npoints
        proj = zeros(Float64,2)
        θ_new = grid.θ[source] + radius_θ * sin(Ω)
        φ_new = grid.φ[source] + radius_φ * cos(Ω)
        gnomonic_proj(proj,θ_new,φ_new,θmid,φmid)
        push!(p,proj)
        for receiver in receivers
            proj = zeros(Float64,2)
            radius_φ = deg2rad(1.0 / (111.319 * cos(grid.φ[receiver])) * radius)
            θ_new = grid.θ[receiver] + radius_θ * sin(Ω)
            φ_new = grid.φ[receiver] + radius_φ * cos(Ω)
            gnomonic_proj(proj,θ_new,φ_new,θmid,φmid)
            push!(p,proj)
        end
        Ω += dΩ
    end

    hull = convex_hull(p)
    vhull = VPolygon(hull)

    gproj = zeros(Float64,2)
    for i in eachindex(grid.θ)
        gnomonic_proj(gproj,grid.θ[i],grid.φ[i],θmid,φmid)
        if !(element(Singleton(gproj)) ∈ vhull) 
            (visited[i] = true)
        end
    end

    # fig,ax = PyPlot.subplots()
    # xs = Float64[]
    # ys = Float64[]
    # gproj = zeros(Float64,2)
    # for i in eachindex(grid.x)[1:1:end]
    #     if !visited[i]
    #         gnomonic_proj(gproj,grid.θ[i],grid.φ[i],θmid,φmid)
    #         push!(xs,gproj[1])
    #         push!(ys,gproj[2])
    #     end
    # end
    # ax.scatter(xs,ys,s=1,color="black")
    # for i in eachindex(p)
    #     if i <= (length(receivers) + 1)
    #         ax.scatter(p[i][1],p[i][2],s=30,color="red")
    #     else
    #         ax.scatter(p[i][1],p[i][2],s=20,color="orange")
    #     end
    # end
    # for i in eachindex(hull)
    #     if i < length(hull)
    #         ax.plot([hull[i][1],hull[i+1][1]],[hull[i][2],hull[i+1][2]])
    #     else
    #         ax.plot([hull[i][1],hull[1][1]],[hull[i][2],hull[1][2]])
    #     end
    # end 
    # PyPlot.show()
end