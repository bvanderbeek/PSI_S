## SETUP ##
using JLD
using PSI_S

# START INPUT GRID PARAMETERS
minimum_latitude, maximum_latitude, n_latitude = -20.2, -10.6, 49
minimum_longitude, maximum_longitude, n_longitude = -121.4, -102.6, 95
minimum_radius, maximum_radius, n_radius = 6371.0-600.0, 6371.0, 31
# STOP INPUT GRID PARAMETERS

name = ARGS[1]
models = load(string(@__DIR__, "/output/", name, "/$(name)_inv.jld"),"MarkovChain")
nmodels = length(models)

function updateT(T,f,𝛙,γ,N)
    i = pi/2 - γ
    ϕ = pi/2 - 𝛙

    T[1,1] += f*sin(i)^2*sin(ϕ)^2/N
    T[1,2] += f*sin(i)^2*sin(ϕ)*cos(ϕ)/N
    T[1,3] += f*sin(i)*cos(i)*sin(ϕ)/N

    T[2,1] += f*sin(i)^2*sin(ϕ)*cos(ϕ)/N
    T[2,2] += f*sin(i)^2*cos(ϕ)^2/N
    T[2,3] += f*sin(i)*cos(i)*cos(ϕ)/N

    T[3,1] += f*sin(i)*cos(i)*sin(ϕ)/N
    T[3,2] += f*sin(i)*cos(i)*cos(ϕ)/N
    T[3,3] += f*cos(i)^2/N
end

latitude = range(start = deg2rad(minimum_latitude), stop = deg2rad(maximum_latitude), length = n_latitude)
longitude = range(start = deg2rad(minimum_longitude), stop = deg2rad(maximum_longitude), length = n_longitude)
radius = range(start = minimum_radius, stop = maximum_radius, length = n_radius)

x_nodes, y_nodes, z_nodes = Float64[], Float64[], Float64[]
θ_nodes, φ_nodes, r_nodes = Float64[], Float64[], Float64[]
for k in eachindex(radius), j in eachindex(longitude), i in eachindex(latitude)
    # x, y, z = @cartesian(latitude[i], longitude[j], radius[k])
    x, y, z = PSI_S.geo_to_cartesian(latitude[i], longitude[j], radius[k])
    push!(x_nodes,x)
    push!(y_nodes,y)
    push!(z_nodes,z)
    push!(θ_nodes,latitude[i])
    push!(φ_nodes,longitude[j])
    push!(r_nodes,radius[k])
end
nnodes = [length(latitude), length(longitude), length(radius)]
grid = GridConst(x_nodes,y_nodes,z_nodes,θ_nodes,φ_nodes,r_nodes,[],[],0,nnodes,nnodes[1]*nnodes[2],0.0,false)
nodes = permutedims(hcat(x_nodes, y_nodes, z_nodes))

dlnVp_mean = zeros(Float64,length(latitude),length(longitude),length(radius))
dlnVp_relstd = zeros(Float64,length(latitude),length(longitude),length(radius))
slowAxes = Array{Vector{Vector{Float64}}}(undef,length(latitude),length(longitude),length(radius))
slowTensors = Array{Matrix{Float64}}(undef,length(latitude),length(longitude),length(radius))
for k in axes(slowTensors,3), j in axes(slowTensors,2), i in axes(slowTensors,1)
    slowTensors[i,j,k] = zeros(Float64,3,3)
    slowAxes[i,j,k] = Vector{Vector{Float64}}()
    [push!(slowAxes[i,j,k],zeros(Float64,3)) for l in 1:3]
end

for imodel in eachindex(models)
    model = models[imodel]
    print("progress: ",round(imodel/nmodels*100; digits=2),"%\r")
    dlnVp = model.fields[1]
    # Fp = model.fields[2] 
    # 𝛙 = model.fields[3]
    # γ = model.fields[4]

    # -- dlnVp
    n = dlnVp.n[1]
    nuclei = permutedims(hcat(dlnVp.c[1,1:n], dlnVp.c[2,1:n], dlnVp.c[3,1:n]))
    inds_dlnVp = PSI_S.NN_interpolation(nodes, nuclei)
    # -- anisotropy
    # nuclei = permutedims(hcat(Fp.x, Fp.y, Fp.z))
    # inds_Fp = NN_interpolation(nodes, nuclei)
    # nuclei = permutedims(hcat(𝛙.x, 𝛙.y, 𝛙.z))
    # inds_𝛙 = NN_interpolation(nodes, nuclei)
    # nuclei = permutedims(hcat(γ.x, γ.y, γ.z))
    # inds_γ = NN_interpolation(nodes, nuclei)

    for n in eachindex(inds_dlnVp)
        i, j, k = PSI_S.CartesianIndex(grid,n)

        if !PSI_S.check_coordinates(models[end].fields[1].slims,grid.θ[n],grid.φ[n],grid.r[n])
            dlnVp.v[inds_dlnVp[n]] = dlnVp.ref_value
        end
        # if !check_coordinates(models[end].fields[2].slims,grid.θ[n],grid.φ[n],grid.r[n])
        #     Fp.v[inds_Fp[n]] = Fp.ref_value
        # end

        dlnVp_mean[i,j,k] += dlnVp.v[inds_dlnVp[n]]
        dlnVp_relstd[i,j,k] += dlnVp.v[inds_dlnVp[n]]^2
        # updateT(slowTensors[i,j,k],Fp.v[inds_Fp[n]],𝛙.v[inds_𝛙[n]],γ.v[inds_γ[n]],nmodels)
    end
end

# Threads.@threads for k in axes(slowTensors,3)
#     for j in axes(slowTensors,2)
#         for i in axes(slowTensors,1)
#             T = slowTensors[i,j,k]
#             vecs = eigvecs(T)
#             vals = abs.(eigvals(T))
#             inds = sortperm(vals)
#             ind_a, ind_b, ind_c = reverse(inds)
#             slowAxes[i,j,k][1] .= [(vals[ind_a])*vecs[3,ind_a],(vals[ind_a])*vecs[1,ind_a],(vals[ind_a])*vecs[2,ind_a]]
#             slowAxes[i,j,k][2] .= [(vals[ind_b])*vecs[3,ind_b],(vals[ind_b])*vecs[1,ind_b],(vals[ind_b])*vecs[2,ind_b]]
#             slowAxes[i,j,k][3] .= [(vals[ind_c])*vecs[3,ind_c],(vals[ind_c])*vecs[1,ind_c],(vals[ind_c])*vecs[2,ind_c]]
#         end
#     end
# end

dlnVp_mean = @. dlnVp_mean / nmodels
dlnVp_relstd = @. dlnVp_relstd / nmodels
dlnVp_relstd = @. sqrt(dlnVp_relstd - dlnVp_mean^2)/abs(dlnVp_mean)
for k in axes(dlnVp_relstd,3), j in axes(dlnVp_relstd,2), i in axes(dlnVp_relstd,1)
    if dlnVp_relstd[i,j,k] > 1
        dlnVp_relstd[i,j,k] = 1.0
    end
end


save(string(@__DIR__, "/output/", name, "/$(name)_grid_model.jld"),"dlnVp_mean",dlnVp_mean,"dlnVp_relstd",dlnVp_relstd,"slowAxes",slowAxes,"slowTensors",slowTensors,"latitude",latitude,"longitude",longitude,"radius",radius)
