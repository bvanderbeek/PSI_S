# POST PROCESSING
# Contains a number of post-processing functions to make extracting chain metrics and building models
# from rj-mcmc output easier and less ad-hoc. Does not include any plotting functions (see plotting_library.jl)
# because loading graphics packagaes can be problematic on clusters.
# - BPV
# include(@__DIR__() * "/dependencies.jl") # I live in the src directory

#################
# CHAIN METRICS #
#################

# MISSING FEATURES!!!
# Multi-observable inversions
# + List of observation types
# + Single RMS value or stored seperately for each observation type?
# + A noise parameter for each observation type
struct ChainMetrics{I, F, S}
    iteration::Vector{I} # Iteration in chain
    num_accepted::Vector{I} # Number of accepted proposals at each iteration
    fobj::Vector{F} # Objective function
    rrms::Vector{F} # Root-mean-squared observation residual; assumed to be in seconds; what about joint inversion of different measurement units?
    pnoise::Vector{F} # Noise parameter; assumed to be in seconds
    fields::Dict{S, I} # Field names (key) and index (val); assumes indexing of fieldslist does not change across each iteration
    num_cells::Array{I,2} # Number of Voronoi cells in each field (columns) at each iteration (rows)
    # ADDITIONAL FIELDS
    # norm_field::Array{F,2} # Norm of each field?
    # norm_event_statics::Array{F,2} # (iteration, data type); ModelConst.dataspace.Obs[i].estatics
end
function ChainMetrics(Chain::Vector{ModelConst}; nsaved = 1)
    L, F = length(Chain), Chain[1].nfields # Length of chain and number of fields (N, P)
    CM = ChainMetrics(zeros(Int, L), zeros(Int, L), zeros(L), zeros(L), zeros(L), Chain[1].fieldslist, zeros(Int, L, F))
    for (i, chain_i) in enumerate(Chain)
        CM.iteration[i] = (i - 1)*nsaved
        CM.num_accepted[i] = chain_i.accepted[1]
        CM.fobj[i] = chain_i.misfit[1]
        CM.rrms[i] = chain_i.rms[1]
        CM.pnoise[i] = chain_i.dataspace.Obs[1].noise[1]
        for j in 1:F
            CM.num_cells[i,j] = chain_i.fields[j].n[1]
        end
    end
    
    return CM
end

function load_chains_metrics(num_chains; chain_directory = pwd, chain_name = "chain", nsaved = 1)
    v = Vector{ChainMetrics{Int, Float64, String}}(undef, num_chains)
    for i in 1:num_chains
        chain_file = chain_directory*"/"*chain_name*"."*string(i)*".jld"
        chain = load(chain_file, "MarkovChain")
        v[i] = ChainMetrics(chain; nsaved = nsaved)
    end
    return v
end


##################
# ENSEMBLE MODEL #
##################

function write_ensemble_model_vtk(vtk_file, dims, grid_points; scalar_fields = ((), ()), vector_fields = ((), ()), tensor_fields = ((), ()))
    # Parse tuples
    sfield_names, sfield = scalar_fields # Scalar fields are Matrices with each column corresponding to a field!
    vfield_names, vfield = vector_fields # Vector fields are Vectors (or Tuples) of Matrices
    tfield_names, tfield = tensor_fields # Tensor fields are Vectors (or Tuples) of Arrays!
    # Reshape the grid
    xg, yg, zg = reshape(grid_points[1,:], dims), reshape(grid_points[2,:], dims), reshape(grid_points[3,:], dims)
    # Write VTK file
    vtk_grid(vtk_file, xg, yg, zg) do vtk
        # Write Scalar fields
        for (i, a_field) in enumerate(sfield_names)
            vtk[a_field] = reshape(sfield[i,:], dims)
        end
        # Write vector fields
        for (i, a_field) in enumerate(vfield_names)
            vtk[a_field] = reshape(vfield[i], (3, dims...)) # Note, 'dims...' uses the splat operator to unpack Tuple
        end
        # Write tensor fields
        for (i, a_field) in enumerate(tfield_names)
            vtk[a_field] = reshape(tfield[i], (3, 3, dims...)) # Note, 'dims...' uses the splat operator to unpack Tuple
        end
    end

    return nothing
end
function write_ensemble_model_vtk(lat_limdim, lon_limdim, rad_limdim, num_chains, num_burn, tf_squeeze;
    ensemble_fields = (), correlation_fields = (), directional_fields = (), f_spherical = (a,b,c) -> (a,b,c),
    chain_directory = pwd, chain_name = "chain", vtk_file = chain_directory * "/EnsembleModel",
    tf_write_major_only = true, tf_write_tensors = true, tf_local_to_global = true)

    # Initialise tuples to collect field names
    names_scalar_fields, names_vector_fields, names_tensor_fields = (), (), ()

    # Means and stadard deviations will be computed for the two correlation fields.
    # Remove them from the ensemble list to prevent double calculations and non-unique field names in vtk file
    if !isempty(ensemble_fields) && !isempty(correlation_fields)
        # Remove correlation field 1 from ensemble list
        irmv = findfirst(ensemble_fields .== correlation_fields[1])
        !isnothing(irmv) && deleteat!(ensemble_fields, irmv)
        # Remove correlation field 2 from ensemble list
        irmv = findfirst(ensemble_fields .== correlation_fields[2])
        !isnothing(irmv) && deleteat!(ensemble_fields, irmv)
    end

    # Unpack regular grid parameters
    minimum_latitude, maximum_latitude, n_latitude = lat_limdim
    minimum_longitude, maximum_longitude, n_longitude = lon_limdim
    minimum_radius, maximum_radius, n_radius = rad_limdim
    num_pts = n_latitude*n_longitude*n_radius

    # Build an array points that defines the ensemble model. Each sampled voronoi model in the posterior
    # will be interpolated to these points. This array should have dimensions (3, number_of_points) where
    # the cartesian x,y,z coordintes are stored along the columns for each point.
    xglobal, _ = build_global_geographic_array(minimum_latitude, maximum_latitude, n_latitude,
        minimum_longitude, maximum_longitude, n_longitude,
        minimum_radius, maximum_radius, n_radius)
    
    # Define ensemble model. Computes mean and standard deviation of the posterior at each point (columns)
    # defined in xglobal for each field (rows). The returned arrays have dimensions (num_fields, num_points).
    if isempty(ensemble_fields)
        mean_model, sdev_model = zeros(0, num_pts), zeros(0, num_pts)
    else
        mean_model, sdev_model = build_ensemble_model(num_chains, num_burn, xglobal, ensemble_fields;
            chain_directory = chain_directory, chain_name = chain_name, tf_squeeze = tf_squeeze)
        # Collect field names
        names_scalar_fields = (names_scalar_fields..., ensemble_fields..., "sdev_" .* ensemble_fields...)
    end

    # Define ensemble model. Computes means, standard deviations, and linear correlation among two parameter fields
    # from the posterior samples at each point (columns) defined in xglobal. The returned mean and standard deviation arrays
    # have dimensions (2, num_points) while the correlation array has dimensions (1, num_pts)
    if isempty(correlation_fields)
        c_mean_model, c_sdev_model, pccf_model = zeros(0, num_pts), zeros(0, num_pts), zeros(0, num_pts)
    else
        c_mean_model, c_sdev_model, pccf_model = build_ensemble_correlation(num_chains, num_burn, xglobal,
            correlation_fields[1], correlation_fields[2]; chain_directory = chain_directory, chain_name = chain_name, tf_squeeze = tf_squeeze)
        # Collect field names
        names_scalar_fields = (names_scalar_fields..., correlation_fields..., "sdev_" .* correlation_fields...,
        "pccf_" * correlation_fields[1] * "_" * correlation_fields[2])
    end

    # Define principal axes model. To average directional parameters (i.e. anisotropy), we construct a directional sampling
    # tensor (se Eq. Appendix B of Munzarova et. al, GJI 2018) that described an ellipsoid. The major axis of this ellipsoid
    # parallels the mean direction while the minor axes reflect uncertianty in this orientation. The returned mean tensor
    # is an array of size (3, 3, num_points) with tensors at each grid point. The principal axes of this tensor (paxes_1,2,3)
    # are also returned and scaled by the eigenvalues (paxes_1 = major). These have dimensions (3, num_points). Lastly, the 
    # directional bias (see Wikipedia entry for fractional anisotropy) is returned as a (1, num_pts) array. A value of 0 implies no
    # preferred orientation while a value of 1 implies a single unique orientation.
    if isempty(directional_fields)
        paxes_1, paxes_2, paxes_3 = zeros(0, num_pts), zeros(0, num_pts), zeros(0, num_pts)
        mean_tensor, directional_bias = zeros(0, 0, num_pts), zeros(0, num_pts)
    else
        mean_tensor, paxes_1, paxes_2, paxes_3, directional_bias = build_ensemble_tensor(num_chains, num_burn, xglobal, directional_fields;
            spherical_conversion = f_spherical, tf_local_to_global = tf_local_to_global,
            chain_directory = chain_directory, chain_name = chain_name)
        # Collect field names
        names_scalar_fields = (names_scalar_fields..., "directional_bias")
        names_vector_fields = (names_vector_fields..., "paxes_1", "paxes_2", "paxes_3")
        names_tensor_fields = (names_tensor_fields..., "directional_tensor")
    end

    # Collect fields (mind the ordering!)
    scalar_fields = (names_scalar_fields, ApplyArray(vcat, mean_model, sdev_model, c_mean_model, c_sdev_model, pccf_model, directional_bias))
    # Option to only write major axis
    if !isempty(directional_fields) && tf_write_major_only
        vector_fields = (names_vector_fields[1:1], [paxes_1])
    else
        vector_fields = (names_vector_fields, [paxes_1, paxes_2, paxes_3])
    end
    # Option to include tensor fields
    tensor_fields = tf_write_tensors ? (names_tensor_fields, [mean_tensor]) : ((), ())

    # Some screen output
    println("Preparing to write the following fields:")
    println("Scalar Fields: ", scalar_fields[1])
    println("Vector Fields: ", vector_fields[1])
    println("Tensor Fields: ", tensor_fields[1])

    # Save model file
    write_ensemble_model_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal;
    scalar_fields = scalar_fields,
    vector_fields = vector_fields,
    tensor_fields = tensor_fields)

    return xglobal, scalar_fields, vector_fields, tensor_fields
end

# Store posterior samples for given set of points and fields
function build_posterior(num_chains, num_burn, num_out, sample_points, field_names;
    chain_directory = pwd, chain_name = "chain")
    # Allocate model arrays
    num_pts = size(sample_points, 2) # Number of sample points
    num_fields = length(field_names) # Number of fields
    posterior_samples = zeros(num_fields, num_chains*(num_out - num_burn), num_pts)

    j = 0 # Sample counter
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        for chain_n_j in chain_n # Loop over iterations in chain
            j += 1 # Increment sample counter
            for (i, field_i) in enumerate(field_names) # Loop over parameter fields
                fid = chain_n_j.fieldslist[field_i]
                index = nearest_voronoi_cell(sample_points, chain_n_j.fields[fid])
                posterior_samples[i,j,:] .= chain_n_j.fields[fid].v[index]
            end
        end

        println("Finished interpolating chain "*string(n)*" of "*string(num_chains)*"...")
    end

    return posterior_samples
end

# ASSUMPTIONS
# + All chains contain same set of fields (order can differ)
# + Chain file naming format is, chain_name.number.jld
function build_ensemble_model(num_chains, num_burn, grid_points, field_names;
    chain_directory = pwd, chain_name = "chain", tf_squeeze = false)
    # Allocate model arrays
    num_pts = size(grid_points, 2) # Number of points in model grid
    num_fields = length(field_names)
    mean_model = zeros(num_fields, num_pts) # Field average
    sdev_model = zeros(num_fields, num_pts) # Field standard deviation

    # Accumulate Chains
    num_models = 0
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        num_models += length(chain_n) # Number of models in chain
        accumulate_chain!(mean_model, sdev_model, grid_points, field_names, tf_squeeze, chain_n)

        println("Finished interpolating chain "*string(n)*" of "*string(num_chains)*"...")
    end

    # Finish statistics calculations
    mean_model ./= num_models
    sdev_model ./= num_models
    sdev_model .-= (mean_model.^2)
    sdev_model .= sqrt.(sdev_model)

    return mean_model, sdev_model
end

function build_ensemble_correlation(num_chains, num_burn, grid_points, field_a, field_b;
    chain_directory = pwd, chain_name = "chain", tf_squeeze = false)
    # Allocate model arrays
    num_pts = size(grid_points, 2) # Number of points in model grid
    mean_model = zeros(2, num_pts) # Field average
    sdev_model = zeros(2, num_pts) # Field standard deviation
    pccf_model = zeros(1, num_pts) # Pearson correlation coefficient

    # Accumulate Chains
    num_models = 0
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        num_models += length(chain_n) # Number of models in chain
        for chain_n_j in chain_n # Loop over iterations in chain
            accumulate_correlation!(mean_model, sdev_model, pccf_model, grid_points, field_a, field_b, tf_squeeze, chain_n_j)
        end

        println("Finished interpolating chain "*string(n)*" of "*string(num_chains)*"...")
    end

    # Average model
    mean_model ./= num_models
    # Standard deviation
    sdev_model ./= num_models
    sdev_model .-= (mean_model.^2)
    sdev_model .= sqrt.(sdev_model)
    # Linear correlation: (μ₁₂ - μ₁μ₂)/(σ₁σ₂)
    pccf_model ./= num_models
    @views @. pccf_model[1,:] -= mean_model[1,:] * mean_model[2,:]
    @views @. pccf_model[1,:] /= sdev_model[1,:] * sdev_model[2,:]

    return mean_model, sdev_model, pccf_model
end

function build_ensemble_tensor(num_chains, num_burn, grid_points, directional_fields;
    spherical_conversion = (a,b,c) -> (a,b,c), tf_local_to_global = false,
    chain_directory = pwd, chain_name = "chain")
    # Allocate model arrays
    num_pts = size(grid_points, 2) # Number of points in model grid
    paxes_1, paxes_2, paxes_3 = zeros(3, num_pts), zeros(3, num_pts), zeros(3, num_pts) # Principal axes of tensor
    mean_tensor, directional_bias = zeros(3, 3, num_pts), zeros(1, num_pts)

    # Accumulate Chains
    num_models = 0
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        num_models += length(chain_n) # Number of models in chain
        for chain_n_i in chain_n # Loop over iterations in chain
            accumulate_tensor!(mean_tensor, grid_points, directional_fields, chain_n_i; spherical_conversion = spherical_conversion)
        end

        println("Finished interpolating chain "*string(n)*" of "*string(num_chains)*"...")
    end
    mean_tensor ./= num_models

    # Compute symmetry axes
    root_half = sqrt(0.5)
    Threads.@threads for i in 1:num_pts
        T = @view mean_tensor[:,:,i]
        F = eigen(T)
        # Sorting indices for increasing *magnitude*
        p1, p2, p3 = sortperm(abs.(F.values), rev = true)
        paxes_1[:,i] .= F.values[p1]*F.vectors[:,p1]
        paxes_2[:,i] .= F.values[p2]*F.vectors[:,p2]
        paxes_3[:,i] .= F.values[p3]*F.vectors[:,p3]
        # Compute directional bias (fractional anisotropy of diffusion; see Wikipedia entry for fractional anisotropy)
        λ1, λ2, λ3 = abs(F.values[p1]), abs(F.values[p2]), abs(F.values[p3])
        directional_bias[1,i] = root_half*sqrt( ((λ1 - λ2)^2) + ((λ1 - λ3)^2) + ((λ2 - λ3)^2) )/sqrt( (λ1^2) + (λ2^2) + (λ3^2) )
    end

    # Rotate to global cartesian coordinates
    if tf_local_to_global
        local_to_global_vectors!(paxes_1, paxes_2, paxes_3, grid_points)
        local_to_global_tensors!(mean_tensor, grid_points)
    end

    return mean_tensor, paxes_1, paxes_2, paxes_3, directional_bias
end

# Accumulate models in single chain
function accumulate_chain!(mean_model, sdev_model, grid_points, field_names, tf_squeeze, chain::Vector{ModelConst})
    for chain_j in chain # Loop over iterations in chain
        for (i, field_i) in enumerate(field_names) # Loop over parameter fields
            fid = chain_j.fieldslist[field_i]
            out_sum = @view mean_model[i,:]
            out_squared_sum = @view sdev_model[i,:]
            accumulate_voronoi_field!(out_sum, out_squared_sum, grid_points, tf_squeeze, chain_j.fields[fid])
        end
    end

    return nothing
end

function accumulate_voronoi_field!(out_sum, out_squared_sum, query_points, tf_squeeze, Field::Voronoi; leafsize = 10)
    nuclei = @view Field.c[:,1:Field.n[1]] # View into valid nuclei
    Tree = KDTree(nuclei; leafsize = leafsize)
    index, _ = nn(Tree, query_points) # Very slow serially with lots of allocations
    tf_in_bounds = true # Assume initially in bounds
    for (jq, icell) = enumerate(index)
        if tf_squeeze
            xj, yj, zj = query_points[:, jq]
            λj, ϕj, rj = atan(yj, xj), atan(zj, sqrt((xj^2) + (yj^2))), sqrt((xj^2) + (yj^2) + (zj^2))
            tf_in_bounds = (ϕj >= Field.slims[1][1]) && (ϕj <= Field.slims[1][2]) &&
                           (λj >= Field.slims[2][1]) && (λj <= Field.slims[2][2]) &&
                           (rj >= Field.slims[3][1]) && (rj <= Field.slims[3][2])
        end
        val = tf_in_bounds ? Field.v[icell] : Field.ref_value
        out_sum[jq] += val
        out_squared_sum[jq] += val^2
    end
    return nothing
end

function accumulate_tensor!(tensor_sum, query_points, directional_fields, Model::ModelConst; spherical_conversion = (a,b,c) -> (a,b,c))
    # Index directional model fields
    j_1 = Model.fieldslist[directional_fields[1]]
    j_2 = Model.fieldslist[directional_fields[2]]
    j_3 = Model.fieldslist[directional_fields[3]]
    # Index nearest voronoi neighbors
    index_1 = nearest_voronoi_cell(query_points,  Model.fields[j_1])
    index_2 = nearest_voronoi_cell(query_points,  Model.fields[j_2])
    index_3 = nearest_voronoi_cell(query_points,  Model.fields[j_3])
    # Accumulate tensor elements
    for k = axes(query_points, 2)
        # Convert directional parameters to spherical
        a, b, c = Model.fields[j_1].v[index_1[k]], Model.fields[j_2].v[index_2[k]], Model.fields[j_3].v[index_3[k]]
        f_j, azm_j, elv_j = spherical_conversion(a,b,c)
        T = @view tensor_sum[:,:,k]
        sum_directional_tensor!(T, f_j, azm_j, elv_j)
    end
    
    return nothing
end

# Assumes only two fields are passed
function accumulate_correlation!(f_sum, f_squared_sum, f_prod, query_points, field_a, field_b, tf_squeeze, Model::ModelConst)
    # Nearest neighbor interpolation fields
    j_a, j_b = Model.fieldslist[field_a], Model.fieldslist[field_b]
    index_a = nearest_voronoi_cell(query_points,  Model.fields[j_a])
    index_b = nearest_voronoi_cell(query_points,  Model.fields[j_b])
    # Sum product of fields
    tf_in_bounds_a, tf_in_bounds_b = true, true # Assume initially in bounds
    for k = axes(query_points, 2)
        if tf_squeeze
            xj, yj, zj = query_points[:, k]
            λj, ϕj, rj = atan(yj, xj), atan(zj, sqrt((xj^2) + (yj^2))), sqrt((xj^2) + (yj^2) + (zj^2))
            tf_in_bounds_a = (ϕj >= Model.fields[j_a].slims[1][1]) && (ϕj <= Model.fields[j_a].slims[1][2]) &&
                           (λj >= Model.fields[j_a].slims[2][1]) && (λj <= Model.fields[j_a].slims[2][2]) &&
                           (rj >= Model.fields[j_a].slims[3][1]) && (rj <= Model.fields[j_a].slims[3][2])
            tf_in_bounds_b = (ϕj >= Model.fields[j_b].slims[1][1]) && (ϕj <= Model.fields[j_b].slims[1][2]) &&
                           (λj >= Model.fields[j_b].slims[2][1]) && (λj <= Model.fields[j_b].slims[2][2]) &&
                           (rj >= Model.fields[j_b].slims[3][1]) && (rj <= Model.fields[j_b].slims[3][2])
        end
        val_a = tf_in_bounds_a ? Model.fields[j_a].v[index_a[k]] : Model.fields[j_a].ref_value
        val_b = tf_in_bounds_b ? Model.fields[j_b].v[index_b[k]] : Model.fields[j_b].ref_value
        # Field A
        f_sum[1,k] += val_a
        f_squared_sum[1,k] += val_a^2
        # Field B
        f_sum[2,k] += val_b
        f_squared_sum[2,k] += val_b^2
        # AB Correlation Product
        if tf_in_bounds_a && tf_in_bounds_b 
            f_prod[k] += val_a * val_b
        end
    end

    return nothing
end

function nearest_voronoi_cell(query_points, Field::Voronoi; leafsize = 10)
    nuclei = @view Field.c[:,1:Field.n[1]] # View into valid nuclei
    Tree = KDTree(nuclei; leafsize = leafsize)
    index, _ = nn(Tree, query_points) # Very slow serially with lots of allocations
    # Field.v[index] # Nearest value
    return index
end

function sum_directional_tensor!(T, magnitude, azimuth, elevation)
    # Compute trigonometric terms
    sinλ, cosλ = sincos(azimuth)
    sinϕ, cosϕ = sincos(elevation)
    cosλ2, cosϕ2 = cosλ^2, cosϕ^2
    sinλ2, sinϕ2 = 1.0 - cosλ2, 1.0 - cosϕ2
    # Update global directional tensor
    T[1,1] += magnitude*cosϕ2*cosλ2
    T[1,2] += magnitude*cosϕ2*cosλ*sinλ
    T[1,3] += magnitude*cosϕ*cosλ*sinϕ
    T[2,2] += magnitude*cosϕ2*sinλ2
    T[2,3] += magnitude*cosϕ*sinϕ*sinλ
    T[3,3] += magnitude*sinϕ2
    # Symmetries
    T[2,1], T[3,1], T[3,2] = T[1,2], T[1,3], T[2,3]

    return nothing
end

function local_to_global_vectors!(paxes_1, paxes_2, paxes_3, xglobal)
    for i in axes(xglobal, 2)
        # Global cartesian coordinates to spherical (i.e. geographic)
        xg, yg, zg = xglobal[1,i], xglobal[2,i], xglobal[3,i]
        λ, ϕ = atan(yg, xg), atan(zg, sqrt((xg^2) + (yg^2)))
        # East-North-Elevation components of vectors
        a = (paxes_1[1,i], paxes_1[2,i], paxes_1[3,i])
        b = (paxes_2[1,i], paxes_2[2,i], paxes_2[3,i])
        c = (paxes_3[1,i], paxes_3[2,i], paxes_3[3,i])
        # Rotate vectors into Earth-Centered Earth-Fixed (i.e. global cartesian) coordinates
        paxes_1[1,i], paxes_1[2,i], paxes_1[3,i] = ecef_vector(a, ϕ, λ; c = 1.0) # Radians (c = 1.0)!
        paxes_2[1,i], paxes_2[2,i], paxes_2[3,i] = ecef_vector(b, ϕ, λ; c = 1.0) # Radians (c = 1.0)!
        paxes_3[1,i], paxes_3[2,i], paxes_3[3,i] = ecef_vector(c, ϕ, λ; c = 1.0) # Radians (c = 1.0)!
    end

    return nothing
end

function local_to_global_vectors!(vec_array, xglobal)
    for i in axes(xglobal, 2)
        # Global cartesian coordinates to spherical (i.e. geographic)
        xg, yg, zg = xglobal[1,i], xglobal[2,i], xglobal[3,i]
        λ, ϕ = atan(yg, xg), atan(zg, sqrt((xg^2) + (yg^2)))
        # East-North-Elevation components of vectors
        a = (vec_array[1,i], vec_array[2,i], vec_array[3,i])
        # Rotate vectors into Earth-Centered Earth-Fixed (i.e. global cartesian) coordinates
        vec_array[1,i], vec_array[2,i], vec_array[3,i] = ecef_vector(a, ϕ, λ; c = 1.0) # Radians (c = 1.0)!
    end

    return nothing
end

function local_to_global_tensors!(local_tensors, xglobal)
    # Local tensor coordinates:
    # (1 -> east -> yglob), (2 -> north -> zglob), (3 -> radial -> xglob)
    for i in axes(xglobal, 2)
        # Global cartesian coordinates to spherical (i.e. geographic)
        xg, yg, zg = xglobal[1,i], xglobal[2,i], xglobal[3,i]
        λ, ϕ = atan(yg, xg), atan(zg, sqrt((xg^2) + (yg^2)))
        # Local tensor components in global cartesian coordinates
        T = @SMatrix [
            local_tensors[3,3,i] local_tensors[3,1,i] local_tensors[3,2,i];
            local_tensors[1,3,i] local_tensors[1,1,i] local_tensors[1,2,i];
            local_tensors[2,3,i] local_tensors[2,1,i] local_tensors[2,2,i];
        ]
        # Rotate from (0N, 0E) to (ϕ, λ)
        R = rotation_matrix((-ϕ, λ), (2, 3))
        local_tensors[:,:,i] .= R*T*transpose(R)
    end

    return nothing
end

function prior_directional_tensor(a, b, c; spherical_conversion = (a,b,c) -> (a,b,c))
    T = zeros(3,3)
    for i in eachindex(a)
        r_i, λ_i, ϕ_i = spherical_conversion(a[i], b[i], c[i])
        sum_directional_tensor!(T, r_i, λ_i, ϕ_i);
    end
    T ./= length(a)
    F = eigen(T)
    p1, p2, p3 = sortperm(abs.(F.values), rev = true)
    λ1, λ2, λ3 = abs(F.values[p1]), abs(F.values[p2]), abs(F.values[p3])
    db = sqrt(0.5)*sqrt( ((λ1 - λ2)^2) + ((λ1 - λ3)^2) + ((λ2 - λ3)^2) )/sqrt( (λ1^2) + (λ2^2) + (λ3^2) )

    return db, T
end

#############
# UTILITIES #
#############

function build_global_geographic_array(minimum_latitude, maximum_latitude, n_latitude,
    minimum_longitude, maximum_longitude, n_longitude,
    minimum_radius, maximum_radius, n_radius)
    # Regular geographic coordinate arrays
    latitude = range(start = deg2rad(minimum_latitude), stop = deg2rad(maximum_latitude), length = n_latitude)
    longitude = range(start = deg2rad(minimum_longitude), stop = deg2rad(maximum_longitude), length = n_longitude)
    radius = range(start = minimum_radius, stop = maximum_radius, length = n_radius)

    # Global cartesian coordinate array -- irregularly spaced
    xglobal = zeros(3, n_longitude*n_latitude*n_radius)
    for k in eachindex(radius), j in eachindex(longitude), i in eachindex(latitude)
        n = subscripts_to_index((n_latitude, n_longitude, n_radius), (i, j, k))
        xglobal[1,n], xglobal[2,n], xglobal[3,n] = geo_to_cartesian(latitude[i], longitude[j], radius[k])
    end

    return xglobal, latitude, longitude, radius
end

function subscripts_to_index(dimsize::NTuple{N, Int}, subscripts::NTuple{N, Int}) where {N}
    # Loop over array dimensions
    nd = length(dimsize)
    ind = min(max(subscripts[1], 1), dimsize[1]) - 1
    stride = 1
    for i in 2:nd
        sub = min(max(subscripts[i], 1), dimsize[i]) - 1
        # Count indices
        stride *= dimsize[i - 1]
        ind += stride*sub
    end
    ind += 1

    return ind
end

function ecef_vector(east_north_elv::NTuple{3, T}, latitude, longitude; c = π/180.0) where {T}
    # ECEF components for vector at (0°N, 0°E)
    w = (east_north_elv[3], east_north_elv[1], east_north_elv[2])
    # Rotate to geographic position
    R = rotation_matrix((-c*latitude, c*longitude), (2, 3))
    sx = R[1,1]*w[1] + R[1,2]*w[2] + R[1,3]*w[3]
    sy = R[2,1]*w[1] + R[2,2]*w[2] + R[2,3]*w[3]
    sz = R[3,1]*w[1] + R[3,2]*w[2] + R[3,3]*w[3]

    return sx, sy, sz
end

function rotation_matrix(α, n)
    sinα, cosα = sincos(α)
    if n == 1
        return @SMatrix [1.0 0.0 0.0; 0.0 cosα -sinα; 0.0 sinα cosα]
    elseif n == 2
        return @SMatrix [cosα 0.0 sinα; 0.0 1.0 0.0; -sinα 0.0 cosα]
    elseif n == 3
        return @SMatrix [cosα -sinα 0.0; sinα cosα 0.0; 0.0 0.0 1.0]
    else
        error("Requested rotation axis out-of-bounds!")
    end
end

function rotation_matrix(α::Tuple, n::Tuple)
    R = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(n)
        R = rotation_matrix(α[i],n[i])*R
    end
    return R
end