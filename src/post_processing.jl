# POST PROCESSING
# Contains a number of post-processing functions to make extracting chain metrics and building models
# from rj-mcmc output easier and less ad-hoc. Does not include any plotting functions (see plotting_library.jl)
# because loading graphics packagaes can be problematic on clusters.
# - BPV
# include(@__DIR__() * "/dependencies.jl") # I live in the src directory
#
# How to parallelize this?
# > Multi-thread loops that do the point-wise interpolation?

#################
# CHAIN METRICS #
#################

struct ChainMetrics{I, F, S}
    obs::Dict{I, S} # Observable index (key) and name (val); assumes indexing of observables does not change across each iteration
    fields::Dict{S, I} # Field names (key) and index (val); assumes indexing of fieldslist does not change across each iteration
    iteration::Vector{I} # Iteration in chain
    num_accepted::Vector{I} # Number of accepted proposals at each iteration
    fobj::Vector{F} # Objective function
    rrms::Array{F,2} # Root-mean-squared observation residual for each observable (columns) at each iteration (rows)
    pnoise::Array{F,2} # Noise parameter for each observable (columns) at each iteration (rows)
    num_cells::Array{I,2} # Number of Voronoi cells in each field (columns) at each iteration (rows)
    norm_field::Array{F,2} # Norm of each field?
    rms_event_statics::Array{F,2} # (iteration, data type); ModelConst.dataspace.Obs[i].estatics
end
function ChainMetrics(Chain::Vector{ModelConst}; nsaved = 1)
    N = length(Chain[1].dataspace.Obs) # Number of observables
    L, F = length(Chain), Chain[1].nfields # Length of chain and number of fields (N, P)
    # Construct observable dictionary
    obs = Dict{Int, String}()
    [obs[i] = B.obsname for (i, B) in enumerate(Chain[1].dataspace.Obs)]
    CM = ChainMetrics(obs, Chain[1].fieldslist, zeros(Int, L), zeros(Int, L), zeros(Float64, L), zeros(Float64, L, N), zeros(Float64, L, N), zeros(Int, L, F), zeros(Float64, L, F), zeros(Float64, L, N))
    for (i, chain_i) in enumerate(Chain) # Chain loop
        CM.iteration[i] = (i - 1)*nsaved
        CM.num_accepted[i] = chain_i.accepted[1]
        CM.fobj[i] = chain_i.misfit[1]
        for j in 1:N # Observable loop
            CM.rrms[i,j] = chain_i.rms[j]
            CM.pnoise[i,j] = chain_i.dataspace.Obs[j].noise[1]
            CM.rms_event_statics[i,j] = sum(x->x^2, chain_i.dataspace.Obs[j].estatics)
            CM.rms_event_statics[i,j] /= length(chain_i.dataspace.Obs[j].estatics)
            CM.rms_event_statics[i,j] = sqrt(CM.rms_event_statics[i,j])
        end
        for j in 1:F
            CM.num_cells[i,j] = chain_i.fields[j].n[1]
            CM.norm_field[i,j] = sum(x->x^2, chain_i.fields[j].v)
            CM.norm_field[i,j] = sqrt(CM.norm_field[i,j])
        end
    end
    
    return CM
end

function load_chains_metrics(num_chains; chain_directory = pwd(), chain_name = "chain", nsaved = 1)
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

# Intended to deprecate write_ensemble_model_vtk. Updated so scalar_fields are indexed as the others.
# (i.e. no longer assumed to be arrays but vectors)
function write_ensemble_vtk(vtk_file, dims, grid_points; scalar_fields = ((), ()), vector_fields = ((), ()), tensor_fields = ((), ()))
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
            vtk[a_field] = reshape(sfield[i], dims)
        end
        # Write vector fields
        for (i, a_field) in enumerate(vfield_names)
            nrow = size(vfield[i], 1)
            vtk[a_field] = reshape(vfield[i], (nrow, dims...)) # Note, 'dims...' uses the splat operator to unpack Tuple
        end
        # Write tensor fields
        for (i, a_field) in enumerate(tfield_names)
            nrow, ncol, _ = size(tfield[i])
            vtk[a_field] = reshape(tfield[i], (3, 3, dims...)) # Note, 'dims...' uses the splat operator to unpack Tuple
        end
    end

    return nothing
end
# Posterior Tendencies
function write_ensemble_vtk(vtk_file, dims, grid_points, PostTend, field_names)

    # Allocations
    num_fields, num_pts = size(PostTend.mean) # Number of moments, grid points, and fields
    num_quantiles = size(PostTend.quantiles, 1)
    # Number of scalar fields. For each parameter field, we have the following:
    # Mean, median, mode, additional quantiles
    num_sca = num_fields*(3 + num_quantiles)
    if haskey(PostTend, :h_distance)
        num_sca += 3*num_fields # Nuclei distances, Hellinger distance, Kolmogorov metric
    end
    vtk_scalar_field = Vector{String}(undef, num_sca)
    vtk_scalar_data = Vector{AbstractVector{Float64}}(undef, num_sca)
    # Number of vector fields (median interval, highest density interval, low-high interval, directional vectors)
    num_vec = 3*num_fields
    if haskey(PostTend, :directional_vectors)
        num_vec += 1
        field_names[2:3] .= "directional_strength", "angular_deviation"
    end
    vtk_vector_field = Vector{String}(undef, num_vec)
    vtk_vector_data = Vector{AbstractArray{Float64}}(undef, num_vec)

    # Define vtk scalar fields
    l = 0
    # Mean
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_scalar_field[l] = "mean_"*field_i
        vtk_scalar_data[l] = view(PostTend.mean, i, :)
    end
    # Median
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_scalar_field[l] = "median_"*field_i
        vtk_scalar_data[l] = view(PostTend.median, i, :)
    end
    # Mode
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_scalar_field[l] = "mode_"*field_i
        vtk_scalar_data[l] = view(PostTend.mode, i, :)
    end
    # Quantiles
    for (i, field_i) in enumerate(field_names)
        for j in 1:num_quantiles
            l += 1
            vtk_scalar_field[l] = "q"*string(j)*"_"*field_i
            vtk_scalar_data[l] = view(PostTend.quantiles, j, i, :)
        end
    end
    # New metrics computed by posterior_tendencies (not yet by posterior_directional_tendencies)
    if haskey(PostTend, :h_distance)
        # Nuclei distances
        for (i, field_i) in enumerate(field_names)
            l += 1
            vtk_scalar_field[l] = "ndist_"*field_i
            vtk_scalar_data[l] = view(PostTend.dist, i, :)
        end
        # Posterior differences
        for (i, field_i) in enumerate(field_names)
            l += 1
            vtk_scalar_field[l] = "hval_"*field_i
            vtk_scalar_data[l] = view(PostTend.h_distance, i, :)
        end
        # Posterior differences
        for (i, field_i) in enumerate(field_names)
            l += 1
            vtk_scalar_field[l] = "kval_"*field_i
            vtk_scalar_data[l] = view(PostTend.k_metric, i, :)
        end
    end


    # Vector Fields: Credible Intervals
    l = 0
    # Median Interval
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_vector_field[l] = "mi_"*field_i
        vtk_vector_data[l] = view(PostTend.median_interval, :, i, :)
    end
    # Highest Density Interval
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_vector_field[l] = "hdi_"*field_i
        vtk_vector_data[l] = view(PostTend.highest_density_interval, :, i, :)
    end
    # Highest Density Interval
    for (i, field_i) in enumerate(field_names)
        l += 1
        vtk_vector_field[l] = "lhi_"*field_i
        vtk_vector_data[l] = view(PostTend.low_high_interval, :, i, :)
    end
    # Directional Vectors
    if haskey(PostTend, :directional_vectors)
        l += 1
        vtk_vector_field[l] = "symmetry_axis"
        vtk_vector_data[l] = view(PostTend.directional_vectors, :, 3, :)
    end

    write_ensemble_vtk(
        vtk_file, dims, grid_points;
        scalar_fields = (vtk_scalar_field, vtk_scalar_data),
        vector_fields = (vtk_vector_field, vtk_vector_data)
    )

    return nothing
end
# Posterior Moments
function write_ensemble_vtk(vtk_file, dims, grid_points, post_moments, post_correlation, field_names, correlation_pairs)

    # Allocations
    num_moments, num_pts, num_fields = size(post_moments) # Number of moments, grid points, and fields
    num_corr = length(correlation_pairs) # Number of correlations
    num_out = num_corr + num_moments*num_fields # Total number of scalar fields
    vtk_scalar_field = Vector{String}(undef, num_out)
    vtk_scalar_data = Vector{AbstractVector{Float64}}(undef, num_out)

    # Define vtk scalar fields
    l = 0
    # Moments
    for (k, field_k) in enumerate(field_names)
        for i in axes(post_moments, 1)
            l += 1
            vtk_scalar_field[l] = "u"*string(i)*"_"*field_k
            vtk_scalar_data[l] = view(post_moments, i, :, k)
        end
    end
    # Correlations
    for (k, pairs_k) in enumerate(correlation_pairs)
        l += 1
        vtk_scalar_field[l] = "pcc_"*pairs_k[1]*"_"*pairs_k[2]
        vtk_scalar_data[l] = view(post_correlation, 1, :, k)
    end

    write_ensemble_vtk(
        vtk_file, dims, grid_points;
        scalar_fields = (vtk_scalar_field, vtk_scalar_data)
    )

    return nothing
end
# Posterior Directional Moments
function write_ensemble_vtk(vtk_file, dims, grid_points, post_moments, post_correlation,
    directional_bias, post_vectors, directional_fields, correlation_fields)

    # Allocations
    num_moments, num_pts, num_fields = size(post_moments) # Number of moments, grid points, and fields
    num_corr = length(correlation_fields) # Number of correlations
    num_out = 1 + num_corr + num_moments*num_fields # Total number of scalar fields
    vtk_scalar_field = Vector{String}(undef, num_out)
    vtk_scalar_data = Vector{AbstractVector{Float64}}(undef, num_out)

    # Define moment fields
    field_names = Vector{String}(undef, num_fields)
    field_names[1], field_names[2], field_names[3] = directional_fields[1], "directional_strength", "angular_deviation"
    [field_names[i+3] = vi for (i, vi) in enumerate(correlation_fields)]
    # Define correlation pairs
    correlation_pairs = Vector{NTuple{2, String}}(undef, num_corr)
    [correlation_pairs[i] = (field_names[2], vi) for (i, vi) in enumerate(correlation_fields)]

    # Define vtk scalar fields
    l = 0
    # Moments
    for (k, field_k) in enumerate(field_names)
        for i in axes(post_moments, 1)
            l += 1
            vtk_scalar_field[l] = "u"*string(i)*"_"*field_k
            vtk_scalar_data[l] = view(post_moments, i, :, k)
        end
    end
    # Correlations
    for (k, pairs_k) in enumerate(correlation_pairs)
        l += 1
        vtk_scalar_field[l] = "pcc_"*pairs_k[1]*"_"*pairs_k[2]
        vtk_scalar_data[l] = view(post_correlation, 1, :, k)
    end
    # Directional Bias
    l += 1
    vtk_scalar_field[l] = "directional_bias"
    vtk_scalar_data[l] = view(directional_bias, 1, :)

    write_ensemble_vtk(
        vtk_file, dims, grid_points;
        scalar_fields = (vtk_scalar_field, vtk_scalar_data),
        vector_fields = (["symmetry_axis"], [view(post_vectors, :, 3, :)])
    )

    return nothing
end

function read_ensemble_vtk(the_vtk; fields = nothing)
    # Get vtk-file info
    vtk = VTKFile(the_vtk)
    point_data = get_point_data(vtk)

    # User-defined fields to load? Default to all fields.
    if isnothing(fields)
        fields = keys(point_data)
    end

    # Read fields
    Model = Dict{String, Array{Float64}}()
    for field_i in fields
        Model[field_i] = get_data_reshaped(point_data[field_i])
    end

    # Read coordinates (global cartesian)
    dims = size(Model[fields[1]])
    dims = length(dims) > 3 ? dims[2:4] : dims
    point_coords = get_points(vtk)
    Model["xgc"] = reshape(point_coords[1,:], dims)
    Model["ygc"] = reshape(point_coords[2,:], dims)
    Model["zgc"] = reshape(point_coords[3,:], dims)

    return Model
end
# Method to write a Model loaded via read_ensemble_vtk
function write_ensemble_vtk(vtk_file, Model::Dict)
    xg = pop!(Model, "xgc")
    yg = pop!(Model, "ygc")
    zg = pop!(Model, "zgc")
    vtk_grid(vtk_file, xg, yg, zg) do vtk
        for (k, v) in Model
            vtk[k] = v
        end
    end
    return nothing
end

# Load all Markov chains into a single array
function load_markov_chains(num_chains, num_burn; chain_directory = pwd(), chain_name = "chain")
    markov_chains = Vector{Any}(undef, num_chains)
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        markov_chains[n] = chain_n
    end

    return markov_chains
end

# Build nearest-neighbor interpolation trees for all diagrams
function build_nn_trees(markov_chains, field_names; leafsize = 10)
    num_chains, num_samp, num_fields = length(markov_chains), length(markov_chains[1]), length(field_names)
    nn_trees = Array{Any}(undef, num_fields, num_samp, num_chains)
    for (k, chain_k) in enumerate(markov_chains) # Loop over Markov chains
        for (j, chain_j_k) in enumerate(chain_k) # Loop over samples in chain
            for (i, field_i) in enumerate(field_names) # Loop over fields in chain
                fid = chain_j_k.fieldslist[field_i]
                Field = chain_j_k.fields[fid]
                nuclei = @view Field.c[:,1:Field.n[1]]
                nn_trees[i,j,k] = KDTree(nuclei; leafsize = leafsize)
            end
        end
    end
    return nn_trees
end

# Collect posterior samples
function collect_posterior_samples(num_chains, num_burn, num_out, sample_points, field_names;
    chain_directory = pwd(), chain_name = "chain")
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
# Collect posterior samples for multiple points
function collect_posterior_samples(sample_points, markov_chains, nn_trees, field_list)
    num_pts = size(sample_points, 2)
    num_fields, num_samples, num_chains = size(nn_trees)
    posterior_samples = Vector{Array{Float64, 2}}(undef, num_pts)
    nuclei_distances = Vector{Array{Float64, 2}}(undef, num_pts)
    [posterior_samples[i] = zeros(num_fields, num_samples*num_chains) for i in 1:num_pts]
    [nuclei_distances[i] = zeros(num_fields, num_samples*num_chains) for i in 1:num_pts]
    for (i, vi) in enumerate(posterior_samples)
        ni = nuclei_distances[i]
        xi = @view sample_points[:,i]
        collect_posterior_samples!(vi, ni, xi, markov_chains, nn_trees, field_list)
    end

    return posterior_samples, nuclei_distances
end
# Collect posterior samples for all fields at a single point and store in pre-allocated array
function collect_posterior_samples!(posterior_samples, nuclei_distances, sample_point, markov_chains, nn_trees, field_list::Vector{String})
    l = 0 # Initialise sample counter
    icell, dist = [0], [0.0] # Arrays for NearestNeighbor indexing
    for (k, chain_k) in enumerate(markov_chains) # Loop over Markov chains
        for (j, chain_k_j) in enumerate(chain_k) # Loop over samples in a Markov chain
            l += 1 # Increment sample counter
            for (i, field_i) in enumerate(field_list) # Loop over parameter fields
                fid = chain_k_j.fieldslist[field_i]
                knn!(icell, dist, nn_trees[i, j, k], sample_point, 1)
                posterior_samples[i, l] = chain_k_j.fields[fid].v[icell[1]]
                nuclei_distances[i, l] = dist[1]
            end
        end
    end

    return nothing
end

# --- Posterior Tendencies --- #

function posterior_tendencies(query_points::Array{T,2}, markov_chains, nn_trees, field_names, priors;
    credible_level = 0.9, q_quantiles = [], mode_intervals = zeros(length(field_names))) where {T}

    # Allocations
    num_fields, num_out, num_chains = size(nn_trees) # Number of fields, samples output per chain, and chains
    num_samples, num_pts = num_out*num_chains, size(query_points, 2) # Number of posterior samples and interpolation points
    num_quantiles = length(q_quantiles) # Number of quantiles
    post_samples = zeros(num_fields, num_samples) # Container for posterior samples
    nuclei_dist = zeros(num_fields, num_samples) # Distances to Voronoi cells
    PostTend = (
        dist = zeros(num_fields, num_pts),
        mean = zeros(num_fields, num_pts),
        median = zeros(num_fields, num_pts),
        mode = zeros(num_fields, num_pts),
        quantiles = zeros(num_quantiles, num_fields, num_pts),
        median_interval = zeros(2, num_fields, num_pts),
        low_high_interval = zeros(2, num_fields, num_pts),
        highest_density_interval = zeros(2, num_fields, num_pts),
        h_distance = zeros(num_fields, num_pts),
        k_metric = zeros(num_fields, num_pts)
    )

    # Compute tendencies
    percent_done, print_interval = 0.1, 0.1 # For printing progress
    for j in 1:num_pts # Loop over interpolation points
        # Collect posterior samples 
        x_j = @view query_points[:, j]
        collect_posterior_samples!(post_samples, nuclei_dist, x_j, markov_chains, nn_trees, field_names)
        # Sort samples for interval caluclations (increasing order)
        sort!(post_samples, dims = 2)
        for i in eachindex(field_names) # Loop over fields to compute metrics
            d_i = @view nuclei_dist[i,:]
            PostTend.dist[i,j] = mean(d_i)

            s_i = @view post_samples[i,:]
            PostTend.mean[i,j] = mean(s_i)
            PostTend.median[i,j] = quantile(s_i, 0.5; sorted = true)
            for (h, q_h) in enumerate(q_quantiles)
                PostTend.quantiles[h,i,j] = quantile(s_i, q_h; sorted = true)
            end
            PostTend.median_interval[:,i,j] .= median_interval(s_i, credible_level; is_sorted = true)
            PostTend.low_high_interval[:,i,j] .= low_high_interval(s_i, credible_level; is_sorted = true)
            PostTend.highest_density_interval[:,i,j] .= highest_density_interval(s_i, credible_level; is_sorted = true)
            # Quantify how much the posterior differs from the prior
            PostTend.h_distance[i,j] = hellinger_distance!(s_i, priors[i]...; is_sorted = true)
            # Needs to updated for non-uniform priors...also I don't use this metric
            # PostTend.k_metric[i,j] = kolmogorov_metric!(s_i, priors[i]...; is_sorted = true)

            # Compute last because sample values are rounded in-place
            PostTend.mode[i,j] = binned_mode!(s_i; mode_interval = mode_intervals[i], is_sorted = true)
        end
        # Print progress
        if (j/num_pts) >= percent_done
            println("Posterior tendency calculation is "*string(round(Int, 100.0*percent_done))*"% complete...")
            percent_done += print_interval
        end
    end

    return PostTend
end

function posterior_directional_tendencies(query_points::Array{T,2}, markov_chains, nn_trees, field_names;
    credible_level = 0.9, q_quantiles = [], mode_intervals = zeros(length(field_names)),
    spherical_conversion = (a,b,c) -> (a,b,c), tf_global = true) where {T}
    # Allocations
    num_fields, num_out, num_chains = size(nn_trees) # Number of fields, samples output per chain, and chains
    num_samples, num_pts = num_out*num_chains, size(query_points, 2) # Number of posterior samples and interpolation points
    num_quantiles = length(q_quantiles) # Number of quantiles
    post_samples = zeros(num_fields, num_samples) # Container for posterior samples
    nuclei_dist = zeros(num_fields, num_samples) # Distances to Voronoi cells
    PostTend = (
        dist = zeros(num_fields, num_pts),
        directional_vectors = zeros(3, 3, num_pts),
        mean = zeros(num_fields, num_pts),
        median = zeros(num_fields, num_pts),
        mode = zeros(num_fields, num_pts),
        quantiles = zeros(num_quantiles, num_fields, num_pts),
        median_interval = zeros(2, num_fields, num_pts),
        low_high_interval = zeros(2, num_fields, num_pts),
        highest_density_interval = zeros(2, num_fields, num_pts)
    )

    # Compute tendencies
    percent_done, print_interval = 0.1, 0.1 # For printing progress
    for j in 1:num_pts # Loop over interpolation points
        # Collect posterior samples 
        x_j = @view query_points[:, j]
        collect_posterior_samples!(post_samples, nuclei_dist, x_j, markov_chains, nn_trees, field_names)
        # Convert to directional posterior
        Tj, dbj, vj, wj = directional_posterior!(post_samples; spherical_conversion = spherical_conversion)
        PostTend.directional_vectors[:,1,j] .= vj[1]*wj[:,1]
        PostTend.directional_vectors[:,2,j] .= vj[2]*wj[:,2]
        PostTend.directional_vectors[:,3,j] .= vj[3]*wj[:,3]
        # Sort samples for interval caluclations (increasing order)
        sort!(post_samples, dims = 2)
        for i in eachindex(field_names) # Loop over fields to compute metrics
            d_i = @view nuclei_dist[i,:]
            PostTend.dist[i,j] = mean(d_i)

            s_i = @view post_samples[i,:]
            PostTend.mean[i,j] = mean(s_i)
            PostTend.median[i,j] = quantile(s_i, 0.5; sorted = true)
            for (h, q_h) in enumerate(q_quantiles)
                PostTend.quantiles[h,i,j] = quantile(s_i, q_h; sorted = true)
            end
            PostTend.median_interval[:,i,j] .= median_interval(s_i, credible_level; is_sorted = true)
            PostTend.low_high_interval[:,i,j] .= low_high_interval(s_i, credible_level; is_sorted = true)
            PostTend.highest_density_interval[:,i,j] .= highest_density_interval(s_i, credible_level; is_sorted = true)
            # Compute last because sample values are rounded in-place
            PostTend.mode[i,j] = binned_mode!(s_i; mode_interval = mode_intervals[i], is_sorted = true)
        end
        # Print progress
        if (j/num_pts) >= percent_done
            println("Posterior tendency calculation is "*string(round(Int, 100.0*percent_done))*"% complete...")
            percent_done += print_interval
        end
    end

    # Convert vectors from local geographic to global cartesian
    if tf_global
        v1 = @view PostTend.directional_vectors[:,1,:]
        v2 = @view PostTend.directional_vectors[:,2,:]
        v3 = @view PostTend.directional_vectors[:,3,:]
        local_to_global_vectors!(v1, v2, v3, query_points)
    end

    return PostTend
end

function directional_posterior!(posterior_samples; spherical_conversion = (a,b,c) -> (a,b,c))
    # Compute directional tensor
    T = zeros(3,3)
    for j in axes(posterior_samples, 2)
        a_j, b_j, c_j = posterior_samples[1,j], posterior_samples[2,j], posterior_samples[3,j]
        f_j, azm_j, elv_j = spherical_conversion(a_j, b_j, c_j)
        sum_directional_tensor!(T, f_j, azm_j, elv_j)
        posterior_samples[1,j], posterior_samples[2,j], posterior_samples[3,j] = f_j, azm_j, elv_j
    end
    # Compute principal axes sorted by *increasing magnitude*
    T ./= size(posterior_samples, 2)
    val, vec = eigen(T; sortby = abs)
    # Compute directional bias (fractional anisotropy of diffusion; see Wikipedia entry for fractional anisotropy)
    λ1, λ2, λ3 = val[3], val[2], val[1]
    directional_bias = sqrt(0.5)*sqrt( ((λ1 - λ2)^2) + ((λ1 - λ3)^2) + ((λ2 - λ3)^2) )/sqrt( (λ1^2) + (λ2^2) + (λ3^2) )

    # Distribution of strength and orientation
    v1, v2, v3 = val[3]*vec[1,3], val[3]*vec[2,3], val[3]*vec[3,3]
    v_norm = abs(val[3])
    for j in axes(posterior_samples, 2)
        # Directional vector
        f_j, azm_j, elv_j = posterior_samples[1,j], posterior_samples[2,j], posterior_samples[3,j]
        sinλ, cosλ = sincos(azm_j)
        sinϕ, cosϕ = sincos(elv_j)
        u1, u2, u3 = f_j*cosϕ*cosλ, f_j*cosϕ*sinλ, f_j*sinϕ
        # Metrics
        u_norm = sqrt((u1^2) + (u2^2) + (u3^2))
        prj, puv = abs(u1*v1 + u2*v2 + u3*v3), u_norm*v_norm # Projection of sampled vectors onto principal orientation <--- Compute dlnVp correlation with this field???
        cos2Δ = puv > 0.0 ? 2.0*((prj/puv)^2) - 1.0 : 1.0 # Angular deviation with principal orientation; cos2Δ = 2.0*(cosΔ^2) - 1.0
        prj = sqrt(prj) # Projected magnitude
        # Store result
        posterior_samples[2,j], posterior_samples[3,j] = prj, 0.5*acos(cos2Δ)
    end

    return T, directional_bias, val, vec
end

# --- Scalar Posterior Moments --- #

# Point-wise posterior moments
function posterior_moments(query_points::Array{T,2}, markov_chains, nn_trees, field_names; correlation_pairs = [], num_moments = 2) where {T}
    # Allocations
    num_pairs = length(correlation_pairs) # Number of correlation coefficients to compute
    num_fields, num_out, num_chains = size(nn_trees) # Number of fields, samples output per chain, and chains
    num_samples, num_pts = num_out*num_chains, size(query_points, 2) # Number of posterior samples and interpolation points
    num_moments = (num_moments < 2) && (num_pairs > 0) ? 2 : num_moments # Number of moments (correlation calculation requires at least 2)
    post_samples = zeros(num_fields, num_samples) # Container for posterior samples
    nuclei_dist = zeros(num_fields, num_samples) # Distances to Voronoi cells...not currently used
    post_moments = zeros(num_moments, num_pts, num_fields) # Container for posterior moments
    post_correlation = zeros(1, num_pts, num_pairs) # Container for correlation products

    # Index field names for correlation coefficient calculation
    D = Dict{String, Int}()
    [D[key] = i for (i, key) in enumerate(field_names)]

    # Compute statistics
    for j in 1:num_pts # Loop over interpolation points
        # Collect posterior samples 
        x_j = @view query_points[:, j]
        collect_posterior_samples!(post_samples, nuclei_dist, x_j, markov_chains, nn_trees, field_names)
        # Moment sums
        for i in 1:num_moments # Loop over moments
            post_moments[i,j,:] .+= sum(x -> x^i, post_samples; dims = 2)
        end
        # Correlation sums
        for (k, pair_k) in enumerate(correlation_pairs) # Loop over correlation pairs
            i1, i2 = D[pair_k[1]], D[pair_k[2]]
            @views v1, v2 = post_samples[i1,:], post_samples[i2,:]
            post_correlation[1,j,k] += dot(v1, v2)
        end
    end
    # Finish moment calculations
    fill_field_moments!(post_moments, num_samples)
    if num_pairs > 0
        fill_field_pcc!(post_correlation, post_moments, field_names, correlation_pairs, num_samples)
    end

    return post_moments, post_correlation
end

# Model-wise posterior moments
function posterior_moments(num_chains::Int, num_burn::Int, query_points::Array{T,2}, field_names; correlation_pairs = [], num_moments = 2,
    chain_directory = pwd(), chain_name = "chain", tf_squeeze = false, tf_cart_to_geo = true) where {T}

    # Allocations
    num_pts = size(query_points, 2) # Number of interpolation points
    num_fields = length(field_names) # Number of parameter fields
    num_pairs = length(correlation_pairs) # Number of correlation coefficients to compute
    num_moments = (num_moments < 2) && (num_pairs > 0) ? 2 : num_moments # Correlation computation requires second moment
    post_moments = zeros(num_moments, num_pts, num_fields) # Sample sums
    post_correlation = zeros(1, num_pts, num_pairs) # Sample correlation products

    # Accumulate Chains
    num_models = 0
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        num_models += length(chain_n) # Number of models in chain
        for chain_n_i in chain_n # Loop over iterations in chain
            accumulate_moments!(post_moments, query_points, field_names, chain_n_i; tf_squeeze = tf_squeeze, tf_cart_to_geo = tf_cart_to_geo)
            if num_pairs > 0 # Some redundant interpolations here
                accumulate_correlation_products!(post_correlation, query_points, correlation_pairs, chain_n_i;
                tf_squeeze = tf_squeeze, tf_cart_to_geo = tf_cart_to_geo)
            end
        end
        println("Finished moment accumulation for chain "*string(n)*" of "*string(num_chains)*"...")
    end
    num_models == 0 && error("Cannot build ensemble. Empty chains!")
    # Finish moment calculations
    fill_field_moments!(post_moments, num_models)
    if num_pairs > 0
        fill_field_pcc!(post_correlation, post_moments, field_names, correlation_pairs, num_models)
    end
    
    return post_moments, post_correlation
end

function accumulate_moments!(post_moments, query_points, field_names, Model;
    tf_squeeze = false, tf_cart_to_geo = true)

    tf_in_bounds = true # Initialize in-bounds flag
    for (k, field_k) in enumerate(field_names) # Loop over parameter fields
        # Interpolate posterior
        fid = Model.fieldslist[field_k]
        Tree = return_voronoi_tree(Model.fields[fid])
        Threads.@threads for j = axes(query_points, 2) # Parallel loop over query points
            if tf_squeeze
                tf_in_bounds = check_field_bounds(query_points[1, j], query_points[2, j], query_points[3, j], 
                Model.fields[fid].slims; tf_cart_to_geo = tf_cart_to_geo)
            end
            if tf_in_bounds
                ind, dist = nn(Tree, @view query_points[:, j])
                val = Model.fields[fid].v[ind]
            else
                val = Model.fields[fid].ref_value
            end
            for i in axes(post_moments, 1) # Loop over moments
                post_moments[i, j, k] += (val^i)
                # coord_moments[i, j, k] += (dist^i) # Also compute positional statistics
            end
        end
    end

    return nothing
end

function accumulate_correlation_products!(post_correlation, query_points, correlation_pairs, Model;
    tf_squeeze = false, tf_cart_to_geo = true)
    # Compute correlation products for pairs of fields
    tf_in_bounds_1, tf_in_bounds_2 = true, true # Initialize in-bounds flag
    cj, rj = [0], [0.0] # Allocate arrays to store nearest neighbor index and distance for small performance boost
    for (k, pair_k) in enumerate(correlation_pairs) # Loop over field pairs
        fid_1, fid_2 = Model.fieldslist[pair_k[1]], Model.fieldslist[pair_k[2]]
        Tree_1, Tree_2 = return_voronoi_tree(Model.fields[fid_1]), return_voronoi_tree(Model.fields[fid_2])
        for j = axes(query_points, 2) # Loop over query points
            xj = @view query_points[:, j]
            if tf_squeeze # Check if query_point is in field domain
                tf_in_bounds_1 = check_field_bounds(xj[1], xj[2], xj[3], Model.fields[fid_1].slims; tf_cart_to_geo = tf_cart_to_geo)
                tf_in_bounds_2 = check_field_bounds(xj[1], xj[2], xj[3], Model.fields[fid_2].slims; tf_cart_to_geo = tf_cart_to_geo)
            end
            val_1 = tf_in_bounds_1 ? nearest_neighbor_value(xj, Model.fields[fid_1], Tree_1; icell = cj, dist = rj) : Model.fields[fid_1].ref_value
            val_2 = tf_in_bounds_2 ? nearest_neighbor_value(xj, Model.fields[fid_2], Tree_2; icell = cj, dist = rj) : Model.fields[fid_2].ref_value
            if tf_in_bounds_1 && tf_in_bounds_2 # Compute correlation where field domains overlap
                post_correlation[1,j,k] += val_1 * val_2
            end
        end
    end

    return nothing
end

# --- Directional Posterior Moments --- #
# Correlations are most sensible between scalar fields and vector components of directional fields <--- Implement this!

function posterior_directional_moments(num_chains, num_burn, query_points, directional_fields; correlation_fields = [], num_moments = 2,
    chain_directory = pwd(), chain_name = "chain", tf_squeeze = false, tf_cart_to_geo = true, tf_global = true, spherical_conversion = (a,b,c) -> (a,b,c))
    
    # Allocations
    num_pts = size(query_points, 2) # Number of interpolation points
    num_pairs = length(correlation_fields) # Number of correlation coefficients to compute
    num_moments = (num_moments < 2) && (num_pairs > 0) ? 2 : num_moments # Correlation computation requires second moments
    post_vectors = zeros(3, 3, num_pts) # Principal vectors (3-components, 3-vectors)
    post_moments = zeros(num_moments, num_pts, 3 + num_pairs) # Moments for (1) fraction, (2) projected magnitude, (3) angular deviation, and correlation fields
    post_correlation = zeros(1, num_pts, num_pairs) # Sample correlation products with projected magnitude
    directional_strength = zeros(1, num_pts) # Container for directional strength (i.e. projected magnitude) calculations...equivalent to post_moments[1,:,2]?
    angular_weights = zeros(num_pts) # For weighted moment calculation

    # Build directional tensor
    num_models = 0
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        num_models += length(chain_n) # Number of models in chain
        for chain_n_i in chain_n # Loop over iterations in chain
            accumulate_directional_tensor!(post_vectors, query_points, directional_fields, chain_n_i; spherical_conversion = spherical_conversion)
        end
        println("Finished tensor accumulation for chain "*string(n)*" of "*string(num_chains)*"...")
    end
    num_models == 0 && error("Cannot build ensemble. Empty chains!")
    post_vectors ./= num_models

    # Convert tensors to principal eigenvectors scaled by eigenvalues
    directional_bias = principal_directions!(post_vectors)

    # Compute directional moments
    major_axis = @view post_vectors[:,3,:]
    @views directional_moments, scalar_moments = post_moments[:,:,1:3], post_moments[:,:,4:end]
    for n in 1:num_chains # Loop over chains
        chain_n = load(chain_directory*"/"*chain_name*"."*string(n)*".jld", "MarkovChain")
        splice!(chain_n, 1:num_burn) # Remove burn-in
        for chain_n_i in chain_n # Loop over iterations in chain
            accumulate_directional_moments!(directional_moments, directional_strength, angular_weights, major_axis, query_points, directional_fields, chain_n_i;
            spherical_conversion = spherical_conversion)
            if num_pairs > 0 # Compute correlations with directional strength
                accumulate_directional_correlation!(scalar_moments, post_correlation, directional_strength, query_points, correlation_fields, chain_n_i)
            end
        end
        println("Finished moment accumulation for chain "*string(n)*" of "*string(num_chains)*"...")
    end
    # Finish moment calculations -- with angular weights
    angular_moments = @view post_moments[:,:,3]
    [angular_moments[:,j] ./= w_j for (j, w_j) in enumerate(angular_weights)]
    fill_field_moments!(angular_moments, 1)
    # Remaining non-weighted moments
    scalar_moments = @view post_moments[:,:,1:2]
    fill_field_moments!(scalar_moments, num_models)
    scalar_moments = @view post_moments[:,:,4:end]
    fill_field_moments!(scalar_moments, num_models)

    # Finish moment calculations
    # fill_field_moments!(post_moments, num_models)
    if num_pairs > 0
        # Update list of field names
        field_names = Vector{String}(undef, 3 + num_pairs)
        field_names[1], field_names[2], field_names[3] = "a", "directional_strength", "angular_deviation"
        [field_names[i+3] = vi for (i, vi) in enumerate(correlation_fields)]
        # Define correlation pairs with directional strength
        correlation_pairs = Vector{NTuple{2, String}}(undef, num_pairs)
        [correlation_pairs[i] = (field_names[2], vi) for (i, vi) in enumerate(correlation_fields)]
        fill_field_pcc!(post_correlation, post_moments, field_names, correlation_pairs, num_models)
    end

    # Convert vectors from local geographic to global cartesian
    if tf_global
        v1 = @view post_vectors[:,1,:]
        v2 = @view post_vectors[:,2,:]
        v3 = @view post_vectors[:,3,:]
        local_to_global_vectors!(v1, v2, v3, query_points)
    end

    return post_moments, post_correlation, post_vectors, directional_bias
end

function accumulate_directional_tensor!(tensor_sum, query_points, directional_fields, Model; spherical_conversion = (a,b,c) -> (a,b,c))
    # Index directional model fields
    k_1 = Model.fieldslist[directional_fields[1]]
    k_2 = Model.fieldslist[directional_fields[2]]
    k_3 = Model.fieldslist[directional_fields[3]]
    # Nearest neighbor trees
    Tree_1 = return_voronoi_tree(Model.fields[k_1])
    Tree_2 = return_voronoi_tree(Model.fields[k_2])
    Tree_3 = return_voronoi_tree(Model.fields[k_3])
    # Accumulate tensor elements
    cj, rj = [0], [0.0] # Allocate arrays to store nearest neighbor index and distance for small performance boost
    for j = axes(query_points, 2)
        xj = @view query_points[:, j]
        # Convert directional parameters to spherical
        a_j = nearest_neighbor_value(xj, Model.fields[k_1], Tree_1; icell = cj, dist = rj)
        b_j = nearest_neighbor_value(xj, Model.fields[k_2], Tree_2; icell = cj, dist = rj)
        c_j = nearest_neighbor_value(xj, Model.fields[k_3], Tree_3; icell = cj, dist = rj)
        f_j, azm_j, elv_j = spherical_conversion(a_j, b_j, c_j)
        T = @view tensor_sum[:, :, j]
        sum_directional_tensor!(T, f_j, azm_j, elv_j)
    end

    return nothing
end

function accumulate_directional_moments!(post_moments, directional_strength, angular_weights, major_axis, query_points, directional_fields, Model;
    spherical_conversion = (a,b,c) -> (a,b,c))

    # Index directional model fields
    k_1 = Model.fieldslist[directional_fields[1]]
    k_2 = Model.fieldslist[directional_fields[2]]
    k_3 = Model.fieldslist[directional_fields[3]]
    # Nearest neighbor trees
    Tree_1 = return_voronoi_tree(Model.fields[k_1])
    Tree_2 = return_voronoi_tree(Model.fields[k_2])
    Tree_3 = return_voronoi_tree(Model.fields[k_3])
    # Accumulate elements
    cj, rj = [0], [0.0] # Allocate arrays to store nearest neighbor index and distance for small performance boost
    for j = axes(query_points, 2)
        xj = @view query_points[:, j]
        # Convert directional parameters to spherical
        a_j = nearest_neighbor_value(xj, Model.fields[k_1], Tree_1; icell = cj, dist = rj)
        b_j = nearest_neighbor_value(xj, Model.fields[k_2], Tree_2; icell = cj, dist = rj)
        c_j = nearest_neighbor_value(xj, Model.fields[k_3], Tree_3; icell = cj, dist = rj)
        f_j, azm_j, elv_j = spherical_conversion(a_j, b_j, c_j)
        # Directional vectors
        sinλ, cosλ = sincos(azm_j)
        sinϕ, cosϕ = sincos(elv_j)
        u1, u2, u3 = f_j*cosϕ*cosλ, f_j*cosϕ*sinλ, f_j*sinϕ
        v1, v2, v3 = major_axis[1,j], major_axis[2,j], major_axis[3,j]
        # Metrics
        u0, v0 = sqrt((u1^2) + (u2^2) + (u3^2)), sqrt((v1^2) + (v2^2) + (v3^2))
        prj, puv0 = abs(u1*v1 + u2*v2 + u3*v3), u0*v0 # Unsigned dot-product
        cos2Δ = puv0 > 0.0 ? 2.0*((prj/puv0)^2) - 1.0 : 1.0 # Angular deviation with principal orientation; cos2Δ = 2.0*(cosΔ^2) - 1.0
        cos2Δ = min(cos2Δ, 1.0) # Avoid floating point errors that cause subsequent acos call to fail
        prj = sqrt(prj) # Projected magnitude
        # Accumulate result
        for i in axes(post_moments, 1)
            post_moments[i,j,1] += (f_j^i)
            post_moments[i,j,2] += (prj^i)
            # post_moments[i,j,3] += ((0.5*acos(cos2Δ))^i)
            # Consider weighting angular deviation by magnitude
            post_moments[i,j,3] += u0*((0.5*acos(cos2Δ))^i)
        end
        angular_weights[j] += u0 # Need also the sum of angular weights for each interpolation point
        directional_strength[1,j] = prj # Directional correlation with this field
    end

    return nothing
end

function accumulate_directional_correlation!(post_moments, post_correlation, directional_strength, query_points, correlation_fields, Model)
    cj, rj = [0], [0.0] # Allocate arrays to store nearest neighbor index and distance for small performance boost
    for (k, field_k) in enumerate(correlation_fields)
        fid = Model.fieldslist[field_k]
        Tree = return_voronoi_tree(Model.fields[fid])
        for j = axes(query_points, 2)
            x_j = @view query_points[:, j]
            v_j = nearest_neighbor_value(x_j, Model.fields[fid], Tree; icell = cj, dist = rj)
            for i in axes(post_moments, 1)
                post_moments[i, j, k] += (v_j^i)
            end
            post_correlation[1, j, k] += v_j*directional_strength[1,j]
        end
    end
    return nothing
end

# --- Statistical Metrics --- #

# Compute moments from sample sums
function fill_field_moments!(post_moments, num_samples)
    num_moments = size(post_moments, 1)
    # Mean
    post_moments ./= num_samples
    # Standard Deviation
    if num_moments > 1
        post_moments[2,:,:] .-= (post_moments[1,:,:].^2)
        post_moments[2,:,:] .= sqrt.(post_moments[2,:,:])
    end
    # Skewness
    if num_moments > 2
        post_moments[3,:,:] .-= (3.0*post_moments[1,:,:].*(post_moments[2,:,:].^2) .+ (post_moments[1,:,:].^3))
        post_moments[3,:,:] ./= (post_moments[2,:,:].^3)
    end
    if num_moments > 3
        @warn "Moments of order > 3 not defined!"
    end
    return nothing
end
function fill_field_pcc!(post_correlation, post_moments, field_names, correlation_pairs, num_samples)
    # Index field names
    D = Dict{String, Int}()
    [D[key] = i for (i, key) in enumerate(field_names)]
    # Linear correlation: (μ₁₂ - μ₁μ₂)/(σ₁σ₂)
    post_correlation ./= num_samples
    for (k, pair_k) in enumerate(correlation_pairs)
        k1, k2 = D[pair_k[1]], D[pair_k[2]]
        @views @. post_correlation[1,:,k] -= post_moments[1,:,k1] * post_moments[1,:,k2]
        @views @. post_correlation[1,:,k] /= post_moments[2,:,k1] * post_moments[2,:,k2]
    end

    return nothing
end

function binned_mode!(xi; mode_interval = 0.0, is_sorted = false)
    # Compute interval using (Freedman-Diaconis rule)
    if mode_interval <= 0.0
        vmin, vmax = median_interval(xi, 0.5; is_sorted = is_sorted)
        mode_interval = 2.0*(vmax - vmin)/(length(xi)^(1/3))
    end
    # Compute binned mode
    xi ./= mode_interval
    xi .= round.(xi)
    return mode_interval*mode(xi)
end
# Value at which some percentage of samples are below
function lowest_interval(samples, credible_level; is_sorted = false)
    return quantile(samples, credible_level; sorted = is_sorted)
end
# Values at which some percentage of samples are above
function highest_interval(samples, credible_level; is_sorted = false)
    return quantile(samples, 1.0 - credible_level; sorted = is_sorted)
end
function low_high_interval(samples, credible_level; is_sorted = false)
    return quantile(samples, (credible_level, 1.0 - credible_level); sorted = is_sorted)
end
# Interval with equal probability of being below as above (aka equal-tailed interval)
function median_interval(samples, credible_level; is_sorted = false)
    α = 1.0 - credible_level
    interval = (0.5*α, 1.0 - 0.5*α)
    return quantile(samples, interval; sorted = is_sorted) # lower, upper interval boundaries
end
# Smallest interval containing some percentage of samples...see PosteriorStats package for fancier options
function highest_density_interval(samples, credible_level; is_sorted = false)
    n_samp = length(samples)
    n_ci = floor(Int, credible_level*n_samp)
    n_intervals = 1 + n_samp - n_ci

    sorted_samples = is_sorted ? samples : sort(samples)
    min_width, hdi_min, hdi_max = Inf, -Inf, Inf
    for i in 1:n_intervals
        low = sorted_samples[i]
        high = sorted_samples[i + n_ci - 1]
        width = high - low
        if width < min_width
            min_width = width
            hdi_min = low
            hdi_max = high
        end
    end

    return hdi_min, hdi_max
end
# Probability inside (or outside) a specified interval
function probability_interval(samples, interval)
    low, high = interval
    w = 1.0/length(samples)
    f_in = w*sum(x -> (x > low) && (x < high), samples)
    return f_in
end
# Equiavlent to above with low = -Inf
function probability_tail(samples, val)
    w = 1.0/length(samples)
    f_in = w*sum(x -> (x < val), samples)
    return f_in
end

# Computes the largest deviation between the CDFs of a sampled distribution
# and a uniform distribution...needs to be updated for other priors (see hellinger_distance!)
function kolmogorov_metric!(samples, a, b; is_sorted = false)
    !is_sorted && sort!(samples)
    n, dba, kmax = length(samples), b - a, 0.0
    for (i, s_i) in enumerate(samples)
        p_il, p_ir = (i-1)/n, i/n
        u_i = (s_i - a)/dba
        kmax = max(kmax, abs(p_il - u_i), abs(p_ir - u_i))
    end
    return kmax
end

# Hellinger distance assuming a uniform prior
function hellinger_distance!(samples, a, b; is_sorted = false)
    # Sort samples if necessary
    !is_sorted && sort!(samples)

    # Compute bin interval using -- simple √n rule
    ns = length(samples)
    dx = (b - a)/round(Int, sqrt(ns))

    # Bin probability -- uniform
    u_0 = dx/(b-a)

    # Sampled probability
    p_sum, bin_max, ni = 0.0, a + dx, 0
    for s_i in samples
        while s_i > bin_max
            p_sum += sqrt(ni/ns)
            bin_max += dx # Shift bin limit
            ni = 0 # Reset bin counter
        end
        ni += 1
    end
    p_sum += sqrt(ni/ns) # Last bin

    H2 = 1.0 - sqrt(u_0)*p_sum
    H2 = max(0.0, min(H2, 1.0)) # Avoid rounding errors
    return sqrt(H2)
end
# Hellinger distance given an arbitrary prior cdf function
# a, b, σ, n = -0.2, 0.2, 0.05, 1000
# samples_u = a .+ (b-a)*rand(n)
# samples_n = σ*randn(n)
# cdf_u = (x; v_min = a, v_max = b) -> min(1.0, max(0.0, (x-v_min)/(v_max-v_min)))
# cdf_n = (x; u = 0.0, sdev = σ) -> 0.5 + 0.5*erf((x-u)/(sqrt(2.0)*sdev))
# dx = (b-a)/round(Int, sqrt(length(samples)))
function hellinger_distance!(samples, dx::Float64, cdf::Function; is_sorted = false)
    # Bin probabilities
    # Uniform: (x; v_min = -0.1, v_max = 0.1) -> min(1.0, max(0.0, (x-v_min)/(v_max-v_min)))
    # Normal: (x; u = 0.0, sdev = 0.1) -> 0.5 + 0.5*erf((x-u)/(sqrt(2.0)*sdev))
    # Half-normal: (x; sdev = 0.1) -> max(0.0, erf(x/(sqrt(2.0)*sdev)))

    # Sort samples if necessary
    !is_sorted && sort!(samples)

    # Determine bin interval (Freedman-Diaconis rule)
    if dx <= 0.0
        xmin, xmax = median_interval(samples, 0.5; is_sorted = true)
        dx = 2.0*(xmax - xmin)/(length(samples)^(1/3))
    end

    # Sampled probability
    a, b = samples[1], samples[1] + dx
    p_sum, ni, ns = 0.0, 0, length(samples)
    for s_i in samples
        while s_i > b # Update counters
            p_ab = cdf(b) - cdf(a)
            p_sum += sqrt((ni/ns)*p_ab)
            a += dx
            b += dx
            ni = 0
        end
        ni += 1
    end
    # Integrate last bin
    p_ab = cdf(b) - cdf(a)
    p_sum += sqrt((ni/ns)*p_ab)

    H2 = 1.0 - p_sum
    H2 = max(0.0, min(H2, 1.0)) # Avoid rounding errors
    return sqrt(H2)
end
# Analytic Hellinger distance between a uniform on the interval (a,b) and a
# normal distribution with mean = 0.0 and standard deviation = σ.
function hellinger_uniform_v_normal(a,b,σ)
   f = (x,a,b,σ) -> ((b-a)^(-0.5))*((2.0*pi*(σ^2))^(-0.25))*sqrt(pi)*σ*erf(0.5*x/σ)
   return sqrt(1.0 - f(b,a,b,σ) + f(a,a,b,σ))
end


######################
# OLD ENSEMBLE MODEL #
######################


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
    chain_directory = pwd(), chain_name = "chain", vtk_file = chain_directory * "/EnsembleModel",
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



# ASSUMPTIONS
# + All chains contain same set of fields (order can differ)
# + Chain file naming format is, chain_name.number.jld
# Re-name...posterior_moments
function build_ensemble_model(num_chains, num_burn, grid_points, field_names;
    chain_directory = pwd(), chain_name = "chain", tf_squeeze = false)
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
    num_models == 0 && error("Cannot build ensemble. Empty chains!")
    # Finish statistics calculations
    mean_model ./= num_models
    sdev_model ./= num_models
    sdev_model .-= (mean_model.^2)
    sdev_model .= sqrt.(sdev_model)
    # # Skewness (where skew_model = sum(v^3))
    # skew_model ./= num_models
    # skew_model .-= (3.0*mean_model.*(sdev_model.^2) .+ (mean_model.^3))
    # skew_model ./= (sdev_model.^3)

    return mean_model, sdev_model
end
# Re-name...posterior_correlations. Update to accept vector of field pairs
function build_ensemble_correlation(num_chains, num_burn, grid_points, field_a, field_b;
    chain_directory = pwd(), chain_name = "chain", tf_squeeze = false)
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
    num_models == 0 && error("Cannot build ensemble. Empty chains!")
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
# Re-name...posterior_directional_moments. Update to extract angular deviations and projected magnitudes
function build_ensemble_tensor(num_chains, num_burn, grid_points, directional_fields;
    spherical_conversion = (a,b,c) -> (a,b,c), tf_local_to_global = false,
    chain_directory = pwd(), chain_name = "chain")
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
    num_models == 0 && error("Cannot build ensemble. Empty chains!")
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

#############
# UTILITIES #
#############

function principal_directions!(T)
    # Allocations
    num_pts, root_half = size(T,3), sqrt(0.5)
    directional_bias = zeros(1, num_pts)
    # Solve eigenvalue problem
    Threads.@threads for j in 1:num_pts
        T_j = @view T[:,:,j]
        λ, v = eigen(T_j; sortby = abs)
        T[:,1,j] .= λ[1]*v[:,1]
        T[:,2,j] .= λ[2]*v[:,2]
        T[:,3,j] .= λ[3]*v[:,3]
        # Compute directional bias (fractional anisotropy of diffusion; see Wikipedia entry for fractional anisotropy)
        directional_bias[1,j] = root_half*sqrt( ((λ[1] - λ[2])^2) + ((λ[1] - λ[3])^2) + ((λ[2] - λ[3])^2) )
        directional_bias[1,j] /= sqrt( (λ[1]^2) + (λ[2]^2) + (λ[3]^2) )
    end

    return directional_bias
end

function nearest_voronoi_cell(query_points, Field::Voronoi; leafsize = 10)
    nuclei = @view Field.c[:,1:Field.n[1]] # View into valid nuclei
    Tree = KDTree(nuclei; leafsize = leafsize)
    index, _ = nn(Tree, query_points) # Very slow serially with lots of allocations
    # Field.v[index] # Nearest value
    return index
end

# Fastest interpolation strategy using NearestNeighbors
# Also avoids allocating large arrays for storing indices and distances
function interpolate_voronoi_field!(v, query_points, Field::Voronoi; leafsize = 10, icell = [0], dist = [0.0])
    Tree = return_voronoi_tree(Field; leafsize = leafsize)
    interpolate_voronoi_field!(v, query_points, Field, Tree; leafsize = leafsize, icell = icell, dist = dist)
    return nothing
end
function interpolate_voronoi_field!(v, query_points, Field::Voronoi, Tree; leafsize = 10, icell = [0], dist = [0.0])
    for j in axes(query_points, 2)
        xj = @view query_points[:,j]
        knn!(icell, dist, Tree, xj, 1)
        v[j] = Field.v[icell[1]]
    end
    return nothing
end
function return_voronoi_tree(Field::Voronoi; leafsize = 10)
    nuclei = @view Field.c[:,1:Field.n[1]] # View into valid nuclei
    return KDTree(nuclei; leafsize = leafsize)
end
function nearest_neighbor_value(query_point, Field::Voronoi, Tree; icell = [0], dist = [0.0])
    knn!(icell, dist, Tree, query_point, 1)
    return Field.v[icell[1]]
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
        R = extrinsic_rotation_matrix((-ϕ, λ), (2, 3))
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

    v = F.vectors[:, p1]
    vl = sqrt((v[1]^2) + (v[2]^2) + (v[3]^2))
    prj, cosΔ = zeros(size(a)), zeros(size(a))
    for i in eachindex(a)
        r_i, λ_i, ϕ_i = spherical_conversion(a[i], b[i], c[i])
        sinλ, cosλ = sincos(λ_i)
        sinϕ, cosϕ = sincos(ϕ_i)
        u1, u2, u3 = r_i*cosϕ*cosλ, r_i*cosϕ*sinλ, r_i*sinϕ
        prj[i] = u1*v[1] + u2*v[2] + u3*v[3]

        ul = sqrt((u1^2) + (u2^2) + (u3^2))
        cosΔ[i] = prj[i]/(ul*vl) # cos2Δ = 2.0*(cosΔ^2) - 1.0
    end

    return db, λ1, λ2, λ3, T, prj, cosΔ
end



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
    R = extrinsic_rotation_matrix((-c*latitude, c*longitude), (2, 3))
    sx = R[1,1]*w[1] + R[1,2]*w[2] + R[1,3]*w[3]
    sy = R[2,1]*w[1] + R[2,2]*w[2] + R[2,3]*w[3]
    sz = R[3,1]*w[1] + R[3,2]*w[2] + R[3,3]*w[3]

    return sx, sy, sz
end

function extrinsic_rotation_matrix(α, n)
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

function extrinsic_rotation_matrix(α::Tuple, n::Tuple)
    R = @SMatrix [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    for i in eachindex(n)
        R = extrinsic_rotation_matrix(α[i],n[i])*R
    end
    return R
end

function check_field_bounds(xj, yj, zj, field_lims; tf_cart_to_geo = true)
    if tf_cart_to_geo
        λj, ϕj, rj = atan(yj, xj), atan(zj, sqrt((xj^2) + (yj^2))), sqrt((xj^2) + (yj^2) + (zj^2))
    else
        λj, ϕj, rj = xj, yj, zj
    end
    tf_in_bounds = (ϕj >= field_lims[1][1]) && (ϕj <= field_lims[1][2]) &&
                   (λj >= field_lims[2][1]) && (λj <= field_lims[2][2]) &&
                   (rj >= field_lims[3][1]) && (rj <= field_lims[3][2])
    return tf_in_bounds
end