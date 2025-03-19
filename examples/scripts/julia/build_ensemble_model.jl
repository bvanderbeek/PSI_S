############################################
### Make Ensemble Model from Source Code ###
############################################
# Julia script to compute and save the ensemble model as a vts-file via the shell call,
# >> julia build_ensemble_model.jl
# Several input parameters need to be specified below. A vts-file that can be visualised
# in Paraview will be saved to the user-defined location specified below.

# Set-Up Environment Variables #
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/plotting_library.jl")
# include(ENV["PSI_S"]*"/src/dependencies.jl")

# START INPUT

# Inversion Parameters
# Burn-in number. Remove iterations 1:num_burn from each Markov chain.
# Note! The num_burn refers to the index of the output model and NOT the absolute iteration.
# It is defined, num_burn = 1 + burn_iteration/saving_interval (+1 because first saved model is starting model).
# For example, if we want to discard the first 5e5 iterations and our saving interval is 1e4, then num_burn = 51.
num_burn = 6 # -- INT
num_chains = 4 # Collect Markov chains 1:num_chains -- INT
chain_directory = "output/iso_ttp" # Where chains are stored -- STRING
chain_name = "chain" # Chain names (usually "chain"). Chain file names follow the format, chain_name.chain_number.jld -- STRING
vtk_file = chain_directory * "/PosteriorMoments_Chains" * string(num_chains) * "_Burn" * string(num_burn) # Name of vtk file that will be written -- STRING

# Ensemble Grid Parameters: Define the redular grid parameters for interpolating the model.
lat_limdim = (37.2, 37.8, 61) # Latitude (DEGREES) sampling Tuple, (minimum latitude, maximum latitude, number of samples) -- (Float, Float, Int)
lon_limdim = (14.7, 15.3, 61) # Longitude (DEGREES) sampling Tuple, (minimum longitude, maximum longitude, number of samples) -- (Float, Float, Int)
rad_limdim = (6371.0 - 12.0, 6371.0 + 1.0, 14) # Radial (KILOMETERS) sampling Tuple, (minimum radius, maximum radius, number of samples) -- (Float, Float, Int)
tf_squeeze = false # If true, will use extrapolation value defined in parameter file when interpolating Voronoi models outside field domain -- BOOL
tf_cart_to_geo = true # If true, domain limits are defined in geographic coordinates but interpolation grid is global geographic

# Ensemble Model Options (depends on use-case)
ensemble_option = 1

# (1) Posterior Moments
if ensemble_option == 1
    num_moments = 3 # Compute moments 1:num_moments; Int < 4
    field_names = ["dlnVp"] # Compute moments for these fields; Vector{String}
    correlation_pairs = [] # Compute correlation coefficients for these field pairs; Empty or Vector{NTuple{2,String}}
end

# (2) Posterior Directional Moments
if ensemble_option == 2
    num_moments = 3 # Compute moments 1:num_moments; Int < 4
    directional_fields = ["fp", "psi", "gamma"] # Compute moments for these directional fields; Vector{String}
    spherical_conversion = (a, b, c) -> (a, b, asin(c)) # Generic function
    tf_global = true # Connvert vectors to global ECEF coordinates
    # Compute moments and correlation coefficient between the following fields and the directional strength
    correlation_fields = ["dlnVp"] # Empty or Vector{String}
end

# (3) Posterior Tendencies
if ensemble_option == 3
    field_names = ["dlnVp"] # Compute moments for these fields; Vector{String}
    credible_level = 0.68
    q_quantiles = []
    mode_intervals = zeros(length(field_names))
end

# (4) Posterior Directional Tendencies
if ensemble_option == 4
    directional_fields = ["fp", "psi", "gamma"] # Compute moments for these directional fields; Vector{String}
    spherical_conversion = (a, b, c) -> (a, b, asin(c)) # Generic function
    tf_global = true # Connvert vectors to global ECEF coordinates
    credible_level = 0.68
    q_quantiles = []
    mode_intervals = zeros(length(directional_fields))
end

# END INPUT


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
    minimum_radius, maximum_radius, n_radius
)

# --- Posterior Moments --- #

# Compute posterior moments
if ensemble_option == 1
    post_moments, post_correlation = posterior_moments(
        num_chains, num_burn, xglobal, field_names;
        correlation_pairs = correlation_pairs, num_moments = num_moments, chain_directory = chain_directory,
        chain_name = chain_name, tf_squeeze = tf_squeeze, tf_cart_to_geo = tf_cart_to_geo
    )

    write_ensemble_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal, post_moments, post_correlation, field_names, correlation_pairs)
end

# --- Posterior Directional Moments --- #

# Compute directional posterior moments
if ensemble_option == 2
    post_moments, post_correlation, post_vectors, directional_bias = posterior_directional_moments(
        num_chains, num_burn, xglobal, directional_fields; correlation_fields = correlation_fields, num_moments = num_moments,
        chain_directory = chain_directory, chain_name = chain_name, tf_squeeze = tf_squeeze, tf_cart_to_geo = tf_cart_to_geo,
        tf_global = tf_global, spherical_conversion = spherical_conversion
    )

    write_ensemble_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal, post_moments, post_correlation,
        directional_bias, post_vectors, directional_fields, correlation_fields
    )
end

# --- Posterior Tendencies --- #

if ensemble_option == 3
    # Load chains and nearest neighbor trees
    markov_chains = load_markov_chains(num_chains, num_burn; chain_directory = chain_directory, chain_name = chain_name)
    nn_trees = build_nn_trees(markov_chains, field_names; leafsize = 10)
    # Compute posterior tendencies
    PostTend = posterior_tendencies(xglobal, markov_chains, nn_trees, field_names;
        credible_level = credible_level, q_quantiles = q_quantiles, mode_intervals = mode_intervals
    )

    # Write vtk file
    write_ensemble_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal, PostTend, field_names)
end

# --- Posterior Directional Tendencies --- #

if ensemble_option == 4
    # Load chains and nearest neighbor trees
    markov_chains = load_markov_chains(num_chains, num_burn; chain_directory = chain_directory, chain_name = chain_name)
    nn_trees = build_nn_trees(markov_chains, directional_fields; leafsize = 10)
    # Compute posterior tendencies
    PostTend = posterior_directional_tendencies(xglobal, markov_chains, nn_trees, directional_fields;
        credible_level = credible_level, q_quantiles = q_quantiles, mode_intervals = mode_intervals,
        spherical_conversion = spherical_conversion, tf_global = tf_global
    )

    # Write vtk file
    write_ensemble_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal, PostTend, directional_fields)
end