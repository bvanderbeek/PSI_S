ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/plotting_library.jl")

###########################
#### PLOT CHAIN SUMMARY ###
###########################

# START INPUT #
chain_directory = "output/cas_ani_ref_NN200_f10" # Where chains are stored
chain_name = "chain" # Chain file name follows format, chain_name.chain_number.jld
num_chains = 28 # Plots metrics for chains 1:num_chains
nsaved = 1e4 # Chain saving interval defined in buildIP.jl  (inconsequential for plots)
istart = 1 # Plot chain metrics starting at iteration istart (zoom in to later iterations by increasing)
tf_display_figures = true # If true, will display figures. Define false for machines with no display.
# END INPUT #

# Load chain metrics. This is an array of ChainMetrics structures. See definition in src/post_processing.jl for field descriptions.
cm = load_chains_metrics(num_chains; chain_directory = chain_directory, chain_name = chain_name, nsaved = nsaved)

# Construct plots. Note!
# To convert residual variance reduction to absolute root-mean_square residual; rms_final = rms_initial*sqrt(1.0 + variance_reduction)
fig_array, fig_names = plot_chain_metrics(cm; istart = istart, tf_display_figures = tf_display_figures)
# Using the returned figure handles, you can customise plot options by adding you own code below (e.g. axes limits etc.)
# Alternatively, access all the chain metrics from the 'cm' vector defined above to make your own plots.

# Print figures (hard-coded figure file names)
for (i, h) in enumerate(fig_array)
    png(h, chain_directory*"/"*fig_names[i]*".png")
end



####################################
#### PLOT POSTERIOR DISTRIBUTION ###
####################################

# START INPUT #
chain_directory = "output/cas_ani_ref_NN200_f10" # Where chains are stored
chain_name = "chain" # Chain file name follows format, chain_name.chain_number.jld
num_chains = 28 # Collect chains 1:num_chains
num_burn = 49 # Remove iterations 1:num_burn from each chain when computing ensemble
num_out  = 101 # Number of models output (1 + total_iterations/saving_interval)
# Note! The num_burn refers to the index of the output model and NOT the absolute iteration.
# It is defined, num_burn = 1 + burn_iteration/saving_interval (+1 because first saved model is starting model).
# For example, if we want to discard the first 5e5 iterations and our saving interval is 1e4, then num_burn = 51.

# Collect posterior samples for these fields
field_names = ("dlnVp", "fp", "psi", "gamma")
# Collect posterior samples at these points
qlat, qlon, qrad = [45.0, 45.0], [-123.0, -123.0], 6371.0 .- [100.0, 200.0]
# END INPUT #


# Build array of sample points at which we want to collect posterior values
xg, yg, zg = geo_to_cartesian(qlat, qlon, qrad)
sample_points = vcat(transpose(xg), transpose(yg), transpose(zg))

# Define posterior for each field. Collects samples in each chain for every field at the requested points.
# Using these samples, you can plot the posterior PDF (i.e. histogram of parameter values). The returned
# posterior array has dimensions (num_fields, num_samples, num_points).
posterior_samples = build_posterior(num_chains, num_burn, num_out, sample_points, field_names;
    chain_directory = chain_directory, chain_name = chain_name)

# Plot posterior PDFs.
ifld, kpt = 1, 2 # Plot i'th field and k'th sample point
histogram(posterior_samples[ifld, :, kpt], label = field_names[ifld])

# Look for correlation between two fields
kpt = 2 # Correlation at k'th sample point
i1, i2 = 1, 2 # Correlation between fields i1 and i2
scatter(posterior_samples[i1, :, kpt], posterior_samples[i2, :, kpt],
    xlabel=field_names[i1], ylabel=field_names[i1], label=nothing)
# Instead of individual points, make 2D histogram
histogram2d(posterior_samples[i1, :, kpt], posterior_samples[i2, :, kpt])



#############################
#### BUILD ENSEMBLE MODEL ###
#############################

# START INPUT #
# Output Information
vtk_name = "EnsembleModel" # Name of vtk file that will be written
chain_directory = "output/cas_ani_ref_NN200_f10" # Where chains are stored
chain_name = "chain" # Chain file name follows format, chain_name.chain_number.jld
num_chains = 5 # Collect chains 1:num_chains
num_burn = 98 # Remove iterations 1:num_burn from each chain when computing ensemble
# Note! The num_burn refers to the index of the output model and NOT the absolute iteration.
# It is defined, num_burn = 1 + burn_iteration/saving_interval (+1 because first saved model is starting model).
# For example, if we want to discard the first 5e5 iterations and our saving interval is 1e4, then num_burn = 51.

# Compute means and standard deviations for these fields
field_names = ("dlnVp", "fp", "psi", "gamma")
tf_squeeze = true # If true, will not interpolate models outside the field domains

# Define output grid
minimum_latitude, maximum_latitude, n_latitude = 36.0, 54.0, 181
minimum_longitude, maximum_longitude, n_longitude = -134.0, -112.0, 221
minimum_radius, maximum_radius, n_radius = 6371.0 - 800.0, 6371.0, 81
# END INPUT #


# Build an array points that defines the ensemble model. Each sampled voronoi model in the posterior
# will be interpolated to these points. This array should have dimensions (3, number_of_points) where
# the cartesian x,y,z coordintes are stored along the columns for each point.
xglobal, _ = build_global_geographic_array(minimum_latitude, maximum_latitude, n_latitude,
    minimum_longitude, maximum_longitude, n_longitude,
    minimum_radius, maximum_radius, n_radius)

# Define ensemble model. Computes mean and standard deviation of the posterior at each point (columns)
# defined in xglobal for each field (rows). The returned arrays have dimensions (num_fields, num_points).
mean_model, sdev_model = build_ensemble_model(num_chains, num_burn, xglobal, field_names;
    chain_directory = chain_directory, chain_name = chain_name, tf_squeeze = tf_squeeze)

# Write ensemble model to VTK
vtk_file = chain_directory*"/"*vtk_name
# Scalar fields. Define list of field names and store each field's value in an array with dimension (num_fields, num_points).
names_scalar_fields = (field_names...,"sdev_".*field_names...) # Note, '...' syntax unpacks tuples so they can be concatonated
scalar_fields_array = vcat(mean_model, sdev_model)
# Save model file
write_ensemble_model_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal;
    scalar_fields = (names_scalar_fields, scalar_fields_array))



#########################################
#### BUILD ANISOTROPIC ENSEMBLE MODEL ###
#########################################

# START INPUT #
# Output Information
vtk_name = "EnsembleModel" # Name of vtk file that will be written
chain_directory = "output/cas_ani_ref_NN200_f10" # Where chains are stored
chain_name = "chain" # Chain file name follows format, chain_name.chain_number.jld
num_chains = 5 # Collect chains 1:num_chains
num_burn = 98 # Remove iterations 1:num_burn for each chain
# Note! The num_burn refers to the index of the output model and NOT the absolute iteration.
# It is defined, num_burn = 1 + burn_iteration/saving_interval (+1 because first saved model is starting model).
# For example, if we want to discard the first 5e5 iterations and our saving interval is 1e4, then num_burn = 51.

# Compute means, standard deviations, and linear correlation coefficient for these fields
field_a, field_b = "dlnVp", "fp"
tf_squeeze = true # If true, will not interpolate models outside the field domains
# Compute tensoral average of the THREE anisotropy orientatin fields (squeezing not yet implemented for this calculation)
directional_fields = ("fp", "psi", "gamma")
tf_local_to_global = true # If true, will rotate symmetry axes from an assumed *local* geographic coordinate system to global cartesian
# Function that describes how to convert the three directional fields into spherical parameters
# Useful when anisotropic inversion uses different parameterisations
f_spherical = (a,b,c) -> (a,b,c)

# Define output grid
minimum_latitude, maximum_latitude, n_latitude = 36.0, 54.0, 181
minimum_longitude, maximum_longitude, n_longitude = -134.0, -112.0, 221
minimum_radius, maximum_radius, n_radius = 6371.0 - 800.0, 6371.0, 81
# END INPUT #

# Build array of interpolation points that defines the ensemble model.
xglobal, _ = build_global_geographic_array(minimum_latitude, maximum_latitude, n_latitude,
    minimum_longitude, maximum_longitude, n_longitude,
    minimum_radius, maximum_radius, n_radius)

# Define ensemble model. Computes means, standard deviations, and linear correlation among two parameter fields
# from the posterior samples at each point (columns) defined in xglobal. The returned mean and standard deviation arrays
# have dimensions (2, num_points) while the correlation array has dimensions (1, num_pts)
mean_model, sdev_model, pccf_model = build_ensemble_correlation(num_chains, num_burn, xglobal, field_a, field_b;
    chain_directory = chain_directory, chain_name = chain_name, tf_squeeze = tf_squeeze)

# Define principal axes model. To average directional parameters (i.e. anisotropy), we construct a directional sampling
# tensor (se Eq. Appendix B of Munzarova et. al, GJI 2018) that described an ellipsoid. The major axis of this ellipsoid
# parallels the mean direction while the minor axes reflect uncertianty in this orientation. The returned mean tensor
# is an array of size (3, 3, num_points) with tensors at each grid point. The principal axes of this tensor (paxes_1,2,3)
# are also returned and scaled by the eigenvalues (paxes_1 = major). These have dimensions (3, num_points). Lastly, the 
# directional bias (see Wikipedia entry for fractional anisotropy) is returned as a (1, num_pts) array. A value of 0 implies no
# preferred orientation while a value of 1 implies a single unique orientation.
mean_tensor, paxes_1, paxes_2, paxes_3, directional_bias = build_ensemble_tensor(num_chains, num_burn, xglobal, directional_fields;
    spherical_conversion = f_spherical, tf_local_to_global = tf_local_to_global,
    chain_directory = chain_directory, chain_name = chain_name)

# Write ensemble model to VTK
vtk_file = chain_directory*"/"*vtk_name
# Scalar fields. Define list of field names and store each field's value in an array with dimension (num_fields, num_points).
names_scalar_fields = (field_a, field_b, "sdev_"*field_a, "sdev_"*field_b, "pccf_"*field_a*"_"*field_b, "directional_bias")
scalar_fields_array = vcat(mean_model, sdev_model, pccf_model, directional_bias)
# Vector fields. Define list of field names and store each vector field with dimensions (3, num_points)
# in a vector (or tuple) of length 3.
names_vector_fields = ("paxes_1", "paxes_2", "paxes_3")
vector_fields_tuple = (paxes_1, paxes_2, paxes_3)
# Tensor fields. Define list of tensor field names; tenosr fields stored in array of size (3,3,num_points)
names_tensor_fields = ["ani_tensor"]
tensor_fields_array = [mean_tensor] # Must be able to index such that tensor_fields_array[1] = Array{Float64}(3,3,num_points)
# Save model file
write_ensemble_model_vtk(vtk_file, (n_latitude, n_longitude, n_radius), xglobal;
    scalar_fields = (names_scalar_fields, scalar_fields_array),
    vector_fields = (names_vector_fields, vector_fields_tuple),
    tensor_fields = (names_tensor_fields, tensor_fields_array))


# TESTING DIRECTIONAL FIELDS #

# Reference local vector (east, north, radial)
u0 = [0.05; 0.0; 0.0]
# Reference local tensor (east, north, radial)
T0 = [0.05 0.0 0.0; 0.0 0.03 0.0; 0.0 0.0 0.01]
# Add some rotation
Rm = rotation_matrix(deg2rad(30.0), 2)
u = Rm*u0
T = Rm*T0*transpose(Rm)
# Define homogeneous directional fields
mean_vector = zeros(3, n_latitude*n_longitude*n_radius);
mean_tensor = zeros(3, 3, n_latitude*n_longitude*n_radius);
[mean_vector[:,i] .= u for i in axes(xglobal, 2)];
[mean_tensor[:,:,i] .= T for i in  axes(xglobal, 2)];
# Convert to global cartesian coordinates
local_to_global_vectors!(mean_vector, xglobal)
local_to_global_tensors!(mean_tensor, xglobal)
# Write file
write_ensemble_model_vtk(chain_directory*"/TEST", (n_latitude, n_longitude, n_radius), xglobal;
    scalar_fields = (["mag"], norm(u0)*ones(1, n_latitude*n_longitude*n_radius)),
    vector_fields = (["vector"], [mean_vector]),
    tensor_fields = (["tensor"], [mean_tensor]));



########################
# CHECK DERIVED PRIORS #
########################
N = 1_000_000
fmax = 0.1 # Needed for reference but does not effect directional bias calculation

# Directional bias for spherical parameterisation approaches ~0.408...
# Note! Directional bias is independent of mean magnitude (which is quite small for uniform prior)
f, azm, elv = fmax*rand(N), 2.0*pi*(2.0*rand(N) .- 1.0), 0.5*pi*rand(N);
db_0, T_0 = prior_directional_tensor(f, azm, elv; spherical_conversion = (a,b,c) -> (a,b,c))
# To account for the magnitude, we could normalise the directional bias by (major_axis)/mean(magnitude)?
λ = eigvals(T_0)
db_0 *= maximum(abs.(λ))/mean(f)

# Directional bias using projection of vertical component approaches 0
prj = 2.0*rand(N) .- 1;
db_1, T_1 = prior_directional_tensor(f, azm, prj; spherical_conversion = (a,b,c) -> (a, b, asin(c)))

# Directional bias using vectoral parameterisation approaches 0
s1, s2, s3 = fmax*(2.0*rand(N) .- 1), fmax*(2.0*rand(N) .- 1), fmax*(2.0*rand(N) .- 1)
db_2, T_2 = prior_directional_tensor(s1, s2, s3; spherical_conversion = (a,b,c) -> (sqrt(a^2 + b^2 + c^2), atan(b, a), atan(c, sqrt(a^2 + b^2))))