############################################
### Make Ensemble Model from Source Code ###
############################################
# Julia script to compute and save the ensemble model as a vts-file via the shell call,
# >> julia make_ensemble_vts.jl
# Several input parameters need to be specified below. A vts-file that can be visualised
# in Paraview will be saved to the user-defined location specified below.

# Set-Up Environment Variables #
# ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/dependencies.jl")

# START INPUT

# Inversion Parameters
# Burn-in number. Remove iterations 1:num_burn from each Markov chain.
# Note! The num_burn refers to the index of the output model and NOT the absolute iteration.
# It is defined, num_burn = 1 + burn_iteration/saving_interval (+1 because first saved model is starting model).
# For example, if we want to discard the first 5e5 iterations and our saving interval is 1e4, then num_burn = 51.
num_burn = 50 # -- INT
num_chains = 28 # Collect Markov chains 1:num_chains -- INT
chain_directory = "output/test_1" # Where chains are stored -- STRING
chain_name = "chain" # Chain names (usually "chain"). Chain file names follow the format, chain_name.chain_number.jld -- STRING
vtk_file = chain_directory * "/EnsembleModel_Burn" * string(num_burn) # Name of vtk file that will be written -- STRING

# Ensemble Grid Parameters: Define the redular grid parameters for interpolating the model.
lat_limdim = (-9.0, 9.0, 181) # Latitude (DEGREES) sampling Tuple, (minimum latitude, maximum latitude, number of samples) -- (Float, Float, Int)
lon_limdim = (-9.0, 9.0, 181) # Longitude (DEGREES) sampling Tuple, (minimum longitude, maximum longitude, number of samples) -- (Float, Float, Int)
rad_limdim = (6371.0 - 800.0, 6371.0, 81) # Radial (KILOMETERS) sampling Tuple, (minimum radius, maximum radius, number of samples) -- (Float, Float, Int)
tf_squeeze = false # If true, will use extrapolation value defined in parameter file when interpolating Voronoi models outside field domain -- BOOL

# Ensemble Fields to Compute
# Compute means and standard deviations for these fields (e.g. ["dlnVp", "dlnVs"] or [] is none)
ensemble_fields = ["fp", "psi", "gamma"] # -- Vector{String}(number of fields)
# Compute means, standard deviations, and linear correlation coefficient for these TWO fields (e.g. ["dlnVp", "dlnVs"])
correlation_fields = ["dlnVp", "fp"] # -- Vector{String}(2)
# Compute tensoral average of these THREE directional fields and the directional bias; Vector{String}(3) (e.g. ["fp", "psi", "gamma"])
directional_fields = ["fp", "psi", "gamma"] # -- Vector{String}(3)
# Anonymous function that describes how to convert the three directional fields into spherical parameters (magnitude, azimuth, and elevation)
# Useful when anisotropic inversion uses different parameterisations. For example, if inverting directly for the components of the anisotropic
# vector, then define f_spherical = (a,b,c) -> (sqrt(a^2 + b^2 + c^2), atan(b, a), atan(c, sqrt(a^2 + b^2)))
f_spherical = (a, b, c) -> (a, b, c)
# If true, will only write the major axis of the directional tensor (i.e. the mean anisotropy orientation).
# If false, all three eigenvector are written.
tf_write_major_only = true # -- BOOL
# If true, write the 3x3 mean directional tensor to the vtk file
tf_write_tensors = true # -- BOOL
# If true, convert the directional tensor (and eigen vectors) from local anisotropic coordinate system
# (east, north, radial) to the global earth-centered earth-fixed coordinates (xg, yg, zg). This should
# generally be true
tf_local_to_global = true # -- BOOL

# END INPUT



# Build Ensemble Model
write_ensemble_model_vtk(
    lat_limdim, lon_limdim, rad_limdim, num_chains, num_burn, tf_squeeze;
    ensemble_fields = ensemble_fields, correlation_fields = correlation_fields,
    directional_fields = directional_fields, f_spherical = f_spherical,
    chain_directory = chain_directory, chain_name = chain_name, vtk_file = vtk_file,
    tf_write_major_only = tf_write_major_only, tf_write_tensors = tf_write_tensors,
    tf_local_to_global = tf_local_to_global
)
