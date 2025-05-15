###################################
### Call PSI_S from Source Code ###
###################################
# Julia script to run rj-mcmc inversion via the shell call,
# >> julia call_psi_s.jl parameter_file num_chains_per_chunk num_chunk
# Can also call with --check-bounds=no flag for potenital performance gain
# To use multi-threading when 3D ray-tracing,
# >> julia --threads 28 call_psi_s.jl parameter_file num_chains_per_chunk num_chunk
# Where parameter_file is a TOML-file specifying the PSI_S inversion parameters, num_chains_per_chunk is the
# number of Markov chains to run in parallel (must match available resources) and num_chunk specifies which
# batch of chains we are running. For example, if we want to run a total of 56 chains on two independent
# compute nodes (i.e. different jobs) with 28 cores each, we would need to run 2 chunks of 28 Markov chains each.
# Note! JavaCall, which is required for using TauP, breaks mulit-threading in Julia. This appears to be the result
# of defining the environment variable JULIA_COPY_STACKS = 1 (https://github.com/JuliaLang/julia/issues/44589).
# For this reason, TauP is no longer loaded via dependencies.jl. If needed for teleseismic inversions, load the
# package as part of this script. TauP does not use any parallel functionality.
#
# Results are saved to the output location defined in parameter_file. Following files are written:
# 1) parameters.toml -- copy of parameter used to run inversion
# 2) IPConst.jld
# 3) chain.x.jld (one for each chain)
# 4) psi_s_utilities.jld

# Set-Up Environment Variables #
# Alternatively, define in HPC job submission script. 
# The "TAUP_JAR" variable is optional (default is that which ships with TauP.jl)
# The following three lines are only needed if using TauP to trace rays (e.g. teleseismic imgaing)
# If not using TauP, these lines should not be active as they will break parallel 3D ray tracing capabilities
# ENV["TAUP_JAR"] = "TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
# ENV["JULIA_COPY_STACKS"] = "yes" # For non-windows users, this environment variable may be required for JavaCall to function properly; may break multi-threading
# using TauP
# PSI_S source code directory
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"

# Load Required Packages #
# Only need to load Distributed.jl to initialise parallel processes
# All other required packages files are specified in the PSI_S dependencies.jl file loaded below
using Distributed

# Parse Command-Line Arguments #
parameter_file, num_chains, ichunk = ARGS[1], parse(Int64, string(ARGS[2])), parse(Int64, string(ARGS[3]))

# Initialise Parallel Workers #
# Starts workers and copies source code to process
addprocs(num_chains)
@everywhere include(ENV["PSI_S"]*"/src/dependencies.jl")

# Call PSI_S #
run_psi_s(parameter_file, ichunk)
