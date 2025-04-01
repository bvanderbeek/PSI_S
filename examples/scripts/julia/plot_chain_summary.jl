###########################################
### Plot Chain Summary from Source Code ###
###########################################
# Julia script to plot the variance reduction, model dimensions, acceptance ratio,
# likelihood function, and data noise parameter per iteration for each Markov chain
# via the shell call,
# >> julia plot_chain_summary.jl parameter_file istart
# Where parameter_file is a TOML-file specifying the PSI_S inversion parameters and istart
# is the iteration at which to begin the plot (istart = 0 corresponds to the starting model)
#
# The resulting figures are saved to the ouput location defined in parameter_file.
#
# Making these plots is not computationally heavy and can be done locally/interactively.


# Set Environment Variables #
ENV["GKSwstype"] = "nul" # Disable display for remote plotting
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S" # Path to code directory
include(ENV["PSI_S"]*"/src/plotting_library.jl") # Load plotting functions

# Parse Command-Line Arguments #
parameter_file, istart = ARGS[1], parse(Int64, string(ARGS[2]))

# Plot Chain Metrics #
plot_chain_metrics(parameter_file; istart = istart, tf_display_figures = false, tf_print = true)
