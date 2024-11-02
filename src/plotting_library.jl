# Plotting Library
# Collection of useful plotting functions
include(@__DIR__() * "/dependencies.jl") # I live in the src directory
using Plots

# WILL NEED TO BE UPDATED FOR MULTI-OBSERVABLE INVERSIONS
# + Assumes single noise parameter
function plot_chain_metrics(parameter_file; istart = 1, tf_display_figures = false, tf_print = false)
    # Load parameter file
    P = TOML.parsefile(parameter_file)
    # Load chains
    chain_directory = P["out_directory"] * "/" * P["run_id"]
    cm = load_chains_metrics(P["MonteCarloSolver"]["num_chains"];
        chain_directory = chain_directory, nsaved = P["MonteCarloSolver"]["save_interval"])
    # Plot chains
    fig_array, fig_names = plot_chain_metrics(cm; istart = istart, tf_display_figures = tf_display_figures)
    # Print figures (hard-coded figure file names)
    tf_print && [png(h, chain_directory * "/" * fig_names[i] * ".png") for (i, h) in enumerate(fig_array)]

    return fig_array, fig_names
end
function plot_chain_metrics(CM::ChainMetrics;
    istart = 1, fig_array = [plot(),plot(),plot(),plot(),plot()], tf_display_figures = true,
    field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink))
    # Number of iterations
    istart += 1 # Offset iteration index because index 1 = iteration 0 (starting model)
    num_its = CM.iteration[end] - 1
    b = floor(Int, log10(num_its))
    a = round(num_its/(10^b), digits = 1)
    xlabel_its = "percent iteration (of "*string(a)*"e"*string(b)*")"
    # Hard-coded figure names
    fig_names = ("AcceptanceRatio", "Objective", "VarianceReduction", "Noise", "Dimensions")
    # Define percent iteration
    rit = 100.0*(CM.iteration[istart:end] .- 1)./num_its
    # Acceptance Ratio
    plot!(fig_array[1], rit, 100.0*CM.num_accepted[istart:end]./CM.iteration[istart:end],
    xlabel = xlabel_its, ylabel = "acceptance ratio (%)", legend = false)
    # Objective function
    f0 = string(round(log10(CM.fobj[1]), digits = 4))
    plot!(fig_array[2], rit, (CM.fobj[istart:end] .- CM.fobj[1])./CM.fobj[1],
    xlabel = xlabel_its, ylabel = "likelihood reduction (%)", title = "Initial Log Likelihood = "*f0, legend = false)
    # Residual variance reduction
    rms0 = string(round(CM.rrms[1], digits = 4))
    plot!(fig_array[3], rit, (CM.rrms[istart:end].^2 .- CM.rrms[1]^2)./(CM.rrms[1]^2),
    xlabel = xlabel_its, ylabel = "variance reduction (%)", title = "Initial RMS = "*rms0, legend = false)
    # Noise parameter
    plot!(fig_array[4], rit, CM.pnoise[istart:end],
    xlabel = xlabel_its, ylabel = "noise parameter", legend = false)
    # Number of cells per field
    tf_add_labels = fig_array[5].n > 0 #  Do not add more line labels when figure already populated
    for (field_name, j) in CM.fields
        f_label = tf_add_labels ? nothing : field_name
        plot!(fig_array[5], rit, CM.num_cells[istart:end,j], linecolor = field_colors[j], labels = f_label)
    end
    xlabel!(fig_array[5],xlabel_its)
    ylabel!(fig_array[5], "number of cells")

    tf_display_figures ? [display(h) for h in fig_array] : nothing
    return fig_array, fig_names
end
function plot_chain_metrics(CM::Array{<:ChainMetrics, N};
    istart = 1, tf_display_figures = true,
    field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink)) where {N}

    fig_array, fig_names = plot_chain_metrics(CM[1]; istart = istart, tf_display_figures = false, field_colors = field_colors)
    for i in 2:length(CM)
        plot_chain_metrics(CM[i]; istart = istart, fig_array = fig_array, tf_display_figures = false, field_colors = field_colors)
    end

    tf_display_figures ? [display(h) for h in fig_array] : nothing
    return fig_array, fig_names
end