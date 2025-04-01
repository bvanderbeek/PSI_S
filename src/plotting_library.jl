# Plotting Library
# Collection of useful plotting functions
include(@__DIR__() * "/dependencies.jl") # I live in the src directory
using Plots

# Plot evolution of all chains in a PSI_S inversion
function plot_chain_metrics(parameter_file; istart = 1, tf_display_figures = false, tf_print = false)
    # Load parameter file
    P = TOML.parsefile(parameter_file)
    # Load chains
    chain_directory = P["out_directory"] * "/" * P["run_id"]
    cm = load_chains_metrics(P["MonteCarloSolver"]["num_chains"] * P["MonteCarloSolver"]["num_chunks"];
        chain_directory = chain_directory, nsaved = P["MonteCarloSolver"]["save_interval"])
    # Plot chains
    fig_array, fig_names = plot_chain_metrics(cm; istart = istart, tf_display_figures = tf_display_figures)
    # Print figures (hard-coded figure file names)
    tf_print && [png(h, chain_directory * "/" * fig_names[i] * "_i" * string(istart) * ".png") for (i, h) in enumerate(fig_array)]

    return fig_array, fig_names
end
# Plot evolution of a collection of Markov chains
function plot_chain_metrics(CM::Array{<:ChainMetrics, N};
    istart = 1, tf_display_figures = true, field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink)) where {N}

    # Average initial residual rms
    ref_rms = zeros(length(CM[1].obs))
    [ref_rms .+= cm_i.rrms[1,:] for cm_i in CM]
    ref_rms ./= length(CM)

    # Plot all chains
    fig_array, fig_names = plot_chain_metrics(CM[1]; istart = istart, ref_rms = ref_rms, tf_display_figures = false, field_colors = field_colors)
    for i in 2:length(CM)
        plot_chain_metrics(CM[i]; istart = istart, ref_rms = ref_rms, fig_array = fig_array, tf_display_figures = false, field_colors = field_colors)
    end

    tf_display_figures ? [display(h) for h in fig_array] : nothing
    return fig_array, fig_names
end
# Plot evolution of a single Markov chain
function plot_chain_metrics(CM::ChainMetrics;
    istart = 1, ref_rms = CM.rrms[1,:], fig_array = [plot(),plot(),plot(),plot(),plot()],
    tf_display_figures = true, field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink))

    # Hard-coded figure names
    fig_names = ("AcceptanceRatio", "Objective", "VarianceReduction", "Noise", "Dimensions")

    # Initialise container for chain labels
    N, F = length(CM.obs), length(CM.fields)
    chain_labels = Vector{Union{String, Nothing}}(undef, max(N,F))

    # Number of (fractional) iterations. More convenient plot axis.
    istart += 1 # Offset iteration index because index 1 = iteration 0 (starting model)
    num_its = CM.iteration[end] # No offset here because CM.iteration[1] = 0
    rit = 100.0*CM.iteration[istart:end]./num_its
    # Construct nice label for x-axis
    b = floor(Int, log10(num_its))
    a = round(num_its/(10^b), digits = 1)
    xlabel_its = "percent iteration (of "*string(a)*"e"*string(b)*")"

    # 1. Acceptance Ratio
    plot!(fig_array[1], rit, 100.0*CM.num_accepted[istart:end]./CM.iteration[istart:end],
    xlabel = xlabel_its, ylabel = "acceptance ratio (%)", legend = false)

    # 2. Objective Function
    plot!(fig_array[2], rit, log10.(CM.fobj[istart:end]), xlabel = xlabel_its, ylabel = "log-likelihood", legend = false)

    # 3. Observable Variance Reduction
    # Compute variance reduction for each observable
    var_red = CM.rrms[istart:end,:].^2
    [var_red[:,j] .= (var_red[:,j] .- rms_j^2)./(rms_j^2) for (j, rms_j) in enumerate(ref_rms)]
    # Construct labels for each observable
    if fig_array[3].n > 0
        fill!(chain_labels, nothing)
    else
        [chain_labels[j] = obs_name for (j, obs_name) in CM.obs]
        [chain_labels[j] *= " ("*string(round(rms_j, digits = 4))*")" for (j, rms_j) in enumerate(ref_rms)]
    end
    fig_array[3] = plot_data_array!(fig_array[3], rit, 100.0*var_red; chain_labels = chain_labels, field_colors = field_colors)
    xlabel!(fig_array[3], xlabel_its)
    ylabel!(fig_array[3], "variance reduction (%)")

    # 4. Noise Parameter
    # Construct labels for each observable
    if fig_array[4].n > 0
        fill!(chain_labels, nothing)
    else
        [chain_labels[j] = obs_name for (j, obs_name) in CM.obs]
        [chain_labels[j] *= " ("*string(round(rms_j, digits = 4))*")" for (j, rms_j) in enumerate(ref_rms)]
    end
    fig_array[4] = plot_data_array!(fig_array[4], rit, CM.pnoise[istart:end,:]; chain_labels = chain_labels, field_colors = field_colors)
    xlabel!(fig_array[4], xlabel_its)
    ylabel!(fig_array[4], "noise parameter")

    # 5. Field Dimensions
    # Construct labels for each field
    if fig_array[5].n > 0
        fill!(chain_labels, nothing)
    else
        [chain_labels[j] = field_name for (field_name, j) in CM.fields]
    end
    fig_array[5] = plot_data_array!(fig_array[5], rit, CM.num_cells[istart:end,:]; chain_labels = chain_labels, field_colors = field_colors)
    xlabel!(fig_array[5], xlabel_its)
    ylabel!(fig_array[5], "number of cells")

    tf_display_figures ? [display(h) for h in fig_array] : nothing
    return fig_array, fig_names
end
# Convenience function to iteratively plot data stored in columns of an array
function plot_data_array!(Fig, iterations, chain_metric;
    chain_labels = Vector{Nothing}(nothing, size(chain_metric, 2)),
    field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink))
    for j in axes(chain_metric, 2)
        plot!(Fig, iterations, chain_metric[:,j], linecolor = field_colors[j], labels = chain_labels[j])
    end
    return Fig
end