# Plotting Library
# Collection of useful plotting functions
include(@__DIR__() * "/dependencies.jl") # I live in the src directory
using Plots

# Plot evolution of all chains in a PSI_S inversion
function plot_chain_metrics(parameter_file; istart = 1, tf_display_figures = false, tf_print = false, num_burn = 1)
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

    # Additional plots
    print_path = chain_directory * "/hist_dimensions" * "_burn" * string(num_burn) * ".png"
    plot_dimension_histograms(cm, num_burn; bi = 25, tf_display_figures = tf_display_figures, print_path = print_path)
    print_path = chain_directory * "/hist_rms" * "_burn" * string(num_burn) * ".png"
    plot_residual_histograms(cm, num_burn; tf_display_figures = tf_display_figures, tf_ms = true, print_path = print_path)

    return cm, fig_array, fig_names
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
    istart = 1, ref_rms = CM.rrms[1,:], fig_array = [plot(),plot(),plot(),plot(),plot(),plot()],
    tf_display_figures = true, field_colors = (:black, :blue, :green, :red, :purple, :cyan, :pink))

    # Hard-coded figure names
    fig_names = ("AcceptanceRatio", "Objective", "VarianceReduction", "Noise", "Dimensions", "EventStatics")

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

    # 6. RMS Event Statics
    # Construct labels for each observable
    if fig_array[6].n > 0
        fill!(chain_labels, nothing)
    else
        [chain_labels[j] = obs_name for (j, obs_name) in CM.obs]
    end
    fig_array[6] = plot_data_array!(fig_array[6], rit, CM.rms_event_statics[istart:end,:]; chain_labels = chain_labels, field_colors = field_colors)
    xlabel!(fig_array[6], xlabel_its)
    ylabel!(fig_array[6], "rms event static")

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

function plot_dimension_histograms(CM::Array{<:ChainMetrics, N}, num_burn; bi = 25, tf_display_figures = true, print_path = nothing) where {N}
    num_chains, num_out, num_flds = length(CM), length(CM[1].iteration), length(CM[1].fields)
    num_samples = num_out - num_burn

    # Collect cell numbers for each field
    v = zeros(Int, num_chains*num_samples, num_flds)
    for (i_chain, chain_i) in enumerate(CM)
        ia = 1 + (i_chain - 1)*num_samples
        ib = ia + num_samples - 1
        for (_, j) in chain_i.fields
            v[ia:ib,j] .= chain_i.num_cells[(num_burn+1):num_out,j]
        end
    end

    pan = num_flds == 1 ? (1,1) : (ceil(Int, num_flds/2), 2)
    h = plot(; layout = pan)
    nmin, nmax = extrema(v)
    nmin = bi*floor(Int, nmin/bi)
    nmax = bi*ceil(Int, nmax/bi)
    for (field_name, j) in CM[1].fields
        histogram!(h,
            v[:,j]; label = nothing, title = field_name, bins = nmin:bi:nmax,
            normalize = :probability, seriesalpha = 0.5, subplot = j
        )
    end
    tf_display_figures && display(h)
    !isnothing(print_path) && png(h, print_path)
    return h
end

function plot_residual_histograms(CM::Array{<:ChainMetrics, N}, num_burn; tf_display_figures = true, tf_ms = true, print_path = nothing) where {N}
    num_chains, num_out, num_obs = length(CM), length(CM[1].iteration), size(CM[1].rrms,2)
    num_samples = num_out - num_burn

    v = zeros(num_chains*num_samples, num_obs)
    for (i_chain, chain_i) in enumerate(CM)
        ia = 1 + (i_chain - 1)*num_samples
        ib = ia + num_samples - 1
        for (j, _) in chain_i.obs
            v[ia:ib,j] .= chain_i.rrms[(num_burn+1):num_out,j]
        end
    end
    tf_ms && (v .*= 1000.0)

    pan = num_obs == 1 ? (1,1) : (ceil(Int, num_obs/2), 2)
    h = plot(; layout = pan)
    for (j, field_name) in CM[1].obs
        histogram!(h,
            v[:,j]; label = nothing, title = field_name,
            normalize = :probability, seriesalpha = 0.5, subplot = j
        )
    end
    tf_display_figures && display(h)
    !isnothing(print_path) && png(h, print_path)

    return h
end