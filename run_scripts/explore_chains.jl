## SETUP ##
### Using PSI_S Package ###
# using Plots
# using JLD
# using PSI_S
### Using Local PSI_S Package ###
# ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
# using Pkg
# Pkg.activate(ENV["PSI_S"])
# using Plots
# using JLD
# using PSI_S
### Using PSI_S Source Code ###
ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
include(ENV["PSI_S"]*"/src/dependencies.jl")
using Plots

## START INPUT ##
niterations = 1e5
nsaved = 1e4
nchains = 4
chain_name = "test"
out_dir = "/Users/bvanderbeek/research/software/GitRepos/PSI_S/example/MELT/output/test"
## END INPUT ##

function extract_chain_metrics(chain_name, n_chains, n_iterations, n_saved; chain_directory = pwd())
    i_chain = Vector{Vector{Int}}(undef, n_chains)
    n_cells = Vector{Vector{Int}}(undef, n_chains)
    n_accepted = Vector{Vector{Float64}}(undef, n_chains)
    n_rms = Vector{Vector{Float64}}(undef, n_chains)
    n_noise = Vector{Vector{Float64}}(undef, n_chains)
    for n in 1:n_chains
        chain_file = chain_directory*"/"*chain_name*"."*string(n)*".jld"
        Chain_n = load(chain_file, "MarkovChain")

        i_chain[n], n_cells[n], n_accepted[n], n_rms[n], n_noise[n] = extract_chain_metrics(Chain_n, n_iterations, n_saved)
    end
    return i_chain, n_cells, n_accepted, n_rms, n_noise
end

function extract_chain_metrics(Chain::Vector{ModelConst}, niterations, nsaved)
    N = length(Chain)
    i_chain = Vector{Int}(undef, N)
    n_cells = Vector{Int}(undef, N)
    n_accepted = Vector{Float64}(undef, N)
    n_rms = Vector{Float64}(undef, N)
    # Compute norm of each Voronoi diagram?
    # n_es_rms = Vector{Float64}(undef, N)
    # n_ss_rms = Vector{Float64}(undef, N)
    n_noise = Vector{Float64}(undef, N)
    for i in eachindex(Chain)
        i_chain[i] = (i - 1)*nsaved
        n_cells[i] = sum(f -> f.n[1], Chain[i].fields) # All fields; field 'n' is always a single-element vector
        n_accepted[i] = Chain[i].accepted[1]
        n_rms[i] = Chain[i].rms[1]
        # The statics and noise parameter are stored in dataspace
        # The 'Obs' field is always a single-element vector
        # n_es_rms[i] = rms(Chain_n[i].dataspace.Obs[1].estatics) # All event statics
        # n_ss_rms[i] = rms(Chain_n[i].dataspace.Obs[1].sstatics) # All station statics
        n_noise[i] = Chain[i].dataspace.Obs[1].noise[1] # Noise parameter; always a single-element vector
    end
    
    return i_chain, n_cells, n_accepted, n_rms, n_noise
end

i_chain, n_cells, n_accepted, n_rms, n_noise = extract_chain_metrics(chain_name, nchains, niterations, nsaved; chain_directory = out_dir);


h = plot();
for n in eachindex(i_chain)
    plot!(h, i_chain[n]./nsaved, n_cells[n], legend = false)
end
xlabel!(h, "n'th model out");
ylabel!(h, "number of cells");
display(h)
png(out_dir*"/num_cells.png")

h = plot();
for n in eachindex(i_chain)
    #plot!(h, i_chain[n][1:end]./nsaved, n_rms[n][1:end], legend = false)
    plot!(h, i_chain[n][2:end]./nsaved, (n_rms[n][2:end].^2 .- n_rms[n][1]^2)./(n_rms[n][1]^2), legend = false)
end
xlabel!(h, "n'th model out");
ylabel!(h, "Var. Red. (%)");
display(h)
png(out_dir*"/var.png")

h = plot();
for n in eachindex(i_chain)
    plot!(h, i_chain[n][2:end]./nsaved, n_noise[n][2:end], legend = false)
end
xlabel!(h, "n'th model out");
ylabel!(h, "Noise (s)");
display(h)
png(out_dir*"/noise.png")

h = plot();
for n in eachindex(i_chain)
    plot!(h, i_chain[n][2:end]./nsaved, 100*n_accepted[n][2:end]./i_chain[n][2:end], legend = false)
end
xlabel!(h, "n'th model out");
ylabel!(h, "percent accepted");
display(h)
png(out_dir*"/acceptance.png")