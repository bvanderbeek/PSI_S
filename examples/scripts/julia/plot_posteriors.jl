# Set-Up Environment Variables #
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/plotting_library.jl")

num_burn = 75 # -- INT
num_chains = 28 # Collect Markov chains 1:num_chains -- INT
chain_directory = "/Users/bvanderbeek/research/CASCADIA/Joint/PSI_S/output/S/TEST_ANI_PARAM/REF_N200" #TEST_ISO_PARAMS/REF" #TEST_ANI_PARAM/REF_N200" # Where chains are stored -- STRING
chain_name = "chain" # Chain names (usually "chain"). Chain file names follow the format, chain_name.chain_number.jld -- STRING
nsaved = 1e4

cm = load_chains_metrics(num_chains; chain_directory = chain_directory, chain_name = chain_name, nsaved = nsaved)
rrms = Vector{Float64}()
num_cells = Vector{Float64}()
[append!(rrms, cm_i.rrms[num_burn:101]) for cm_i in cm]
[append!(num_cells, sum(cm_i.num_cells[num_burn:101,:]; dims = 2)) for cm_i in cm]

nbins = 18
h_rms = histogram(rrms; label = "iso", bins = nbins, color = :lightgray)
h_num = histogram(num_cells; label = "iso", bins = nbins, color = :lightgray)

histogram!(h_rms, rrms; label = "ani", bins = nbins, color = :coral)
histogram!(h_num, num_cells; label = "ani", bins = nbins, color = :coral)

# Posterior

scalar_fields = ["dlnVs"]
directional_fields = ["fsh", "psi", "gamma"] # Compute moments for these directional fields; Vector{String}
spherical_conversion = (a, b, c) -> (a, b, asin(c)) # Generic function
credible_level = 0.68
q_quantiles = [0.25, 0.5, 0.75]

# (1) Southern Slab, (2) Northwestern Gorda, (3) Northeastern Cobb, (4) hole, (5) NE null
sample_points = [
    6205.07 6202.5 6190.33 6069.18 6142.25;
    97.4771 -379.361 -585.158 148.328 656.418;
    -433.954 -293.049 195.406 10.5959 736.59
]

markov_chains = load_markov_chains(num_chains, num_burn; chain_directory = chain_directory, chain_name = chain_name)
iso_nn_trees = build_nn_trees(markov_chains, scalar_fields; leafsize = 10)
ani_nn_trees = build_nn_trees(markov_chains, directional_fields; leafsize = 10)

# Define Priors
nsamp = length(markov_chains[1])*num_chains
iso_prior = 0.08*(2.0*rand(1,nsamp) .- 1.0)
ani_prior = rand(3, nsamp)
ani_prior[1,:] .*= 0.08
ani_prior[2,:] .= pi*(2.0*ani_prior[2,:] .- 1.0)
ani_prior[3,:] .= (2.0*ani_prior[3,:] .- 1.0)
directional_posterior!(ani_prior; spherical_conversion = spherical_conversion)
ani_prior[3,:] .*= (180.0/pi)

iso_samples = collect_posterior_samples(sample_points, markov_chains, iso_nn_trees, scalar_fields)
ani_samples = collect_posterior_samples(sample_points, markov_chains, ani_nn_trees, directional_fields)
for xi in ani_samples
    directional_posterior!(xi; spherical_conversion = spherical_conversion)
    xi .*= (180.0/pi)
end

# Extract posterior to plot
ifld, jpt, fld_label = 1, 1, "dlnVs"
@views x, x0 = iso_samples[jpt][ifld, :], iso_prior[ifld, :]

# Plot posterior

nbins = 15
h = histogram(x0; label = "prior", bins = nbins, color = :lightgray);
histogram!(h, x; label = fld_label, bins = nbins, color = :slateblue) # :firebrick, :seagreen, :slateblue
histogram!(h, x0; label = nothing, bins = nbins, color = :lightgray, seriesalpha = 0.5)

# Compute credible intervals
iqt = 1
qi = quantile(x, q_quantiles[iqt])
mi = median_interval(x, credible_level; is_sorted = false)
lhi = low_high_interval(x, credible_level; is_sorted = false)
hdi = highest_density_interval(x, credible_level; is_sorted = false)

int = lhi
bheight = 125 #2.0*nsamp/nbins
h = histogram(x; label = fld_label, bins = nbins, color = :slateblue);
plot!(h, [int[1], int[1]], [0, bheight], linewidth = 2, color = :red, label = nothing);
plot!(h, [int[2], int[2]], [0, bheight], linewidth = 2, color = :blue, label = nothing)