
# Set-up
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/dependencies.jl")

cd("/Users/bvanderbeek/research/software/GitRepos/PSI_S/examples/iso_local_eq")
parameter_file = "input/psi_s_parameter_file.toml"

P, IP_obs, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager = build_inputs(parameter_file; tf_save = false, tf_serial = true);