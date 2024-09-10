## SETUP ##
### Define TauP jar-file path (if necessary) ###
# ENV["TAUP_JAR"] = "/u/unipd/vanderbeek/TauP_Toolkit/TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
### Using PSI_S Package ###
# using JLD
# using Distributed
# addprocs(parse(Int64, ARGS[2]))
# @everywhere using PSI_S
### Using Local PSI_S Package ###
# ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
# using Pkg
# using JLD
# using Distributed
# addprocs(parse(Int64, ARGS[2]))
# Pkg.activate(ENV["PSI_S"])
# @everywhere using PSI_S
### Using PSI_S Source Code ###
ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
using Distributed
addprocs(parse(Int64, ARGS[2]))
@everywhere include(ENV["PSI_S"]*"/src/dependencies.jl")

const wrk_dir = @__DIR__
const name = string(ARGS[1])
const chain_id = parse(Int64,string(ARGS[3]))
run_RJMCMC(wrk_dir, name, chain_id)
