module PSI_S
println("PSI_S: Plateform for Seismic Imaging - Stochastic")

include("dependencies.jl")

# Structures
export IPConst, DomainGeoBoundaries, MonteCarloSolver, ParallelTempering, StaticCorrections, RayTracing, EarthquakesRelocation
export Bayesian4Dimaging, fieldinfo, ModelConst, GridConst
# Functions
export initialize_IP, add_obs, obsinfo, add_field, run_RJMCMC, tt_iso_P_dlnVp

end
