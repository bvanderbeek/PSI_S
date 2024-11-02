module PSI_S

include("dependencies.jl")
println("PSI_S: Plateform for Seismic Imaging - Stochastic")

export psi_s, run_psi_s
export load_chains_metrics, build_ensemble_model, build_ensemble_correlation, build_ensemble_tensor, build_posterior, write_ensemble_model_vtk

# No idea how to get this package to load with parallelisation from Distributed.jl

# @everywhere begin
#     using Pkg
#     pkg_path, _ = splitdir(@__DIR__)
#     Pkg.activate(pkg_path)
#     Pkg.instantiate()
#     Pkg.precompile()
# end

# @everywhere begin
#     include("dependencies.jl")
#     export psi_s, run_psi_s
#     export load_chains_metrics, build_ensemble_model, build_ensemble_correlation, build_ensemble_tensor, build_posterior, write_ensemble_model_vtk
# end

# Main function to run inversion...idea was to have processes added from within the function...but not working
function psi_s(parameter_file::String, ichunk; tf_serial = false, tf_save = true)
    # Make input structures
    P, _, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager =
        build_inputs(parameter_file; tf_save = tf_save, tf_serial = tf_serial)
    # Initialise parallel processes...cannot do from inside the package...how?
    # num_chains = P["MonteCarloSolver"]["num_chains"]
    # addprocs(num_chains)
    # @everywhere using PSI_S...this will fail
    # @everywhere include(@__DIR__() * "/dependencies.jl")
    # @eval @everywhere using $(ModuleName(@__MODULE__))
    # Run serial or parallel Markov chains
    if tf_serial
        println("Serial rj-McMC inversion.")
        println("Chunks: ",P["MonteCarloSolver"]["num_chunks"],". Chains per chunk: ", P["MonteCarloSolver"]["num_chains"])
        for nchunk in 1:P["MonteCarloSolver"]["num_chunks"]
            run_psi_s(nchunk, P, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager)
        end
    else
        println("Parallel rj-McMC inversion.")
        println("Chunk: ",ichunk,". Parallel chains in chunk: ", P["MonteCarloSolver"]["num_chains"])
        run_psi_s(ichunk, P, IP_fields, IP, rnodes, rays, evtsta, observables, LocalRaysManager)
    end

    return nothing
end

end
