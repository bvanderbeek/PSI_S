## SETUP ##
### Define TauP jar-file path (if necessary) ###
# ENV["TAUP_JAR"] = "/u/unipd/vanderbeek/TauP_Toolkit/TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
### Using PSI_S Package ###
# using JLD
# using PSI_S
### Using Local PSI_S Package ###
# ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
# using JLD
# using Pkg
# Pkg.activate(ENV["PSI_S"])
# using PSI_S
### Using PSI_S Source Code ###
ENV["PSI_S"] = "/storage2/unipd/vanderbeek/PSI_S"
include(ENV["PSI_S"]*"/src/dependencies.jl")
### Always include datafields.jl ###
wrk_dir = @__DIR__
include(wrk_dir*"/datafieldsIP.jl")

################################
## -- parameters selection -- ##
################################

# -- DEFINE GENERAL PARAMETERS FOR SEISMIC IMAGING
IP = IPConst(
    "test",                                 # -- name of the run -> Create output/name directory to store results
	wrk_dir*"/input/sources_melt.dat",      # -- events data file [InputData/ implicit]
	wrk_dir*"/input/receivers_melt.dat",    # -- stations data file [InputData/ implicit]
	wrk_dir*"/input/iasp91.txt",            # -- 1D velocity model [VelocityModels/ implicit]                              
	DomainGeoBoundaries(                    # -- DEFINE GEOGRAPHIC BOUNDARIES FOR THE INVERSION DOMAIN
        [-20.2,-10.6],                      # -- latitude limits [°]
        [-121.4,-102.6],                    # -- longitude limits [°]
        [-600.0,0.0]                        # -- depth limits [km]
    ),      
    MonteCarloSolver(           # -- DEFINE PARAMETERS OF RJMCMC SOLVER
        2,                      # -- misfit function L-n norm (L1, L2...)
        1e5,                    # -- total number of iterations
        4,                      # -- number of gianMarkov chains
        "linear",               # -- Nearest-Neighbour(NN) interpolation algorithm ("linear" -> linear search,"kdtree" -> NearestNeighbors.jl [deprecated])
        0.01,                   # -- perturb.size (fraction of lims provided for each voronoi diagram in datafieldsIP.jl)
        0,                      # -- initial step-1 range
        "uniform",              # -- prior choice for the number of nuclei ("uniform","log-uniform")   
        100,                    # -- maximum number of nuclei [memory usage depends on this limit]
        10,                     # -- initial number of nuclei
        false,                  # -- hierarchical Bayesian inversion (include noise random variables in posterior definition for each obs.)
        false,                  # -- initial random models in chain (avoid for teleseismic imaging...)
        false,                  # -- squeezing (outside Voronoi diagrams the reference values are interpolated)
        1e2,                    # -- interval to print out on progress' state 
        1e4                     # -- saving interval (every n iterations the current model is saved)
    ),   
    DelayedRejection(           # -- DEFINE PARAMETERS OF DELAYED REJECTION ALGORITHM (PERTURB VALUE/POSITION STEPS)
        false,                  # -- status
        4                       # -- proposal st.dev rescaling
    ),
    ParallelTempering(          # -- DEFINE PARAMETERS OF PARALLEL TEMPERING
        false,                  # -- enable parallel tempering
        1e3,                    # -- first P.T. swaps after n iterations
        1.0,                    # -- cold chains temperature
        4,                      # -- number of cold chains (T=1)
        10,                     # -- hottest chain temperature
        1e3,                    # -- swaps every n iterations
        true                    # -- increase prosal perturbation size as sqrt(T) [check maximum...]
    ),  
    StaticCorrections(          # -- DEFINE PARAMETERS OF EVENT/STATION STATIC CORRECTIONS INVERSION
        5e4,                    # -- first static inversions after n iteration [set 1 to run at the beginning]
        1e4,                    # -- static inversion runs every n iteration
        1.0e6                   # -- expected ratio data/station-statics uncertainties [damping parameter for station statics]
    ),
    RayTracing(                 # -- DEFINE PARAMETERS OF RAY-TRACING (TELESEISMIC AND LOCAL EARTHQUAKES DATA)
        10,                     # -- ray-paths discretization length [km]
        "linear",               # -- interpolation along ray-paths ["linear","spline"]
        "iasp91",               # -- TauP velocity model used for teleseismic ray-tracing [if used]
        true,                   # -- trace all the rays using TauP
        6371.0,                 # -- maximum Earth's radius according to TauP extended velocity model
        [250,250,150],          # -- nodes along each coordinate (latitude, longitude, radius)
        5,                      # -- forward star level
        false,                  # -- noise perturbation to grid's nodes 
        false,                  # -- grid-carving, might reduce significantly Dijkstra's execution time
        1e5                     # -- ray-tracing running every n iterations
    ),                                        
    EarthquakesRelocation(      # -- DEFINE PARAMETERS OF EARTHQUAKES RELOCATION (ONLY FOR 3D LET) [EXPERIMENTAL]
        false,                  # -- enable earthquakes relocation
        3                       # -- relocation runs every n ray-tracing
    ),
    Bayesian4Dimaging(          # -- DEFINE PARAMETERS OF 4D IMAGING [EXPERIMENTAL]
        false                   # -- execute full-4D Bayesian Inversion
    )
)     


##########################
## -- initialization -- ##
##########################

print("\ninitializing the inverse problem...\n")
IP_fields = build_fieldslist()
IP_obs = build_obslist()
rnodes, rays, evtsta, observables, LocalRaysManager = initialize_IP(IP, IP_obs, IP_fields)

mkpath(wrk_dir*"/output/"*IP.name)
save(
    string(wrk_dir, "/output/", IP.name, "/", IP.name, ".jld"), "IP", IP,
    "IP_fields", IP_fields, 
    "rnodes", rnodes, 
    "rays", rays, 
    "evtsta", evtsta, 
    "observables", observables,
    "LocalRaysManager", LocalRaysManager
)

# Save copy of input parameters
cp(wrk_dir*"/buildIP.jl", wrk_dir*"/output/"*IP.name*"/buildIP.jl")
cp(wrk_dir*"/datafieldsIP.jl", wrk_dir*"/output/"*IP.name*"/datafieldsIP.jl")
