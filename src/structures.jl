struct DomainGeoBoundaries
    lat::Vector{Float64}
    lon::Vector{Float64}
    depth::Vector{Float64}
end
function DomainGeoBoundaries(D::Dict)
    return DomainGeoBoundaries(D["lat_lims"], D["lon_lims"], D["elv_lims"])
end

struct MonteCarloSolver
    Lnorm::Int64
    niterations::Int64
    nchains::Int64
    NN_algorithm::String
    pert_size::Float64
    pertv_initrange::Int64
    nuclei_prior::String
    max_nuclei::Int64
    ini_nuclei::Int64
    hierarchical::Bool
    rand_init_values::Bool
    squeezing::Bool
    printon::Int64
    saveint::Int64
end
function MonteCarloSolver(D::Dict)
    return MonteCarloSolver(D["l_norm"], D["num_iterations"], D["num_chains"], D["nn_algorithm"], D["pert_size"],
        D["vpert_iterations"], D["dim_prior"], D["max_dim"], D["init_dim"], D["tf_hierarchical"], D["tf_rand_init"],
        D["tf_squeeze"], D["print_interval"], D["save_interval"])
end

struct DelayedRejection
    status::Bool
    rescale::Float64
end
function DelayedRejection(D::Dict)
    return DelayedRejection(D["tf_delayed_rejection"], D["rescale_value"])
end

struct ParallelTempering
    status::Bool
    pt_pause::Int64
    Tcold::Float64
    ncold::Int64
    maxT::Float64
    intsw_its::Int64
    temp_pert::Bool
end
function ParallelTempering(D::Dict)
    return ParallelTempering(D["tf_parallel_tempering"], D["swap_start"], D["min_temp"], D["num_cold"],
        D["max_temp"], D["swap_interval"], D["tf_increase_pert"])
end

struct StaticCorrections
    statics_pause::Int64
    statics_its::Int64
    ss_uncertainty::Float64
end
function StaticCorrections(D::Dict)
    return StaticCorrections(D["static_start"], D["static_interval"], D["static_damping"])
end

struct RayTracing
    dray::Float64
    itp::String
    TauP_Model::String
    allTauP::Bool
    Rmax::Float64
    nnodes::Vector{Int64}
    fw_level::Int64
    perturb::Bool
    carving::Bool
    sub_its::Int64
    topography_file::String
    velmod_shift::Bool
end
function RayTracing(D::Dict)
    if D["tf_1D_rays"]
        V = load_taup_model(D["taup_model"])
        Rmodel = V.r[end]
    else
        Rmodel = 6371.0
    end
    return RayTracing(D["sampling_interval"], D["interp_method"], D["taup_model"], D["tf_1D_rays"], Rmodel,
        D["grid_dims"], D["forward_star"], D["tf_grid_noise"], D["tf_grid_carve"], D["tracing_interval"],
        D["topography_file"], D["tf_hang_model"])
end

struct EarthquakesRelocation
    relocate::Bool
    relocations_init::Bool
    relocations_pause::Int64
    relocations_its::Int64
    fw_reloc::Int64
    wadati_init::Bool
end
function EarthquakesRelocation(D::Dict)
    return EarthquakesRelocation(D["tf_run_relocation"], D["tf_initial_1D_relocation"], D["first_relocation"],
    D["relocation_interval"], D["search_width"], D["tf_wadati_init"])
end

struct Bayesian4Dimaging
    status::Bool
end
function Bayesian4Dimaging(D::Dict)
    return Bayesian4Dimaging(D["tf_4D_imaging"])
end

# -- IP parameters
struct IPConst
    name::String
	InputEvt::String
	InputSta::String
    velocitymodel::String
    lims::DomainGeoBoundaries
    MCS::MonteCarloSolver
    DR::DelayedRejection
    PT::ParallelTempering
    SC::StaticCorrections
    RayTracingInit::RayTracing
    EQRLOC::EarthquakesRelocation
    B4DI::Bayesian4Dimaging
end
function IPConst(P::Dict)
    return IPConst(P["run_id"], P["AquisitionFiles"]["events_file"], P["AquisitionFiles"]["stations_file"], P["AquisitionFiles"]["starting_model"],
        DomainGeoBoundaries(P["DomainBoundaries"]), MonteCarloSolver(P["MonteCarloSolver"]), DelayedRejection(P["DelayedRejection"]),
        ParallelTempering(P["ParallelTempering"]), StaticCorrections(P["StaticCorrections"]), RayTracing(P["RayTracing"]),
        EarthquakesRelocation(P["EarthquakeRelocation"]), Bayesian4Dimaging(P["Bayesian4D"]))
end

# -- model structures
struct Voronoi
    fieldname::String 
    prior::String       # -- prior choice     
    n::Vector{Int64}    # -- number of nuclei in the diagram
    c::Matrix{Float64}  # -- coords of the nuclei
    r::Matrix{Float64}  # -- radial coords on the nuclei
    v::Vector{Float64}  # -- values assigned to the nuclei
    nuclei2rays::Vector{Set{Int64}}
    nodes2nuclei::Vector{Int64}         # -- associates every ray-node the the influecing nucleus
    vlims::Vector{Float64}              # -- limits of field's values
    slims::Vector{Vector{Float64}}      # -- space-limits of diagram
    ref_value::Float64
end

struct DataPConst
    obsname::String
    estatics::Vector{Float64}
    sstatics::Vector{Float64}
    noise::Vector{Float64}
end

struct DataSpaceConst
    Obs::Vector{DataPConst}
end

struct ModelConst                      # -- House of the Voronoi Diagrams
    nfields::Int64
    fieldslist::Dict{String, Int64}
    fields::Vector{Voronoi}
    misfit::Vector{Float64}            # -- updated every time a new model is accepted
    rms::Vector{Float64}               # -- updated every time a new model is accepted
    T::Vector{Float64}                 # -- model temperature (PT)
    accepted::Vector{Int64}
    dataspace::DataSpaceConst
end




# -- data structures

struct ObsConst                         # -- Observable
    obsname::String                     # -- name
    forward_name::Symbol                # -- name of the function that computes predictions
    obsval::Vector{Float64}             # -- array of measured values
    prdval::Vector{Float64}             # -- array of predicted values (updated during the inversion)
    ref_t::Vector{Float64}
    evtids::Vector{Int64}               # -- true events' ids
    staids::Vector{Int64}               # -- true stations' ids
    obs2evt::Vector{Int64}              # -- maps the specific observed values to the event the belong to in estatics
    obs2sta::Vector{Int64}              # -- maps the specific observed values to the station the belong to in sstatics
    noise_guess::Float64                # -- noise guess (used to scale the hierarchical proposal)
    ray2obs::Dict{Int64, Int64}         # -- maps the ray-number of the observation in this structure
    demean::Bool
    solve_evt::Int64                    # -- it only saves if you are inverting for event statics for this obs
    solve_sta::Int64                    # -- it only saves if you are inverting for station statics for this obs
end

struct ObservablesConst                 # -- House of the Observables
    nobs::Int64
    obslist::Dict{String, Int64}
    Obs::Vector{ObsConst}
end


struct NodesConst                 # -- collection of ray-nodes
    c::Matrix{Float64}
    r::Matrix{Float64}            # -- radius of the ray-nodes
    r2n::Vector{Vector{Int64}}    # -- maps each ray to the indeces of the nodes in the arrays defined above
    n2r::Vector{Int64}            # -- maps each node to the ray it belong to
end

struct RayConst                         # -- structure for saving some utilities for predictions computation
    evt::Int64
    sta::Int64
    phase::String
    L::Vector{Float64}          # -- vector of lengths of ray-path segments
    v_1D_P::Vector{Float64}
    v_1D_S::Vector{Float64}
    ϕ::Vector{Float64}          # -- ray azimuth
    θ::Vector{Float64}          # -- ray elevation
    ζ::Vector{Float64}          # -- ray polarization
    t_1D::Vector{Float64}       # -- for tele phases this is the 1D travel time inside the inversion domain
    outdomain::Vector{Vector{Int64}}
end

struct RayTmpConst
    fields::Matrix{Float64}
    u::Vector{Float64}
    u_path::Vector{Float64}
    dt::Vector{Float64}
end

struct InfluenceConst
    rays::Vector{Bool}
    nodes::Vector{Bool}
    nuclei::Vector{Bool}
    neighbours::Set{Int64}
    new_neighborhood::Set{Int64}
end

struct TrialSampleConst
    n::Vector{Int64}
    dv1::Vector{Float64}
    dv2::Vector{Float64}
    ds1::Vector{Float64}
    ds2::Vector{Float64}
    misfit::Vector{Float64}
    term::Vector{Float64}
end

# -- carries some information needed by the local Voronoi updating scheme
struct UpdateConst              
    field::Vector{Int64}
    action::Vector{Int64}
    ind::Vector{Int64}
    obs_id::Vector{Int64}
    term::Vector{Float64}
    buffer::Vector{Vector{Float64}}
    pred::Vector{Vector{Float64}}
    noise::Vector{Vector{Float64}}
    statics_residuals::Vector{Float64}
    statics_values::Vector{Float64}
    Gg::SparseArrays.SPQR.QRSparse{Float64, Int64}
    Gt::SparseMatrixCSC
    influence::InfluenceConst
    old_voronoi::Voronoi
    misfit::Vector{Float64}
    rms::Vector{Float64}
    raytmp::RayTmpConst
    trial_sample::TrialSampleConst
end

struct fieldinfo
    name::String
    prior::String
    vlims::Vector{Float64}
    θlims::Vector{Float64}
    φlims::Vector{Float64}
    depthlims::Vector{Float64}
    timelims::Vector{Float64}
    init_val::Float64
end

struct obsinfo
    name::String
    file::String
    demean::Bool
    eventstatics::Int64
    stationstatics::Int64
    noise_guess::Float64
    forward_name::String
    # ... ?
end

struct fields_list
    fields::Vector{fieldinfo}
end

struct obs_list
    obs::Vector{obsinfo}
end

struct LocalRaysManagerConst
    status::Bool
    nnodes::Vector{Int64}
    fw_level::Int64
    perturb::Bool
    carving::Bool
    ray2source_receiver::Dict{Int64,Vector{Int64}}
    source2receivers::Dict{Int64,Vector{Int64}}
    pair2ray::Dict{Vector{Int64},Int64}
    local_evts::Set{Int64}
    local_stats::Set{Int64}
    source_nodes::Dict{Int64,Int64}
    receiv_nodes::Dict{Int64,Int64}
    sub_its::Int64
    relocations::Bool
    topography_status::Bool
    topography::Matrix{Float64}
end

mutable struct ObjsInChainConst
    model::ModelConst
    observables::ObservablesConst
    update::UpdateConst
end


# -- needed structures (something may change...!)

struct evt
    id::Int64
    lat::Float64
    lon::Float64
    depth::Float64
    T0::Float64
end

struct sta
    id::Int64
    lat::Float64
    lon::Float64
    elevation::Float64
end

struct evtstaConst
    evts::Vector{evt}
    stas::Vector{sta}
end

struct pathConst
    evt::Int64
    sta::Int64
    phase::String
    lat::Vector{Float64}
    lon::Vector{Float64}
    rad::Vector{Float64}
    t::Vector{Float64}
    sp_t::Vector{Float64}
end

function copy_model(model,update)
    old_voronoi = update.old_voronoi
    (; n, c, r, v, nuclei2rays, nodes2nuclei) = model.fields[update.field[1]]
    v .= 0.0
    r .= 0.0
    c .= 0.0
    n[1] = old_voronoi.n[1]
    for i in axes(old_voronoi.c,2)[1:n[1]]
        v[i] = old_voronoi.v[i]
        r[1,i] = old_voronoi.r[1,i]
        for j in axes(c,1)
            c[j,i] = old_voronoi.c[j,i]
        end
    end
    for (i,j) in enumerate(findall(update.influence.nuclei))
        nuclei2rays[j] = old_voronoi.nuclei2rays[i]
    end
    nodes2nuclei .= old_voronoi.nodes2nuclei
end

function output_model(model)
    vorovec = Vector{Voronoi}()
    for field in eachindex(model.fields)
        voronoi = model.fields[field]
        n = copy(voronoi.n)
        c = copy(voronoi.c[:,1:voronoi.n[1]])
        r = copy(voronoi.r[:,1:voronoi.n[1]])
        v = copy(voronoi.v[1:voronoi.n[1]])
        nodes2nuclei = zeros(Float64,1)
        vlims = voronoi.vlims
        slims = voronoi.slims
        ref_value = voronoi.ref_value
        nuclei2rays = Vector{Set{Int64}}()
        out_voronoi = Voronoi(
            voronoi.fieldname,
            voronoi.prior,
            n,
            c,
            r,
            v,
            nuclei2rays,
            nodes2nuclei,
            vlims,
            slims,
            ref_value
        )
        push!(vorovec,out_voronoi)
    end
    out_model = ModelConst(
        model.nfields,
        model.fieldslist,
        vorovec,
        copy(model.misfit),
        copy(model.rms),
        copy(model.T),
        copy(model.accepted),
        deepcopy(model.dataspace)
    )
    return out_model
end

