# Pre-Processing functions for PSI_S
# using TauP
using Printf

# Compares up and downgoing phases and selects first arriving
function check_taup_predictions(Events, Stations, Data, ref_model, phase_up, phase_down; dR = 0.0, Rref = 6371.0)
    tf_extended = dR > 0.0
    nobs = length(Data.observation)
    kbad, sld = falses(nobs), zeros(nobs)
    tt, dist, baz, phase_out = zeros(nobs), zeros(nobs), zeros(nobs), Vector{String}(undef, nobs)

    TimeObj = buildTimeObj(ref_model)
    for k in eachindex(Data.observation)
        # Indexing
        evt_k, sta_k = Data.event_id[k], Data.station_id[k] # Event and station ID
        j_evt, j_sta = Events.position[evt_k], Stations.position[sta_k] # Event and station linear index
        # Event-station distance
        evt_lat, evt_lon, evt_elv = Events.latitude[j_evt], Events.longitude[j_evt], Events.elevation[j_evt]
        sta_lat, sta_lon, sta_elv = Stations.latitude[j_sta], Stations.longitude[j_sta], Stations.elevation[j_sta]
        dist[k], baz[k] = TauP.taup_geoinv(sta_lat, sta_lon, evt_lat, evt_lon)
        # Adjusted source and receiver depths
        evt_dpt, sta_dpt = tf_extended ? (dR - evt_elv, dR - sta_elv) : (-evt_elv, 0.0)

        # Compare phases
        tt_u, _ = taup_time!(TimeObj, phase_up, dist[k], evt_dpt, sta_dpt)
        tt_u = isnan(tt_u) ? Inf : tt_u
        tt_d, _ = taup_time!(TimeObj, phase_down, dist[k], evt_dpt, sta_dpt)
        tt_d = isnan(tt_d) ? Inf : tt_d
        # Select first arriving phase
        tt[k], phase_out[k] = tt_d < tt_u ? (tt_d, phase_down) : (tt_u, phase_up)
        kbad[k] = isinf(tt[k]) || isnan(tt[k]) ? true : false

        # Compute straight-line-distance
        r_1, r_2 = Rref + evt_elv, Rref + sta_elv
        sinϕ_1, cosϕ_1 = sincosd(evt_lat)
        sinλ_1, cosλ_1 = sincosd(evt_lon)
        sinϕ_2, cosϕ_2 = sincosd(sta_lat)
        sinλ_2, cosλ_2 = sincosd(sta_lon)
        x_1, y_1, z_1 = r_1*cosϕ_1*cosλ_1, r_1*cosϕ_1*sinλ_1, r_1*sinϕ_1
        x_2, y_2, z_2 = r_2*cosϕ_2*cosλ_2, r_2*cosϕ_2*sinλ_2, r_2*sinϕ_2
        sld[k] = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2)
    end
    nbad = sum(kbad)
    nbad > 0 && (@warn "There are "*string(nbad)*" bad events!")

    return kbad, tt, dist, baz, phase_out, sld
end

# 'check_data' will loop over the arrivals and check that TauP can predict a travel-time
# for every event-station-phase pair. It will return a boolean vector (kbad) that is true
# at the indices of arrivals that could not be predicted, a modified list of phase names
# (phase_out), and the distance (dist) and back-azimuth (baz) of the arrivals.
function check_data(Events, Stations, Data, ref_model; dR = 0.0, AltPhase = return_alt_phases())
    tf_extended = dR > 0.0 # If true, event depths, station elevations, and surface reflection phase names will be modified
    itag = dR - floor(Int, dR) == 0.0 ? "^"*string(floor(Int, dR)) : "^"*string(dR) # Modifier for surface reflection phase names
    
    nobs = length(Data.observation)
    kbad = falses(nobs)
    tt, dist, baz, phase_out = zeros(nobs), zeros(nobs), zeros(nobs), Vector{String}(undef, nobs)

    TimeObj = buildTimeObj(ref_model)
    for (k, phase_k) in enumerate(Data.phase)
        # If using extended model, modify surface reflected phase names accordingly (e.g. pP becomes p^9P)
        if tf_extended
            tf_surf_refl = (length(phase_k) > 1) && (phase_k[1:2] == "pP" || phase_k[1:2] == "sS")
            tf_surf_refl && (phase_k = phase_k[1]*itag*phase_k[2:end])
            tf_surf_refl && @warn("Surface reflection phase name modified!")
        end
        phase_out[k] = phase_k
        # Indexing
        evt_k, sta_k = Data.event_id[k], Data.station_id[k] # Event and station ID
        j_evt, j_sta = Events.position[evt_k], Stations.position[sta_k] # Event and station linear index
        # Event-station distance
        evt_lat, evt_lon, evt_elv = Events.latitude[j_evt], Events.longitude[j_evt], Events.elevation[j_evt]
        sta_lat, sta_lon, sta_elv = Stations.latitude[j_sta], Stations.longitude[j_sta], Stations.elevation[j_sta]
        dist[k], baz[k] = TauP.taup_geoinv(sta_lat, sta_lon, evt_lat, evt_lon)
        # Predict traveltime
        evt_dpt, sta_dpt = tf_extended ? (dR - evt_elv, dR - sta_elv) : (-evt_elv, 0.0)
        tt_k, _ = taup_time!(TimeObj, phase_k, dist[k], evt_dpt, sta_dpt)
        if isnan(tt_k) # Try an alternative phase
            phase_alt = (dist[k] >= 160.0) && (phase_k == "P") ? "PKIKP" : AltPhase[phase_k] # AD-HOC ASSUMPTION!!!
            println("Trying phase: ", phase_alt)
            tt_k, _ = taup_time!(TimeObj, phase_alt, dist[k], evt_dpt, sta_dpt)
            if isnan(tt_k)
                kbad[k] = true
                println("Phase, range, event depth, station depth: ",phase_k,", ",dist[k],", ",evt_dpt,", ",sta_dpt)
            else
                phase_out[k] = phase_alt
            end
        end
        tt[k] = tt_k
    end
    nbad = sum(kbad)
    nbad > 0 && (@warn "There are "*string(nbad)*" bad events!")

    return kbad, tt, dist, baz, phase_out
end

# 'return_alt_phases' creates a dictionary of alternative phase names to try when TauP fails to predict an arrival.
function return_alt_phases()
    AltPhase = Dict{String, String}()
    AltPhase["P"] = "Pdiff"
    AltPhase["S"] = "Sdiff"
    AltPhase["Pdiff"] = "P"
    AltPhase["Sdiff"] = "S"

    AltPhase["pP"] = "pPdiff"
    AltPhase["sS"] = "sSdiff"
    AltPhase["pPdiff"] = "pP"
    AltPhase["sSdiff"] = "sS"

    AltPhase["p^9P"] = "p^9Pdiff"
    AltPhase["s^9S"] = "s^9Sdiff"
    AltPhase["p^9Pdiff"] = "p^9P"
    AltPhase["s^9Sdiff"] = "s^9S"
    return AltPhase
end

# 'read_aquisition_file' will read an event or station file formated for PSI_S/D.
function read_aquisition_file(aquisition_file; id_type = Int, order = (1,2,3,4), tf_depth = false, dlm = isspace)
    tf_string_id = id_type == String # Treat IDs as Strings?
    position = Dict{id_type, Int}() # Dictionary to index of specific ID
    id, lat, lon, elv = Vector{id_type}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
    j_id, j_lat, j_lon, j_elv = order # Column indexing order

    # Loop over lines in data file
    k, nline, num_col = 0, 0, maximum(order)
    for line in readlines(aquisition_file)
        nline += 1
        line = split(line, dlm; keepempty = false)
        if length(line) >= num_col
            k += 1
            id_k = tf_string_id ? string(strip(line[j_id])) : parse(id_type, line[j_id])
            position[id_k] = k
            push!(id, id_k)
            push!(lat, parse(Float64, line[j_lat]))
            push!(lon, parse(Float64, line[j_lon]))
            push!(elv, parse(Float64, line[j_elv]))
        else
            @warn "Skipping line " * string(nline) * ". Expected " * string(num_col) * " lines but found " * string(length(line)) * "!"
        end
    end
    if tf_depth
        elv .*= -1.0
    end

    return (position = position, id = id, latitude = lat, longitude = lon, elevation = elv) # Returns NamedTuple
end
# 'read_observation_file' will read a data file formated for PSI_S/D
function read_observation_file(observation_file; eid_type = Int, sid_type = String, order = (1,2,3,4), dlm = isspace)
    tf_string_event_id = eid_type == String # Treat event IDs as Strings?
    tf_string_station_id = sid_type == String # Treat station IDs as Strings?
    evt_id, sta_id, b, phs = Vector{eid_type}(), Vector{sid_type}(), Vector{Float64}(), Vector{String}()
    j_evt, j_sta, j_obs, j_phs = order

    # Loop over lines in data file
    k, nline, num_col = 0, 0, maximum(order)
    for line in readlines(observation_file)
        nline += 1
        line = split(line, dlm; keepempty = false)
        if length(line) >= num_col
            k += 1
            eid_k = tf_string_event_id ? string(strip(line[j_evt])) : parse(eid_type, line[j_evt])
            sid_k = tf_string_station_id ? string(strip(line[j_sta])) : parse(sid_type, line[j_sta])
            push!(evt_id, eid_k)
            push!(sta_id, sid_k)
            push!(b, parse(Float64, line[j_obs]))
            push!(phs, string(strip(line[j_phs])))
        else
            @warn "Skipping line " * string(nline) * ". Expected " * string(num_col) * " lines but found " * string(length(line)) * "!"
        end
    end

    return (event_id = evt_id, station_id = sta_id, observation = b, phase = phs) # Returns NamedTuple
end
# 'read_observation_file' will read a data file formated for PSI_D.
function read_psi_d_observation_file(observation_file; eid_type = Int, sid_type = String, dlm = isspace, num_col = 7)
    tf_string_event_id = eid_type == String # Treat event IDs as Strings?
    tf_string_station_id = sid_type == String # Treat station IDs as Strings?
    obs, unc, prd = Vector{Float64}(),  Vector{Float64}(), Vector{Float64}()
    phs, evt_id, sta_id, chan = Vector{String}(), Vector{eid_type}(), Vector{sid_type}(), Vector{String}()

    # Loop over lines in data file
    k, nline = 0, 0
    for line in readlines(observation_file)
        nline += 1
        line = split(line, dlm; keepempty = false)
        if length(line) == num_col
            k += 1
            eid_k = tf_string_event_id ? string(strip(line[5])) : parse(eid_type, line[5])
            sid_k = tf_string_station_id ? string(strip(line[6])) : parse(sid_type, line[6])
            push!(obs, parse(Float64, line[1]))
            push!(unc, parse(Float64, line[2]))
            push!(prd, parse(Float64, line[3]))
            push!(phs, string(strip(line[4])))
            push!(evt_id, eid_k)
            push!(sta_id, sid_k)
            push!(chan, string(strip(line[7])))
        else
            @warn "Skipping line " * string(nline) * ". Expected " * string(num_col) * " lines but found " * string(length(line)) * "!"
        end
    end

    return (observation=obs, uncertainty=unc, period=prd, phase=phs,
        event_id=evt_id, station_id=sta_id, channel=chan) # Returns NamedTuple
end
# Convenience function to delete a set of IDs missing from a reference list
function subset_aquisition_structure!(Aqui, id_list)
    # Flag events and stations missing from dataset
    tf_out = trues(size(Aqui.id))
    for id in id_list
        pos = Aqui.position[id]
        tf_out[pos] = false
    end
    # Delete position maps
    [out_i && delete!(Aqui.position, Aqui.id[i]) for (i, out_i) in enumerate(tf_out)]
    # Remove missing ids
    deleteat!(Aqui.id, tf_out)
    deleteat!(Aqui.latitude, tf_out)
    deleteat!(Aqui.longitude, tf_out)
    deleteat!(Aqui.elevation, tf_out)
    # Re-map positions
    [Aqui.position[id] = k for (k, id) in enumerate(Aqui.id)]

    return nothing
end



# These functions write properly formated PSI_S data files.
function write_psi_s_events(out_file, Events; tf_overwrite = false)
    !tf_overwrite && isfile(out_file) && error("Output file "*out_file*" already exists!")
    io = open(out_file, "w")
    for k in eachindex(Events.id)
        @printf(io, "%5i %12.6f %12.6f %12.6f %12.6f \n", k, Events.latitude[k], Events.longitude[k], -Events.elevation[k], 0.0) # Depth and a null origin time
    end
    close(io)
    return nothing
end
function write_psi_s_stations(out_file, Stations; tf_overwrite = false)
    !tf_overwrite && isfile(out_file) && error("Output file "*out_file*" already exists!")
    io = open(out_file, "w")
    for k in eachindex(Stations.id)
        @printf(io, "%5i %12.6f %12.6f %12.6f \n", k, Stations.latitude[k], Stations.longitude[k], Stations.elevation[k]) # Elevation!
    end
    close(io)
    return nothing
end
function write_psi_s_observations(out_file, Data, Events, Stations; tf_overwrite = false)
    !tf_overwrite && isfile(out_file) && error("Output file "*out_file*" already exists!")
    io = open(out_file, "w")
    for (k, b) in enumerate(Data.observation)
        evt_k, sta_k = Data.event_id[k], Data.station_id[k]
        i_evt, j_sta = Events.position[evt_k], Stations.position[sta_k]
        @printf(io, "%5i %5i %10.4f %5s \n", i_evt, j_sta, b, Data.phase[k])
    end
    close(io)
    return nothing
end
function write_psi_s_observations(out_file, Data; tf_overwrite = false)
    !tf_overwrite && isfile(out_file) && error("Output file "*out_file*" already exists!")
    io = open(out_file, "w")
    for (k, b) in enumerate(Data.observation)
        evt_k, sta_k = Data.event_id[k], Data.station_id[k]
        @printf(io, "%5i %5i %10.4f %5s \n", evt_k, sta_k, b, Data.phase[k])
    end
    close(io)
    return nothing
end


function compute_event_demeaned_delays(Data, Events; delays = Data.observation)
    rdt = zeros(size(delays))
    med, num = compute_event_means(Data, Events; delays = delays)
    for (k, v) in enumerate(delays)
        evt_id = Data.event_id[k]
        ievt = Events.position[evt_id]
        rdt[k] = v - med[ievt]
    end
    return rdt, num
end

function compute_event_means(Data, Events; delays = Data.observation)
    med = zeros(length(Events.id))
    num = zeros(length(Events.id))
    for (k, v) in enumerate(delays)
        if !isnan(v)
            evt_id = Data.event_id[k]
            ievt = Events.position[evt_id]
            med[ievt] += v
            num[ievt] += 1.0
        end
    end
    med ./= num
    
    return med, num
end

function compute_station_means(Data, Stations; delays = Data.observation)
    msd = zeros(length(Stations.id))
    num = zeros(length(Stations.id))
    for (k, v) in enumerate(delays)
        if !isnan(v)
            sta_id = Data.station_id[k]
            jsta = Stations.position[sta_id]
            msd[jsta] += v
            num[jsta] += 1.0
        end
    end
    msd ./= num

    return msd, num
end