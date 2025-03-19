######################
# Convert Data Files #
######################
# This routine converts source, receiver, and observation files (e.g., those from PSI_D) into
# the correct format for PSI_S. Specifically, stations are re-assigned sequential integer IDs
# (corresponding station names in data file are also updated). Proper phase naming (in the case
# TauP is used for travel-times/ray tracing) is also checked.

# Set TauP related environment variables if necessary
# ENV["TAUP_JAR"] = "/Users/bvanderbeek/research/software/TauP_Toolkit/TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
# ENV["JULIA_COPY_STACKS"] = 1 # For non-windows users
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/pre_processing.jl") # Include the pre-processing functions
using Plots

# Extract path to SinkingSlab example directory
cd(@__DIR__)

# START INPUT #

# Input Data Files
dlm = "," # What is used to deliminate columns in the below data files (isspace or ",")
event_file = "../../../../PSI_D/examples/iso_local_eq/input/Sources.dat"
station_file = "../../../../PSI_D/examples/iso_local_eq/input/Receivers.dat"
data_file = "../../../../PSI_D/examples/iso_local_eq/input/fast_block_50ms_TravelTime_CompressionalWave.dat"
# Does this data file contain travel-times that need to be converted to delay times?
tf_convert_to_delays = false
# Need to know column indexing in the data files
j_id, j_lat, j_lon, j_elv = 1, 3, 2, 4 # At which columns are the IDs, latitude, longitude, and elevation/depth found in the event_, station_file?
j_evt, j_sta, j_obs, j_phs = 5, 6, 1, 4 # At which columns are the event, station, data, and phase found in the data_file?
# Define true if the input file stores depths; false if it stores elevation. Defined for both event and station files
tf_event_depth, tf_station_depth = false, false
# What type of data, Int or String, represents the event (eid_type) and station (sid_type) IDs?
# Often events are identified by integers (Int) while station IDs are characters (String).
eid_type, sid_type = Int, String
# TauP formatted 1D model (or recognized TauP model name, e.g., "ak135")
taup_model = "../../../../PSI_D/examples/iso_local_eq/input/etna_ak135_ext3500m.nd"
# If this TauP model is extended above sealevel, define dR > 0.0; otherwise dR = 0.0
# This is used to define the appropriate source/receiver depths in the TauP model
dR = 3.5;
# This script will write the following new input files
out_event = "../input/evt.dat"
out_station = "../input/sta.dat"
out_data = "../input/ttp_iso_fast_block_50ms.dat"

# END INPUT #

# Load data into structures
Events = read_aquisition_file(event_file; id_type = eid_type, order = (j_id, j_lat, j_lon, j_elv), tf_depth = tf_event_depth, dlm = dlm)
Stations = read_aquisition_file(station_file; id_type = sid_type, order = (j_id, j_lat, j_lon, j_elv), tf_depth = tf_station_depth, dlm = dlm)
Data = read_observation_file(data_file; eid_type = eid_type, sid_type = sid_type, order = (j_evt, j_sta, j_obs, j_phs), dlm = dlm)

# Check for bad travel-time predictions. This function will throw a warning if there are any bad predictions.
# Returns the logical index (kbad) of bad data and the predicted travel-times (tt), range (dist), back-azimuth (baz),
# and output phase (phase_out) for all arrivals.
kbad, tt1D, dist, baz, phase_out = check_data(Events, Stations, Data, taup_model; dR = dR)

if any(isinf.(tt1D)) || any(isnan.(tt1D))
    error("Missing travel-times!")
end

# Update phase names
if !all(Data.phase .== phase_out)
    @warn "Modifying phase names."
    Data.phase .= phase_out
end

# Convert to delays for teleseismic data
if tf_convert_to_delays
    @warn "Converting travel-times to delays."
    Data.observation .-= tt1D # Subtract 1D predictions from data
end

histogram(Data.observation .- tt1D, title = "Absolute Delays")

# Write new input files
write_psi_s_events(out_event, Events)
write_psi_s_stations(out_station, Stations)
write_psi_s_observations(out_data, Data, Events, Stations)