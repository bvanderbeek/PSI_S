##############
# Check Data #
##############
# This routine will check delay times to make sure they are predicted by TauP. If TauP fails,
# an attempt will be made to correct the phase name by selecting the first arriving P or S phase
# for the given event-station parameters. The modified data will then be written to file.
#
# Note! Unlike PSI_S, this routine can read event and station files in which the IDs are INTEGERS
# or STRING values in ANY order. When the new data files are written, these IDs will be
# modified according to PSI_S requirements (i.e. integers and sequential).
#
# Common Issues with TauP
# 1. Bad phase names
# 1.1. At long ranges, "P" becomes "Pdiff" (and "S" is "Sdiff")
# 1.2. For extended models, surface reflections should be renamed to account for the new reference sea-level
#      For example, "sS" becomes "s^9S" if the TauP model is radially extended by 9 km from sea-level.
# 2. Stations have non-zero elevation, but the TauP model referenced has not been radially extended to allow elevation
# 3. Strange TauP errors can occur when using extended models and TauP-2.6.1 or earlier versions.
# 3.1. Make sure your TAUP_JAR environment variable points to the correct jar-file version before using TauP.jl
# 4. Other strange TauP errors related to Java may occur for the following reasons
# 4.1. Mac users, Julia not started with flag --handle-signals=no 
# 4.2. For non-windows users, environment variable JULIA_COPY_STACKS = 1 not defined.
#
# Last important note! How were your delays computed? With respect to TauP travel-time predictions with or without elevation?
# For PSI_S, delays should be computed accounting for elevation! This is a bit of a pain because often times delays are measured
# with respect to a 1D TauP prediction without elevation. Additionally, if you want to change the 1D starting model, you also
# should recompute the delays. For these reasons it would be better to store travel-times and have PSI_S internally compute the
# delays with respect to the provided reference model...I will hopefully implement this soon.

# Set TauP related environment variables if necessary
# ENV["TAUP_JAR"] = "/Users/bvanderbeek/research/software/TauP_Toolkit/TauP-2.7.0-SNAPSHOT5/lib/TauP-2.7.0-SNAPSHOT5.jar"
# ENV["JULIA_COPY_STACKS"] = 1 # For non-windows users
ENV["PSI_S"] = "/Users/bvanderbeek/research/software/GitRepos/PSI_S"
include(ENV["PSI_S"]*"/src/pre_processing.jl") # Include the pre-processing functions
using Plots

# START INPUT #

# PSI_S input files
dlm = isspace # What is used to deliminate columns in the below data files (isspace or ",")
event_file = "Sources_Local_mIGUANA.dat"
station_file = "Receivers_Local_all.dat"
data_file = "SYN_TravelTime_CompressionalWave.dat"
# Does this data file contain travel-times that need to be converted to delay times?
tf_convert_to_delays = true
# Need to know column indexing in the data files
j_id, j_lat, j_lon, j_elv = 1, 3, 2, 4 # At which columns are the IDs, latitude, longitude, and elevation/depth found in the event_, station_file?
j_evt, j_sta, j_obs, j_phs = 5, 6, 1, 4 # At which columns are the event, station, data, and phase found in the data_file?
# Define true if the input file stores depths; false if it stores elevation. Defined for both event and station files
tf_event_depth, tf_station_depth = false, false
# What type of data, Int or String, represents the event (eid_type) and station (sid_type) IDs?
# Often events are identified by integers (Int) while station IDs are characters (String).
eid_type, sid_type = Int, String
# TauP formatted 1D model (or recognized TauP model name, e.g., "ak135")
taup_model = "/Users/bvanderbeek/research/software/TauP_Toolkit/custom_models/ak135_ext9km.tvel"
# If true, the above TauP model is assumed to be a custom model that has been radially extended by 9 km to include elevation in calculations
tf_extended_9km = true
# This script will write the following new input files
out_event = "0_evt.dat"
out_station = "0_sta.dat"
out_data = "0_delays.dat"

# END INPUT #

# Load data into structures
Events = read_aquisition_file(event_file; id_type = eid_type, order = (j_id, j_lat, j_lon, j_elv), tf_depth = tf_event_depth, dlm = dlm)
Stations = read_aquisition_file(station_file; id_type = sid_type, order = (j_id, j_lat, j_lon, j_elv), tf_depth = tf_station_depth, dlm = dlm)
Data = read_observation_file(data_file; eid_type = eid_type, sid_type = sid_type, order = (j_evt, j_sta, j_obs, j_phs), dlm = dlm)

# Check for bad travel-time predictions. This function will throw a warning if there are any bad predictions.
# Returns the logical index (kbad) of bad data and the predicted travel-times (tt), range (dist), back-azimuth (baz),
# and output phase (phase_out) for all arrivals.
kbad, tt, dist, baz, phase_out = check_data(Events, Stations, Data, taup_model; tf_extended_9km = tf_extended_9km)

# Update phase names
Data.phase .= phase_out

# Compute delays (if travel-time data was read)
if tf_convert_to_delays
    Data.observation .-= tt # Subtract 1D predictions from data
end

# Relative delay times (rdt) and number of stations that recorded each event (nsta)
rdt, nsta = compute_event_demeaned_delays(Data, Events)
# Relative mean station delay times (rmsd) and number of events recorded by each station (nevt)
rmsd, nevt = compute_station_means(Data, Stations; delays = rdt)
# Mean station delay times (msd); uses absolute delay times
msd, _ = compute_station_means(Data, Stations)

# Some quality control plots
color_limits, msize = (-0.4, 0.4), 5
scatter(Stations.longitude, Stations.latitude; markersize=msize, zcolor=msd,
    color=:balance, clims=color_limits, label=nothing, title="Absolute Delays")
scatter(Stations.longitude, Stations.latitude; markersize=msize, zcolor=rmsd,
    color=:balance, clims=color_limits, label=nothing, title="Relative Delays")
histogram(Data.observation, title = "Absolute Delays")
histogram(rdt, title = "Relative Delays")

# Write new input files
write_psi_s_events(out_event, Events)
write_psi_s_stations(out_station, Stations)
# If you want to write the relative delays to the data file then define,
# >> Data.observation .= rdt
write_psi_s_observations(out_data, Data, Events, Stations)