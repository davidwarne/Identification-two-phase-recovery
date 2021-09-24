##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_filter_recovery_trajectories.R
# Summary: script to filter recovery periods from the Australian Institute
#  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marin 
#  Monitoring Program (MMP) data.
#
# Details: Applies filter conditions to the recovery trajectory data set that
# is derived using the script_extract_recovery_trajectories_transect.R 
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021

# User functions for filtering
#-------------------------------------------------------------------------------

##
# Summary: append derived data function
# Detail: This function adds the number of observations and the number of days 
# since the last disturbance for each visit in the recovery trajectory
#
# param: traj the data for a single recovery trajectory
# returns: a modified data.frame with appended variables
#
# note: this user function can be customised the include other derived variables that could
# be useful for filtering.
append_num_obs <- function(traj) {
    num_obs <- length(unique(traj$VISIT_NO))
    visits <- traj %>% 
                select(VISIT_NO,Date) %>% 
                unique() %>% 
                arrange(VISIT_NO)
    diff_in_days = as.numeric(visits$Date[num_obs]-visits$Date[1],units = "days")
    return(mutate(traj,NUM_OBS = num_obs,DURATION = diff_in_days))
}

##
# Summary: Trajectory filter function
# Detail: User defined filter function implementing the condition
#   HC_0 < max_init and n >= min_obs
#
# param: traj the data for a single recovery trajectory with appended data
#
# returns: TRUE if traj should be included
#
# note: this user function can be customised to implement any filtering rule based on data
trajectory_filter_condition_no_end <- function(traj) {
    # Short-cut if not enough observations
    if (traj$NUM_OBS[1] < min_obs) {
        return(FALSE)
    } else {
        # condition total coral cover % = total HC % + total SC % < thresh_cover %
        first_visit <- min(traj$VISIT_NO)
        last_visit <- max(traj$VISIT_NO)
        init_state <- traj %>% 
                      filter(VISIT_NO == first_visit,GROUP_CODE == 'HC')
        final_state <- traj %>% 
                       filter(VISIT_NO == last_visit,GROUP_CODE == 'HC')
        # site level cover is derived from the mean of transect level cover
        tf <- (mean(init_state$COVER) < max_init 
               && mean(final_state$COVER) > mean(init_state$COVER))
        # check no complete zeros
        for (visit in traj$VISIT_NO) {
            state <- traj %>% 
                      filter(VISIT_NO == visit,GROUP_CODE == 'HC')
            tf <- tf && (mean(state$COVER) > 0.0)
        }
        return(tf)
    }
}


# process filter
#-------------------------------------------------------------------------------

# load recovery trajectory data, cover data, and spatial data
load(paste(PROC_DATA_DIR,REC_TRAJ_DAT_FILE,sep=""))
load(paste(DATA_DIR,GRP_LEV_TRANS_DAT_FILE,sep=""))

# perform filtering 
filt.rec.traj <- filter_recovery_trajectories_transect(recovery.trajectories,
                                     groups.transect,
                                     derive_func = append_num_obs,
                                     filter_func = trajectory_filter_condition_no_end)
# save the trajectory data
save(filt.rec.traj,file=paste(PROC_DATA_DIR,
                        sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
                        ,".RData",sep=""))
