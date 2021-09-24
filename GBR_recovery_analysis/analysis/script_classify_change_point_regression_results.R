##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_classify_change_point_regression_results.R
# Summary: Classifies trajectories accourding to the slopes of the change-point regression
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021


# load filter results with per-capita data 
load(paste(PROC_DATA_DIR,"pcHC.",sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
# load regression summaries
load(paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))

# can arrange by F-statistic/p-Value or Adj-R^2
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% arrange(desc(ADJ_R2)) 
# 1st phase duration based on the change point
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% mutate(DELAY_DURATION = 0)
# 1st phase duration based on regression line intersections
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% mutate(DELAY_DURATION_LB = 0) 
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% mutate(DELAY_DURATION_UB = 0) 
# classification as per the change point test
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% mutate(CLASSIF = "") 
for (j in 1:length(filt.rec.traj.proc)) {
    
    k <- filt.rec.traj.proc.reg$INDEX[j]
    # get the trajectory for this regression result
    traj <- filt.rec.traj.proc[[k]]
   
    # look at slope combinations to classify
    if (filt.rec.traj.proc.reg$SLOPE1[j] > 0 && 
        filt.rec.traj.proc.reg$SLOPE2[j] <= 0 ) {
        filt.rec.traj.proc.reg$CLASSIF[j] <- 'D'
    } else if (filt.rec.traj.proc.reg$SLOPE1[j] <= 0 && 
               filt.rec.traj.proc.reg$SLOPE2[j] > 0) {
        filt.rec.traj.proc.reg$CLASSIF[j] <- '?'
    } else if (filt.rec.traj.proc.reg$SLOPE1[j] > 0 && 
               filt.rec.traj.proc.reg$SLOPE2[j] > 0 ) {
        filt.rec.traj.proc.reg$CLASSIF[j] <- 'D'
    } else {
        filt.rec.traj.proc.reg$CLASSIF[j] <- 'L/G'
    }

    # in the case of delay, estimate the 1st phase duration (direct approach)
    if (filt.rec.traj.proc.reg$CLASSIF[j] %in% c('D')) {
        # direct method (use data point identified as change-point)
        filt.rec.traj.proc.reg$DELAY_DURATION[j] <- traj$T[filt.rec.traj.proc.reg$CHANGE_POINT[j]]
        filt.rec.traj.proc.reg$DELAY_DURATION_LB[j] <- traj$T[filt.rec.traj.proc.reg$CHANGE_POINT_LOWER[j]]
        filt.rec.traj.proc.reg$DELAY_DURATION_UB[j] <- traj$T[filt.rec.traj.proc.reg$CHANGE_POINT_UPPER[j]]
        filt.rec.traj.proc.reg$DELAY_HC_COVER[j] <- filt.rec.traj.proc[[k]]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[j]]
    }
}

# save summaries with classifications appended
save(filt.rec.traj.proc.reg,file=paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
