##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_process_filter_results.R
# Summary: script to transform recovery trajectories into a more suitable form for
# individual analysis.
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021

# load filter results 
load(paste(PROC_DATA_DIR,sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))

filt.rec.traj.proc <- list()

j <- 1
# reformat all trajectories and store in a list
for (rp_id in unique(filt.rec.traj$RP_ID)) {
   cover.dat <- filt.rec.traj %>% filter(RP_ID == rp_id)
   filt.rec.traj.proc[[j]] <- reformat_recovery_trajectories(cover.dat)
   j <- j + 1
}

# save the cleaned-up trajectory data
save(filt.rec.traj.proc,file=paste(PROC_DATA_DIR,
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
