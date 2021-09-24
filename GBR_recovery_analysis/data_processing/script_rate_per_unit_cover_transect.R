##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_rate_per_unit_cover_transect.R
# Summary: script to convert recovery HC % cover data into per-capita (in % cover 
# sense) using interpolation and numerical differentiation.
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021

# load the cleaned-up trajectory data
load(paste(PROC_DATA_DIR,sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))

# plot ALL recoveries just to get a feel for the variation
summary.stats.rec.traj.proc <- data.frame('HC_START'=double(),'HC_END'=double(),'DUR'=double(),'NUM_OBS'=double())

for (j in 1:length(filt.rec.traj.proc)) {

    # collect info for summary stats
    n <- length(filt.rec.traj.proc[[j]]$RP_ID)
    endm1_hc <- filt.rec.traj.proc[[j]]$HC[n-1]
    start_hc <- filt.rec.traj.proc[[j]]$HC[1]
    end_hc <- filt.rec.traj.proc[[j]]$HC[n]
    dur <- filt.rec.traj.proc[[j]]$T[n]
    
    summary.stats.rec.traj.proc[j,] <- c(start_hc,end_hc,dur/365,n) 
    

    # Computing derivatives and standard errors
    # linear interp 
    lf <- approxfun(x = filt.rec.traj.proc[[j]]$T,
                    y = filt.rec.traj.proc[[j]]$HC,
                    method = "linear")
    lf_sig <- approxfun(x = filt.rec.traj.proc[[j]]$T,
                    y = filt.rec.traj.proc[[j]]$HC_se,
                    method = "linear")
    filt.rec.traj.proc[[j]] <- mutate(filt.rec.traj.proc[[j]], 
                                      pcHC_lin_fd = 0, 
                                      pcHC_lin_qlb = 0, 
                                      pcHC_lin_qub = 0)
    h <- 10
    reef <- filt.rec.traj.proc[[j]]$REEF[1]
    
    # finite differences (2nd order internal nodes, 1st order boundaries)
    pc_lin_uq <- finite_diff_pc_uq(lf,lf_sig,filt.rec.traj.proc[[j]]$T,h)
    filt.rec.traj.proc[[j]] <- mutate(filt.rec.traj.proc[[j]],pcHC_lin_fd = pc_lin_uq[[1]])
    filt.rec.traj.proc[[j]] <- mutate(filt.rec.traj.proc[[j]],pcHC_lin_qlb = pc_lin_uq[[2]])
    filt.rec.traj.proc[[j]] <- mutate(filt.rec.traj.proc[[j]],pcHC_lin_qub = pc_lin_uq[[4]])
}
#print summary stats for dataset
print(summary(summary.stats.rec.traj.proc) )

# save trajectory data with rate per unit cover data
save(filt.rec.traj.proc,file=paste(PROC_DATA_DIR,"pcHC.",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
