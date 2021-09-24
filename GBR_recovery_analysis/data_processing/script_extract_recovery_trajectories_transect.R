##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_extract_recovery_trajectories_transect.R
# Summary: script to extract out recovery periods from the Australian Institute
#  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marine 
#  Monitoring Program (MMP) data.
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
# 
# Last Modified: 17 September 2021
#

# import the samples and disturbance data
load(paste(DATA_DIR,DIST_DAT_FILE,sep=""))
load(paste(DATA_DIR,SAMPLE_DAT_FILE,sep=""))
load(paste(DATA_DIR,GRP_LEV_TRANS_DAT_FILE,sep=""))

# extract the recovery data based on LTMP disturbance records
recovery.trajectories <- extract_recovery_trajectories_transect(disturbance,
                                                                samples,
                                                                groups.transect)

# save resulting data table
save(recovery.trajectories, file=paste(PROC_DATA_DIR,REC_TRAJ_DAT_FILE,sep=""))
