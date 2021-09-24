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

# load cover data at benthos code level 
load(paste(DATA_DIR,BNTHS_LEV_TRANS_DAT_FILE,sep=""))

# load regression summaries
load(paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))

# compute Benthic code version of trajectories
filt.rec.traj.benthos <- right_join(benthos.transect,
                                    select(filt.rec.traj,-c("GROUP_CODE","COVER")))

filt.rec.traj.proc.reg.delay <- filt.rec.traj.proc.reg %>% 
                                    filter(CLASSIF == 'D')

FinalCover <- data.frame(COVER = -1,DELAY = 'Two-phase',YEARS = -1,stringsAsFactors = FALSE)

# compute final cover under different scenarios
for (T_int in c(1,5,7.5,10)) {
    for (j in 1:nrow(filt.rec.traj.proc.reg.delay)) {
        group.code.traj <- reformat_recovery_trajectories(filt.rec.traj %>% 
                            filter(RP_ID == filt.rec.traj.proc.reg.delay$RP_ID[j]),
                            codes = "GROUP_CODE",
                            annotate = "_G")
        benthos.code.traj <- reformat_recovery_trajectories(filt.rec.traj.benthos %>% 
                            filter(RP_ID == filt.rec.traj.proc.reg.delay$RP_ID[j]),
                            codes = "BENTHOS_CODE",
                            annotate = "_B")
        group.benthos.code.traj <- right_join(group.code.traj,benthos.code.traj)

        res <- quantify_delay_impact(group.benthos.code.traj,
                                         filt.rec.traj.proc.reg.delay$CHANGE_POINT[j],
                                         T_int,30) 
        FinalCover <- rbind(FinalCover, list(res[[1]][1],'Single-phase',T_int))
        FinalCover <- rbind(FinalCover, list(res[[3]][1],'Two-phase',T_int))
    }
}


Dist_freq <- 5
YearsAboveCover <- data.frame(COVER = -1,DELAY = 'Single-phase', YEARS_FREQ1 = -1,
                              YEARS_FREQ5 = -1,YEARS_FREQ10 = -1,
                              stringsAsFactors = FALSE)
# compute years above a given threshold cover
for (thresh_cover in c(10,15,20)) {
    for (j in 1:nrow(filt.rec.traj.proc.reg.delay)) {
        group.code.traj <- reformat_recovery_trajectories(filt.rec.traj %>% 
                            filter(RP_ID == filt.rec.traj.proc.reg.delay$RP_ID[j]),
                            codes = "GROUP_CODE",
                            annotate = "_G")
        benthos.code.traj <- reformat_recovery_trajectories(filt.rec.traj.benthos %>% 
                            filter(RP_ID == filt.rec.traj.proc.reg.delay$RP_ID[j]),
                            codes = "BENTHOS_CODE",
                            annotate = "_B")
        group.benthos.code.traj <- right_join(group.code.traj,benthos.code.traj)

        res <- quantify_delay_impact(group.benthos.code.traj,
                                         filt.rec.traj.proc.reg.delay$CHANGE_POINT[j],
                                         10,thresh_cover) 
        YearsAboveCover <- rbind(YearsAboveCover,
                                 list(thresh_cover,'Single-phase', 
                                      10*max(1 - res[[2]][1],0), 
                                      2*max(5 - res[[2]][1],0),
                                      max(10-res[[2]][1],0)))
        YearsAboveCover <- rbind(YearsAboveCover,
                                 list(thresh_cover,'Two-phase',
                                      10*max(1 - res[[4]][1],0),
                                      2*max(5 - res[[4]][1],0),
                                      max(10-res[[4]][1],0)))
    }
}

# results tables
FinalCover <- FinalCover %>% filter(COVER >= 0)
FinalCover$DELAY <- as.factor(FinalCover$DELAY)
FinalCover$YEARS <- as.factor(FinalCover$YEARS)

YearsAboveCover <- YearsAboveCover %>% filter(COVER >= 0)
YearsAboveCover$DELAY <- as.factor(YearsAboveCover$DELAY)
YearsAboveCover$COVER <- as.factor(YearsAboveCover$COVER)

sumFinalCover <- FinalCover %>% group_by(DELAY,YEARS) %>%
    summarise(mu = mean(COVER),
              sd = sd(COVER), 
              se = sd(COVER)/sqrt(length(COVER)),
              CIl = mean(COVER) - 1.96*sd(COVER)/sqrt(length(COVER)),
              CIU = mean(COVER) + 1.96*sd(COVER)/sqrt(length(COVER)),
              q1 = quantile(COVER,0.25),
              q3 = quantile(COVER,0.75),
              q2 = median(COVER))
sumFinalCover$YEARS <- as.numeric(as.character(sumFinalCover$YEARS)) 

sumYearsAbove <- YearsAboveCover %>% group_by(DELAY,COVER) %>%
    summarise(mu10 = mean(YEARS_FREQ10),
              sd10 = sd(YEARS_FREQ10), 
              se10 = sd(YEARS_FREQ10)/sqrt(length(YEARS_FREQ10)),
              CIl10 = mean(YEARS_FREQ10) - 1.96*sd(YEARS_FREQ10)/sqrt(length(YEARS_FREQ10)),
              CIU10 = mean(YEARS_FREQ10) + 1.96*sd(YEARS_FREQ10)/sqrt(length(YEARS_FREQ10)),
              q110 = quantile(YEARS_FREQ10,0.25),
              q310 = quantile(YEARS_FREQ10,0.75),
              q210 = median(YEARS_FREQ10),
              mu5 = mean(YEARS_FREQ5),
              se5 = sd(YEARS_FREQ5), 
              se5 = sd(YEARS_FREQ5)/sqrt(length(YEARS_FREQ5)),
              CIl5 = mean(YEARS_FREQ5) - 1.96*sd(YEARS_FREQ5)/sqrt(length(YEARS_FREQ5)),
              CIU5 = mean(YEARS_FREQ5) + 1.96*sd(YEARS_FREQ5)/sqrt(length(YEARS_FREQ5)),
              q15 = quantile(YEARS_FREQ5,0.25),
              q35 = quantile(YEARS_FREQ5,0.75),
              q25 = median(YEARS_FREQ5),
              mu1 = mean(YEARS_FREQ1),
              se1 = sd(YEARS_FREQ1), 
              se1 = sd(YEARS_FREQ1)/sqrt(length(YEARS_FREQ1)),
              CIl1 = mean(YEARS_FREQ1) - 1.96*sd(YEARS_FREQ1)/sqrt(length(YEARS_FREQ1)),
              CIU1 = mean(YEARS_FREQ1) + 1.96*sd(YEARS_FREQ1)/sqrt(length(YEARS_FREQ1)),
              q11 = quantile(YEARS_FREQ1,0.25),
              q31 = quantile(YEARS_FREQ1,0.75),
              q21 = median(YEARS_FREQ1))

sumYearsAbove$COVER <- as.numeric(as.character(sumYearsAbove$COVER))

