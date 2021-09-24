##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_apply_regression_rank_fit.R
# Summary: script applies standard linear regression to eahc per-capita recovery,
# trajectory. Cases with good fit and negative slope indicate consitency with
# logistic-growth family of models. Those with poor fit or good fite with 
# positive slope indicate deviations from logistic growth functions.
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
load(paste(PROC_DATA_DIR,"pcHC.",sprintf(FILTER_REFMT_OUT_DATA_FMT,
                            max_init,min_final,min_obs)
           ,".RData",sep=""))

# initial analysis, just use linear 
linregs = list()
cp_regs = list()
cp_nlregs = list()
errs = c()
filt.rec.traj.proc.reg = data.frame(INDEX = 1:length(filt.rec.traj.proc),
                                    REEF = "",
                                    RP_ID = 0, 
                                    SITE_NO = 0,
                                    DISTURBANCE = "",
                                    DURATION = 0,
                                    DIFF = 0,
                                    INIT_Date = as.Date(0,origin="1970-01-01"),
                                    INIT_HC_COVER = 0,
                                    SLOPE = 0, 
                                    INTERCEPT = 0, 
                                    ADJ_R2 = 0, 
                                    FSTAT = 0,
                                    PVAL = 0,
                                    RHO = 0,
                                    PVAL_RHO = 0,
                                    CHANGE_POINT = 0,
                                    CHANGE_POINT_LOWER = 0,
                                    CHANGE_POINT_UPPER = 0,
                                    FSTAT_CP = 0,
                                    PVAL_CP = 0,
                                    SLOPE1 = 0,
                                    INTERCEPT1 = 0,
                                    SLOPE2 = 0,
                                    INTERCEPT2 = 0,
                                    ADJ_R2_CP = 0, 
                                    stringsAsFactors = FALSE)

for (j in 1:length(filt.rec.traj.proc)) {
    
    linregs <- lm(pcHC_lin_fd ~ HC,filt.rec.traj.proc[[j]])
    res <- summary(linregs)
    # correlation coefficient for simple linear regression
    spear_rho <- cor.test(filt.rec.traj.proc[[j]]$HC, 
                          filt.rec.traj.proc[[j]]$pcHC_lin_fd,
                          method="spearman")
    
    cp_regs <- sup_Like_Ratio(filt.rec.traj.proc[[j]])
    cp <- cp_regs[[1]]
    res1 <- summary(cp_regs[[2]])
    res2 <- summary(cp_regs[[3]])
    cp_res <- cp_regs[[4]]
    adj_R2_cp <- cp_regs[[6]] 
    cp_l <- cp_regs[[8]]
    cp_u <- cp_regs[[9]]
    
    # store model output and statistics
    filt.rec.traj.proc.reg$INDEX[j] <- j
    filt.rec.traj.proc.reg$REEF[j] <-  as.character(filt.rec.traj.proc[[j]]$REEF[1])
    filt.rec.traj.proc.reg$RP_ID[j] <-  filt.rec.traj.proc[[j]]$RP_ID[1]
    filt.rec.traj.proc.reg$SITE_NO[j] <-  filt.rec.traj.proc[[j]]$SITE_NO[1]
    filt.rec.traj.proc.reg$DISTURBANCE[j] <-  as.character(filt.rec.traj.proc[[j]]$DISTURBANCE[1])
    filt.rec.traj.proc.reg$DURATION[j] <-  filt.rec.traj.proc[[j]]$DURATION[1]
    filt.rec.traj.proc.reg$DIFF[j] <- filt.rec.traj.proc[[j]]$DIFF[1]
    filt.rec.traj.proc.reg$INIT_Date[j] <-  as.Date(filt.rec.traj.proc[[j]]$Date[1],origin="1970-01-01")
    filt.rec.traj.proc.reg$INIT_HC_COVER[j] <-  filt.rec.traj.proc[[j]]$HC[1]
    filt.rec.traj.proc.reg$SLOPE[j] <-  res$coefficients[2,1]
    filt.rec.traj.proc.reg$INTERCEPT[j] <-  res$coefficients[1,1]
    filt.rec.traj.proc.reg$ADJ_R2[j] <-  res$adj.r.squared
    filt.rec.traj.proc.reg$FSTAT[j] <-  res$fstatistic[1]
    filt.rec.traj.proc.reg$PVAL[j] <-  res$coefficients[2,4]
    filt.rec.traj.proc.reg$RHO[j] <- spear_rho$estimate
    filt.rec.traj.proc.reg$PVAL_RHO[j] <- spear_rho$p.value
    filt.rec.traj.proc.reg$CHANGE_POINT[j] <- cp
    filt.rec.traj.proc.reg$CHANGE_POINT_LOWER[j] <- cp_l
    filt.rec.traj.proc.reg$CHANGE_POINT_UPPER[j] <- cp_u
    filt.rec.traj.proc.reg$FSTAT_CP[j] <- as.numeric(cp_res[1])
    filt.rec.traj.proc.reg$PVAL_CP[j] <- as.numeric(cp_res[4])
    filt.rec.traj.proc.reg$SLOPE1[j] <-  res1$coefficients[2,1]
    filt.rec.traj.proc.reg$INTERCEPT1[j] <-  res1$coefficients[1,1]
    filt.rec.traj.proc.reg$SLOPE2[j] <-  res2$coefficients[2,1]
    filt.rec.traj.proc.reg$INTERCEPT2[j] <-  res2$coefficients[1,1]
    filt.rec.traj.proc.reg$ADJ_R2_CP[j] <-  adj_R2_cp
}

filt.rec.traj.proc.reg <- mutate(filt.rec.traj.proc.reg,REEF = as.factor(REEF),
                                    DISTURBANCE = as.factor(DISTURBANCE))

filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% arrange(FSTAT)

# save regression summaries
save(filt.rec.traj.proc.reg,file=paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))

