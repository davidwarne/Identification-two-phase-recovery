##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_sensistivity_analysis.R
# Summary: checks the detection rate of delays using a simulation study
#
# Author: David J. Warne (1,2,3)
#         1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#         2. Centre for Data Science, Queensland University of Technology
#         3. ARC Centre of Excellence for Mathematical and Statistical Frontiers
# 
# Email: david.warne@qut.edu.au
#

library(tidyverse)
library(ggmap)
library(gridExtra)
library(functional)
library(tictoc)
library(gap)

source("../LTMPTools/LTMPDataTools.R")
source("../LTMPTools/LTMPModellingTools.R")

# to use the scripts
PROC_DATA_DIR <- "./"
FILTER_REFMT_OUT_DATA_FMT <- "temp_file"
max_init <- 5.0
min_final <- 0
min_obs <- 5

# model solutions for simulation study
#
# Logistic
logistic_sol <- function(t,alpha,K,c0) {
    return((K*c0)/((K-c0)*exp(-alpha*t) + c0))
}
# Gompertz
gompertz_sol <- function(t,alpha,K,c0) {
    return(K*exp(log(c0/K)*exp(-alpha*t)))
}
delay_sol <- function(t,alpha,c0) {
    return(c0*exp(alpha*t))
}
# delay logistic
delay_logistic_sol <- function(t,T,alpha1,alpha2,K,c0) {
    if (t <= T) {
        return(delay_sol(t,alpha1,c0))
    } else {
        return(logistic_sol(t-T,alpha2,K,delay_sol(T,alpha1,c0)))
    }
}
# delay Gompertz
delay_gompertz_sol <- function(t,T,alpha1,alpha2,K,c0) {
    if (t <= T) {
        return(delay_sol(t,alpha1,c0))
    } else {
        return(gompertz_sol(t-T,alpha2,K,delay_sol(T,alpha1,c0)))
    }
}

# hold these as constant
K <- 80
c0 <- 2.5
n_years <- 12
visit_std <- 1/36 # i.e., 2std = 1 month

N_samples <- 400

# run simulation study over these configutations
# alpha from Thompson and Dolman (2010 Arc = 0.933, OthrHC = 0.349, SC = 0.282)
# obs_error based on LTMP and MMP data standard errors for observations < 20%
results_table = expand.grid(alpha = c(0.349,0.641,0.933),
                     alpha_scale = c(0.25,0.5,0.75),
                     DTime = c(3.5,5.5,7.5),
                     sampling = c('annual','biannual'),
                     model = c('logistic','Gompertz'),
                     obs_error = c(0.5,1.5,2.5),
                     stringsAsFactors = TRUE) %>% 
                mutate(TPR = -1,FPR = -1,TP=-1,TN=-1,FP=-1,FN=-1)
# range of rate parameters
for (i in 1:nrow(results_table)) {
   
    # extract the configuration
    alpha <- results_table$alpha[i]
    alpha_d <- alpha*(results_table$alpha_scale[i])
    Ti <- results_table$DTime[i]
    se_mag <- results_table$obs_error[i]
    model <- results_table$model[i]
    sampling <- results_table$sampling[i]
    t <- c()
    if (sampling =='annual') {
        t <- seq(0,n_years,by=1)
    } else {
        t <- seq(0,n_years,by=2)
    }
    t <- t + rnorm(length(t),0,visit_std)
    t[t<0] <- 0 
    t[1] <- 0
    # do sampling to construct simulated data set
    filt.rec.traj.proc <- list()
    true_class <-  data.frame(RP_ID = seq(1,2*N_samples,by=1),
                              TRUE_CLASS = '',stringsAsFactors = FALSE)
    k <- 1
    print(sprintf("Sampling simulated data for config %d %f %f %f %f %s %s...",i, 
                  alpha,alpha_d,Ti,se_mag,model,sampling))
    tic()
    for (j in 1:N_samples) {
        hc <- seq(0,0,length=length(t))
        for (class in c('D','L/G')) {
            if (model == 'logistic' && class == 'L/G') { 
                for (ii in 1:length(t)) { 
                    hc[ii] <- logistic_sol(t[ii],alpha,K,c0) + 
                                rnorm(1,0,se_mag)
                } 
            } else if (model == 'logistic' && class == 'D') {
                for (ii in 1:length(t)) { 
                    hc[ii] <- delay_logistic_sol(t[ii],Ti,alpha_d,alpha,K,c0) + 
                                rnorm(1,0,se_mag)
                }
            } else if (model == 'Gompertz' && class == 'L/G') {
                for (ii in 1:length(t)) {
                    hc[ii] <- gompertz_sol(t[ii],alpha,K,c0) + 
                                rnorm(1,0,se_mag)
                }
            } else if (model == 'Gompertz' && class == 'D') {
                for (ii in 1:length(t)) {
                    hc[ii] <- delay_gompertz_sol(t[ii],Ti,alpha_d,alpha,K,c0) + 
                                rnorm(1,0,se_mag)
                }
            }
            true_class$RP_ID[k] <- k
            true_class$TRUE_CLASS[k] <- class
            filt.rec.traj.proc[[k]] <- data.frame(P_CODE = "none",
                                          REEF = sprintf("sim %s",class),
                                          SITE_NO = 1, DEPTH = 6,RP_ID = k,
                                          Date = as.Date("2005-12-05") + t*365.0,
                                          VISIT_NO = 1:length(t),
                                          DURATION = max(t)*365.0,
                                          DISTURBANCE = "n",T = t*365.0,
                                          HC = hc,HC_sd = 0.0,HC_se = se_mag)
            k <- k + 1
        }
    }
    toc()
    print("done.")

    i_old <- i
    
    # process analysis pipeline for simulated data set
    tic()
    save(filt.rec.traj.proc,file=paste(PROC_DATA_DIR,
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
    source('./data_processing/script_rate_per_unit_cover_transect.R')
    source('./model_fitting/script_apply_cp_regression_rank_fit.R')
    source('./analysis/script_classify_change_point_regression_results.R')
    load(paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
             sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
             ,".RData",sep=""))
    toc()
    # estimate true positive and false positive rates
    out_comes <- right_join(true_class,
                     select(filt.rec.traj.proc.reg,
                            .data$RP_ID,.data$CLASSIF),by = c('RP_ID'))
    TPos <- out_comes %>% filter(TRUE_CLASS == "D" & CLASSIF == "D")
    TNeg <- out_comes %>% filter(TRUE_CLASS == "L/G" & 
                                 (CLASSIF == "L/G" | CLASSIF == "?"))
    FPos <- out_comes %>% filter(TRUE_CLASS == "L/G" & CLASSIF == "D")
    FNeg <- out_comes %>% filter(TRUE_CLASS == "D" & 
                                 (CLASSIF == "L/G" | CLASSIF == "?"))
    i <- i_old
    # power or sensitivity
    results_table$TPR[i] <- nrow(TPos)/(nrow(TPos) + nrow(FNeg)) 
    results_table$FPR[i] <- nrow(FPos)/(nrow(FPos) + nrow(TNeg)) 
    # 1 - specifity
    results_table$TP[i] <- nrow(TPos)
    results_table$TN[i] <- nrow(TNeg)
    results_table$FP[i] <- nrow(FPos)
    results_table$FN[i] <- nrow(FNeg)
}

