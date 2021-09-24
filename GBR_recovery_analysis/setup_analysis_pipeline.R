##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: setup_analysis_pipeline.R
# Summary: Load packages and utility functions.
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

library(tidyverse)
library(ggmap)
library(gridExtra)
library(functional)
library(tictoc)
library(gap)
library(stats4)
library(scales)
library(ggExtra)
library(cowplot)

source("../LTMPTools/LTMPDataTools.R")
source("../LTMPTools/LTMPModellingTools.R")

# input dir and data
DATA_DIR                  <- "../data/primary/"
DIST_DAT_FILE             <- "disturbance.RData"
SAMPLE_DAT_FILE           <- "samples.RData"
VISIT_DAT_FILE            <- "ker.code.visit.RData"
SPAT_DAT_FILE             <- "spatial.dat.RData"
GRP_LEV_TRANS_DAT_FILE    <- "groups.transect.RData"
BNTHS_LEV_TRANS_DAT_FILE  <- "benthos.transect.RData"


# output dir and data
PROC_DATA_DIR           <- "../data/processed/"
REC_TRAJ_DAT_FILE       <- "time.series.site.transect.RData"
TIME_SERIES_DAT_FILE    <- "recovery.trajectories.site.transect.RData"

FILTER_OUT_DATA_FMT         <- "rec.traj.trans.init%f.final%f.obs%d"
FILTER_REFMT_OUT_DATA_FMT   <- "rec.traj.trans.proc.init%f.final%f.obs%d"

# filter criteria (initial cover <= 10% and number of observations >= 5) 
max_init    <- 10.0
min_final   <- 0.0
min_obs     <- 5.0