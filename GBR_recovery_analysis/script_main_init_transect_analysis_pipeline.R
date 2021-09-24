##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: script_main_analysis_pipeline.R
# 
# Summary: 
# Complete analysis pipeline for identification of two-phase recovery patterns 
# in the Great Barrier Reef.
#
# Description: Performs all of the analysis produces figures from the journal article:
#
# DJ Warne, KA Crossman, W Jin, K Mengersen, K Osborne, MJ Simpson, AA Thompson, P Wu, J-C Ortiz (2021).
#     Identification of two-phase recovery for interpretation of coral reef monitoring data. Journal of
#     Applied Ecology (in press).
#   
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

# Set up analysis pipeline
source('./setup_analysis_pipeline.R')

# Extract all recovery trajectories 
tic()
print(" Extracting recovery trajectories...")
source('./data_processing/script_extract_recovery_trajectories_transect.R')
print(" done!")
toc()

# Filter on selected criteria
tic()
print(" Applying filter ...")
source('./data_processing/script_filter_recovery_trajectories_transect.R')
print(" done!")
toc()

# Reformat filter results
tic()
print(" Reformatting filter results ...")
source('./data_processing/script_process_filter_results_transect.R')
print(" done!")
toc()

# Compute growth rate per unit cover
tic()
print(" Computing growth rate per unit cover ...")
source('./data_processing/script_rate_per_unit_cover_transect.R')
print(" done!")
toc()

# Apply change-point model and rank results by fit
tic()
print(" applying regression and ranking ...")
source('./model_fitting/script_apply_cp_regression_rank_fit.R')
print(" done!")
toc()

# Assign delay classification based on trends before/after change points
tic()
print(" classification  ...")
source('./analysis/script_classify_change_point_regression_results.R')
print(" done!")
toc()

# optional power analysis (warning takes > 24 hours)
#tic()
#print("Simulations for power analysis ...")
#source("./analysis/script_power_analysis.R")
#toc()


# estimate effect of delay on % cover and years of reef services
tic()
print("estimating delay effects ...")
source("./analysis/script_estimate_delay_effect.R")
toc()

print('plotting figures')
source("./plotting/Figures.R")

