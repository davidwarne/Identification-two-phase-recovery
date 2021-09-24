##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# File: LTMPDataTools.R
# Summary: Data processing functions to deal with the LTMP (Long Term Monitoring 
# Program) and MMP (marine Monitoring Program) Data from AIMS (Australian Institute
# of Marine Sciences).
#
# Functions are specifically focused on the filtering of AIMS coral cover 
# data at a recovery trajectory level (sequences between two disturbance periods
# for a given site). 
#
# Current features:
#   - data processing for splitting data into recovery trajectories from 
#     site-level using LTMP disturbance classifications.
#   - Filtering of recovery trajectories based on upper-bound for initial 
#     % hard coral cover, lower-bound for final % hard coral cover, and lower 
#     on number of observations.
#
#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------


## ****
# Summary: extraction of recovery trajectories using transect level data.
# Details: Uses differences in successive transect observations at each site over
# time to extract to detect start and end of site level recovery trajectories. 
# Assigns a recovery id to each trajectory at a site level.
#
# param: disturbance        data frame of disturbance records
# param: samples            data frame of visit time/place records
# param: transect.cover     data frame of % cover at transect level, can be at
#                           major group level or benthic code (ker code) level
#                           by major group ("GROUP_CODE") or ker code ("BENTHOS_CODE")
# param: code               code base trajectory tests on
# param: alpha              significance level for paired t-test
#
extract_recovery_trajectories_transect <- function(disturbance,samples,
                                                   transect.cover,
                                                   code = "HC",alpha = 0.95) {

    # to avoid factor warnings
    combined <- sort(union(union(levels(disturbance$REEF),levels(samples$REEF)),
                           levels(transect.cover$REEF)))
    disturbance <- disturbance %>% 
        mutate(REEF=factor(REEF,levels=combined))
    samples <- samples %>% 
        mutate(REEF=factor(REEF,levels=combined))
    transect.cover <- transect.cover %>% 
        mutate(REEF=factor(REEF,levels=combined))
   
    # extend disturbance factor levels to distinguish between recorded disturbances
    # and our own entries
    combined <- sort(union(levels(disturbance$DISTURBANCE), c('U','T','S')))
    disturbance <- disturbance %>% 
        mutate(DISTURBANCE=factor(DISTURBANCE,levels=combined))


    # join disturbances and samples to transects
    samples.disturbance.transect.cover <- right_join(
                      right_join(disturbance,samples,by = c("REEF","VISIT_NO")),
                      select(transect.cover, P_CODE,REEF,DEPTH,VISIT_NO,SITE_NO),
                      by = c("REEF","DEPTH","VISIT_NO")) %>%
                      mutate(RP_ID = 0,DIFF = 0, Pval = 0)
    recovery.trajectories.site <- data.frame()
    
    # for every REEF/SITE/DEPTH triple by visit
    unique_sites <- unique(select(samples.disturbance.transect.cover,
                                   REEF,DEPTH,SITE_NO))
    rp_id <- 1
    for (i in 1:nrow(unique_sites)){
        visit.sequence.site <- samples.disturbance.transect.cover %>% 
                                  filter(REEF == unique_sites$REEF[i],
                                         DEPTH == unique_sites$DEPTH[i],
                                         SITE_NO == unique_sites$SITE_NO[i]) %>%
                                  unique() %>%
                                  arrange(VISIT_NO)
        # step through sequence applying one-sided paired t-test to subsequent pairs
        num_dist <- 1
        visit.seq <- visit.sequence.site$VISIT_NO
        for (j in 1:length(visit.seq)) {
            # first visit is always a new recovery trajectory
            if (j == 1) { 
                rp_id <- rp_id + 1
                num_dist <- 1
                visit.sequence.site$DISTURBANCE[j] <- 'S' 
            } else { 
                # Get transect cover for level code provided for time t and t-1
                transect.cover.pair <- transect.cover %>% 
                          filter(REEF == visit.sequence.site$REEF[1],
                                 DEPTH == visit.sequence.site$DEPTH[1],
                                 SITE_NO == visit.sequence.site$SITE_NO[1],
                                 VISIT_NO == visit.seq[j] |
                                 VISIT_NO == visit.seq[j-1],
                                 GROUP_CODE == code) %>%
                          select(TRANSECT_NO,VISIT_NO,COVER)
                transect.cover.t1 <- transect.cover.pair %>%
                          filter(VISIT_NO == visit.seq[j-1])
                transect.cover.t2 <- transect.cover.pair %>%
                          filter(VISIT_NO == visit.seq[j])

                #compute drop/growth in cover of more than (drop > 0 and growth < 0)
                diff <- mean(transect.cover.t1$COVER) - mean(transect.cover.t2$COVER) 
                visit.sequence.site$DIFF[j] <- diff # store difference

                if (visit.sequence.site$DISTURBANCE[j] %in% 
                       c('d','s','b','c','m','u','f')) { # recorded disturbance
                    rp_id <- rp_id + 1
                    num_dist <- 1
                } else if (diff >= 5.0) { # unrecorded absolute drop >= 5%
                    rp_id <- rp_id + 1
                    num_dist <- num_dist + 1
                    # use U for unknown disturbance that was not recorded 
                    visit.sequence.site$DISTURBANCE[j] <- 'U' 
                } else { # otherwise perform paired t-test to check for small but 
                         # statistically significant drops

                    # catch cases of missing transect observation, ensure we keep
                    # only the pairs that match to ensure the paired t-test is
                    # validly applied
                    if (length(transect.cover.t1$VISIT_NO) < 
                        length(transect.cover.t2$VISIT_NO)) {
                        transect.cover.t2 <- transect.cover.t2 %>%
                               filter(TRANSECT_NO %in% transect.cover.t1$TRANSECT_NO)
                    } else if (length(transect.cover.t1$VISIT_NO) > 
                               length(transect.cover.t2$VISIT_NO)) {
                        transect.cover.t1 <- transect.cover.t1 %>%
                               filter(TRANSECT_NO %in% transect.cover.t2$TRANSECT_NO)
                    }
                    transect.cover.t1 <- transect.cover.t1 %>% arrange(TRANSECT_NO)
                    transect.cover.t2 <- transect.cover.t2 %>% arrange(TRANSECT_NO)
                
                    # perform paired t-test if there are enough transects
                    if (length(transect.cover.t1$TRANSECT_NO) > 1 &&
                        length(transect.cover.t2$TRANSECT_NO) > 1) {
                        # run paired t-test
                        res <- t.test(x = transect.cover.t1$COVER, y = transect.cover.t2$COVER,
                                      alternative = "greater", paired = TRUE, conf.level = alpha)
                        # any statistically significant pair update the RP_ID
                        if (res$p.value < 1.0 - alpha) {
                            rp_id <- rp_id + 1
                            num_dist <- num_dist + 1
                            # use T for small unknown disturbance as per t-test 
                            visit.sequence.site$DISTURBANCE[j] <- 'T'
                        }
                        # store p-value and diff estimate regardless of test result 
                        visit.sequence.site$Pval[j] <- res$p.value
                        visit.sequence.site$DIFF[j] <- res$estimate
                    }
                } 
            }
            visit.sequence.site$RP_ID[j] <- rp_id
            visit.sequence.site$NUM_DIST[j] <- num_dist
        }    
        recovery.trajectories.site <- bind_rows(recovery.trajectories.site,
                                          visit.sequence.site)
    }
    return(recovery.trajectories.site)
}


## ***
# Summary: script to filter recovery periods from the Australian Institute
#  for Marine Science (AIMS) Long Term Monitoring Program (LTMP) data and Marin 
#  Monitoring Program (MMP) data.
#
# Detail: applies filter conditions to the recovery trajectory data set that
#  is derived using the extract_recoveries.R script. filter conditions are specified
#  by a generic user function.
#
# param: recovery.trajectory data frame as produced by 
#                            extract_recovery_trajectories()
# param: transect.cover      transect level % cover data grouped by group codes
#                             or benthic codes
# param: derive_func         user prescribed function to append derived quantities
#                            to trajectories data. This function receives all cover 
#                            cover data for a single recovery trajectory with RP_ID 
#                            already appended
# param: filter_func         user prescribed function to allow a trajectory through
#                            trhe filter. This function receives the output from
#                            derive_func() as input
#
# return: a data frame of trajectories with coral cover data that satisfied the 
#         filter conditions
#
filter_recovery_trajectories_transect <- function(recovery.trajectories,
                                    transect.cover,
                                    derive_func = function(x){return(x)},
                                    filter_func = function(x){return(FALSE)}) {

    # create new data set containing only recovery trajectories with the 
    # required conditions
    filt.rec.traj <- data.frame() 
    # to avoid factor warnings/errors
    combined <- sort(union(levels(recovery.trajectories$REEF),
                           levels(transect.cover$REEF)))
    recovery.trajectories <- recovery.trajectories %>% 
                           mutate(REEF=factor(REEF,levels=combined))
    transect.cover <- transect.cover %>% 
                      mutate(REEF=factor(REEF,levels=combined))

    for (rp_id in unique(recovery.trajectories$RP_ID)) {
        # extract the trajectory
        visit.sequence <- recovery.trajectories %>%
                          filter(RP_ID == rp_id) %>%
                          arrange(VISIT_NO)
        
        # And coral cover data associated with this trajectory
        cover.dat <- transect.cover %>%
                     filter(REEF == visit.sequence$REEF[1],
                            DEPTH == visit.sequence$DEPTH[1],
                            SITE_NO == visit.sequence$SITE_NO[1],
                            VISIT_NO %in% visit.sequence$VISIT_NO) %>%
                     mutate(RP_ID = rp_id)
        
        ## join to get the dates for each visit to be available for the user
        cover.dat <- right_join(select(visit.sequence,VISIT_NO,Date,DIFF,Pval),
                                         cover.dat, by = "VISIT_NO")

        cover.dat <- derive_func(cover.dat)
            
        if (filter_func(cover.dat) == TRUE) {
             filt.rec.traj <- bind_rows(filt.rec.traj,cover.dat)
        }

    }
    filt.rec.traj <- right_join(select(recovery.trajectories,
                                       RP_ID,VISIT_NO,DISTURBANCE),
                                filt.rec.traj, by = c("RP_ID","VISIT_NO") )
    return(filt.rec.traj)
}

## ***
# Summary: reformat a single recovery trajectory as produced by the function
# filter_recover_trajectories for easier analysis.
# Details: the multiple coral cover records are aggregated by group code and 
# included as columns, so that each visit is a single row. Also the variable T is
# included as the number of days since the start of recovery.
#
# param: traj a single recovery trajectory (RP_ID/SITE_NO is unique)
# param: groups the Group codes to include in the final output.
#
# return: a data table such that result[i,] is the data for the ith visit in this
# recovery sequence.
reformat_recovery_trajectories <- function(traj, codes = 'GROUP_CODE', annotate = "") {

    # extract the constant data and unique visits (date is unique for each visit
    # so the number of rows is the same as the number of visits
    output.seq <- traj %>%
                  select(P_CODE,REEF,SITE_NO,DEPTH,RP_ID,Date,VISIT_NO,
                         DURATION,DISTURBANCE,DIFF,Pval) %>%
                  unique() %>%
                  arrange(VISIT_NO) %>%
                  mutate('T' =  as.double(Date - Date[1]))
    groups <- levels(traj[[codes]])
    # append a column for each group and populate
    for (g in groups){
        res <- traj %>% 
                filter(!!rlang::sym(codes) == g) %>%
                group_by(VISIT_NO) %>%
                summarise(mu = mean(COVER),
                          sd = sd(COVER), 
                          se = sd(COVER)/sqrt(length(COVER))) %>%
                arrange(VISIT_NO)
        output.seq <- mutate(output.seq,!!paste(g,annotate,sep="") := res$mu, 
                             !!paste(g,"_sd",annotate,sep="") := res$sd,
                             !!paste(g,"_se",annotate,sep="") := res$se)
    }
    return(output.seq)
}

