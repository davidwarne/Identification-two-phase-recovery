##
# Identification-two-phase-recovery  Copyright (C) 2021  David J. Warne
#    This program comes with ABSOLUTELY NO WARRANTY.
#    This is free software, and you are welcome to redistribute it
#    under certain conditions; See the LICENSE file for details.
#
# Summary:produces figures from the journal article (and Supplementary Material):
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

# basic libraries needed
library("ggplot2")
library("ggspatial")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

# load filter results with rate per unit cover data 
load(paste(PROC_DATA_DIR,"pcHC.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
# load regression summaries
load(paste(PROC_DATA_DIR,"cp.reg.sum.pcHC.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))


# load spatial data 
load(paste(DATA_DIR,SPAT_DAT_FILE,sep=""))
spatial.dat <- spatial.dat %>% filter(!is.na(LATITUDE) & !is.na(LONGITUDE))

# test removal of low adj-R^2
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% filter(ADJ_R2_CP >= 0.2)

# combine regression results with spatial data
filt.rec.traj.proc.reg.spat <- right_join(spatial.dat,
                                          filt.rec.traj.proc.reg)
filt.rec.traj.proc.reg.spat$CLASSIF <- as.factor(filt.rec.traj.proc.reg.spat$CLASSIF) 

# dark-on-light theme appropriate for maps
theme_set(theme_bw())

# load GBR boundary provided by the Great Barrier Reef Marine Park Authority (GBRMPA)
gbr = st_read(paste(DATA_DIR,
              "/gbr-spatial-data/zipfolder/Great_Barrier_Reef_Features.shp",
              sep = ""))
reefs <- subset(gbr,FEAT_NAME == "Reef")

AUS_states <- ne_states(country = "australia", returnclass = "sf")
AUS_states <- cbind(AUS_states, st_coordinates(st_centroid(AUS_states)))

# Figure 2
fig2 <- ggplot(data = AUS_states) +
            geom_sf(fill ="darkgray") +
            geom_sf(data = reefs, fill ="antiquewhite", color = gray(0.8))  + 
            geom_point(data = spatial.dat, mapping = aes(LONGITUDE,LATITUDE),
               shape = 21,color="black",fill="white",size = 2,stroke = 0.5) + 
            geom_point(filt.rec.traj.proc.reg.spat, 
                       mapping = aes(LONGITUDE,LATITUDE),fill="black",
                       size=2,alpha=1,shape = 21,color = "black",stroke=1) +
            annotate(geom = "text", x = 144.55, y = -19, label = "Queensland", 
                     color = "white", size = 6) +
            annotate(geom = "text", x = 148.5, y = -15.5, 
                     label = "Great Barrier Reef", color = "grey22", size = 6) +
            annotation_scale(location = "bl", width_hint = 0.5) +
            annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering()) +
            coord_sf(xlim = c(min(spatial.dat$LONGITUDE)-0.5,
                              max(spatial.dat$LONGITUDE)+0.5), 
                     ylim = c(min(spatial.dat$LATITUDE)-0.5, 
                              max(spatial.dat$LATITUDE)+0.5), expand = FALSE) +
            xlab("Longitude") + ylab("Latitude") + 
            theme(panel.grid.major = element_line(color=gray(0.5), 
                                                  linetype = "dashed", 
                                                  size = 0.5),
                  panel.background = element_rect(fill="aliceblue"))

fig2
# Figure 6A
fig6a <- ggplot(data = AUS_states) +
            geom_sf(fill ="darkgray") +
            geom_sf(data = reefs, fill ="antiquewhite", color = gray(0.8))  + 
            geom_point(filt.rec.traj.proc.reg.spat %>% filter(CLASSIF == "D"), 
                       mapping = aes(LONGITUDE,LATITUDE,size = DELAY_DURATION/365),
                       fill = "red",shape = 21,alpha=0.3,color = "black",stroke=1) + scale_radius() +
            annotation_scale(location = "bl", width_hint = 0.5) +
            annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering()) +
            coord_sf(xlim = c(145.5,max(spatial.dat$LONGITUDE)+0.5), 
                     ylim = c(min(spatial.dat$LATITUDE)-0.5, -16), 
                     expand = FALSE) +
            labs(x = "Longitude", y = "Latitude", size = "1st Phase (years)") +
            theme(panel.grid.major = element_line(color=gray(0.5), 
                                                  linetype = "dashed", 
                                                  size = 0.5),
                  panel.background = element_rect(fill="aliceblue"))
fig6a
# Figure 6B
fig6b <- ggplot(data = AUS_states) +
            geom_sf(fill ="darkgray") +
            geom_sf(data = reefs, fill ="antiquewhite", color = gray(0.8))  + 
            geom_point(filt.rec.traj.proc.reg.spat %>% filter(CLASSIF!="D"), 
                       mapping = aes(LONGITUDE,LATITUDE,fill = CLASSIF),size = 3, 
                       alpha=0.5,shape = 21,color = "black",stroke=1) +
            scale_fill_manual(breaks = c("?","L/G"), 
                              labels = c("Uncomfirmed","Single-phase"),
                              values = c("blue","green")) +
            annotation_scale(location = "bl", width_hint = 0.5) +
            annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering()) +
            coord_sf(xlim = c(145.5,max(spatial.dat$LONGITUDE)+0.5), 
                     ylim = c(min(spatial.dat$LATITUDE)-0.5, -16), 
                     expand = FALSE) +
            labs(x = "Longitude", y = "Latitude", fill = "Class") + 
            theme(panel.grid.major = element_line(color=gray(0.5), 
                                                  linetype = "dashed", 
                                                  size = 0.5),
                  panel.background = element_rect(fill="aliceblue"))
fig6b
# Figure 6
fig6 <- grid.arrange(grobs = list(fig6a,fig6b),nrow = 1,ncol = 2)

# Figure 5A
fig5a <- ggplot(filt.rec.traj.proc.reg %>% filter(DELAY_DURATION > 0),
               aes(x = DELAY_DURATION/365)) +
        geom_histogram(color="darkblue",bins=10,fill="lightblue") +
        labs(x="1st phase duration (years)",y="Number of recovery trajectories") + 
        theme_bw()
fig5a

# Figure 5B
fig5b <- ggplot(filt.rec.traj.proc.reg %>% filter(DELAY_DURATION > 0),
               aes(x = DELAY_HC_COVER)) +
  geom_histogram(color="darkblue",bins=10,fill="lightblue") +
  labs(x="1st phase final cover (%)",y="Number of recovery trajectories") + 
  theme_bw()
fig5b
# Figure 5
fig5 <- grid.arrange(grobs = list(fig5a,fig5b),nrow = 1,ncol = 2)

# load and format power analysis data
power_analysis_table <- read_csv("power_analysis_table.csv")

power_analysis_table_annual <- power_analysis_table %>% 
            select(Power = Power_annual,Specificity = Specificity_annual,
                   alpha,alpha_scale,model,DTime,obs_error) %>% 
            mutate(Sampling = "Annual")

            
power_analysis_table_biennial <- power_analysis_table %>% 
            select(Power = Power_biennial,Specificity = Specificity_biennial,
                   alpha,alpha_scale,model,DTime,obs_error) %>% 
            mutate(Sampling = "Biennial")

power_analysis_table_grouped <- rbind(power_analysis_table_annual,
                                      power_analysis_table_biennial)
# Figure A4
fig4a <- ggplot(power_analysis_table_grouped,aes(x=Power))  + 
            geom_histogram(position = "identity",binwidth = 0.050,
                           fill = 'lightblue',color = 'blue') + 
            labs(x = 'Power',y = 'Number of simulations') + 
            theme_bw() 
fig4a

#Figure 4B
fig4b <- ggplot(power_analysis_table_grouped,aes(x=Specificity))  + 
            geom_histogram(position = "identity",binwidth = 0.050,
                           fill = 'lightblue',color = 'blue') + 
            labs(x = 'Specificity',y = 'Number of simulations') +
            xlim(c(0,1))+
            theme_bw()
fig4b
fig4 <- grid.arrange(grobs = list(fig4a,fig4b),nrow = 1,ncol = 2)

# Figure C2 (Supporting Material)
figC2 <- ggplot(power_analysis_table, aes(y = Power_annual, x = Power_biannual)) +
            geom_point() +
            geom_abline(slope = 1, intercept = 0) +
            xlim(c(0,1)) +
            ylim(c(0,1)) +
            labs(x = "Power (biennual sampling)",
                 y = "Power (annual sampling)")  +
            theme_bw()
figC2

# plot figure 7 errorbar style
d <- 0.5
fig7a <- ggplot(sumFinalCover,aes(x=YEARS,y=q2/100.0, color=DELAY))+ ggtitle("Expected maximum cover") +
        geom_errorbar(aes(ymin=q1/100,ymax = q3/100),width = d , 
                  position = position_dodge(d)) + 
        geom_line(linetype = "dashed",position = position_dodge(d)) +
        geom_point(position = position_dodge(d)) +
        scale_color_brewer(palette="Dark2") +
        ylim(c(0,1)) +
        scale_x_continuous(breaks = c(5,7.5,10), labels = c("5","7.5","10")) +
        scale_y_continuous(labels = percent) +
        labs(x = 'Major disturbance frequency (years)',
             y = 'Cover', color = 'Recovery pattern') +
        theme_bw() + theme(
            legend.position = c(0.5, .99),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(1, 1, 1, 1),
            legend.box.background = element_rect(color="black", size=1)
        )
d <- 1
fig7b <- ggplot(sumYearsAbove,aes(x=COVER,y=q210,color=DELAY)) + ggtitle("10-year disturbance frequency") +  
        geom_errorbar(aes(ymin=q110,ymax = q310),width =d , 
                    position = position_dodge(d)) + 
        geom_line(linetype = "dashed",position = position_dodge(d)) +
        geom_point(position = position_dodge(d)) +
        scale_color_brewer(palette="Dark2") +
        scale_x_continuous(breaks = c(10,15,20), labels = c("10%","15%","20%")) +
        ylim(c(0,10)) +
        labs(x = 'Cover threshold',
             y = 'Time above threshold (years/decade)', color = 'Recovery pattern') + 
        theme_bw() +
        theme(legend.position = "none")

fig7c <- ggplot(sumYearsAbove,aes(x=COVER,y=q25,color=DELAY))  + ggtitle("5-year disturbance frequency") + 
        geom_errorbar(aes(ymin=q15,ymax = q35),width = d , 
                      position = position_dodge(d)) + 
        geom_line(linetype = "dashed",position = position_dodge(d)) +
        geom_point(position = position_dodge(d)) +
        scale_color_brewer(palette="Dark2") +
        scale_x_continuous(breaks = c(10,15,20), labels = c("10%","15%","20%")) +
        ylim(c(0,10)) +
        labs(x = 'Cover threshold',
             y = 'Time above threshold (years/decade)', color ='Recovery pattern') + 
        theme_bw() +
        theme(legend.position = "none")
    
fig7 <- grid.arrange(grobs = list(fig7a,fig7b,fig7c),nrow = 1,ncol = 3)

# SI version plot figure 7 errorbar style with 1 year disturbances
d <- 0.5
figD7a <- ggplot(sumFinalCover,aes(x=YEARS,y=q2/100.0, color=DELAY)) + ggtitle("Expected maximum cover")+
  geom_errorbar(aes(ymin=q1/100,ymax = q3/100),width = d , 
                position = position_dodge(d)) + 
  geom_line(linetype = "dashed",position = position_dodge(d)) +
  geom_point(position = position_dodge(d)) +
  scale_color_brewer(palette="Dark2") +
  ylim(c(0,1)) +
  scale_x_continuous(breaks = c(1,5,7.5,10), labels = c("1","5","7.5","10")) +
  scale_y_continuous(labels = percent) +
  labs(x = 'Major disturbance frequency (years)',
       y = 'Cover', color = 'Recovery pattern') +
  theme_bw() + theme(
    legend.position = c(0.5, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(1, 1, 1, 1),
    legend.box.background = element_rect(color="black", size=1)
  )
d <- 1
figD7b <- ggplot(sumYearsAbove,aes(x=COVER,y=q210,color=DELAY)) + ggtitle("10-year disturbance frequency") + 
  geom_errorbar(aes(ymin=q110,ymax = q310),width =d , 
                position = position_dodge(d)) + 
  geom_line(linetype = "dashed",position = position_dodge(d)) +
  geom_point(position = position_dodge(d)) +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks = c(10,15,20), labels = c("10%","15%","20%")) +
  ylim(c(0,10)) +
  labs(x = 'Cover threshold',
       y = 'Time above threshold (years/decade)', color = 'Recovery pattern') + 
  theme_bw() +
  theme(legend.position = "none")

figD7c <- ggplot(sumYearsAbove,aes(x=COVER,y=q25,color=DELAY)) + ggtitle("5-year disturbance frequency")  + 
  geom_errorbar(aes(ymin=q15,ymax = q35),width = d , 
                position = position_dodge(d)) + 
  geom_line(linetype = "dashed",position = position_dodge(d)) +
  geom_point(position = position_dodge(d)) +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks = c(10,15,20), labels = c("10%","15%","20%")) +
  ylim(c(0,10)) +
  labs(x = 'Cover threshold',
       y = 'Time above threshold (years/decade)', color ='Recovery pattern') + 
  theme_bw() +
  theme(legend.position = "none")

figD7d <- ggplot(sumYearsAbove,aes(x=COVER,y=q21,color=DELAY)) + ggtitle("1-year disturbance frequency")  + 
  geom_errorbar(aes(ymin=q11,ymax = q31),width = d , 
                position = position_dodge(d)) + 
  geom_line(linetype = "dashed",position = position_dodge(d)) +
  geom_point(position = position_dodge(d)) +
  scale_color_brewer(palette="Dark2") +
  scale_x_continuous(breaks = c(10,15,20), labels = c("10%","15%","20%")) +
  ylim(c(0,10)) +
  labs(x = 'Cover threshold',
       y = 'Time above threshold (years/decade)', color ='Recovery pattern') + 
  theme_bw() +
  theme(legend.position = "none")
figD7 <- grid.arrange(grobs = list(figD7a,figD7b,figD7c,figD7d),nrow = 2,ncol = 2)

# SI disturbance magnitude plot
figD6 <- ggplot(filt.rec.traj.proc.reg %>% group_by(CLASSIF),
                aes(x = DIFF)) +
  geom_histogram(aes(fill=CLASSIF),color="black",bins=8,alpha=0.6) +
  scale_fill_manual(breaks = c("?","D","L/G"), 
                    labels = c("Uncomfirmed","Two-phase","Single-phase"),
                    values = c("blue","red","green")) +
  labs(x="Disturbance magnitude (% cover decrease)",y="Number of recovery trajectories",fill="Recovery pattern") + theme_bw()
figD6

# can arrange by F-statistic/p-Value or Adj-R^2
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% arrange(desc(ADJ_R2_CP)) 
filt.rec.traj.proc.reg <- filt.rec.traj.proc.reg %>% mutate(REL_ADJ_R2 = (ADJ_R2_CP - ADJ_R2) / abs(ADJ_R2))
figs <- list()
nD   <- 1
figsL1 <- list()
nL     <- 1
figsD1 <- list()
nD     <- 1
figsU1 <- list()
nU    <- 1
numDelay <- 0
numLog   <- 0
numUnk   <- 0
for (j in 1:length(filt.rec.traj.proc)) {
  
  k     <- filt.rec.traj.proc.reg$INDEX[j]
  reef  <- filt.rec.traj.proc[[k]]$REEF[1]
  p_code <- filt.rec.traj.proc[[k]]$P_CODE[1]
  depth <- filt.rec.traj.proc[[k]]$DEPTH[1]
  rp_id <- filt.rec.traj.proc[[k]]$RP_ID[1]
  site  <- filt.rec.traj.proc[[k]]$SITE_NO[1]
  
  # sanity check !
  if (rp_id != filt.rec.traj.proc.reg$RP_ID[j] || 
      site != filt.rec.traj.proc.reg$SITE_NO[j]) {
    warning('mismatch in regression/per-capita data!')
  }
  
  
  if (filt.rec.traj.proc.reg$CLASSIF[j] %in% c('L/G')) {
    numLog <- numLog + 1
    figsL1[[nL]] <- ggplot(data=filt.rec.traj.proc[[k]],aes(x=Date,y=HC)) + 
      ggtitle(sprintf("%s Site %d %d[m]",reef,site,depth)) +
      geom_line(linetype = "solid", color="blue",size=1) + 
      geom_point() +
      geom_errorbar(aes(ymin=HC-HC_se,ymax=HC+HC_se)) +
      geom_hline(yintercept = filt.rec.traj.proc[[k]]$HC[1] + filt.rec.traj.proc[[k]]$DIFF[1],
                 color = "blue",linetype="dashed") + 
      labs(x="Date",y="Coral cover") + theme_bw()
    nL <- nL + 1
    figsL1[[nL]] <- ggplot() + 
      ggtitle(sprintf("Adj-R^2 %0.2f, Rel Adj-R^2 %0.2f",
                      filt.rec.traj.proc.reg$ADJ_R2_CP[j],filt.rec.traj.proc.reg$REL_ADJ_R2[j])) + 
      geom_errorbar(aes(x=HC,ymin=pcHC_lin_qlb,ymax=pcHC_lin_qub),
                    filt.rec.traj.proc[[k]])
    sl1 <- filt.rec.traj.proc.reg$SLOPE1[j]
    int1 <- filt.rec.traj.proc.reg$INTERCEPT1[j]
    sl2 <- filt.rec.traj.proc.reg$SLOPE2[j]
    int2 <- filt.rec.traj.proc.reg$INTERCEPT2[j]
    min_HC <- min(filt.rec.traj.proc[[k]]$HC)
    max_HC <- max(filt.rec.traj.proc[[k]]$HC)
    cp_HC <- filt.rec.traj.proc[[k]]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[j]]
    x1 <- min_HC
    y1 <- min_HC*sl1 + int1
    xend1  <- cp_HC
    yend1 <- cp_HC*sl1 + int1
    x2 <- cp_HC
    y2 <- cp_HC*sl2 + int2
    xend2 <- max_HC
    yend2 <- max_HC*sl2 + int2
    eval(parse(text = paste0(
      'figsL1[[nL]] <- figsL1[[nL]] +
          geom_vline(xintercept = filt.rec.traj.proc[[',k,']]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[',j,']],color = "red",linetype="dashed") +
          geom_point(aes(HC,pcHC_lin_fd),filt.rec.traj.proc[[',k,']]) + 
          geom_segment(aes(x = ',x1,',y =', y1,', xend =', xend1,', yend = ',yend1,'),size=1, color = "red",linetype = "solid") +
          geom_segment(aes(x = ',x2,', y = ',y2,', xend =', xend2,', yend = ',yend2,'),size=1, color = "green",linetype = "solid") +
          labs(x="Cover",y="Rate per unit cover")  + theme_bw()')))
    nL <- nL + 1
  } else if (filt.rec.traj.proc.reg$CLASSIF[j] %in% c('D')) {
    numDelay <- numDelay + 1
    figsD1[[nD]] <- ggplot(data=filt.rec.traj.proc[[k]],aes(x=Date,y=HC)) + 
      ggtitle(sprintf("%s Site %d %d[m]",reef,site,depth)) +
      geom_line(linetype = "solid", color="blue",size=1) + 
      geom_point() +
      geom_errorbar(aes(ymin=HC-HC_se,ymax=HC+HC_se)) +
      geom_vline(xintercept = filt.rec.traj.proc[[k]]$Date[filt.rec.traj.proc.reg$CHANGE_POINT[j]],
                 color = "red",linetype="dashed") +
      geom_hline(yintercept = filt.rec.traj.proc[[k]]$HC[1] + filt.rec.traj.proc[[k]]$DIFF[1],
                 color = "blue",linetype="dashed") +
      labs(x="Date",y="Coral cover") + theme_bw()
    nD <- nD + 1
    figsD1[[nD]] <- ggplot() + ggtitle(sprintf("Adj-R^2 %0.2f, Rel Adj-R^2 %0.2f",
                                               filt.rec.traj.proc.reg$ADJ_R2_CP[j],filt.rec.traj.proc.reg$REL_ADJ_R2[j])) + 
      geom_errorbar(aes(x=HC,ymin=pcHC_lin_qlb,ymax=pcHC_lin_qub),
                    filt.rec.traj.proc[[k]])
    sl1 <- filt.rec.traj.proc.reg$SLOPE1[j]
    int1 <- filt.rec.traj.proc.reg$INTERCEPT1[j]
    sl2 <- filt.rec.traj.proc.reg$SLOPE2[j]
    int2 <- filt.rec.traj.proc.reg$INTERCEPT2[j]
    min_HC <- min(filt.rec.traj.proc[[k]]$HC)
    max_HC <- max(filt.rec.traj.proc[[k]]$HC)
    cp_HC <- filt.rec.traj.proc[[k]]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[j]]
    x1 <- min_HC
    y1 <- min_HC*sl1 + int1
    xend1  <- cp_HC
    yend1 <- cp_HC*sl1 + int1
    x2 <- cp_HC
    y2 <- cp_HC*sl2 + int2
    xend2 <- max_HC
    yend2 <- max_HC*sl2 + int2
    eval(parse(text = paste0(
      'figsD1[[nD]] <- figsD1[[nD]] +
          geom_vline(xintercept = filt.rec.traj.proc[[',k,']]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[',j,']],color = "red",linetype="dashed") +
          geom_point(aes(HC,pcHC_lin_fd),filt.rec.traj.proc[[',k,']]) + 
          geom_segment(aes(x = ',x1,',y =', y1,', xend =', xend1,', yend = ',yend1,'),size=1, color = "red",linetype = "solid") +
          geom_segment(aes(x = ',x2,', y = ',y2,', xend =', xend2,', yend = ',yend2,'),size=1, color = "green",linetype = "solid") +
          labs(x="Cover",y="Rate per unit cover")  + theme_bw()'))) 
    nD <- nD + 1
    #       }
  } else if (filt.rec.traj.proc.reg$CLASSIF[j] %in% c('?')) {
    numUnk <- numUnk + 1
    figsU1[[nU]] <- ggplot(data=filt.rec.traj.proc[[k]],aes(x=Date,y=HC)) + 
      ggtitle(sprintf("%s Site %d %d[m]",reef,site,depth)) +
      
      geom_line(linetype = "solid", color="blue",size=1) + 
      geom_point() +
      geom_errorbar(aes(ymin=HC-HC_se,ymax=HC+HC_se)) +
      geom_hline(yintercept = filt.rec.traj.proc[[k]]$HC[1] + filt.rec.traj.proc[[k]]$DIFF[1],
                 color = "blue",linetype="dashed") +
      labs(x="Date",y="Cover") + theme_bw()
    nU <- nU + 1
    figsU1[[nU]] <- ggplot() + ggtitle(sprintf("Adj-R^2 %0.2f, Rel Adj-R^2 %0.2f",
                                               filt.rec.traj.proc.reg$ADJ_R2_CP[j],filt.rec.traj.proc.reg$REL_ADJ_R2[j])) + 
      geom_errorbar(aes(x=HC,ymin=pcHC_lin_qlb,ymax=pcHC_lin_qub),
                    filt.rec.traj.proc[[k]])
    sl1 <- filt.rec.traj.proc.reg$SLOPE1[j]
    int1 <- filt.rec.traj.proc.reg$INTERCEPT1[j]
    sl2 <- filt.rec.traj.proc.reg$SLOPE2[j]
    int2 <- filt.rec.traj.proc.reg$INTERCEPT2[j]
    min_HC <- min(filt.rec.traj.proc[[k]]$HC)
    max_HC <- max(filt.rec.traj.proc[[k]]$HC)
    cp_HC <- filt.rec.traj.proc[[k]]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[j]]
    x1 <- min_HC
    y1 <- min_HC*sl1 + int1
    xend1  <- cp_HC
    yend1 <- cp_HC*sl1 + int1
    x2 <- cp_HC
    y2 <- cp_HC*sl2 + int2
    xend2 <- max_HC
    yend2 <- max_HC*sl2 + int2
    eval(parse(text = paste0(
      'figsU1[[nU]] <- figsU1[[nU]] +
          geom_vline(xintercept = filt.rec.traj.proc[[',k,']]$HC[filt.rec.traj.proc.reg$CHANGE_POINT[',j,']],color = "red",linetype="dashed") +
          geom_point(aes(HC,pcHC_lin_fd),filt.rec.traj.proc[[',k,']]) + 
          geom_segment(aes(x = ',x1,',y =', y1,', xend =', xend1,', yend = ',yend1,'), size=1, color = "red",linetype = "solid") +
          geom_segment(aes(x = ',x2,', y = ',y2,', xend =', xend2,', yend = ',yend2,'), size=1, color = "green",linetype = "solid") +
          labs(x="Cover",y="Rate per unit cover")  + theme_bw()')))
    nU <- nU + 1
  }
}

sd <- summary(filt.rec.traj.proc.reg$DELAY_DURATION[filt.rec.traj.proc.reg$DELAY_DURATION > 0 & filt.rec.traj.proc.reg$DELAY_DURATION <= 3650.0]/365)

figD1 <- grid.arrange(grobs = figsL1,nrow=6,ncol=16)
figD2 <- grid.arrange(grobs = figsD1,nrow=8,ncol=16)
figD3 <- grid.arrange(grobs = figsU1,nrow=3,ncol=16)
