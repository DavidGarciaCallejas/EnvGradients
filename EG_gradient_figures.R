
source("./EG_population_structure.R")

############################
############################

grid.rows <- 50
timesteps <- 500
# transient.phase <- 100
saturated <- TRUE
facilitation <- TRUE

my.file.name <- paste(grid.rows,"rows_",timesteps,"t",sep="")
if(!facilitation){
  my.file.name <- paste("NO_FACILITATION_",my.file.name,sep="")
}

average.interaction.biomass <- readr::read_delim(paste("./results/average_interaction_biomass_",my.file.name,".csv",sep=""),delim = ";")
average.diversity <- readr::read_delim(paste("./results/average_diversity_",my.file.name,".csv",sep=""),delim = ";")
average.persistence <- readr::read_delim(paste("./results/average_persistence_",my.file.name,".csv",sep=""),delim = ";")

##################################
##################################
# gradient figures

######################
# average persistence
avg.persistence.plot <- ggplot(average.persistence) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = mean.persistence),interpolate = FALSE) + #scale_x_reverse() +
  #facet_wrap(~interaction) + 
  scale_fill_viridis_c(name = "") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("Average persistence time") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white"))
# theme_Publication()
# avg.persistence.plot

######################
# normalised persistence
norm.persistence.plot <- ggplot(average.persistence) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = norm.persistence),interpolate = FALSE) + #scale_x_reverse() +
  #facet_wrap(~interaction) + 
  scale_fill_viridis_c(name = "") + #,direction = -1) + 
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  # scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("Average persistence time") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) + 
  guides(fill = FALSE)
# theme_Publication()
# avg.persistence.plot

######################
# intensity of interaction types (regarding biomass involved), averaged across all timesteps
avg.interaction.plot <- ggplot(average.interaction.biomass) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = average.biomass),interpolate = FALSE) + #scale_x_reverse() +
  facet_wrap(~interaction) + 
  scale_fill_viridis_c(name = "") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("Interaction intensity") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white"))
# theme_Publication()
# avg.interaction.plot

######################
# normalised intensity of interaction types
# it is actually the same as previous plot, without legend
norm.interaction.plot <- ggplot(average.interaction.biomass) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = average.biomass),interpolate = FALSE) + #scale_x_reverse() +
  facet_wrap(~interaction) + 
  scale_fill_viridis_c(name = "") + #,direction = -1) + 
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  #scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("Interaction intensity") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) + 
  guides(fill = FALSE)
# theme_Publication()
# avg.interaction.plot

######################
# average hill diversity
avg.hill.plot <- ggplot(average.diversity) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = average.hill),interpolate = FALSE) + #scale_x_reverse() +
  # facet_wrap(~interaction) +
  scale_fill_viridis_c(name = "") + #, option = "plasma") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("effective number of species") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white"))
# theme_Publication()
# avg.hill.plot

######################
# normalised hill diversity
norm.hill.plot <- ggplot(average.diversity) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = norm.diversity),interpolate = FALSE) + #scale_x_reverse() +
  # facet_wrap(~interaction) +
  scale_fill_viridis_c(name = "") + #, option = "plasma") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle("effective number of species") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) + 
  guides(fill = FALSE)
# theme_Publication()
# avg.hill.plot

######################

# plots in their own scale
legend.plot <- avg.interaction.plot + {avg.hill.plot + avg.persistence.plot} + plot_layout(ncol=1)

# normalised plots
normalised.plot <- norm.interaction.plot + {norm.hill.plot + norm.persistence.plot} + plot_layout(ncol=1)

# tiff(paste("/home/david/CREAF/FPU/plant_interactions/results/images/Fig_2_",my.file.name,".tiff",sep=""), res=600, compression = "lzw", width = 4000, height = 4000, units = "px")
# normalised.plot
# legend.plot
# dev.off()

##############################
# facilitation versus no-facilitation plots

my.file.name <- paste(grid.rows,"rows_",timesteps,"t",sep="")
facilitation.file.name <- paste("NO_FACILITATION_",my.file.name,sep="")

average.interaction.biomass <- readr::read_delim(paste("./results/average_interaction_biomass_",my.file.name,".csv",sep=""),delim = ";")
average.diversity <- readr::read_delim(paste("./results/average_diversity_",my.file.name,".csv",sep=""),delim = ";")
average.persistence <- readr::read_delim(paste("./results/average_persistence_",my.file.name,".csv",sep=""),delim = ";")

average.interaction.biomass.nof <- readr::read_delim(paste("./results/average_interaction_biomass_",facilitation.file.name,".csv",sep=""),delim = ";")
average.diversity.nof <- readr::read_delim(paste("./results/average_diversity_",facilitation.file.name,".csv",sep=""),delim = ";")
average.persistence.nof <- readr::read_delim(paste("./results/average_persistence_",facilitation.file.name,".csv",sep=""),delim = ";")

competition.fac <- subset(average.interaction.biomass,interaction == "competition")
competition.nof <- subset(average.interaction.biomass.nof,interaction == "competition")
# small hack for adding zero to competition intensities not listed in the simulations without facilitation
missing.nrsf <- unique(competition.fac$NRSF)
missing.nrsf <- missing.nrsf[which(!(missing.nrsf %in% unique(competition.nof$NRSF)))]
missing.nrsf <- expand.grid(missing.nrsf,unique(competition.fac$RSF))#data.frame(NRSF = missing.nrsf,RSF = unique(competition.fac$RSF))
names(missing.nrsf) <- c("NRSF","RSF")

missing.nrsf$interaction <- "competition"
missing.nrsf$average.biomass <- 0
competition.nof = bind_rows(competition.nof,missing.nrsf)

competition.fac = arrange(competition.fac,RSF,NRSF)
competition.nof = arrange(competition.nof,RSF,NRSF)

average.diff = data.frame(RSF = competition.fac$RSF,NRSF = competition.fac$NRSF, 
                          diversity = 0, persistence = 0, competition = 0)

full.grid = expand(average.diff,RSF,NRSF)

average.diversity = left_join(full.grid,average.diversity)
average.diversity.nof = left_join(full.grid,average.diversity.nof)
average.persistence = left_join(full.grid,average.persistence)
average.persistence.nof = left_join(full.grid,average.persistence.nof)

average.diversity$average.hill[is.na(average.diversity$average.hill)] = 0
average.diversity.nof$average.hill[is.na(average.diversity.nof$average.hill)] = 0

average.persistence$mean.persistence[is.na(average.persistence$mean.persistence)] = 0
average.persistence.nof$mean.persistence[is.na(average.persistence.nof$mean.persistence)] = 0

average.persistence = arrange(average.persistence,RSF,NRSF)
average.persistence.nof = arrange(average.persistence.nof,RSF,NRSF)

average.diversity = arrange(average.diversity,RSF,NRSF)
average.diversity.nof = arrange(average.diversity.nof,RSF,NRSF)

######

average.diff = arrange(average.diff,RSF,NRSF)

average.diff$competition = competition.fac$average.biomass - competition.nof$average.biomass
# average.diff$competition = LinMap(average.diff$competition,-1,1)

average.diff$diversity = average.diversity$average.hill - average.diversity.nof$average.hill
# average.diff$diversity = LinMap(average.diff$diversity,-1,1)

average.diff$persistence = average.persistence$mean.persistence - average.persistence.nof$mean.persistence
# average.diff$persistence = LinMap(average.diff$persistence,-1,1)

# average.diff = gather(average.diff,key = "metric",value = "value",diversity,persistence,competition)

##############################
# plot

diff.competition <- ggplot(average.diff) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = competition),interpolate = FALSE) + #scale_x_reverse() +
  # facet_wrap(~metric) +
  scale_fill_viridis_c(name = "",option = "inferno") + #, option = "plasma") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle(label = "",subtitle = "competition") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) +
  NULL

diff.diversity <- ggplot(average.diff) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = diversity),interpolate = FALSE) + #scale_x_reverse() +
  # facet_wrap(~metric) +
  scale_fill_viridis_c(name = "",option = "inferno") + #, option = "plasma") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle(label = "",subtitle = "Hill diversity") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) +
  NULL

diff.persistence <- ggplot(average.diff) + 
  geom_raster(aes(x = RSF, y = NRSF, fill = persistence),interpolate = FALSE) + #scale_x_reverse() +
  # facet_wrap(~metric) +
  scale_fill_viridis_c(name = "",option = "inferno") + #, option = "plasma") + #,direction = -1) + 
  # scale_x_discrete(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  scale_x_reverse(breaks = NULL,expand = c(0,0)) + scale_y_discrete(breaks = NULL,expand = c(0,0)) +
  xlab("Resource stress") + ylab("Non-resource stress") + 
  ggtitle(label = "",subtitle = "persistence time") +
  # geom_hline(aes(yintercept = 10),color = "lightgrey") + geom_vline(aes(xintercept = 10), color = "lightgrey") + 
  theme(strip.text.x = element_text(size=10, face = "bold"),
        strip.background = element_rect(colour="white", fill="white")) +
  NULL
# theme_Publication()

tiff(paste("./results/images/metric_differences_",my.file.name,".tiff",sep=""), res=600, compression = "lzw", width = 7000, height = 2000, units = "px")
diff.competition + diff.diversity + diff.persistence
dev.off()




