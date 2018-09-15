
# statistical models. GAM

############################
############################

grid.rows <- 50
timesteps <- 500
# transient.phase <- 100
# nevermind this flag
saturated <- TRUE

facilitation <- TRUE

my.file.name <- paste(grid.rows,"rows_",timesteps,"t",sep="")
facilitation.file.name <- paste("NO_FACILITATION_",my.file.name,sep="")

average.interaction.biomass <- readr::read_delim(paste("./Results/average_interaction_biomass_",my.file.name,".csv",sep=""),delim = ";")
average.diversity <- readr::read_delim(paste("./Results/average_diversity_",my.file.name,".csv",sep=""),delim = ";")
average.persistence <- readr::read_delim(paste("./Results/average_persistence_",my.file.name,".csv",sep=""),delim = ";")

average.interaction.biomass.nof <- readr::read_delim(paste("./Results/average_interaction_biomass_",facilitation.file.name,".csv",sep=""),delim = ";")
average.diversity.nof <- readr::read_delim(paste("./Results/average_diversity_",facilitation.file.name,".csv",sep=""),delim = ";")
average.persistence.nof <- readr::read_delim(paste("./Results/average_persistence_",facilitation.file.name,".csv",sep=""),delim = ";")


######################
# statistical analyses

# 0 - effect of facilitation on community properties

# perform wilcoxon signed-rank tests on the pair of values (with-without facilitation) from each cell
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

competition.test <- wilcox.test(x = competition.fac$average.biomass,y = competition.nof$average.biomass,
                                paired = TRUE)
median(competition.fac$average.biomass)
median(competition.nof$average.biomass)

# persistence test
persistence.test <- wilcox.test(x = average.persistence$mean.persistence, y = average.persistence.nof$mean.persistence,
                                paired = TRUE,alternative = "greater")
median(average.persistence$mean.persistence)
median(average.persistence.nof$mean.persistence)

# diversity test
diversity.test <- wilcox.test(x = average.diversity$average.hill, y = average.diversity.nof$average.hill,
                              paired = TRUE,alternative = "greater")
median(average.diversity$average.hill)
median(average.diversity.nof$average.hill)

#########################
#########################
# 1 - data visualization, uncomment etc

# hist(average.interaction.biomass$average.biomass[average.interaction.biomass$interaction == "facilitation"],100,col = "black")
# hist(average.interaction.biomass$average.biomass[average.interaction.biomass$interaction == "competition"],100,col = "black")
# 
# hist(average.diversity$log.diversity,100,col = "black")
# hist(log(average.diversity$log.diversity),100,col = "black")
# 
# hist(average.persistence$mean.persistence,100,col = "black")
# # log does not make distribution normal
# hist(log(average.persistence$mean.persistence),100,col = "black")

# 2 - first models

interaction.spread <- spread(average.interaction.biomass,key = interaction,value = average.biomass)
interaction.spread$facilitation[interaction.spread$facilitation > 1] <- 1
#####################
# GAM approach

# 1 - plot the variation of each response to each predictor

# competition
competition.RSF.plot <- ggplot(interaction.spread,aes(x = RSF,y = competition, group = RSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("competition intensity") +
  scale_x_reverse()+
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "cs"),aes(x = RSF,y = competition, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "a") +
  NULL

competition.NRSF.plot <- ggplot(interaction.spread,aes(x = NRSF,y = competition, group = NRSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("competition intensity") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = competition, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "b") +
  NULL

# facilitation
facilitation.RSF.plot <- ggplot(interaction.spread,aes(x = RSF,y = facilitation, group = RSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("facilitation intensity") +
  scale_x_reverse()+
  # geom_smooth(method = "lm",aes(x = RSF,y = facilitation,group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "c") +
  NULL

facilitation.NRSF.plot <- ggplot(interaction.spread,aes(x = NRSF,y = facilitation, group = NRSF)) + 
  # geom_point()+
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("facilitation intensity") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = facilitation, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "d") +
  NULL
# facilitation.NRSF.plot

# diversity
diversity.RSF.plot <- ggplot(average.diversity,aes(x = RSF,y = average.hill, group = RSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("effective richness") +
  scale_x_reverse()+
  geom_smooth(method = "lm",aes(x = RSF,y = average.hill, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "e") +
  NULL

diversity.NRSF.plot <- ggplot(average.diversity,aes(x = NRSF,y = average.hill, group = NRSF)) + 
  # geom_point() +
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("effective richness") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = average.hill, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "f") +
  NULL
# diversity.NRSF.plot

log.diversity.RSF.plot <- ggplot(average.diversity,aes(x = RSF,y = log(average.hill), group = RSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("log(effective richness)") +
  scale_x_reverse()+
  geom_smooth(method = "lm",aes(x = RSF,y = log(average.hill), group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "e") +
  NULL

log.diversity.NRSF.plot <- ggplot(average.diversity,aes(x = NRSF,y = log(average.hill), group = NRSF)) + 
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("log(effective richness)") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = log(average.hill), group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "f") +
  NULL

# persistence
persistence.RSF.plot <- ggplot(average.persistence,aes(x = RSF,y = mean.persistence, group = RSF)) +
  # geom_point()+
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("persistence times") +
  scale_x_reverse()+
  # geom_smooth(method = "lm",aes(x = RSF,y = mean.persistence,group = 1)) +
  # geom_smooth(method = "lm",aes(x = RSF,y = norm.persistence,group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "g") +
  NULL
# persistence.RSF.plot

persistence.NRSF.plot <- ggplot(average.persistence,aes(x = NRSF,y = mean.persistence, group = NRSF)) + 
  # geom_point()+
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("persistence times") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = mean.persistence, group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "h") +
  NULL

log.persistence.RSF.plot <- ggplot(average.persistence,aes(x = RSF,y = log(mean.persistence), group = RSF)) + 
  # geom_point()+
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("resource capacity") + ylab("log(persistence times)") +
  scale_x_reverse()+
  geom_smooth(method = "lm",aes(x = RSF,y = log(mean.persistence),group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "g") +
  NULL
# log.persistence.RSF.plot

log.persistence.NRSF.plot <- ggplot(average.persistence,aes(x = NRSF,y = log(mean.persistence), group = NRSF)) + 
  # geom_point()+
  geom_boxplot(fill = cbPalette()[2],outlier.size = .8) +
  xlab("non-resource stress") + ylab("log(persistence times)") +
  geom_smooth(method = "gam",formula = y ~ s(x,bs = "ad"),aes(x = NRSF,y = log(mean.persistence), group = 1)) +
  theme_Publication() +
  ggtitle(label = "",subtitle = "h") +
  NULL
# log.persistence.NRSF.plot
######
# generate a single plot (load patchwork package)

# tiff(paste("./results/images/individual_gradients.tiff",sep=""), res=600, compression = "lzw", width = 6000, height = 8000, units = "px")
# competition.RSF.plot + 
#   competition.NRSF.plot + 
#   facilitation.RSF.plot + 
#   facilitation.NRSF.plot +
#   diversity.RSF.plot + 
#   diversity.NRSF.plot + 
#   persistence.RSF.plot + 
#   persistence.NRSF.plot + 
#   # log.diversity.RSF.plot + 
#   # log.diversity.NRSF.plot + 
#   # log.persistence.RSF.plot + 
#   # log.persistence.NRSF.plot + 
#   plot_layout(ncol = 2)
# dev.off()

#########################
# models

#####
# facilitation
faci.gam <- gam(facilitation ~ RSF + s(NRSF,bs = "ad"),data = interaction.spread,method = "REML")

# par(mfrow = c(2,2))
# gam.check(faci.gam)

# par(mfrow = c(1,1))
# plot(faci.gam,scheme = 1,residuals = TRUE,select = 1,xlab = "non-resource factor",ylab = "facilitation intensity")

#####
# competition
comp.gam <- gam(competition ~ s(RSF) + s(NRSF,bs = "ad") + s(RSF,NRSF),data = interaction.spread,method = "REML")

# par(mfrow = c(2,2))
# gam.check(comp.gam)

# par(mfrow = c(1,1))
# plot(comp.gam,scheme = 1,residuals = TRUE,select = 1,xlab = "non-resource factor",ylab = "competition intensity")
# plot(comp.gam,scheme = 1,residuals = TRUE,select = 2,xlab = "non-resource factor",ylab = "competition intensity")

#####
# persistence
persistence.gam <- gam(mean.persistence ~ RSF + s(NRSF,bs = "ad"),data = average.persistence,method = "REML")

# par(mfrow = c(2,2))
# gam.check(persistence.gam)

# par(mfrow = c(1,1))
# plot(persistence.gam,scheme = 1,residuals = TRUE,select = 1,xlab = "non-resource factor",ylab = "species persistence")

#####
# diversity
diversity.gam <- gam(average.hill ~ RSF + s(NRSF,bs = "ad"),data = average.diversity,method = "REML")

# par(mfrow = c(2,2))
# gam.check(diversity.gam)

# par(mfrow = c(1,1))
# plot(diversity.gam,scheme = 1,residuals = TRUE,select = 1,xlab = "non-resource factor",ylab = "species persistence")
