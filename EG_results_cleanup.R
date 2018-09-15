
# clean up the results from the model

########################
# auxiliary
LinMap <- function(x, from, to) {
  # Shifting the vector so that min(x) == 0
  x <- x - min(x)
  # Scaling to the range of [0, 1]
  x <- x / max(x)
  # Scaling to the needed amplitude
  x <- x * (to - from)
  # Shifting to the needed level
  x + from
}

############################
hill.diversity <- function(abundances,q = 1){
  abundances <- abundances[abundances != 0]
  abundances <- abundances/sum(abundances)
  R <- length(abundances)
  # hill diversity is not defined for q = 1,
  # but its limit exists and equals
  # the exponential of shannon entropy 
  # (Jost 2006,2007,Tuomisto 2012)
  if(q == 1){
    D <- exp(-sum(abundances*log(abundances)))
  }else{
    mean.p <- (sum(abundances*abundances^(q-1)))^(1/(q-1))
    D <- 1/mean.p
  }
  return(D)
}

############################

grid.rows <- 50
timesteps <- 500
transient.phase <- 100

# nevermind this flag
saturated <- TRUE

facilitation <- TRUE

my.file.name <- paste(grid.rows,"rows_",timesteps,"t",sep="")
if(!facilitation){
  my.file.name <- paste("NO_FACILITATION_",my.file.name,sep="")
}

spatial.grid <- readr::read_delim(file = paste("./results/spatial_grid_",my.file.name,".csv",sep=""),delim = ";")

if(saturated){
  interaction.records <- readr::read_delim(file = paste("./results/interactions_",my.file.name,"_final.csv",sep=""),delim = ";")
  population.dynamics <- readr::read_delim(file = paste("./results/population_dynamics_",my.file.name,".csv",sep=""),delim = ";")
  full.population <- readr::read_delim(file = paste("./results/population_",my.file.name,".csv",sep=""),delim = ";")
}else{
  interaction.records <- readr::read_delim(file = paste("./results/interactions_unsaturated_",my.file.name,"_final.csv",sep=""),delim = ";")
  population.dynamics <- readr::read_delim(file = paste("./results/population_dynamics_unsaturated_",my.file.name,"_final.csv",sep=""),delim = ";")
}

#################
full.population <- full.population[complete.cases(full.population),]

interaction.records <- subset(interaction.records,timestep >= transient.phase)
population.dynamics <- subset(population.dynamics,timestep >= transient.phase)
full.population <- subset(full.population,timestep >= transient.phase)

#################

NRSF.division <- unique(spatial.grid$NRSF)
RSF.division <- sort(unique(spatial.grid$RSF))

interaction.records <- subset(interaction.records, interaction != "null")
interaction.records$RSF <- spatial.grid$RSF[match(interaction.records$site,spatial.grid$grid.ID)]
interaction.records$NRSF <- spatial.grid$NRSF[match(interaction.records$site,spatial.grid$grid.ID)]

population.dynamics$RSF <- spatial.grid$RSF[match(population.dynamics$site,spatial.grid$grid.ID)]
population.dynamics$NRSF <- spatial.grid$NRSF[match(population.dynamics$site,spatial.grid$grid.ID)]

full.population$RSF <- spatial.grid$RSF[match(full.population$site,spatial.grid$grid.ID)]
full.population$NRSF <- spatial.grid$NRSF[match(full.population$site,spatial.grid$grid.ID)]

############################

# # test how many cells in each division
num.cells <- spatial.grid
# num.cells$NRSF <- cut(num.cells$NRSF,breaks = NRSF.division)#,right = FALSE)
# num.cells$RSF <- cut(num.cells$RSF,breaks = RSF.division)#,right = FALSE)
num.cells <- num.cells %>% group_by(RSF,NRSF) %>% summarize(cell.number = n())

############################
# data transformations

net.biomass <- full.population %>% group_by(timestep,RSF,NRSF) %>% summarize(site.biomass = sum(size))

interaction.biomass.data <- interaction.records %>% 
  group_by(timestep,interaction,RSF,NRSF) %>% 
  summarize(interaction.number = n(),interaction.biomass = sum(biomass.involved))

interaction.biomass.data <- left_join(interaction.biomass.data,net.biomass)
interaction.biomass.data$biomass.ratio <- interaction.biomass.data$interaction.biomass/interaction.biomass.data$site.biomass
interaction.biomass.data$interaction[interaction.biomass.data$interaction == "impeded.growth" | interaction.biomass.data$interaction == "impeded.colonization"] <- "competition"

average.interaction.biomass <- interaction.biomass.data %>% group_by(interaction,RSF,NRSF) %>% summarize(average.biomass = mean(biomass.ratio))
average.interaction.biomass <- average.interaction.biomass[complete.cases(average.interaction.biomass),]

# hack to include non-present values

average.interaction.biomass$RSF <- as.factor(average.interaction.biomass$RSF)
average.interaction.biomass$NRSF <- as.factor(average.interaction.biomass$NRSF)

average.interaction.biomass <- complete(average.interaction.biomass,RSF,NRSF,interaction)
average.interaction.biomass$average.biomass[is.na(average.interaction.biomass$average.biomass)] <- 0

average.interaction.biomass$RSF <- as.numeric(as.character(average.interaction.biomass$RSF))
average.interaction.biomass$NRSF <- as.numeric(as.character(average.interaction.biomass$NRSF))

average.interaction.biomass <- unique(average.interaction.biomass)

# persistence times of every species in every cell
species.persistence <- full.population %>% group_by(site,species,ind.ID) %>% summarise(min.timestep = min(timestep), max.timestep = max(timestep))
species.persistence$persistence.timesteps <- species.persistence$max.timestep - species.persistence$min.timestep

species.persistence <- subset(species.persistence,site != 0)

species.persistence$RSF <- spatial.grid$RSF[match(species.persistence$site,spatial.grid$grid.ID)]
species.persistence$NRSF <- spatial.grid$NRSF[match(species.persistence$site,spatial.grid$grid.ID)]
species.persistence$species <- as.factor(species.persistence$species)

average.persistence <- species.persistence %>% group_by(RSF,NRSF) %>% summarize(mean.persistence = mean(persistence.timesteps))

# richness
hill.diversity.data <- full.population %>% group_by(RSF,NRSF,species,timestep) %>% summarize(sp.biomass = sum(size))
hill.diversity.data <- hill.diversity.data %>% group_by(RSF,NRSF,timestep) %>% summarize(hill.diversity = hill.diversity(sp.biomass))

hill.avg.data <- hill.diversity.data %>% group_by(NRSF,RSF) %>% summarize(average.hill = mean(hill.diversity))

##################################
# scale to 0-1
average.interaction.biomass <- average.interaction.biomass[which(!is.na(average.interaction.biomass$average.biomass)),]

average.persistence$norm.persistence <- LinMap(average.persistence$mean.persistence,0,1)

hill.avg.data$norm.diversity <- LinMap(hill.avg.data$average.hill,0,1)
hill.avg.data$log.diversity <- log(hill.avg.data$average.hill,2)
hill.avg.data$norm.log.diversity <- LinMap(hill.avg.data$log.diversity,0,1)

##################################
# write the cleaned up data
readr::write_delim(average.interaction.biomass,path = paste("./results/average_interaction_biomass_",my.file.name,".csv",sep=""),delim = ";")
readr::write_delim(average.persistence,path = paste("./results/average_persistence_",my.file.name,".csv",sep=""),delim = ";")
readr::write_delim(hill.avg.data,path = paste("./results/average_diversity_",my.file.name,".csv",sep=""),delim = ";")
