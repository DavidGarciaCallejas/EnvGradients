######################
# main community model
#######################

source("./EG_population_structure.R")

###########################
# non-resource response function

PlantLogistic <- function(x,p0,r1,r2,c1,k){
  my.result <- (k * p0 * exp(r1*x)) / (k + p0 * (exp(r1*x)-1)) - c1 * (exp(r2*x)-1)
  return(ifelse(my.result < 0, 0, my.result))
}

###########################
# verbose output (print when an individual dies, etc)
VERBOSE <- FALSE
# store console output on a log file
STORE_LOG <- FALSE
# store the results from the simulation
STORE_RESULTS <- TRUE
# store temporary results
STORE_TEMP <- FALSE
# store every ... timesteps
timestep.store <- 100

# how many species
# note that species parameters file has to be already in ./results/
NUM.SP <- 20

# how many timesteps?
timesteps <- 500

# how many individuals at the beginning?
num.ind <- 10

# proportion of colonized cells at each step
colonization.rate <- 0.2

# lambda parameter of the Poisson distribution governing 
# the number of colonizers to a given cell
colonization.lambda <- 1

# do we consider stress amelioration? 
# remember: amelioration only affects survival
amelioration <- FALSE

# initial size of colonizing individuals
# can be either a single value or a (min,max) interval
# from which random draws are taken
initial.size <- c(0.95,1.05)

# load spatial grid with given number of cells

# spatial.grid <- read.table(file = "./results/grid_2_2factors.csv",header = TRUE,sep = ";",dec = ".")
# spatial.grid <- read.table(file = "./results/grid_10_2factors.csv",header = TRUE,sep = ";",dec = ".")
# spatial.grid <- read.table(file = "./results/grid_20_2factors.csv",header = TRUE,sep = ";",dec = ".")
spatial.grid <- read.table(file = "./results/grid_50_2factors.csv",header = TRUE,sep = ";",dec = ".")
# spatial.grid <- read.table(file = "./results/grid_100_2factors.csv",header = TRUE,sep = ";",dec = ".")

NUM.CELLS <- nrow(spatial.grid)
grid.rows <- max(spatial.grid$X)

# load species parameters
species.parameters <- read.table(file = paste("./results/species_parameters_",NUM.SP,"sp.csv",sep=""),header = TRUE,sep = ";",dec = ".")

# generate population data structure
local.population <- GeneratePopulation(NUM.CELLS,NUM.SP)
ind.count <- 1

################
# for storing interactions
interaction.records <- data.frame(site = rep(1,1e6),
                                  timestep = rep(1,1e6),
                                  interaction = rep("null",1e6),
                                  biomass.involved = rep(0,1e6),stringsAsFactors = FALSE)
count.int <- 1

# also store population dynamics: number of individuals,species and deaths
population.dynamics <- data.frame(site = rep(1:NUM.CELLS,timesteps),
                                  timestep = sort(rep(1:timesteps,NUM.CELLS)),
                                  individuals = 0,
                                  species = 0,
                                  deaths = 0,
                                  colonizations = 0)

# ADDITION: store the full population dynamics, not just the last timestep
# in order to get persistence and diversity estimates
full.population <- NULL

################
my.file.name <- paste(grid.rows,"rows_",timesteps,"t",sep = "")
if(amelioration==FALSE){
  my.file.name <- paste("NO_FACILITATION_",my.file.name,sep="")
}

################
# generate random population of N individuals
min.size <- 0.1
max.size <- 1
spatial.grid$community.size <- 0

for(i.ind in 1:num.ind){
  # which species
  my.sp <- sample(1:NUM.SP,1)
  # how big
  my.size <- runif(1,min.size,max.size)
  # which cell
  my.cell <- sample(1:NUM.CELLS,1)
  # is there capacity in that cell? otherwise, change size and cell
  # keep a count in order not to run forever
  count <- 0
  while(spatial.grid$RSF[my.cell] - spatial.grid$community.size[my.cell] < my.size & count < 1000){
    my.size <- runif(1,min.size,max.size)
    my.cell <- sample(1:NUM.CELLS,1)
    count <- count + 1
  }# while
  # if found, add it to the data structure
  if(count < 1000){
    local.population[[my.cell]][[my.sp]][[length(local.population[[my.cell]][[my.sp]])+1]] <- c(ID = spatial.grid$grid.ID[my.cell],
                                                                                                species = species.parameters$species[my.sp],
                                                                                                size = my.size,
                                                                                                ind.ID = ind.count)
    ind.count <- ind.count + 1
  }# if
}# for i.ind

# my.plot <- PlotPopulation(spatial.grid,local.population)

################################
################################
# start dynamics
if(STORE_LOG) sink(paste("./log_",timesteps,"t.txt",sep=""))

for(i.timestep in 1:timesteps){
  
  ################################
  ################################
  # survival
  
  for(i.cell in 1:NUM.CELLS){
    my.pop.index <- which(population.dynamics$site == i.cell & population.dynamics$timestep == i.timestep)
    
    cell.RSF <- spatial.grid$RSF[i.cell]
    cell.NRSF <- spatial.grid$NRSF[i.cell]
    
    cell.amelioration <- 0
    sp.ameliorating <- 0
    
    # check whether in this cell there are plants with amelioration capacity
    # and store the level and species
    for(i.sp in 1:length(local.population[[i.cell]])){
      
      # if there are individuals of i.sp
      if(length(local.population[[i.cell]][[i.sp]])>0){
        my.sp.amelioration <- species.parameters$amelioration.capacity[species.parameters$species == i.sp]
        if(my.sp.amelioration > cell.amelioration){
          cell.amelioration <- my.sp.amelioration
          sp.ameliorating <- i.sp
        }
      }# if individuals
    }# check each species
    
    # check every species in the target cell
    for(i.sp in 1:length(local.population[[i.cell]])){
      
      # if there are individuals of i.sp
      if(length(local.population[[i.cell]][[i.sp]])>0){
        
        # calculate the survival probability of each individual
        num.plants <- length(local.population[[i.cell]][[i.sp]])
        
        # loop backwards in case an element is removed, 
        # because the length of the object will be altered
        for(i.plant in num.plants:1){
          
          my.sp <- local.population[[i.cell]][[i.sp]][[i.plant]]["species"]
          
          # get the parameters of the ith species
          k.param <- species.parameters$k.survival[species.parameters$species == my.sp]
          p0.param <- species.parameters$p0.survival[species.parameters$species == my.sp]
          r1.param <- species.parameters$r1.survival[species.parameters$species == my.sp]
          r2.param <- species.parameters$r2.survival[species.parameters$species == my.sp]
          c1.param <- species.parameters$c1.survival[species.parameters$species == my.sp]
          
          # calculate survival probability
          my.survival <- PlantLogistic(x = cell.NRSF, 
                                       k = k.param, 
                                       p0 = p0.param, 
                                       r1 = r1.param, 
                                       r2 = r2.param, 
                                       c1 = c1.param)
          
          # is there amelioration?
          net.survival <- my.survival
          if(amelioration & cell.amelioration > 0 & my.survival < 1){
            # a species can't facilitate itself
            if(sp.ameliorating != my.sp){
              
              net.amelioration <- cell.amelioration * (1 - my.survival)
              net.survival <- my.survival + net.amelioration
              
            }# if different species
          }# if amelioration
          
          # draw a number to check if it dies, and if so, remove its record
          if(net.survival < runif(1,0,1)){
            
            if(VERBOSE){
              print(paste("T:",i.timestep,
                          "**** DEATH **** cell:",i.cell,
                          ", species:",i.sp,
                          ", size:", round(local.population[[i.cell]][[i.sp]][[i.plant]]["size"],3),
                          ", survival prob:",round(net.survival,3),sep=""))
            }# if verbose
            
            local.population[[i.cell]][[i.sp]][[i.plant]] <- NULL
            
            ################
            population.dynamics$deaths[my.pop.index] <- population.dynamics$deaths[my.pop.index] + 1
            ################
            
          }else{
            if(my.survival < 1 & my.survival != net.survival){
              
              # the biomass involved in this interaction is the size of the surviving plant, that was
              # facilitated
              # this is a realized facilitation, as opposed to the potential one above
              interaction.biomass <- local.population[[i.cell]][[i.sp]][[i.plant]]["size"]
              
              ################
              interaction.records[count.int,] <- list(i.cell,i.timestep,"facilitation",interaction.biomass)
              count.int <- count.int + 1
              ################
              
            }# if amelioration helped
          }# if-else death
          
        }# for i.plant
      }# if i.sp
    }# for i.sp
    
  }# for i.cell
  
  ################################
  ################################
  # colonization
  # I assume homogeneous dispersal, so that an open space is colonized randomly
  # by any species. No differences in dispersal traits.
  # This is because I am more interested in cell-level responses than in metacommunity spatial patterns, for now
  
  community.size <- GetCommunitySize(spatial.grid,local.population)[2]
  spatial.grid$community.size <- community.size$population.size
  
  num.cells.colonized <- nrow(spatial.grid) * colonization.rate
  
  cells.colonized <- sample(1:NUM.CELLS,num.cells.colonized,replace = F)
  
  for(i.cell in cells.colonized){
    my.pop.index <- which(population.dynamics$site == i.cell & population.dynamics$timestep == i.timestep)
    
    # size of the available space in each cell
    available.capacity <- spatial.grid$RSF[i.cell] - spatial.grid$community.size[i.cell]
    
    # how many individuals colonize this cell
    num.colonizers <- 0
    while(num.colonizers == 0){
      num.colonizers <- rpois(1,colonization.lambda)
    }
    # compare to the available capacity of the cell
    # num.colonizers <- ifelse(available.capacity > num.colonizers, num.colonizers, round(available.capacity))
    if(available.capacity < num.colonizers & available.capacity > 0){
      # num.colonizers <- round(available.capacity)
      
      # the biomass involved in this interaction is the biomass of the colonizers
      # that were not allowed to establish
      impeded.colonizers <- round(num.colonizers - available.capacity)
      if(impeded.colonizers == 0){
        impeded.colonizers <- 1
      }
      biomass.impeded <- runif(impeded.colonizers,initial.size[1],initial.size[2])
      
      ################
      interaction.records[count.int,] <- list(i.cell,i.timestep,"impeded.colonization",biomass.impeded)
      count.int <- count.int + 1
      ################
    }

        
    # add each colonizer
    if(num.colonizers > 0){
      for(i.colonizer in 1:num.colonizers){
        i.sp <- round(runif(1,1,NUM.SP))
        index <- length(local.population[[i.cell]][[i.sp]])
        
        # what is the initial size of the colonizing individual?
        if(length(initial.size)>1){
          my.size <- runif(1,initial.size[1],initial.size[2])
        }else{
          my.size <- initial.size
        }
        
        local.population[[i.cell]][[i.sp]][[index + 1]] <- c(ID = spatial.grid$grid.ID[i.cell],
                                                             species = species.parameters$species[i.sp],
                                                             size = my.size,
                                                             ind.ID = ind.count)
        ind.count <- ind.count + 1
        
        if(VERBOSE){
          print(paste("T:",i.timestep,
                      "**** COLONIZATION **** cell:",i.cell,
                      ", available capacity:",round(available.capacity,3),
                      ", filling up to:", num.colonizers,
                      # ", species:",i.sp,
                      # ", size:", local.population[[i.cell]][[i.sp]][[index+1]]["size"],
                      ", cell RSF:",round(spatial.grid$RSF[i.cell]),sep=""))
        }# if verbose
        
        ################
        population.dynamics$colonizations[my.pop.index] <- population.dynamics$colonizations[my.pop.index] + 1
        ################
        
      }# for i.colonizer
    }# if num.colonizers>0
    
  }# for i.cell
  
  # my.plot <- PlotPopulation(spatial.grid,local.population)
  
  ################################
  ################################
  # growth
  # fastest-growing species will grow more, and overall, no growth is possible when community.size equals carrying capacity
  
  community.size <- GetCommunitySize(spatial.grid,local.population)[2]
  spatial.grid$community.size <- community.size$population.size
  
  for(i.cell in 1:NUM.CELLS){
    
    cell.NRSF <- spatial.grid$NRSF[i.cell]
    
    # size of the available space in each cell
    available.capacity <- spatial.grid$RSF[i.cell] - spatial.grid$community.size[i.cell]
    individual.growth <- data.frame(species = integer(),initial.size = numeric(),growth = numeric(), completed = logical())
    
    # check every species in the target cell
    for(i.sp in 1:length(local.population[[i.cell]])){
      
      # if there are individuals of i.sp
      if(length(local.population[[i.cell]][[i.sp]])>0){
        
        # each individual
        num.plants <- length(local.population[[i.cell]][[i.sp]])
        
        for(i.plant in 1:num.plants){
          my.sp <- local.population[[i.cell]][[i.sp]][[i.plant]]["species"]
          
          # get the parameters of the ith species
          k.param <- species.parameters$k.growth[species.parameters$species == my.sp]
          p0.param <- species.parameters$p0.growth[species.parameters$species == my.sp]
          r1.param <- species.parameters$r1.growth[species.parameters$species == my.sp]
          r2.param <- species.parameters$r2.growth[species.parameters$species == my.sp]
          c1.param <- species.parameters$c1.growth[species.parameters$species == my.sp]
          
          # calculate growth 
          individual.growth[nrow(individual.growth)+1,] <- c(i.sp,
                                                             local.population[[i.cell]][[i.sp]][[i.plant]]["size"],
                                                             PlantLogistic(x = cell.NRSF, 
                                                                           k = k.param, 
                                                                           p0 = p0.param, 
                                                                           r1 = r1.param, 
                                                                           r2 = r2.param, 
                                                                           c1 = c1.param),FALSE)
        }# for i.plant
        
      }# if individuals
    }# for i.sp
    
    # discard those species that theoretically grow but the result of the logistic function is 0,
    # so that in practice do not grow
    individual.growth <- subset(individual.growth, growth > 0)
    
    if(nrow(individual.growth)>0){
      
      # of all the species present, the ones that grow more have prevalence, i.e.
      # are more competitive in occupying the available space
      
      individual.growth$relative.growth <- individual.growth$initial.size * individual.growth$growth
      individual.growth <- dplyr::arrange(individual.growth,desc(relative.growth))
      
      # remaining <- nrow(individual.growth)
      # individuals that grow
      completed <- 0
      
      while(available.capacity > 0 & completed < nrow(individual.growth)){
        
        completed <- completed + 1
        
        # search for the ith individual, that with a given size and species in i.cell
        my.ind <- FALSE
        i.ind <- 1
        
        while(my.ind == FALSE){
          if(local.population[[i.cell]][[individual.growth$species[completed]]][[i.ind]]["size"] == individual.growth$initial.size[completed]){
            local.population[[i.cell]][[individual.growth$species[completed]]][[i.ind]]["size"] <- (individual.growth$initial.size[completed] + individual.growth$growth[completed])
            available.capacity <- available.capacity - individual.growth$growth[completed]
            
            individual.growth$completed[completed] <- TRUE
            my.ind <- TRUE
            
            if(VERBOSE){
              print(paste("T:",i.timestep,
                          "**** GROWTH **** cell:",i.cell,
                          ", species:",individual.growth$species[completed],
                          ", initial size:", round(individual.growth$initial.size[completed],3),
                          ", growth:", round(individual.growth$growth[completed],3),
                          ", new size:",round(local.population[[i.cell]][[individual.growth$species[completed]]][[i.ind]]["size"],3),sep=""))
            }# if verbose
            
          }else{
            i.ind <- i.ind + 1
          }# if/else size found
        }# while individual not yet found
        
      }# while space available for growing
      
      #if not all species could grow, record competition
      if(sum(individual.growth$completed == FALSE)>0){
        impeded.sp <- subset(individual.growth, completed == FALSE)
        
        for(i.non.completed in 1:nrow(impeded.sp)){
          
          # the biomass involved in this case is the growth that was not realized
          ################
          interaction.records[count.int,] <- list(i.cell,i.timestep,"impeded.growth",impeded.sp$growth[i.non.completed])
          count.int <- count.int + 1
          ################
        }# for each non completed growth
      }# if not all species could grow
      
    }# if growth
  }# for i.cell
  
  for(i.cell in 1:NUM.CELLS){
    my.pop.index <- which(population.dynamics$site == i.cell & population.dynamics$timestep == i.timestep)
    
    # fill up population dynamics dataframe
    for(i.sp in 1:length(local.population[[i.cell]])){
      # if there are individuals of i.sp
      if(length(local.population[[i.cell]][[i.sp]])>0){
        population.dynamics$species[my.pop.index] <- population.dynamics$species[my.pop.index] + 1
        population.dynamics$individuals[my.pop.index] <- population.dynamics$individuals[my.pop.index] + length(local.population[[i.cell]][[i.sp]])
      }
    }
  }# for i.cell
  
  print(paste(date()," - t:",i.timestep,sep=""))
  
  # beware reaching the limits of interaction.records
  # store it if so and create a new dataframe
  if(count.int>= 1e6){
    readr::write_delim(x = interaction.records,path = paste("./results/TEMP_interactions_",my.file.name,"_timestep",i.timestep,".csv",sep=""),delim = ";")
    count.int <- 1
    interaction.records <- data.frame(site = rep(1,1e6),
                                      timestep = rep(1,1e6),
                                      interaction = rep("null",1e6),
                                      biomass.involved = rep(0,1e6),stringsAsFactors = FALSE)
  }
  
  
  my.population <- ListIndividuals(local.population)
  my.population$timestep <- i.timestep
  
  full.population <- bind_rows(full.population,my.population)
  
  if(STORE_TEMP & i.timestep%%timestep.store == 0){
   save.image(paste("TEMP_workspace_",my.file.name,"_timestep",i.timestep,".RData",sep="")) 
  }
  
}# for i.timestep

if(STORE_LOG){
  sink()
}

if(STORE_RESULTS){
  readr::write_delim(x = spatial.grid,path = paste("./results/spatial_grid_",my.file.name,".csv",sep=""),delim = ";")
  # save(local.population,file = paste("./plant_population_",grid.rows,"rows_",timesteps,"t",sep=""))
  readr::write_delim(x = full.population,path = paste("./results/population_",my.file.name,".csv",sep=""),delim = ";")
  
  # gather temporary interaction files
  interaction.files <- list.files(path = "./results/",pattern = paste("TEMP_interactions_",my.file.name,"_timestep",sep=""))

  # bind them
  full.interactions <- NULL
  
  if(length(interaction.files)>0){
    for(i.file in 1:length(interaction.files)){
      my.interactions <- readr::read_delim(file = paste("./results/",interaction.files[i.file],sep=""),delim = ";")
      full.interactions <- bind_rows(full.interactions,my.interactions)
    }
    
    interaction.records <- subset(interaction.records,interaction != "null")
    interaction.records$site <- as.integer(interaction.records$site)
    interaction.records$timestep <- as.integer(interaction.records$timestep)
    interaction.records$biomass.involved <- as.numeric(interaction.records$biomass.involved)
    
    full.interactions <- bind_rows(full.interactions,interaction.records)
    
    # remove them
    file.remove(paste("./results/",interaction.files,sep=""))
  }else{
    interaction.records <- subset(interaction.records,interaction != "null")
    full.interactions <- interaction.records
  }
  readr::write_delim(x = full.interactions,path = paste("./results/interactions_",my.file.name,"_final.csv",sep=""),delim = ";")
  readr::write_delim(x = population.dynamics,path = paste("./results/population_dynamics_",my.file.name,".csv",sep=""),delim = ";")
}
