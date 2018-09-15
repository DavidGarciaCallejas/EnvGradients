#################
# functions for handling the data structure of horizontal communities in a grid of cells
#################

# 
#' Generate a data structure for storing individuals in a number of sites
#'
#' The structure is a nested list by cells and species
#' so that at each cell, each species can have an indeterminate number of individuals with a given size (or other traits)
#' e.g. plant.population[[cell]][[species]][[individual]](["size"],etc)

#'
#' @param NUM.CELLS number of sites
#' @param NUM.SP number of potential species
#'
#' @return nested list
#' @export
#'
GeneratePopulation <- function(NUM.CELLS,NUM.SP){
  
  population <- list()
  
  # first level: each cell has a list of species
  for(i.cell in 1:NUM.CELLS){
    population[[i.cell]] <- list(NUM.SP)
    
    # second level: each combination of species and cell is itself a list, containing the sizes (or other traits) of
    # the different individuals of that species in that cell
    
    for(i.sp in 1:NUM.SP){
      population[[i.cell]][[i.sp]] <- list()
    }# for i.sp
  }# for i.cell
  
  population
  
}

#' List all the individuals in a plant.population object
#'
#' @param plant.population object of type plant.population, as given by the function GeneratePopulation
#' @param cell either 0, a single number or a vector of integers with the desired cells
#' @param species either 0, a single number or a vector of integers with the desired species
#'
#' @return a dataframe with all the individuals in the specified cells of the specified species
#' @export
#'
ListIndividuals <- function(plant.population,cell = 0,species = 0){
  
  max.sp <- 1
  max.ind <- 1
  
  # I preallocate a big dataframe with the maximum possible
  # number of species and individuals
  # these loops are for finding max.sp and max.ind
  # and also for formatting the number of species and cells desired
  
  if(cell[1] == 0){
    num.cell <- 1:length(plant.population)
    # which is the maximum sp? for preallocating data
    for(i.cell in 1:length(num.cell)){
      if(length(plant.population[[num.cell[i.cell]]]) > max.sp){
        max.sp <- length(plant.population[[num.cell[i.cell]]])
        
        # maximum number of individuals of a given sp?
        for(i.sp in 1:length(plant.population[[num.cell[i.cell]]])){
          if(length(plant.population[[num.cell[i.cell]]][[i.sp]]) > max.ind){
            max.ind <- length(plant.population[[num.cell[i.cell]]][[i.sp]])
          }# if
        }# for
        
      }#if
    }#for
    
  }else{
    num.cell <- cell
    for(i.cell in 1:length(num.cell)){
      if(length(plant.population[[num.cell[i.cell]]]) > max.sp){
        max.sp <- length(plant.population[[num.cell[i.cell]]])
        
        # maximum number of individuals of a given sp?
        for(i.sp in 1:length(plant.population[[num.cell[i.cell]]])){
          if(length(plant.population[[num.cell[i.cell]]][[i.sp]]) > max.ind){
            max.ind <- length(plant.population[[num.cell[i.cell]]][[i.sp]])
          }# if
        }# for
        
      }#if
    }# for
  }
  df.size <- max.ind*max.sp*max(num.cell)
  return.df <- data.frame(site = integer(df.size), 
                          species = integer(df.size),
                          size = 0,
                          ind.ID = integer(df.size))
  df.index <- 1
  for(i.cell in 1:length(num.cell)){

    # format the desired species
    if(species[1] == 0){
      num.species <- 1:length(plant.population[[num.cell[i.cell]]])
    }else{
      num.species <- species
    }
    
    # convert the list to a dataframe and store it in return.df
    for(i.sp in 1:length(num.species)){
      if(length(plant.population[[num.cell[i.cell]]][[num.species[i.sp]]])>0){
        temp.df <- data.frame(t(sapply(plant.population[[num.cell[i.cell]]][[num.species[i.sp]]], `[`)))
        
        return.df[df.index:(df.index+nrow(temp.df)-1),] <- temp.df
        df.index <- df.index + nrow(temp.df)
      }# if length>0
    }# for each sp
  }# for each cell
  
  # trim and return the dataframe
  return.df <- subset(return.df,site != 0 & species != 0 & size != 0)
  return.df
}

### NOTE: NOT REALLY NUMBER OF INDIVIDUALS. AMEND *********
#' Number of individuals in each cell of a spatial grid. All species are summed up together,
#' for more detailed patterns it's better to obtain the full account of individuals 
#' with the function ListIndividuals
#'
#' @param spatial.grid dataframe with, at least, X and Y coordinates and a  grid.ID parameter
#' @param plant.population object of type plant.population, as given by the function GeneratePopulation
#' @param species.set species to account for, a numeric vector. Default behaviour is accounting for all species
#' @param cell.set set of cells to account for. Default, all cells
#' 
#' @return integer vector of length equal to the number of sites in the spatial.grid parameter, containing the number of individuals in each site
#' @export
#'
GetCommunitySize <- function(spatial.grid,plant.population, species.set = NULL, cell.set = NULL){
  
  if(!is.null(cell.set)){
    my.spatial.grid <- subset(spatial.grid, grid.ID %in% cell.set)
  }else{
    my.spatial.grid <- spatial.grid
  }
  
  community.size <- data.frame(grid.ID = my.spatial.grid$grid.ID, population.size = integer(nrow(my.spatial.grid)))
  
  for(i.ID in my.spatial.grid$grid.ID){
    if(is.null(species.set)){
      my.species <- 1:length(plant.population[[i.ID]])
    }else{
      my.species <- species.set
    }
    
    for(i.sp in my.species){
      indiv.length <- length(plant.population[[i.ID]][[i.sp]])
      if(indiv.length > 0){
        for(i.ind in 1:indiv.length){
          community.size$population.size[community.size$grid.ID == i.ID] <- community.size$population.size[community.size$grid.ID == i.ID] + plant.population[[i.ID]][[i.sp]][[i.ind]]["size"]
        }# for i.ind
      }# if indiv.length>0
    }# for i.sp
  }# for i.ID
  
  community.size
  
}

#' Plot the populations of a number of cells
#'
#' @param spatial.grid the spatial arrangement
#' @param plant.population population data structure
#' @param cell cells to be plotted. If 0, all cells
#'
#' @return ggplot object
#' @export
#'
PlotPopulation <- function(spatial.grid, plant.population, cell = 0){
  require(dplyr)
  
  if(cell == 0){
    cell <- 1:nrow(spatial.grid)
  }
  
  individuals <- ListIndividuals(plant.population,cell)
  empty.space <- data.frame(site = cell, species = 0, size = 0, ind.ID = 0)
  for(i.cell in cell){
    empty.space$size[i.cell] <- spatial.grid$RSF[spatial.grid$grid.ID == i.cell] - sum(individuals$size[individuals$site == i.cell])
  }
  empty.space <- subset(empty.space, size > 0)
  
  individuals <- rbind(individuals, empty.space)
  size.list <- individuals %>% group_by(site,species) %>% summarize(species.size = sum(size))
  size.list <- arrange(size.list,site,species)
  size.list$site <- as.factor(size.list$site)
  size.list$species <- as.factor(size.list$species)
  
  stacked.plot <- ggplot(size.list)
  stacked.plot <- stacked.plot + geom_col(aes(x = site, y = species.size, fill = species),color = "black") + coord_flip()
  # stacked.plot <- stacked.plot + geom_text(aes(x = site, y = species.size, label = species, group = species),position = position_stack(vjust = .5))
  stacked.plot <- stacked.plot + facet_wrap(~site,ncol = 20, scales = "free_y") 
  stacked.plot <- stacked.plot + theme_Publication() + labs(y = NULL, x = NULL) + scale_fill_grey(start = 1,end = 0.2) 
  stacked.plot <- stacked.plot + theme(axis.ticks = element_blank(), 
                                       axis.text = element_blank(), 
                                       legend.position = "none",
                                       strip.background = element_blank(),
                                       strip.text.y = element_blank(),
                                       strip.text.x = element_blank())
  stacked.plot
}






