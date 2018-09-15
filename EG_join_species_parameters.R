
# join survival and growth parameterizations, and add amelioration capacity

# how many "levels" of amelioration do I want to model?
# specify them in ratios of facilitation
# e.g. 0.1 is a ten percent of facilitation relative to the basal survival probability
amelioration.levels <- c(0,0.1,0.3)

num.sp <- 20

survival.param <- read.table(file = paste("./results/survival_parameters_",num.sp,"sp.csv",sep=""),header = TRUE,sep = ";",dec = ".")
growth.param <- read.table(file = paste("./results/growth_parameters_",num.sp,"sp.csv",sep=""),header = TRUE,sep = ";",dec = ".")

# final number of species, each parameterization is repeated for each amelioration level
num.sp.final <- length(amelioration.levels)*nrow(survival.param)

# amelioration is defined according to the baseline survival probability of the facilitated species,
# and the amelioration level of the benefactor
# lower survival probabilities are more facilitated (1 - new.prob)
Amelioration <- function(baseline.prob = 0,amelioration.level = 0){
  new.prob <- baseline.prob
  if(new.prob < 1){
    new.prob <- new.prob + ((1-new.prob)*amelioration.level)
  }
  ifelse(new.prob>1,1,new.prob)
}
# Amelioration(0.1,0.1)
# Amelioration(0.5,0.1)
# Amelioration(0.7,0.1)
# Amelioration(0.1,0.3)
# Amelioration(0.5,0.3)
# Amelioration(0.7,0.3)

species.parameters <- data.frame(species = 1:num.sp.final,
                                 k.survival = rep(survival.param$k.survival,length(amelioration.levels)),
                                 p0.survival = rep(survival.param$p0.survival,length(amelioration.levels)),
                                 r1.survival = rep(survival.param$r1.survival,length(amelioration.levels)),
                                 r2.survival = rep(survival.param$r2.survival,length(amelioration.levels)),
                                 c1.survival = rep(survival.param$c1.survival,length(amelioration.levels)),
                                 k.growth = rep(growth.param$k.growth,length(amelioration.levels)),
                                 p0.growth = rep(growth.param$p0.growth,length(amelioration.levels)),
                                 r1.growth = rep(growth.param$r1.growth,length(amelioration.levels)),
                                 r2.growth = rep(growth.param$r2.growth,length(amelioration.levels)),
                                 c1.growth = rep(growth.param$c1.growth,length(amelioration.levels)),
                                 amelioration.capacity = sort(rep(amelioration.levels,nrow(survival.param))))

write.table(x = species.parameters,file = paste("./results/species_parameters_",num.sp,"sp.csv",sep=""),append = F,sep = ";",dec = ".")

##############################
##############################
# test each species' behaviour
PlantLogistic <- function(x,p0,r1,r2,c1,k){
  (k * p0 * exp(r1*x)) / (k + p0 * (exp(r1*x)-1)) - c1 * (exp(r2*x)-1)
}
###########################
###########################


outcome.data <- NULL
for(i.sp in 1:num.sp){
  for(x.value in seq(0,10,0.1)){
    result <- c(i.sp,
                species.parameters$p0.survival[i.sp],
                species.parameters$r1.survival[i.sp],
                species.parameters$r2.survival[i.sp],
                species.parameters$c1.survival[i.sp],
                species.parameters$k.survival[i.sp],
                species.parameters$p0.growth[i.sp],
                species.parameters$r1.growth[i.sp],
                species.parameters$r2.growth[i.sp],
                species.parameters$c1.growth[i.sp],
                species.parameters$k.growth[i.sp],
                x.value,
                PlantLogistic(x = x.value,
                              p0 = species.parameters$p0.survival[i.sp],
                              r1 = species.parameters$r1.survival[i.sp],
                              r2 = species.parameters$r2.survival[i.sp],
                              c1 = species.parameters$c1.survival[i.sp],
                              k = species.parameters$k.survival[i.sp]),
                PlantLogistic(x = x.value,
                              p0 = species.parameters$p0.growth[i.sp],
                              r1 = species.parameters$r1.growth[i.sp],
                              r2 = species.parameters$r2.growth[i.sp],
                              c1 = species.parameters$c1.growth[i.sp],
                              k = species.parameters$k.growth[i.sp])
                )
    outcome.data <- rbind(outcome.data,result)
  }# for x
}# for i.sp

outcome.data <- as.data.frame(outcome.data,row.names = F)
names(outcome.data) <- c("species",
                         "p0.survival","r1.survival","r2.survival","c1.survival","k.survival",
                         "p0.growth","r1.growth","r2.growth","c1.growth","k.growth",
                         "x","survival.value","growth.value")

outcome.data$species <- as.factor(outcome.data$species)

plot.data <- droplevels(subset(outcome.data,species %in% c(1,5,10,20)))

my.survival.plot <- ggplot(plot.data)
my.survival.plot <- my.survival.plot + geom_line(aes(x = x, y = survival.value, color = species))
my.survival.plot <- my.survival.plot + ylim(0,1)
# my.survival.plot <- my.survival.plot + facet_grid(b~.)

my.growth.plot <- ggplot(plot.data)
my.growth.plot <- my.growth.plot + geom_line(aes(x = x, y = growth.value, color = species))
my.growth.plot <- my.growth.plot + ylim(0,1)
# my.growth.plot <- my.growth.plot + facet_grid(b~.)
# DGC::multiplot(my.survival.plot,my.growth.plot,cols = 2)

