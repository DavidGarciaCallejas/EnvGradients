
# check the behaviour of the logistic function for modelling survival and growth
# also plot it as in maestre 2009
##################

# auxiliary function for mapping to 0-1 range
LinMap <- function (x, from, to) {
  x <- x - min(x)
  x <- x/max(x)
  x <- x * (to - from)
  x + from
}

# this is the function
PlantLogistic <- function(x,p0,r1,r2,c1,k){
  (k * p0 * exp(r1*x)) / (k + p0 * (exp(r1*x)-1)) - c1 * (exp(r2*x)-1)
}

#########################
#########################
# survival (Fig 1a of Maestre2009)

p0.survival <- 1

r1.survival.min <- 0.1 #0.3 
r1.survival.max <- 0.01 #0.05

r2.survival.min <- 1.0 #1.5
r2.survival.max <- 0.85 #0.8

c.survival.min <- 3e-5
c.survival.max <- 3e-4 #5e-4

k.survival.min <- 1 #1
k.survival.max <- 0.1 #0.1

# number of species
num.sp <- 20

# create dataframe
survival.param <- data.frame(species = 1:num.sp,
                             k.survival = seq(k.survival.min,k.survival.max,length.out = num.sp),
                             p0.survival = p0.survival,
                             r1.survival = seq(r1.survival.min,r1.survival.max,length.out = num.sp),
                             c1.survival = seq(c.survival.min,c.survival.max,length.out = num.sp),
                             r2.survival = seq(r2.survival.min,r2.survival.max,length.out = num.sp))

outcome.data <- NULL
for(i.sp in 1:num.sp){
  for(x.value in seq(0,10,0.1)){
    result <- c(i.sp,
                survival.param$p0.survival[i.sp],
                survival.param$r1.survival[i.sp],
                survival.param$r2.survival[i.sp],
                survival.param$c1.survival[i.sp],
                survival.param$k.survival[i.sp],
                x.value,
                PlantLogistic(x = x.value,
                              p0 = survival.param$p0.survival[i.sp],
                              r1 = survival.param$r1.survival[i.sp],
                              r2 = survival.param$r2.survival[i.sp],
                              c1 = survival.param$c1.survival[i.sp],
                              k = survival.param$k.survival[i.sp]))
    outcome.data <- rbind(outcome.data,result)
  }# for x
}# for i.sp

outcome.data <- as.data.frame(outcome.data,row.names = F)
names(outcome.data) <- c("species","p0","r1","r2","c1","k","x","value")

outcome.data$species <- as.factor(outcome.data$species)
# outcome.data$rescaled.x <- LinMap(outcome.data$x,0,1)

my.survival.plot <- ggplot(outcome.data)
my.survival.plot <- my.survival.plot + geom_line(aes(x = x, y = value, color = species),size = 1.2)
# my.survival.plot <- my.survival.plot + geom_line(aes(x = x, y = value, group = species), color = "grey")
#
my.survival.plot <- my.survival.plot + theme_Publication()
my.survival.plot <- my.survival.plot + theme(legend.position = "none") + ylim(0,1)
my.survival.plot <- my.survival.plot + xlab("Non-resource stress") + ylab("P(survival)")
my.survival.plot <- my.survival.plot + scale_color_viridis(discrete = TRUE)
# my.survival.plot

#########################
#########################
# growth (Fig 1b of Maestre2009)
# 
p0.growth.min <- 0.7
p0.growth.max <- 1

r1.growth.min <- 0.1 #0.2 #1.5
r1.growth.max <- 0.05 #0.1 #0.2

r2.growth.min <- 0.5 #0.7
r2.growth.max <- 0.5 #0.7

c.growth.min <- 0.005
c.growth.max <- 0.05 #0.1

k.growth.min <- 0.7
k.growth.max <- 0.1

# number of species
num.sp <- 20

# create dataframe
growth.param <- data.frame(species = 1:num.sp,
                             k.growth = seq(k.growth.min,k.growth.max,length.out = num.sp),
                             p0.growth = seq(p0.growth.min,p0.growth.max,length.out = num.sp),
                             r1.growth = seq(r1.growth.min,r1.growth.max,length.out = num.sp),
                             c1.growth = seq(c.growth.min,c.growth.max,length.out = num.sp),
                             r2.growth = seq(r2.growth.min,r2.growth.max,length.out = num.sp))

growth.data <- NULL
for(i.sp in 1:num.sp){
  for(x.value in seq(0,10,0.1)){
    result <- c(i.sp,
                growth.param$p0.growth[i.sp],
                growth.param$r1.growth[i.sp],
                growth.param$r2.growth[i.sp],
                growth.param$c1.growth[i.sp],
                growth.param$k.growth[i.sp],
                x.value,
                PlantLogistic(x = x.value,
                              p0 = growth.param$p0.growth[i.sp],
                              r1 = growth.param$r1.growth[i.sp],
                              r2 = growth.param$r2.growth[i.sp],
                              c1 = growth.param$c1.growth[i.sp],
                              k = growth.param$k.growth[i.sp]))
    growth.data <- rbind(growth.data,result)
  }# for x
}# for i.sp

growth.data <- as.data.frame(growth.data,row.names = F)
names(growth.data) <- c("species","p0","r1","r2","c1","k","x","value")

growth.data$species <- as.factor(growth.data$species)
# growth.data$rescaled.x <- LinMap(growth.data$x,0,1)

my.growth.plot <- ggplot(growth.data)
my.growth.plot <- my.growth.plot + geom_line(aes(x = x, y = value, color = species),size = 1.2)
my.growth.plot <- my.growth.plot + theme_Publication()
my.growth.plot <- my.growth.plot + theme(legend.position = "none") + ylim(0,1)
my.growth.plot <- my.growth.plot + xlab("Non-resource stress") + ylab("P(growth)")
my.growth.plot <- my.growth.plot + scale_color_viridis(discrete = TRUE)
# my.growth.plot
# 
# #########################
# # tune the plot
# 
# tiff(paste("/home/david/CREAF/FPU/plant_interactions/Results/images/Logistic_function_v4.tiff",sep=""), res=600, compression = "lzw", width = 5000, height = 2500, units = "px")
# multiplot(my.survival.plot,my.growth.plot,cols = 2)
# dev.off()
 


# test the amelioration on the survival values
# 
# # current
# amelioration <- 0.5
# 
# k.survival = 0.91818 
# p0.survival = 1 
# r1.survival = 1.3818 
# c1.survival = 0.009136 
# r2.survival = 1.436
# 
# survival.prob <- PlantLogistic(x = seq(0.1,10,0.1),
#                                p0 = p0.survival,
#                                r1 = r1.survival, 
#                                r2 = r2.survival, 
#                                c1 = c1.survival,
#                                k = k.survival)
# 
# 
# net.amelioration <- amelioration * (1 - survival.prob)
# net.prob <- survival.prob + net.amelioration
# 
# png("/home/david/CREAF/FPU/plant_interactions/Results/images/facilitation.png")
# plot(survival.prob,ylim = c(0,1))
# # points(net.amelioration, col = "darkgreen")
# points(net.prob,col = "red")
# dev.off()
# 
# # old
# amelioration <- 0.3
# 
# survival.prob <- PlantLogistic.old(x = seq(0,10,0.1),a = 1,b = 0.0005, c = 0.3)
# plot(survival.prob)
# 
# net.amelioration <- amelioration * (1 - survival.prob)
# points(net.amelioration)
# 
# net.prob <- survival.prob + net.amelioration
# points(net.prob,col = "red")
# 
