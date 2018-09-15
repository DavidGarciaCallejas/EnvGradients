# growth parameterization for all the species

### parameter space for growth:
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

write.table(x = growth.param,file = paste("./results/growth_parameters_",num.sp,"sp.csv",sep=""),append = F,sep = ";",dec = ".",row.names = FALSE)

