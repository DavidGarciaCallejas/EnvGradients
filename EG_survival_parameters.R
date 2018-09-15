# survival parameterization for all the species

### parameter space for survival probability:

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

write.table(x = survival.param,file = paste("./results/survival_parameters_",num.sp,"sp.csv",sep=""),append = F,sep = ";",dec = ".",row.names = FALSE)

