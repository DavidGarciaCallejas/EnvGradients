##################
# generate environmental 2-d gradient
##################

# grid's number of rows
# grid is square

row.number <- 100

# generate grid

grid <- matrix(data = 1:row.number^2,nrow = row.number,ncol = row.number,byrow = T)

# set grid values

# columns: non-resource stress factor (NRSF), quantified in stress level
# rows: resource stress factor (RSF), quantified in carrying capacity

min.NRSF <- 0 # minimum stress level
max.NRSF <- 10 # maximum stress level
NRSF.gradient <- seq(min.NRSF,max.NRSF,length.out = row.number)

min.RSF <- 10 # minimum carrying capacity
max.RSF <- 100 # maximum carrying capacity
RSF.gradient <- seq(max.RSF,min.RSF,length.out = row.number)

# the environmental values are stored in a separate data frame, thus allowing an indeterminate
# number of factors. In principle, I will use the two already defined

grid.values <- data.frame(grid.ID = 1:row.number^2,
                          X = sort(rep(1:row.number,row.number)),
                          Y = rep(row.number:1,row.number),
                          expand.grid(RSF=RSF.gradient,NRSF=NRSF.gradient))

#####################################################
#####################################################
# library(ggplot2)
# grid.plot <- ggplot(grid.values)
# grid.plot <- grid.plot + geom_point(aes(x = X,y = Y,color = NRSF, size = RSF))
# grid.plot <- grid.plot + scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = NULL)
# grid.plot <- grid.plot + DGC::theme_Publication() + xlab("") + ylab("")
# grid.plot
# 
# #####################################################
# #####################################################
# # save the plot
# 
# tiff(paste("/home/david/CREAF/FPU/plant_interactions/Results/images/Fig_1.tiff",sep=""), res=600, compression = "lzw", width = 5000, height = 5000, units = "px")
# grid.plot
# dev.off()

write.table(x = grid.values,file = paste("./results/grid_",row.number,"_2factors.csv",sep=""),append = FALSE,sep = ";",dec = ".",row.names = FALSE)

