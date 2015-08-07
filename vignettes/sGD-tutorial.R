## ----eval=F,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
#  install.packages("devtools")

## ----eval=F,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
#  library("devtools")
#  install_github("Andrew-Shirk/sGD")
#  library("sGD")

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
library(maptools)
library(raster)

# read in location data and landscape
xy_data <- read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
proj <- "+proj=utm +zone=10 +datum=NAD83"   
xy_points <- SpatialPoints(xy_data[,c(2,3)],proj4string=CRS(proj))
IBD_landscape <- raster(nrows=256, ncols=256, crs=proj, xmn = 0, ymn = 0, xmx = 25600, ymx = 25600, vals=1)
plot(IBD_landscape)
plot(xy_points,add=T,pch=20,col="black")

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
pt1 <- xy_points[255]
pt2 <- xy_points[324]
pt1_NH <- buffer(pt1,6000)
pt2_NH <- buffer(pt2,6000)
pt1_NH <- spChFIDs(pt1_NH,"pt1") 
pt2_NH <- spChFIDs(pt2_NH,"pt2") 
pts_NH <- spRbind(pt1_NH,pt2_NH)
plot(IBD_landscape)
plot(pts_NH,add=T)

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
NH1_N <- nrow(crop(xy_points,pt1_NH)@coords)
NH2_N <- nrow(crop(xy_points,pt2_NH)@coords)
NH1_N
NH2_N

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
pt3 <- xy_points[12]
pt3_NH <- buffer(pt3,6000)
plot(IBD_landscape)
plot(xy_points,add=T,pch=20,col="black")
plot(pt3_NH,add=T)
NH3_N <- nrow(crop(xy_points,pt3_NH)@coords)
NH3_N

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
library(sGD)
IBD_dist <- distmat(xy_points,method="ed")
IBD_dist [1:10,1:10]

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
IBD_genepop_file <- system.file("extdata","sGD_demo_IBD.gen",package="sGD") 
library(adegenet)
IBD_genind_obj = read.genepop(IBD_genepop_file,ncode = 3L,quiet=T)
radius <- 6000 # for this demo, units are in meters 6000
IBD_GD <- sGD(IBD_genind_obj,xy_data,IBD_dist,radius,min_N=20,GD_ans=TRUE,NS_ans=F)

library(ggplot2)

ggplot()  +
geom_point(data=IBD_GD, aes(x=X, y=Y,color=N),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
ggplot()  +
geom_point(data=IBD_GD, aes(x=X, y=Y,color=Ap),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
NeEstimator_dir <- "C:/NeEstimator"
IBD_NS <- sGD(IBD_genind_obj,xy_data,IBD_dist,radius,min_N=20,GD_ans=FALSE,NS_ans=TRUE,NeEstimator_dir=NeEstimator_dir)

ggplot()  +
geom_point(data=IBD_NS, aes(x=X, y=Y,color=NS_ex2pct),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
IBR_landscape_ascii <- system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD") 
IBR_landscape <- raster(IBR_landscape_ascii,crs=CRS(proj))
IBR_landscape.p <- rasterToPoints(IBR_landscape)
IBR_landscape.df <- data.frame(IBR_landscape.p)
colnames(IBR_landscape.df) <- c("X", "Y", "Resistance")

ggplot()  +
geom_raster(data=IBR_landscape.df,aes(x=X,y=Y,fill=Resistance),alpha=I(0.5)) +
  scale_fill_gradient(low="black", high="lightgrey") +
  theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) 

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 5, fig.height = 5----
IBR_dist <- distmat(xy_points,method="cd",landscape=IBR_landscape)
IBR_dist [1:10,1:10]

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
# calculate transition surface and do geocorrection
library(gdistance)
tr <- transition(IBR_landscape,transitionFunction = function(x) {1/mean(x)},directions=8) 
trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
cd_pt1 <- accCost(trCorrC, pt1)
plot(cd_pt1)
plot(pt1,col="red",pch=20,add=T)

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
reclass_mat <- matrix(c(0,16000,1,16000,100000,NA),byrow=T,nrow=2)
cd_pt1_bin <- reclassify(cd_pt1,reclass_mat)
plot(cd_pt1_bin,col="black")
plot(pt1,col="red",pch=20,add=T)

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
IBR_genepop_file <- system.file("extdata","sGD_demo_IBR.gen",package="sGD") 
IBR_genind_obj <- read.genepop(IBR_genepop_file,ncode = 3L,quiet=T)
radius <- 16000 # for this demo, units are in meters 6000
IBR_sGD = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,GD_ans=TRUE,NS_ans=TRUE,NeEstimator_dir=NeEstimator_dir)

ggplot()  +
geom_point(data=IBR_sGD, aes(x=X, y=Y,color=Ap),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
boxplot(IBR_sGD$NS_ex2pct)

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
# get rid of an outlier
IBR_sGD$NS_ex2pct[which(IBR_sGD$NS_ex2pct>1000)] = NA

ggplot()  +
geom_point(data=IBR_sGD, aes(x=X, y=Y,color=NS_ex2pct),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
hist(IBR_sGD$FIS)
median(na.omit(IBR_sGD$FIS))

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
radius <- 24000
IBR_sGD_24k = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,GD_ans=TRUE,NS_ans=F)
hist(IBR_sGD_24k$FIS)
median(na.omit(IBR_sGD_24k$FIS))

## ----eval=T,message=FALSE,warning=FALSE,fig.width = 6, fig.height = 6----
radius <- 10000
IBR_sGD_10k = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,GD_ans=TRUE,NS_ans=F)
hist(IBR_sGD_10k$FIS)
median(na.omit(IBR_sGD_10k$FIS))

