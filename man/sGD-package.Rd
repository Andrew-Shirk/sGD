\name{sGD-package}
\alias{sGD-package}
\docType{package}
\title{
Spatially explicit estimation of genetic diversity indices and Wright's
neighborhood size (NS)
}
\description{
sGD provides spatially explicit estimates of genetic diversity indices
and Wright's neighborhood size in continuous populations isolated by distance or
resistance as described in Shirk & Cushman 2014. Inputs include a set of 
spatially referenced codominant marker genotypes, a matrix of pairwise
distances between all sampled individuals (either in Euclidean units or as
effective distances quantified by methods such as least-cost-path or
circuit theory), and several parameters, including the minimum neighborhood
size and the neighborhood radius (specified in Euclidean or effective
distances). sGD defines a local neighborhood extent (based on a user-specified
radius) around each sample location, and then calculates genetic diversity
indices and/or Wright's neighborhood size (a local measure of effective
population size that accounts for the continuous structure of populations
isolated by distance). Requires the NeEstimator program (>= 2.0) available at
http://molecularfisherieslaboratory.com.au/neestimator-software.
}
\details{
\tabular{ll}{
Package: \tab sGD\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2014-09-26\cr
License: \tab GPL (>= 2)\cr
}}
\author{
Andrew Shirk <ashirk@uw.edu>
}
\references{
Shirk, A.J. and Cushman, S.A. 2014. Spatially-explicit estimation of Wright's 
neighborhood size in continuous populations. Frontiers in Ecology and Evolution 2(62).

Shirk, A.J. and Cushman, S. A. 2011. sGD: software for estimating spatially explicit indices of genetic diversity. Molecular Ecology Resources, 11(5), 922-934.
}
\examples{
library(sGD)
setwd("C:/Temp") # the output file will be written to the working directory
genepop_file <- system.file("extdata","sGD_demo_genepop.gen",package="sGD") 
genind_obj <- read.genepop(genepop_file,quiet=T)
xy = read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
dist_mat <- as.matrix(read.csv(system.file("extdata","sGD_demo_cdmat.csv",package="sGD") ,header=F)) 
radius <- 16238 # for this demo, units are in meters
min_N <- 20
NeEstimator_dir <- "C:/NeEstimator"

sGD_output = sGD(genind_obj,xy,dist_mat,radius,min_N,NS_ans=TRUE,GD_ans=F,NHmat_ans=T,genout_ans=T,file_name=NULL,NeEstimator_dir=NULL)

# specify landscape and sGD output 
landscape <- raster(system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD"))
sGD_output <- read.csv(system.file("extdata","sGD_demo_output_sGD.csv",package="sGD"))
# Convert raster to dataframe for ggplot 
landscape.p <- rasterToPoints(landscape)
landscape.df <- data.frame(landscape.p)
colnames(landscape.df) <- c("X", "Y", "Resistance")

# Plot sGD output (Ap is shown here, but explore all sGD outputs) atop the resistance model
library(ggplot2)
ggplot()  +
geom_raster(data=landscape.df,aes(x=X,y=Y,fill=Resistance),alpha=I(0.5)) +
scale_fill_gradient(low="black", high="lightgrey") + 
geom_point(data=sGD_output, aes(x=X, y=Y,color=Ap),size=5) + 
scale_color_gradient(low="red", high="green",na.value = "white") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()) 

# calculate landscape distance matrices based on Euclidean of cost-weighted distance
library(sGD)
setwd("C:/Temp")
file_name <- "sGD_demo" 
xy_file <- system.file("extdata","sGD_demo_xy.csv",package="sGD")
proj <- "+proj=utm +zone=10 +datum=NAD83"   
sp_points <- SpatialPoints(read.csv(xy_file)[,c(2,3)],proj4string=CRS(proj)) 
landscape_ascii <- system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD") 
landscape <- raster(landscape_ascii,crs=CRS(proj))
distmat(sp_points,file_name,"ed")
distmat(sp_points,file_name,"cd",landscape)
distmat(sp_points,file_name,c("ed","cd"),landscape)
}