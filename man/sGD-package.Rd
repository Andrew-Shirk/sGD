\name{sGD-package}
\alias{sGD-package}
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
}
~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Andrew Shirk <ashirk@uw.edu>
}
\references{
Shirk, A.J. and Cushman, S.A. 2014. Spatially-explicit estimation of Wright's 
neighborhood size in continuous populations. Frontiers in Ecology and Evolution. 
(in press)
}
\examples{
# calculate genetic diversity indices and NS within genetic neighborhoods
setwd("C:/Temp")
genepop_file = system.file("extdata","demo_IBR.gen",package="sGD")
output_name = "demo_IBR"
xy = system.file("extdata","demo_xy.csv",package="sGD")
dist_mat = system.file("extdata","demo_IBR_dist_mat.csv",package="sGD")
radius = 20000
min_N = 20
NS_ans = TRUE # calculate Wright's NS spatially
GD_ans = TRUE # calculate genetic diversity indices spatially
NeEstimator_dir = "C:/NeEstimator"

sGD(genepop_file,output_name,xy,dist_mat,radius,min_N,NS_ans,GD_ans,NeEstimator_dir)


# calculate landscape distance matrices based on Euclidean of cost-weighted distance
setwd("C:/Temp")
raster_file = system.file("extdata","demo_IBR_landscape.asc",package="sGD")
output_name = "demo_IBR"
xy = system.file("extdata","demo_XY.csv",package="sGD")
proj = "+proj=utm +zone=10 +datum=NAD83"

distmat(xy,output_name,"ed")
distmat(xy,output_name,"cd",proj,raster_file)
distmat(xy,output_name,c("ed","cd"),proj,raster_file)
}
