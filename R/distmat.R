#' Calculate a pairwise landscape distance matrix (Euclidian or cost-distance).  
#' 
#' @param method Specify the type of distance matrix to be produced, using "ed" for Euclidean distance and "cd" for cost-weighted (i.e. effective) distance. Accurate distance calculations require a projected coordinate system (e.g. UTM), so do not use geographical coordinates (e.g. longlat). If you calculate cost-weighted distances, make sure the \code{points} projection is the same as the \code{landscape} raster. 
#' @param file_name (optional) A character string that will be appended to the beginning of the output filename. If no name is specified, no file will be written to the working directory.
#' @param sp_points An object of class SpatialPoints (see sp package for details).
#' @param landscape An object of class RasterLayer (see raster package for details)
#' 
#' @return An NxN (N= sample size: i.e. nrow(xy)) matrix of pairwise Euclidean and/or effective landscape distances written to .csv comma delimited files with edmat or cdmat appended to the end of the \code{filename}.
#' 
#' @examples 
#' library(sGD)
#' library(raster)
#' library(sp)
#' 
#' # read in locations and landscape data 
#' xy_file <- system.file("extdata","sGD_demo_xy.csv",package="sGD")
#' landscape_ascii <- system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD") 
#' 
#' # convert locations to SpatialPoints (sp package)
#' proj <- "+proj=utm +zone=10 +datum=NAD83"   
#' sp_points <- SpatialPoints(read.csv(xy_file)[,c(2,3)],proj4string=CRS(proj)) 
#' 
#' # convert landscape_ascii to raster object (raster package)
#' landscape <- raster(landscape_ascii,crs=CRS(proj))
#' 
#' # specify output file_name
#' file_name <- "sGD_demo" 
#' 
#' # run ed and cd matrix calculations
#' ed <- distmat(sp_points,method="ed",file_name = file_name)
#' cd <- distmat(sp_points,method="cd",file_name = file_name,landscape=landscape)
#' @importFrom gdistance transition
#' @importFrom gdistance geoCorrection
#' @importFrom gdistance costDistance
#' @importFrom ecodist full
#' @importFrom sp SpatialPoints
#' @importFrom raster crs 
#' @export
distmat <- function(sp_points,method,file_name=NULL,landscape=NULL)
{
  # check if the points are projected
  if (is.na(crs(sp_points)) == TRUE)
  {
    print("Warning: the input points have no projection defined.")
  }  
  
  for (type in method)
  {
    if (type=="ed")
    {
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      ed = round(full(dist(sp_points@coords,method="euclidean")),2)
      
      # write matrices to csv files
      if(is.null(file_name)==F)
      {
        write.table(ed,paste0(file_name,"_edmat.csv",sep=""),row.names=F,col.names=F,sep=",")         
      }
    }
    
    if (type=="cd")
    {
      # check to see if sp_points and landscape are in the same projection
      if (is.na(crs(landscape)) == TRUE)
      {
        print("Warning: the input raster has no projection defined.")
      } else if(as.character(sp_points@proj4string) != as.character(landscape@crs))
      {
        print("Warning: the projection of the points and landscape is not the same.")
      }
      
      # calculate transition surface, and geocorrect it in E-W dimension
      tr <- transition(landscape,transitionFunction = function(x) {1/mean(x)},directions=8) 
      trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
      
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      cd <- round(full(costDistance(trCorrC, sp_points)),2)
      
      # write matrix to csv files
      if(is.null(file_name)==F)
      {
        write.table(cd,paste(file_name,"_cdmat.csv",sep=""),row.names=F,col.names=F,sep=",")         
      }
    }
  }
  
  if(method=="ed") return(ed)
  if(method=="cd") return(cd)
}