#' Infer the radius of Wright's genetic neighborhood from codominant marker genotypes. The correct radius is equal to 2 sigma, where sigma is the mean parent-offspring dispersal distance.
#' 
#' @param genind_obj A genind object (created by the adegenet package function import2genind or related methods) containing individual genotypes. The order of the individuals must be the same as the order in the xy and dist.mat inputs below.
#' @param radii A vector of neighborhood radii at which to evaluate the evidence for the correct neighborhood radius.
#' @param xy A dataframe containing 3 columns in the following order: individual IDs, X coordinates, and Y coordinates. The order of the rows must match the order in the genind_obj and dist.mat inputs.
#' @param dist.mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective). The \code{distmat} function in the sGD package may be used to produce Euclidean and cost-weighted distance matrices. The order of the rows and columns in the matrix must match the order in the xy and genind_obj inputs.
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param max_N Optional. The maximum sample size per neighborhood for indices to be calculated. If the number of individuals in the neighborhood exceeds \code{max_N}, a sample of size \code{max_N} will be used from the neighborhood to compute the metrics and output files specified by the user. Note that if \code{max_N} is specified, and the value is too small to be representative of the neighobrhood, the results could differ significantly compared to if all individuals in the neighborhood were used.
#' @return 
#'
#' @examples 
#' library(sGD)
#' library(adegenet)
#' 
#' # read in genotypes, locations, and distance matrix
#' genepop.file <- system.file("extdata","sGD_demo_IBR.gen",package="sGD") 
#' xy = read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
#' dist.mat <- as.matrix(read.csv(system.file("extdata","sGD_demo_cdmat.csv",package="sGD"),
#'                                header=FALSE))
#' 
#' # convert genepop to genind (make sure you specify the correct allele code digits - ncode)
#' genind_obj <- read.genepop(genepop.file,ncode=3L,quiet=TRUE) 
#' 
#' # specify radii to evaluate
#' radii = c(8000,12000,16000,20000,24000)
#' 
#' # run infer2sigma
#' est2sigma <- infer2sigma(genind_obj,xy,dist.mat,radii,min_N=20)
#' radii_summary = aggregate(FIS~radius,est2sigma,median)
#' 
#' @export
#' 
infer2sigma <- function(genind_obj,xy,dist.mat,radii,min_N,max_N=NULL){
  
  sGD_df <- data.frame(NH_Index=numeric(0),NH_ID=character(0),N=numeric(0),FIS=numeric(0),radius=numeric(0), stringsAsFactors = F)
  
  N <- dim(genind_obj@tab)[1]
  
  for (radius in radii){
    # define neighborhoods
    cat(paste("\t Evaluating FIS at radius:",radius,"\n"))
    capture.output(sGD_output <- sGD(genind_obj,xy,dist.mat,radius,min_N,max_N,metrics=c("GD")))
    sGD_output$radius <- radius
    sGD_df <- rbind(sGD_df,sGD_output[,c("NH_Index","NH_ID","N","FIS","radius")])
  }
  
  # plot FIS of neighborhoods at all radii
  boxplot(FIS~radius,sGD_df,notch=T,xlab="radius",ylab="FIS",main="FIS at varying neigborhood radii")
  abline(h=0,col="red")
  
  return(sGD_df)
}