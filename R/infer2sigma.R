#' Infer the radius of Wright's genetic neighborhood from codominant marker genotypes. The correct radius is equal to 2 sigma, where sigma is the mean parent-offspring dispersal distance.
#' 
#' @param genind_obj A genind object (created by the adegenet package function import2genind or related methods) containing individual genotypes. The order of the individuals must be the same as the order in the xy and dist.mat inputs below.
#' @param distances A vector if distances at which to evaluate the evidence for the correct neighborhood radius.
#' @param xy A dataframe containing 3 columns in the following order: individual IDs, X coordinates, and Y coordinates. The order of the rows must match the order in the genind_obj and dist.mat inputs.
#' @param dist.mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective). The \code{distmat} function in the sGD package may be used to produce Euclidean and cost-weighted distance matrices. The order of the rows and columns in the matrix must match the order in the xy and genind_obj inputs.
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param max_N Optional. The maximum sample size per neighborhood for indices to be calculated. If the number of individuals in the neighborhood exceeds \code{max_N}, a sample of size \code{max_N} will be used from the neighborhood to compute the metrics and output files specified by the user. Note that if \code{max_N} is specified, and the value is too small to be representative of the neighobrhood, the results could differ significantly compared to if all individuals in the neighborhood were used.
#' @param min_distances The minimum number of distances allowed.
#' @param max_absFIS The maximum absolute FIS value allowed. If NULL, all FIS values are considered valid. 
#' @param raw_out Boolean. Do you want to write the neighborhood FIS values for each radius to a csv file?
#'
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
#' # specify distances to evaluate
#' distances = c(8000,12000,16000,20000,24000)
#' 
#' # run infer2sigma
#' est2sigma <- infer2sigma(genind_obj,xy,dist.mat,distances,min_N=20)
#' 
#' NH_radii = est2sigma$opt_radius
#' NH_radii[which(est2sigma$num_distances < 4)] = NA
#' NH_radii[which(abs(est2sigma$FIS_error) > 0.02)] = NA
#' 
#' @export
#' 
infer2sigma = function(genind_obj,xy,dist.mat,distances,min_N,max_N=NULL,min_distances=1,max_absFIS=NULL,raw_out=FALSE)
{
  sGD_df = data.frame(NH_Index=numeric(0),NH_ID=character(0),N=numeric(0),FIS=numeric(0),distance=numeric(0), stringsAsFactors = F)
  
  N = dim(genind_obj@tab)[1]
  
  for (distance in distances)
  {
    # define neighborhoods
    cat(paste("\t Evaluating FIS at distance:",distance,"\n"))
    capture.output(sGD_output <- sGD(genind_obj,xy,dist.mat,distance,min_N,max_N,metrics=c("GD")))
    sGD_output$distance = distance
    sGD_df = rbind(sGD_df,sGD_output[,c("NH_Index","NH_ID","N","FIS","distance")])
  }
  
  # plot FIS of neighborhoods at all distances
  boxplot(FIS~distance,sGD_df,notch=T,xlab="distance",ylab="FIS",main="FIS at varying neigborhood radii")
  abline(h=0,col="red")
  
  opt_dists = c()
  FIS_vec = c()
  nonNA_totals = c()
  NH_Ns = c()
  n_distances = c()

  for(i in c(1:N))
  {
    # return subset of FIS calculations for the ith neighborhood
    subout = subset(sGD_df,NH_Index == i)
    
    # if the sample size of the neighborhood was < min_N or > max_absFIS, convert FIS estimates to NAs
    subout$FIS[which(subout$N<min_N)] = NA
    subout$FIS[which(abs(subout$FIS)>max_absFIS)] = NA
    n_distance = length(which(!is.na(subout$FIS)))
    
    if(n_distance<min_distances)
    {
      subout$FIS = NA
    }

    # determine the neighborhood N and FIS at the optimal distance
    # (i.e. the distance at which FIS is closest to zero)
    i_opt_dist = subout$distance[which.min(abs(subout$FIS))]
    i_N = subout$N[which.min(abs(subout$FIS))]
    FIS = subout$FIS[which.min(abs(subout$FIS))]
    
    # track the number of distances with sufficient sample size (> min_N) per neighborhood
    if(!is.null(min_distances) && n_distance < min_distances)
    {
      i_opt_dist = NA
      FIS = NA
      i_N = NA
    }
    
    opt_dists = c(opt_dists,i_opt_dist)
    FIS_vec = c(FIS_vec,FIS)
    n_distances = c(n_distances,n_distance)
    NH_Ns = c(NH_Ns,i_N)
  }
  
  NH_stats = data.frame(NH_Index = c(1:N), NH_ID = dimnames(genind_obj@tab)[[1]], N = NH_Ns,
                        opt_radius = opt_dists, FIS = FIS_vec,
                        n_distances = n_distances, stringsAsFactors = F)
  
  # plot FIS of optimized neighborhoods
  boxplot(FIS~opt_radius,NH_stats,notch=T,xlab="distance",ylab="FIS",main="FIS at Optimized Neigborhood Radii")
  abline(h=0,col="red")
  
  if(raw_out==TRUE)
  {
    write.csv(sGD_df,"infer2sigma_raw.csv",row.names = FALSE)
  }

  return(NH_stats)
}