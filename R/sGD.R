#' Calculate spatially explicit indicies of genetic diversity and Wright's neighborhood size (NS).
#' 
#' @param genepop_file The path to a genotype file in genepop format (with a .gen extension).
#' @param output_name A character string that will be appended to the front of the output filename (will end with "_sGD.csv").
#' @param xy A dataframe containing 3 columns in the following order: individual IDs, X coordinates, and Y coordinates.
#' @param dist_mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective).
#' @param radius The radius of the genetic neighborhood in the same units as \code{dist_mat}.
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param NS_ans Boolean (T or F) answer to whether you want sGD to calculate Wright's neighborhood size.
#' @param GD_ans Boolean (T or F) answer to whether you want sGD to calculate genetic diversity indices. This calculation can take a long time depending on how many individuals are in your sample and the \code{radius} of the neighborhood.
#' @param NeEstimator_dir Path to the NeEstimator 2.0 directory. NeEstimator 2.0 is required only if NS_ans = T. It can be downloaded from \url{http://molecularfisherieslaboratory.com.au/neestimator-software}.
#' @param NHmat_ans Boolean (T or F) answer to whether you want sGD to write a matrix defining neighborhood membership. For each row/column of the matrix, a value of 1 occurs at the indices of all individuals in the neighborhood and a value of 0 occurs for all individuals outside the neighborhood.
#' @param genout_ans Boolean (T or F) answer to whether you want sGD to write a genepop file containing the genotypes for all neighborhoods to the workind directory. 
#' 
#' @return sGD produces a comma delimited text file containing estimates of genetic diversity and neighborhood size for neighborhoods surrounding each sample location.
#' @return Columns in the output text file include:
#' @return \code{Indiv_ID} - the ID of the individual at the neighborhood center, taken from individual's ID in the \code{xy_file}.
#' @return \code{X} - the X coordinate of the neighborhood center.
#' @return \code{Y} - the Y coordinate of the neighborhood center.
#' @return \code{N} - the number of individuals within the neighborhood.
#' @return \code{A} - the average number of alleles across all loci/individuals within the neighborhood.
#' @return \code{Ap} - the proportion of alleles from the entire population that area actually present in the neighborhood.
#' @return \code{Ar} - the allelic richness across all loci/individuals within the neighborhood.
#' @return \code{He} - the average expected heterozygosity across all loci/individuals within the neighborhood.
#' @return \code{Ho} - the average observed heterozygosity across all loci/individuals within the neighborhood.
#' @return \code{FIS} - the average inbreeding coefficient across all loci/individuals within the neighborhood.
#' @return \code{NS_ex0pct} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, not exluding rare alleles that could bias the estimate.
#' @return \code{NS_ex2pct} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.02 or less.
#' @return \code{NS_ex5pct} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.05 or less.
#' @return \code{NS_ex10pct} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.10 or less.
#' 
#' @examples 
#' #Make sure your paths are correct for your operating system (e.g. in linux, it might be "/home/yourname/Temp")
#' library(sGD)
#' setwd("C:/Temp") # the output file will be written to the working directory
#' genepop_file <- system.file("extdata","sGD_demo_genepop.gen",package="sGD") 
#' genind_obj = read.genepop(genepop_file,quiet=T)
#' output_name <- "demo"     
#' xy = read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
#' dist_mat = as.matrix(read.csv(system.file("extdata","sGD_demo_cdmat.csv",package="sGD") ,header=F)) 
#' radius <- 16238 # for this demo, units are in meters
#' min_N <- 20
#' NS_ans <- TRUE
#' GD_ans <- TRUE
#' NeEstimator_dir <- "C:/NeEstimator"
#' 
#' sGD(genind_obj,output_name,xy,dist_mat,radius,min_N,NS_ans=TRUE,GD_ans=TRUE,NeEstimator_dir,NHmat_ans=FALSE,genout_ans=FALSE)
#' 
#' # specify landscape and sGD output 
#' landscape <- raster(system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD"))
#' sGD_output <- read.csv(system.file("extdata","sGD_demo_output_sGD.csv",package="sGD"))
#'
#' # Convert raster to dataframe for ggplot 
#' landscape.p <- rasterToPoints(landscape)
#' landscape.df <- data.frame(landscape.p)
#' colnames(landscape.df) <- c("X", "Y", "Resistance")
#' 
#' # Plot sGD output (Ap is shown here, but explore all sGD outputs) atop the resistance model
#' library(ggplot2)
#'
#' ggplot()  +
#'  geom_raster(data=landscape.df,aes(x=X,y=Y,fill=Resistance),alpha=I(0.5)) +
#'  scale_fill_gradient(low="black", high="lightgrey") + 
#'  geom_point(data=sGD_output, aes(x=X, y=Y,color=Ap),size=5) + 
#'  scale_color_gradient(low="red", high="green",na.value = "white") + 
#'  theme(panel.grid.major = element_blank(),
#'        panel.grid.minor = element_blank(),
#'        panel.background = element_blank()) 


sGD <- function(genind_obj,output_name,xy,dist_mat,radius,min_N,NS_ans,GD_ans,NeEstimator_dir=NULL,NHmat_ans=FALSE,genout_ans=FALSE)
{
  # do initial error checking
  if (NS_ans == F & GD_ans == F)
  {
    stop("At least one of the following must be TRUE: NS_ans or GD_ans")
  }
  
  # read input files and convert to required format
  cat ("Reading input files...\n")
  hierf_dat = genind2hierfstat(genind_obj)
  df_dat = genind2df(genind_obj)
  df_dat$pop = as.character(df_dat$pop)
  
  # set parameters
  numloci = length(genind_obj@loc.names)
  numindivs = length(genind_obj@ind.names)
  OS = as.character(Sys.info()['sysname']) # Needs to be a character for Linux
  NH_genepop_filename = paste(output_name,"_Genepop.gen",sep="")
  NH_summary_filename = paste(output_name,"_sGD.csv",sep="")

  # do additional error checking
  if(is.numeric(min_N)==F){stop("min_N must be an integer")}
  if(is.numeric(xy[,2])==F){stop("The second column of the xy file must be a numeric x coordinate (e.g. longitude)")}  
  if(is.numeric(xy[,3])==F){stop("The third column of the xy file must be a numeric y coordinate (e.g. latitude)")}  
  if(nrow(xy)!=numindivs){stop("The number of rows in the xy file does not match the number of individuals in the genepop file")}
  if(nrow(dist_mat)!=numindivs){stop("The number of rows in the dist_mat does not match the number of individuals in the genepop file")}
  if(ncol(dist_mat)!=numindivs){stop("The number of columns in the dist_mat does not match the number of individuals in the genepop file")}
  if(OS=="Windows")
  {
    if(file.exists(file.path(NeEstimator_dir,"Ne2.exe"))==F){stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
  }
  if(OS=="Linux")
  {
    if(file.exists(file.path(NeEstimator_dir,"Ne2L"))==F){stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")}
    if(file.exists(file.path(NeEstimator_dir,"NeEstimator.jar"))==F){stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
  }  
  if(OS=="Darwin")
  {
    if(file.exists(file.path(NeEstimator_dir,"Ne2M"))==F){stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")}
    if(file.exists(file.path(NeEstimator_dir,"NeEstimator.jar"))==F){stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
  }  
  
  # output a summary of the inputs:
  cat("Input summary:\n")
  cat(paste("\t individuals:",numindivs,"\n"))
  cat(paste("\t loci:",numloci,"\n"))  
  cat(paste("\t neighborhood radius:",radius,"\n"))  
  cat(paste("\t minimum sample size:",min_N,"\n"))
  cat(paste("\t output file:",NH_summary_filename,"in",getwd(),"\n"))
  
  #write genepop header and loci
  header_length = 1+numloci
  genepop_header = output_name
  write.table(genepop_header, NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F)
  
  for (locus in 1:numloci)
  {
    write.table(genind_obj@loc.names[locus], NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
  
  # define neighborhoods and write to genepop file
  cat("Determining neighborhood membership from dist_mat and radius...\n")

  neighborhood_mat = ifelse(dist_mat < radius,1,0)
  neighborhood_N = colSums(neighborhood_mat)
  
  # npops tracks the number of neighborhoods with min_N individuals
  valid_pops = which(neighborhood_N >= min_N)
  npops = length(valid_pops)
                           
  # write neigborhood genepop file
  for (pop in valid_pops)
  {
    write.table("POP", NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
    output = df_dat[which(neighborhood_mat[pop,]==1),]
    indivIDs = paste("P",pop,"I",c(1:nrow(output)),sep="")
    output = data.frame(indivIDs,",",output[2:ncol(output)])
    write.table(output, NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
  
  NH_summary = data.frame(genind_obj@ind.names,xy[,c(2,3)],neighborhood_N,stringsAsFactors=F)
  names(NH_summary) = c("Indiv_ID","X","Y","N")

  if (GD_ans == T)
  {
    # Calculate genetic diversity indices
    cat("Calculating genetic diversity indices for neighborhoods...\n")
  
    NH_dat = import2genind(NH_genepop_filename,quiet=T)
    NH_hierf_dat = genind2hierfstat(NH_dat)

    NH_A = nb.alleles(NH_hierf_dat)
    NH_Ar = round(colMeans(allelic.richness(NH_hierf_dat)$Ar),4)
    NH_Ap = round(colSums(NH_A)/sum(NH_dat@loc.nall),4)
    NH_A = round(colMeans(NH_A),4)
    NH_stats = basic.stats(NH_hierf_dat,digits=3)
    NH_Ho = round(colMeans(NH_stats$Ho),4)
    NH_Hs = round(colMeans(NH_stats$Hs),4)
    NH_FIS = round(1 - NH_Ho/NH_Hs,4)
    
    NH_summary$A = NA
    NH_summary$Ap = NA
    NH_summary$Ar = NA
    NH_summary$He = NA
    NH_summary$Ho = NA
    NH_summary$FIS = NA
    
    NH_summary$A[valid_pops] = NH_A
    NH_summary$Ap[valid_pops] = NH_Ap
    NH_summary$Ar[valid_pops] = NH_Ar
    NH_summary$He[valid_pops] = NH_Hs
    NH_summary$Ho[valid_pops] = NH_Ho
    NH_summary$FIS[valid_pops] = NH_FIS
    
    # show histograms of FIS/HWE
    hist(NH_summary$FIS,xlab=paste("FIS (median =",
          sprintf("%f",median(na.omit(NH_summary$FIS))),")"),main=paste("Histogram of FIS (radius =",radius,")"))   
  }
    
  if (NS_ans == T)
  {  
    
    write.table(1,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F)
    write.table(3,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table(cbind(0.1,0.05,0.02),"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table(NH_genepop_filename,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table("Ne2_output.txt","Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    
    # Calculate Ne
    cat("Calculating NS for neighborhoods...\n")

    # run NeEstimator using OS-specific executable
    
    if (OS=="Windows")
    {
      file.copy(file.path(NeEstimator_dir,"Ne2.exe"),getwd())
      system("Ne2.exe  m:Ne2_input.txt", show.output.on.console = F,ignore.stdout = T, ignore.stderr = T)
      file.remove("Ne2.exe")
    }
    
    if (OS=="Darwin")
    {
      file.copy(file.path(NeEstimator_dir,"Ne2M"),getwd())
      file.copy(file.path(NeEstimator_dir,"NeEstimator.jar"),getwd())
      system("./Ne2M  m:Ne2_input.txt", ignore.stdout = T, ignore.stderr = T)
      file.remove("Ne2M")
      file.remove("NeEstimator.jar")
    }
        
    if (OS=="Linux")
    {
      file.copy(file.path(NeEstimator_dir,"Ne2L"),getwd())
      file.copy(file.path(NeEstimator_dir,"NeEstimator.jar"),getwd())
      system("./Ne2L  m:Ne2_input.txt", ignore.stdout = T, ignore.stderr = T)
      file.remove("Ne2L")
      file.remove("NeEstimator.jar")
    }
    
    # read the Ne results file
    LDNe_output = readLines("Ne2_output.txt")
    LDNe_datalines = grep("Estimated Ne",LDNe_output)
    LDNe_data = unlist(strsplit(LDNe_output[LDNe_datalines], "\\s+")) 
    LDNe_estimates = data.frame(matrix(LDNe_data,nrow=npops,ncol=7,byrow=T)[,4:7])
    names(LDNe_estimates) = c("NS_ex10pct","NS_ex5pct","NS_ex2pct","NS_ex0pct")
      
    # create columns to hold the Ne estimates and assign default value of N
    NH_summary$NS_ex0pct = NA
    NH_summary$NS_ex2pct = NA
    NH_summary$NS_ex5pct = NA
    NH_summary$NS_ex10pct = NA
    
    NH_summary$NS_ex0pct[valid_pops] = LDNe_estimates$NS_ex0pct
    NH_summary$NS_ex2pct[valid_pops] = LDNe_estimates$NS_ex2pct
    NH_summary$NS_ex5pct[valid_pops] = LDNe_estimates$NS_ex5pct
    NH_summary$NS_ex10pct[valid_pops] = LDNe_estimates$NS_ex10pct

    # remove temporary files
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
    if(genout_ans==FALSE)
    {
      file.remove(NH_genepop_filename)
    }
  }
  
  # write NH_summary to .csv file
  cat("Appending results to neighborhood summary file...\n")
  write.table (NH_summary,NH_summary_filename,row.names=F,sep=",",na="")
  
  if(NHmat_ans==TRUE)
  {
    cat("Writing neighborhood membership matrix to file...\n")
    write.table (neighborhood_mat,paste(output_name,"_neighborhood_mat.csv",sep=""),row.names=NH_summary$Indiv_ID,col.names=NH_summary$Indiv_ID,sep=",",na="")
  }
  cat("Processing complete.\n")
}


#' Calculate a pairwise landscape distance matrix (Euclidian or cost-distance). Make sure the projection of the inputs is specified correctly by the \code{proj} argument. If you use a landscape raster to calculate cost-weighted distances, it must match the projection of the xy_file points.
#' 
#' @param method Specify the type of distance matrix to be produced, using "ed" for Euclidean distance and "cd" for cost-weighted (i.e. effective) distance or c("ed","cd") for both. 
#' @param output_name A character string that will be appended to the front of the output filenames.
#' @param points A comma delimited text file containing columns for individual ID, X coordinate, and Y coordinate. If you calculate cost-weighted distances using a landscape raster, make sure the coordinates in the xy file are in the same projection as the raster. Specify the projection in proj4 format using the \code{proj} argument.
#' @param landscape A raster file in a format supported by GDAL listed here: \url{http://www.gdal.org/formats_list.html}.
#' @param proj A projection in proj4 format. These can be found at \url{http://www.spatialreference.org/} by typing keywords for the projection of interest into the search box, clicking on the result that matches your criteria, and then selecting the proj4 format syntax. A projection is required for calculating cost-weighted distances. The distance units in the output matrix will be in the same units as the projection (e.g. meters for UTM projections). 
#' 
#' @return An NxN (N= sample size: i.e. nrow(xy)) matrix of pairwise Euclidean and/or effective landscape distances written to .csv comma delimited files with edmat or cdmat appended to the end of the \code{output_name}.
#' 
#' @examples 
#' library(sGD)
#' setwd("C:/Temp")
#' output_name <- "sGD_demo" 
#' xy_file <- system.file("extdata","sGD_demo_xy.csv",package="sGD")
#' proj <- "+proj=utm +zone=10 +datum=NAD83"   
#' points <- SpatialPoints(read.csv(xy_file)[,c(2,3)],proj4string=CRS(proj)) 
#' landscape_ascii <- system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD") 
#' landscape <- raster(landscape_ascii,crs=CRS(proj))
#' distmat(xy_file,output_name,"ed")
#' distmat(xy_file,output_name,"cd",landscape)
#' distmat(xy_file,output_name,c("ed","cd"),landscape)
distmat <- function(points,output_name,method=c("ed","cd"),inraster=NULL)
{
  # check if the points are projected
  if (is.na(crs(points)) == TRUE)
    {
      print("Warning: the input points have no projection defined.")
    }  
  
  # determine appropriate number of decimal places depending on the input units
  xrange = max(points@coords[,1]) - min(points@coords[,1]) 
  yrange = max(points@coords[,2]) - min(points@coords[,2])
  maxrange = max(xrange,yrange)
  
  if(maxrange>180) # if units are degrees, use many decimal places, otherwise, units of feet/meters only need 1 decimal place.
  {
    dec=1
  } else {
    dec=6
  }
    
  for (type in method)
  {
    if (type=="ed")
    {
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      ed = round(full(dist(points@coords,method="euclidean")),dec)
      
      # write matrices to csv files
      write.table(ed,paste(output_name,"_edmat.csv"),row.names=F,col.names=F,sep=",")
    }
    
    if (type=="cd")
    {
      # check to see if points and landscape are in the same projection
      if (is.na(crs(inraster)) == TRUE)
      {
        print("Warning: the input raster has no projection defined.")
      } else if(as.character(points@proj4string) != as.character(inraster@crs))
      {
        print("The projection of the points and landscape are not the same - please correct and rerun...")
      }
      
      # calculate transition surface, and geocorrect it in E-W dimension
      tr <- transition(inraster,transitionFunction = function(x) {1/mean(x)},directions=8) 
      trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
      
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      cd <- round(full(costDistance(trCorrC, points)),dec)
      
      # write matrix to csv files
      write.table(cd,paste(output_name,"_cdmat.csv"),row.names=F,col.names=F,sep=",")      
    }
  }
  cat(paste("Output distance matrix file(s) written to: ",getwd(),sep=""))
}



