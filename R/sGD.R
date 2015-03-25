#' Calculate spatially explicit indicies of genetic diversity and Wright's neighborhood size (NS).
#' 
#' @param genind_obj A genind object (created by the adegenet package function import2genind or similar methods) containing individual genotypes.
#' @param file_name (optional) A character string that will be appended to the front of the output filename (will end with "_sGD.csv"). If none specified, no output file will be written.
#' @param xy A dataframe containing 3 columns in the following order: individual IDs, X coordinates, and Y coordinates.
#' @param dist_mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective).
#' @param radius The radius of the genetic neighborhood in the same units as \code{dist_mat}.
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param NS_ans Boolean (T or F) answer to whether you want sGD to calculate Wright's neighborhood size.
#' @param GD_ans Boolean (T or F) answer to whether you want sGD to calculate genetic diversity indices. This calculation can take a long time depending on how many individuals are in your sample and the \code{radius} of the neighborhood.
#' @param NeEstimator_dir Path to the NeEstimator 2.0 directory. NeEstimator 2.0 is required only if NS_ans = T. It can be downloaded from \url{http://molecularfisherieslaboratory.com.au/neestimator-software}.
#' @param NHmat_ans Boolean (T or F) answer to whether you want sGD to write a matrix defining neighborhood membership. For each row/column of the matrix, a value of 1 occurs at the indices of all individuals inside the neighborhood and a value of 0 occurs for all individuals outside the neighborhood.
#' @param genout_ans Boolean (T or F) answer to whether you want sGD to write a genepop file containing the genotypes for all neighborhoods to the working directory. 
#' 
#' @return sGD returns a data frame containing estimates of genetic diversity and/or neighborhood size for neighborhoods surrounding each sample location.
#' @return Variables in the output data frame include:
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
#' genind_obj <- read.genepop(genepop_file,quiet=T)
#' xy = read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
#' dist_mat <- as.matrix(read.csv(system.file("extdata","sGD_demo_cdmat.csv",package="sGD") ,header=F)) 
#' radius <- 18000 # for this demo, units are in meters
#' min_N <- 20
#' NeEstimator_dir <- "C:/NeEstimator"
#' 
#' sGD_output = sGD(genind_obj,xy,dist_mat,radius,min_N,NS_ans=TRUE,GD_ans=F,NHmat_ans=T,genout_ans=T,file_name=NULL,NeEstimator_dir=NULL)

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


sGD <- function(genind_obj,xy,dist_mat,radius,min_N,NS_ans=F,GD_ans=T,NHmat_ans=F,genout_ans=F,file_name=NULL,NeEstimator_dir=NULL)
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

  NH_genepop_filename = "TempGenepop.gen"
  
  if(is.null(file_name)==FALSE)
  {
    NH_summary_filename = paste(file_name,"_sGD.csv",sep="") 
  }


  # do additional error checking
  if(is.numeric(min_N)==F){stop("min_N must be an integer")}
  if(is.numeric(xy[,2])==F){stop("The second column of the xy file must be a numeric x coordinate (e.g. longitude)")}  
  if(is.numeric(xy[,3])==F){stop("The third column of the xy file must be a numeric y coordinate (e.g. latitude)")}  
  if(nrow(xy)!=numindivs){stop("The number of rows in the xy file does not match the number of individuals in the genepop file")}
  if(nrow(dist_mat)!=numindivs){stop("The number of rows in the dist_mat does not match the number of individuals in the genepop file")}
  if(ncol(dist_mat)!=numindivs){stop("The number of columns in the dist_mat does not match the number of individuals in the genepop file")}
  
  if(is.null(NeEstimator_dir)==F)
  {
    if(OS=="Windows")
    {
      if(file.exists(file.path(NeEstimator_dir,"Ne2.exe"))==F)
      {stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
    }
    if(OS=="Linux")
    {
      if(file.exists(file.path(NeEstimator_dir,"Ne2L"))==F)
      {stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")}
      if(file.exists(file.path(NeEstimator_dir,"NeEstimator.jar"))==F)
      {stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
    }  
    if(OS=="Darwin")
    {
      if(file.exists(file.path(NeEstimator_dir,"Ne2M"))==F)
      {stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")}
      if(file.exists(file.path(NeEstimator_dir,"NeEstimator.jar"))==F)
      {stop("Cannot find NeEstimator executable. Is the path to NeEstimator_dir correct?")} 
    }  
  }
  
  if(is.null(NeEstimator_dir)==T & NS_ans==T)
  {
    stop("NS_ans is TRUE, however, you have not specified the location of NeEstimator_dir")
  }
  
  # output a summary of the inputs:
  cat("Input summary:\n")
  cat(paste("\t individuals:",numindivs,"\n"))
  cat(paste("\t loci:",numloci,"\n"))  
  cat(paste("\t neighborhood radius:",radius,"\n"))  
  cat(paste("\t minimum sample size:",min_N,"\n"))
  if(is.null(file_name)==FALSE)
  {
  cat(paste("\t output file:",NH_summary_filename,"in",getwd(),"\n"))
  }
  
  #write genepop header and loci
  header_length = 1+numloci
  genepop_header = "TempGenepop"
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
  valid_pops = as.numeric(which(neighborhood_N >= min_N))
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
  NH_summary$index = c(1:numindivs)

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
    
    GD_output = data.frame(NH_summary$index[valid_pops],NH_A,NH_Ap,NH_Ar,NH_Hs,NH_Ho,NH_FIS,stringsAsFactors=F)
    names(GD_output) = c("index","A","Ap","Ar","Hs","Ho","FIS")
    
    NH_summary = merge(NH_summary,GD_output,by="index",all=T)
    
    # show histograms of FIS/HWE
    #hist(NH_summary$FIS,xlab=paste("FIS (median =",
    #      sprintf("%f",median(na.omit(NH_summary$FIS))),")"),main=paste("Histogram of FIS (radius =",radius,")"))   
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
    LDNe_estimates = data.frame(matrix(LDNe_data,nrow=npops,ncol=7,byrow=T)[,4:7],stringsAsFactors=F)
    LDNe_estimates = data.frame(NH_summary$index[valid_pops],LDNe_estimates,stringsAsFactors = F)
    names(LDNe_estimates) = c("index","NS_ex10pct","NS_ex5pct","NS_ex2pct","NS_ex0pct")
    
    # create columns to hold the Ne estimates and assign default value of N

    NH_summary = merge(NH_summary,LDNe_estimates,by="index",all=T)

    # remove temporary files
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
    if(genout_ans==FALSE)
    {
      file.remove(NH_genepop_filename)
    }
  }
  
  for(col in c(3:ncol(NH_summary)))
  {
    NH_summary[,col]=as.numeric(NH_summary[,col])
  }
  
  if(is.null(file_name)==FALSE)
  {
    # write NH_summary to .csv file
    cat("Appending results to neighborhood summary file...\n")
    write.table (NH_summary,NH_summary_filename,row.names=F,sep=",",na="")
  }
  
  if(NHmat_ans==TRUE)
  {
    cat("Writing neighborhood membership matrix to file...\n")
    write.table (neighborhood_mat,paste(file_name,"_neighborhood_mat.csv",sep=""),row.names=NH_summary$Indiv_ID,col.names=NH_summary$Indiv_ID,sep=",",na="")
  }
  
  return(NH_summary)
  cat("Processing complete.\n")
}


#' Calculate a pairwise landscape distance matrix (Euclidian or cost-distance).  
#' 
#' @param method Specify the type of distance matrix to be produced, using "ed" for Euclidean distance and "cd" for cost-weighted (i.e. effective) distance or c("ed","cd") for both. If you calculate cost-weighted distances, make sure the \code{points} projection is the same as the \code{landscape} raster.
#' @param file_name (optional) A character string that will be appended to the beginning of the output filename. If no name is specified, no file will be written to the working directory.
#' @param sp_points An object of class SpatialPoints (see raster package).
#' @param landscape An object of class RasterLayer (see raster package)
#' 
#' @return An NxN (N= sample size: i.e. nrow(xy)) matrix of pairwise Euclidean and/or effective landscape distances written to .csv comma delimited files with edmat or cdmat appended to the end of the \code{filename}.
#' 
#' @examples 
#' library(sGD)
#' setwd("C:/Temp")
#' file_name <- "sGD_demo" 
#' xy_file <- system.file("extdata","sGD_demo_xy.csv",package="sGD")
#' proj <- "+proj=utm +zone=10 +datum=NAD83"   
#' sp_points <- SpatialPoints(read.csv(xy_file)[,c(2,3)],proj4string=CRS(proj)) 
#' landscape_ascii <- system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD") 
#' landscape <- raster(landscape_ascii,crs=CRS(proj))
#' distmat(sp_points,file_name,"ed")
#' distmat(sp_points,file_name,"cd",landscape)
#' distmat(sp_points,file_name,c("ed","cd"),landscape)
distmat <- function(sp_points,method=c("ed","cd"),file_name=NULL,landscape=NULL)
{
  # check if the points are projected
  if (is.na(crs(sp_points)) == TRUE)
    {
      print("Warning: the input points have no projection defined.")
    }  
  
  # determine appropriate number of decimal places depending on the input units
  xrange = max(sp_points@coords[,1]) - min(sp_points@coords[,1]) 
  yrange = max(sp_points@coords[,2]) - min(sp_points@coords[,2])
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
      ed = round(full(dist(sp_points@coords,method="euclidean")),dec)
      
      # write matrices to csv files
      if(is.null(file_name==F))
      {
        write.table(ed,paste(file_name,"_edmat.csv"),row.names=F,col.names=F,sep=",")         
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
        print("The projection of the points and landscape are not the same - please correct and rerun...")
      }
      
      # calculate transition surface, and geocorrect it in E-W dimension
      tr <- transition(landscape,transitionFunction = function(x) {1/mean(x)},directions=8) 
      trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
      
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      cd <- round(full(costDistance(trCorrC, sp_points)),dec)
      
      # write matrix to csv files
      if(is.null(file_name==F))
      {
        write.table(cd,paste(file_name,"_cdmat.csv"),row.names=F,col.names=F,sep=",")         
      }
    }
  }

  if(method=="ed") return(ed)
  if(method=="cd") return(cd)
  if(length(method)==2) return(list(ed,cd))
}



