#' Calculate spatially explicit indicies of genetic diversity and Wright's neighborhood size (NS).
#' 
#' @param genepop_file The path to a genotype file in genepop format (with a .gen extension).
#' @param output_name A character string that will be appended to the front of the output filename (will end with "_sGD.csv").
#' @param xy_file A comma delimited file containing columns for individual IDs, X coordinates, and Y coordinates.
#' @param dist_mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective).
#' @param radius The radius of the genetic neighborhood in the same units as \code{dist_mat}.
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param NS_ans Boolean (T or F) answer to whether you want to calculate Wright's neighborhood size.
#' @param GD_ans Boolean (T or F) answer to whether you want to calculate genetic diversity indices. This calculation can take a long time depending on how many individuals are in your sample and the \code{radius} of the neighborhood..
#' @param NeEstimator_dir Path to the NeEstimator 2.0 directory. NeEstimator 2.0 is required only if NS_ans = T. It can be downloaded from \url{http://molecularfisherieslaboratory.com.au/neestimator-software}.
#' @return sGD produces a comma delimited text file containing estimates of genetic diversity and neighborhood size for neighborhoods surrounding each sample location.
#' @return Columns in the output text file include:
#' @return \code{ID} - the ID of the individual at the neighborhood center, corresponding to the ID of the individual in the \code{xy_file}.
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
#' genepop_file = system.file("extdata","sGD_demo_genepop.gen",package="sGD") 
#' output_name = "demo"     
#' xy_file = system.file("extdata","sGD_demo_xy.csv",package="sGD")                  
#' dist_mat = system.file("extdata","sGD_demo_cdmat.csv",package="sGD")  
#' radius = 16238 # for this demo, units are in meters
#' min_N = 20
#' NS_ans = TRUE
#' GD_ans = TRUE
#' NeEstimator_dir = "C:/NeEstimator"
#' 
#' sGD(genepop_file,output_name,xy_file,dist_mat,radius,min_N,NS_ans,GD_ans,NeEstimator_dir)
#' 
#' # specify landscape and sGD output 
#' landscape = raster(system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD"))
#' sGD_output = read.csv(system.file("extdata","sGD_demo_output_sGD.csv",package="sGD"))
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


sGD <- function(genepop_file,output_name,xy_file,dist_mat,radius,min_N,NS_ans,GD_ans,NeEstimator_dir=NULL)
{
  # do initial error checking
  if(file.exists(genepop_file)==F){stop("Cannot read genepop_file. Is the path and filename correct?")}  
  if(file.exists(xy_file)==F){stop("Cannot read xy file. Is the path and filename correct?")} 
  if(file.exists(dist_mat)==F){stop("Cannot read dist_mat. Is the path and filename correct?")} 
  if (NS_ans == F & GD_ans == F)
  {
    stop("At least one of the following must be TRUE: NS_ans or GD_ans")
  }
  
  # read input files and convert to required format
  cat ("Reading input files...\n")
  
  adegen_dat = read.genepop(genepop_file,quiet=T)
  
  sink(file=file()) ## silence genind output 
  hierf_dat = genind2hierfstat(adegen_dat)
  df_dat = genind2df(adegen_dat)
  sink() ## undo silencing
  close(con=file())
  
  xy = read.csv(xy_file)
  dist_mat = read.csv(dist_mat,header=F)
  
  # set parameters
  numloci = length(adegen_dat@loc.names)
  numindivs = length(adegen_dat@ind.names)
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
    write.table(adegen_dat@loc.names[locus], NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
  
  # define neighborhoods and write to genepop file
  cat("Determining neighborhood membership from dist_mat and radius...\n")
  neighborhoods = list() # create an empty list to fill with the indices of neigbhorhood members
  neighborhood_N = c()  # create an empty vector to fill with neighborhood sizes
  
  counter = 1 # counter needed to track neighborhoods with min_N individuals (i.e. npops in gen file)
  
  for (i in 1:numindivs)
  {
    neighborhood_N[i] = length(which(dist_mat[i,] < radius))
    if (neighborhood_N[i]>=min_N)
    {
      neighborhoods[[counter]] = which(dist_mat[i,] < radius)
      counter = counter+1
    }
  }
  
  # npops tracks the number of neighborhoods with min_N individuals
  npops = length(subset(neighborhood_N,neighborhood_N >= min_N))
  
  # write neigborhood files
  for (pop in 1:npops)
  {
    write.table("POP", NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
    output = df_dat[neighborhoods[[pop]],]
    indivIDs = paste("P",pop,"I",c(1:nrow(output)),sep="")
    output = data.frame(indivIDs,",",output[2:ncol(output)])
    write.table(output, NH_genepop_filename,sep="\t",quote=F,col.names=F,row.names=F,append=T)
  }
  
  NH_summary = data.frame(adegen_dat@ind.names,xy[,c(2,3)],neighborhood_N,stringsAsFactors=F)
  names(NH_summary) = c("ID","X","Y","N")

  if (GD_ans == T)
  {
    # Calculate genetic diversity indices and HWE
    cat("Calculating genetic diversity indices for neighborhoods...\n")
  
    NH_adegen_dat = import2genind(NH_genepop_filename,quiet=T)
    NH_hierf_dat = genind2hierfstat(NH_adegen_dat)
    #NH_df_dat = genind2df(NH_adegen_dat)

    NH_A = nb.alleles(NH_hierf_dat)
    NH_Ar = allelic.richness(NH_hierf_dat)$Ar
    NH_Ap = c()
    for (pop in 1:npops)
    {
      NH_Ap = c(NH_Ap,sum(NH_A[,pop])/sum(NH_adegen_dat@loc.nall))
    }
    
    NH_stats = basic.stats(NH_hierf_dat,digits=3)
    NH_Ho = NH_stats$Ho
    NH_Hs = NH_stats$Hs
    NH_FIS = NH_stats$Fis
    
    #sink(file=file()) ## silence genind output 
    #NH_HWE = HWE.test(NH_adegen_dat)
    #sink() ## undo silencing
    
    A = c()
    Ap = c()
    Ar = c()
    He = c()
    Ho = c()
    FIS = c()
    #HWE = c()
    
    counter=1 # counter is needed because if NAs exist, nrow(NH_summary) is > length(GD metrics)
        
    for (i in 1:nrow(NH_summary))
    {
      
      if(NH_summary$N[i]<min_N)
      {
        A[i] = NA
        Ap[i] = NA
        Ar[i] = NA
        He[i] = NA
        Ho[i] = NA
        FIS[i] = NA
        #HWE[i] = NA
      }
      
      else
      {
        A[i] = mean(NH_A[,counter])
        Ap[i] = NH_Ap[counter]
        Ar[i] = mean(NH_Ar[,counter])
        He[i] = mean(NH_Hs[,counter])
        Ho[i] = mean(NH_Ho[,counter])
        FIS[i] = mean(NH_FIS[,counter])
        #HWE[i] = mean(NH_HWE[counter,])
        counter = counter + 1  
      }
    } 
    
    # append GD stats to NH_summary
    NH_summary = data.frame(data.frame(NH_summary,A,Ap,Ar,He,Ho,FIS))
    
    # show histograms of FIS/HWE
    hist(NH_summary$FIS,xlab=paste("FIS (median =",
          median(na.omit(NH_summary$FIS)),")"),main=paste("Histogram of FIS (radius =",radius,")"))   
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
    Ne_results = readLines("Ne2_output.txt")
    LD_Ne_lines = grep("Estimated Ne",Ne_results)
    
    # create empty vectors to hold the Ne estimates
    NS_ex5pct = c()
    NS_ex2pct = c()
    NS_ex10pct = c()
    NS_ex0pct = c()
    
    counter = 1
    
    # population empty vectors with estimates if N > min_N
    for (i in 1:nrow(NH_summary))
    {
      if(NH_summary$N[i]<min_N)
      {
        NS_ex5pct[i] = NA
        NS_ex2pct[i] = NA
        NS_ex10pct[i] = NA
        NS_ex0pct[i] = NA    
      }
      
      else
      {
        LD_Ne_line = Ne_results[LD_Ne_lines[counter]]
        NS_ex5pct[i] = strsplit(LD_Ne_line, "\\s+")[[1]][4]
        NS_ex2pct[i] = strsplit(LD_Ne_line, "\\s+")[[1]][5] 
        NS_ex10pct[i] = strsplit(LD_Ne_line, "\\s+")[[1]][6]
        NS_ex0pct[i] = strsplit(LD_Ne_line, "\\s+")[[1]][7]    
        counter = counter + 1  
      }
    }
    
    # convert text to numeric
    NS_ex5pct = as.numeric(NS_ex5pct)
    NS_ex2pct = as.numeric(NS_ex2pct)
    NS_ex10pct = as.numeric(NS_ex10pct)
    NS_ex0pct = as.numeric(NS_ex0pct)
    
    # append Ne estimates to NH_summary    
    NH_summary = data.frame(NH_summary,NS_ex0pct,NS_ex2pct,NS_ex5pct,NS_ex10pct)
    
    # remove temporary files
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
    file.remove(NH_genepop_filename)
  }
  
  # write NH_summary to .csv file
  cat("Appending results to neighborhood summary file...\n")
  write.table (NH_summary,NH_summary_filename,row.names=F,sep=",",na="")
  cat("Processing complete.\n")
}


#' Calculate a pairwise landscape distance matrix (Euclidian or cost-distance).
#' 
#' @param method Specify the type of distance matrix to be produced, using "ed" for Euclidean distance and "cd" for cost-weighted (i.e. effective) distance or c("ed","cd") for both. 
#' @param output_name A character string that will be appended to the front of the output filenames.
#' @param xy_file A comma delimited text file containing columns for individual IDs, X coordinate, and Y coordinates. The coordinates in the xy file must be in the same projection as the coordinates of the raster if you calculate cost-weighted distances. Specify the projection in proj4 format using the \code{proj} parameter.
#' @param landscape A raster file in a format supported by GDAL listed here: \url{http://www.gdal.org/formats_list.html}.
#' @param proj A projection in proj4 format. These can be easily found at \url{http://www.spatialreference.org/} by typing keywords for the projection of interest into the search box, clicking on the result that matches your criteria, and then selecting the proj4 format syntax. A projection is required for calculating cost-weighted distances, so if projection information is not included with your raster file and no projection is specified, the default of UTM Zone 10 NAD83 will be used. The cost-distance units in the output matrix will be in the same units as the projection (e.g. meters for UTM projections). 
#' @return An NxN (N= sample size: i.e. nrow(xy)) matrix of pairwise Euclidean and/or effective landscape distances written to .csv comma delimited files with edmat or cdmat appended to the end of the \code{output_name}.
#' @examples 
#' library(sGD)
#' setwd("C:/Temp")
#' landscape = system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD")
#' output_name = "sGD_demo"     
#' xy_file = system.file("extdata","sGD_demo_xy.csv",package="sGD") 
#' proj = "+proj=utm +zone=10 +datum=NAD83"                    
#' 
#' distmat(xy_file,output_name,"ed")
#' distmat(xy_file,output_name,"cd",proj,landscape)
#' distmat(xy_file,output_name,c("ed","cd"),proj,landscape)
distmat <- function(xy_file,output_name,method=c("ed","cd"),proj=NULL,landscape=NULL)
{
  for (type in method)
  {
    if (type=="ed")
    {
      # read points csv file and reformat to work with dist function
      locs = read.csv(xy_file)[,c(2,3)]
      
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      ed = round(full(dist(locs,method="euclidean")),0)
      
      # write matrices to csv files
      write.table(ed,paste(output_name,"_edmat.csv"),row.names=F,col.names=F,sep=",")
    }
    
    if (type=="cd")
    {
      # read resistance raster and project it
      resistance <- raster(landscape)
      if (is.null(proj) == TRUE)
      {
        proj = "+proj=utm +zone=10 +datum=NAD83"       
      }
      
      if (is.na(crs(resistance)) == TRUE)
      {
      projection(resistance) <- CRS(proj)
      }
      
      # read points csv file and reformat to work with costDistance function
      locs = SpatialPoints(read.csv(xy_file)[,c(2,3)],proj4string=crs(resistance))
      
      # calculate transition surface, and geocorrect it in E-W dimension
      tr <- transition(resistance,transitionFunction = function(x) {1/mean(x)},directions=8) 
      trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
      
      # calculate costdistance and euclidean distance matrices - be careful with rounding if map units aren't meters
      cd <- round(full(costDistance(trCorrC, locs)),0)
      
      # write matrix to csv files
      write.table(cd,paste(output_name,"_cdmat.csv"),row.names=F,col.names=F,sep=",")      
    }
  }
  cat(paste("Output distance matrix file(s) written to: ",getwd(),sep=""))
}



