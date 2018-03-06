#' Calculate spatially explicit indicies of genetic diversity and Wright's neighborhood size (NS).
#' 
#' @param genind_obj A genind object (created by the adegenet package function import2genind or related methods) containing individual genotypes. The order of the individuals must be the same as the order in the xy and dist.mat inputs below.
#' @param file_name (optional) A character string that will be appended to the front of the output filename (will end with "_sGD.csv"). If none specified, no output file will be written.
#' @param xy A dataframe containing 3 columns in the following order: individual IDs, X coordinates, and Y coordinates. The order of the rows must match the order in the genind_obj and dist.mat inputs.
#' @param dist.mat An NxN (N= sample size) matrix of pairwise landscape distances (Euclidean or effective). The \code{distmat} function in the sGD package may be used to produce Euclidean and cost-weighted distance matrices. The order of the rows and columns in the matrix must match the order in the xy and genind_obj inputs.
#' @param NH_radius A single value to be used for all neighborhoods, or a vector of genetic neighborhood radii, optionally obtained from using the \code{infer2sigma} function. Note that if you specify a vector of radii, sGD will not calculate metrics when the radius value is NA.
#' @param metrics Optional. Provide a vector of the metrics you would like sGD to produce. Options include "GD" (genetic diversity indices), "NS" (Wright's genetic neighobrhood size), "HWE" (tests for Hardy-Weinberg equilibrium, heterozygote excess, and homozygote excess), and "pFST" (a matrix of pairwise FST values for all neighborhoods)". Note that calculating pFST takes considerable time (several hours using the sGD demo data).
#' @param min_N The minimum sample size per neighborhood for indices to be calculated. NA is returned for neighborhoods < \code{min_N}.
#' @param max_N Optional. The maximum sample size per neighborhood for indices to be calculated. If the number of individuals in the neighborhood exceeds \code{max_N}, a sample of size \code{max_N} will be used from the neighborhood to compute the metrics and output files specified by the user. Note that if \code{max_N} is specified, and the value is too small to be representative of the neighobrhood, the results could differ significantly compared to if all individuals in the neighborhood were used.
#' @param NeEstimator_dir Optional. Path to the NeEstimator directory. NeEstimator 2.0 or > is required only if you include the "NS" metric. It can be downloaded from \url{http://molecularfisherieslaboratory.com.au/neestimator-software}.
#' @param NHmat_ans Logical (Default = FALSE). If TRUE, a matrix defining neighborhood membership is written to the working directory. For each row in the matrix, a value of 1 occurs at the indices of all individuals inside the neighborhood and a value of 0 occurs for all individuals outside the neighborhood. 
#' @param genout_ans Logical (Default = FALSE). If TRUE, a genepop file containing the genotypes for all neighborhoods is written to the working directory. 
#' 
#' @return sGD returns a data frame containing estimates of genetic diversity and/or neighborhood size for neighborhoods surrounding each sample location. The order of the rows in the output matches the order of the samples in the inputs.
#' @return Variables in the output data frame include (depending on the metrics selected):
#' @return \code{NH_Index} - an index of the neighborhoods, from 1 to the total number of neighborhoods. 
#' @return \code{NH_ID} - the ID of the individual at the neighborhood center, taken from individual's ID in the \code{xy_file}.
#' @return \code{X} - the X coordinate of the neighborhood center.
#' @return \code{Y} - the Y coordinate of the neighborhood center.
#' @return \code{N} - the number of individuals within the neighborhood.
#' @return \code{A} - the average number of alleles across all loci/individuals within the neighborhood.
#' @return \code{Ap} - the proportion of alleles from the entire population that area actually present in the neighborhood.
#' @return \code{Ar} - the allelic richness across all loci/individuals within the neighborhood.
#' @return \code{He} - the average expected heterozygosity across all loci/individuals within the neighborhood.
#' @return \code{Ho} - the average observed heterozygosity across all loci/individuals within the neighborhood.
#' @return \code{FIS} - the average inbreeding coefficient across all loci/individuals within the neighborhood.
#' @return \code{NS_ex0} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, not exluding rare alleles that could bias the estimate.
#' @return \code{NS_ex0.02} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.02 or less.
#' @return \code{NS_ex0.05} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.05 or less.
#' @return \code{NS_ex0.10} - an estimate of the effective number of breeding indviduals (Wright's neighborhood size) present within the neighborhood, exluding alleles with a frequency of 0.10 or less.
#' @return If specified in the sGD arguments, the following output files will also be written to the working directory:
#' @return \code{NHmat} - if NHmat_and = TRUE, sGD writes the NH membership matrix described above to a .csv file in the working directory. The row and column names match the individual ID's in the input files, and are in the same order as the input files. 
#' @return \code{genout} - if genout_ans = TRUE, sGD writes the NH genepop file to the working directory.
#' 
#' @examples 
#' library(sGD)
#' library(adegenet)
#' library(raster)
#' 
#' # read in genotypes, locations, and distance matrix
#' genepop.file <- system.file("extdata","sGD_demo_IBR.gen",package="sGD") 
#' xy = read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
#' dist.mat <- as.matrix(read.csv(system.file("extdata","sGD_demo_cdmat.csv",package="sGD"),
#'                                header=FALSE))
#' 
#' # convert genepop to genind (make sure you specify the correct allele code digits - ncode)
#' genind_obj <- read.genepop(genepop.file,ncode=3L,quiet=TRUE)
#' pop(genind_obj) = xy$Indiv_ID # give each location a unique population ID
#' 
#' # run sGD
#' sGD_output <- sGD(genind_obj,xy,dist.mat,NH_radius=16000,min_N=20,max_N=NULL,
#'                   metrics=c("GD","NS","HWE"), NHmat_ans=TRUE,genout_ans=TRUE,
#'                   file_name="sGD_demo", NeEstimator_dir="C:/NeEstimator")
#'
#' # read in the landscape raster to use in plots
#' landscape <- raster(system.file("extdata","sGD_demo_IBR_landscape.asc",package="sGD"))
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
#'        
#' @importFrom adegenet genind2df
#' @importFrom adegenet read.genepop
#' @importFrom hierfstat allelic.richness
#' @importFrom hierfstat basic.stats
#' @importFrom hierfstat nb.alleles 
#' @importFrom hierfstat genind2hierfstat 
#' @importFrom diveRsity basicStats
#' 
#' @export
#' 
sGD <- function(genind_obj,xy,dist.mat,NH_radius,min_N,max_N=NULL,metrics=NULL,NHmat_ans=FALSE,genout_ans=FALSE,file_name=NULL,NeEstimator_dir=NULL)
{
  
  #### check to make sure correct package versions are installed ####
  if(packageVersion("adegenet") < "2.0.0") {
    stop("Please install the latest version of the adegenet package (>= 2.0.1)")
  }
  
  if(packageVersion("hierfstat") < "0.04.15") {
    stop("Please install the development version of the hierfstat package (>= 0.04.22) 
       First, install the devtools package, and then run:
               library(devtools)
               install_github(\"jgx65/hierfstat\")
       After installing hierfstat, please restart R before running sGD.")
  }
  
  #### read input files and convert to required format ####
  cat ("Reading input files...\n")
  df_dat = genind2df(genind_obj)
  
  allele.digits = nchar(genind_obj@all.names[[1]][1])
  
  #### set parameters ####
  numloci = length(names(genind_obj@all.names))
  N = dim(genind_obj@tab)[1]
  OS = as.character(Sys.info()['sysname']) # Needs to be a character for Linux
  
  genout_file = paste(file_name,"_genepop.gen",sep="")
  
  #### perform error checking on inputs ####
  if (is.null(metrics) == T & NHmat_ans == F & genout_ans == F)
  {
    stop("At least one metric must be specified or one of the following must be TRUE: NHmat_ans,
         or genout_ans")
  }
  
  if (is.null(file_name) & NHmat_ans == T)
  {
    stop("If NHmat_ans = TRUE, you must specify a file_name")
  }
  
  if (is.null(file_name) & genout_ans == T)
  {
    stop("If genout_ans = TRUE, you must specify a file_name")
  }
  
  if (is.null(file_name) & "pFST" %in% metrics)
  {
    stop("If you include pFST as a metric, you must specify a file_name")
  }
  
  if(is.numeric(min_N)==F){stop("min_N must be an integer")}
  if(is.numeric(xy[,2])==F){stop("The second column of the xy file must be a numeric x coordinate (e.g. longitude)")}  
  if(is.numeric(xy[,3])==F){stop("The third column of the xy file must be a numeric y coordinate (e.g. latitude)")}  
  if(nrow(xy)!=N){stop("The number of rows in the xy file does not match the number of individuals in the genepop file")}
  if(nrow(dist.mat)!=N){stop("The number of rows in the dist.mat does not match the number of individuals in the genepop file")}
  if(ncol(dist.mat)!=N){stop("The number of columns in the dist.mat does not match the number of individuals in the genepop file")}
  
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
  
  if(is.null(NeEstimator_dir)==T & "NS" %in% metrics)
  {
    stop("You have specified NS as a metric, but you have not specified the location of NeEstimator_dir")
  }
  
  #### output a summary of the inputs ####
  cat("Input summary:\n")
  cat(paste("\t individuals:",N,"\n"))
  cat(paste("\t loci:",numloci,"\n"))  
  cat(paste("\t minimum sample size:",min_N,"\n"))
  if(is.numeric(max_N))
  {
    cat(paste("\t maximum sample size:",max_N,"\n"))
  } else {
    cat(paste("\t maximum sample size: NA","\n"))
  }
  
  if(length(NH_radius)>1)
  {
    min_R = min(NH_radius,na.rm = T)
    max_R = max(NH_radius,na.rm = T)
    cat(paste0("\t variable neighborhood radius (min = ",min_R,";max = ",max_R,")","\n"))
  } else {
    cat(paste("\t neighborhood radius:",NH_radius,"\n"))
  }
  
  
  if(!is.null(file_name))
  {
    sGD_outfile = paste0(file_name,"_sGD.csv") 
    cat(paste("\t sGD output file:",sGD_outfile,"in",getwd(),"\n"))
  } 
  
  if(NHmat_ans==TRUE)
  {
    NHmat_file = paste0(file_name,"_neighborhood_matrix.csv")
    cat(paste("\t NH matrix file:",NHmat_file,"in",getwd(),"\n"))
  }  
  
  if(genout_ans==TRUE)
  {
    cat(paste("\t NH genepop file:",genout_file,"in",getwd(),"\n"))
  } 
  
  #### define neighborhoods and write to genepop file ####
  cat("Determining neighborhood membership from dist.mat and NH_radius...\n")
  
  # create NH matrix
  NHmat = c()
  
  for(r in c(1:N))
  {
    if(length(NH_radius==1))
    {
      radius = NH_radius
    } else {
      radius = NH_radius[r]
    }
    
    rvals = ifelse(dist.mat[r,] < radius,1,0)
    NHmat = rbind(NHmat,rvals)
  }
  
  #### calculate neighborhood parameters ####
  NH_N = as.numeric(rowSums(NHmat,na.rm = T)) # number of individuals per neighborhood
  
  if(max(NH_N)<min_N)
  {
    stop("There are no neighborhoods with at least min_N individuals at the radius selected. Aborting sGD.")
  }
  
  if(length(NH_radius)==1)
  {
    valid_pops = as.numeric(which(NH_N >= min_N))
  } else {
    valid_pops = as.numeric(which(!is.na(NH_radius) & NH_N >= min_N)) # which pops have > min_N individuals
  }
  
  npops = length(valid_pops) # tracks the number of neighborhoods with min_N individuals
  NH_Index = c(1:N) # a numerical means to consistently refer to neighborhoods (all, including those < min_N)
  
  # create NH summary file
  NH_summary = data.frame(NH_Index,as.character(genind_obj@pop),xy[,c(2,3)],NH_N,stringsAsFactors=F)
  names(NH_summary) = c("NH_Index","NH_ID","X","Y","N")
  
  #### write genepop file if required ####
  if(!is.null(metrics) | genout_ans == T)
  {
    # write header and loci
    header_length = 1+numloci
    genepop_header = "sGD neighboorhood genepop file (each POP is a genetic neighborhood)"
    write.table(genepop_header, genout_file,sep="\t",quote=F,col.names=F,row.names=F)
    
    for (locus in 1:numloci)
    {
      write.table(names(genind_obj@all.names)[locus], genout_file,sep="\t",quote=F,col.names=F,row.names=F,append=T)
    }
    
    # write genotypes for each POP
    for (pop in valid_pops)
    {
      # write POP to start each population
      write.table("POP", genout_file,sep="\t",quote=F,col.names=F,row.names=F,append=T)
      
      # check if maxN is specified and, if so, randomly subsample NH to max_N indivs
      if(!is.null(max_N) && NH_summary$N[pop] > max_N){
        pop_ind = sample(which(NHmat[pop,]==1),max_N)
      } else {
        pop_ind = which(NHmat[pop,]==1) 
      }
      
      # get genotypes for NH indivs and append to genepop file
      NH_IDs = paste0("POP",pop,"_",NH_summary$NH_ID[pop_ind])
      output = df_dat[pop_ind,]
      output = data.frame(NH_IDs,",",output[2:ncol(output)])
      write.table(output, genout_file,sep="\t",quote=F,col.names=F,row.names=F,append=T)
    }
  }
  
  #### if selected, run GD routine ####
  if ("GD" %in% metrics)
  {
    # Calculate genetic diversity indices
    cat("Calculating genetic diversity indices...\n")
    
    # read in NH genepop file
    NH_genind = read.genepop(genout_file,ncode=allele.digits,quiet=T)
    
    # fix issue with population names being taken from NH_IDs
    #NH_genind@pop = as.factor(rep(c(1:npops),NH_summary$N[valid_pops]))
    
    # convert to hierfstat format 
    NH_hierfstat = genind2hierfstat(NH_genind)
    
    # calculate GD indices
    NH_A = nb.alleles(NH_hierfstat)
    NH_Ar = round(colMeans(allelic.richness(NH_hierfstat)$Ar),4)
    NH_Ap = round(colSums(NH_A)/sum(NH_genind@loc.n.all),4)
    NH_A = round(colMeans(NH_A),4)
    NH_stats = basic.stats(NH_hierfstat,digits=4)
    NH_Ho = round(colMeans(NH_stats$Ho),4)
    NH_Hs = round(colMeans(NH_stats$Hs),4)
    NH_FIS = round(colMeans(NH_stats$Fis),4)
    
    # assemble GD indices into dataframe and merge with NH_summary
    GD_output = data.frame(NH_Index[valid_pops],NH_A,NH_Ap,NH_Ar,NH_Hs,NH_Ho,NH_FIS,stringsAsFactors=F)
    names(GD_output) = c("NH_Index","A","Ap","Ar","Hs","Ho","FIS")
    
    NH_summary = merge(NH_summary,GD_output,by="NH_Index",all=T)
  }
  
    #### if selected, run HWE routine ####
  if ("HWE" %in% metrics)
  {
    # Calculate genetic diversity indices
    cat("Testing Hardy-Weinberg equilibrium for neighborhoods...\n")
    
    NH_diveRsity = basicStats(genout_file)
    NH_HWE_p = as.numeric(round(colMeans(NH_diveRsity$hwe_llr_p[,-1],na.rm = T),4))
    NH_HetEx = as.numeric(round(colMeans(NH_diveRsity$hwe_het[,-1],na.rm = T),4))
    NH_HomEx = as.numeric(round(colMeans(NH_diveRsity$hwe_hom[,-1],na.rm = T),4))
    
    # assemble HWE indices into dataframe and merge with NH_summary
    HWE_output = data.frame(NH_Index[valid_pops],NH_HWE_p,NH_HetEx,NH_HomEx,stringsAsFactors=F)
    names(HWE_output) = c("NH_Index","HWE_p","HetEx_p","HomEx_p")
    NH_summary = merge(NH_summary,HWE_output,by="NH_Index",all=T)
  }
  
  #### if selected, run NS routine ####
  if ("NS" %in% metrics)
  {  
    # Calculate Ne
    cat("Calculating NS...\n")
    
    write.table(1,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F)
    write.table(3,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table(cbind(0.1,0.05,0.02),"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table(genout_file,"Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    write.table("Ne2_output.txt","Ne2_input.txt",sep="\t",col.names=F,row.names=F,quote=F,append=T)
    
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
    LDNe_data = suppressWarnings(as.numeric(unlist(strsplit(LDNe_output[LDNe_datalines], "\\s+")))) 
    LDNe_estimates = data.frame(matrix(LDNe_data,nrow=npops,ncol=7,byrow=T)[,4:7],stringsAsFactors=F)
    LDNe_estimates = data.frame(NH_Index[valid_pops],LDNe_estimates,stringsAsFactors = F)
    names(LDNe_estimates) = c("NH_Index","NS_ex0.10","NS_ex0.05","NS_ex0.02","NS_ex0.00")
    
    # append NS results to NH_summary
    NH_summary = merge(NH_summary,LDNe_estimates,by="NH_Index",all=T)
    
    # remove temporary files
    file.remove("Ne2_input.txt")
    file.remove("Ne2_output.txt")
  }
  
  #### write selected files to disk ####
  if(!is.null(file_name))
  {
    # write NH_summary to .csv file
    cat("Writing sGD output file...\n")
    write.table (NH_summary,sGD_outfile,row.names=F,sep=",",na="")
  }
  
  if(NHmat_ans==TRUE)
  {
    cat("Writing neighborhood membership matrix to file...\n")
    dimnames(NHmat) = list(NH_summary$NH_ID,NH_summary$NH_ID)
    write.table(matrix(c(NA,NH_summary$NH_ID),nrow=1),NHmat_file,row.names=F,col.names=F,sep=",",na="")
    write.table(NHmat,NHmat_file,row.names=NH_summary$NH_ID,col.names=F,sep=",",na="",append=T)
  }
  
  if ("Dr" %in% metrics)
  {
    # Calculate genetic diversity indices
    cat("Writing Roger's r genetic distance matrix to file...\n")
    
    # read in NH genepop file
    if(!exists("NH_genind"))
    {
      NH_genind = read.genepop(genout_file,ncode=allele.digits,quiet=T)
    }
    
    # fix issue with population names being taken from NH_IDs
    #NH_genind@pop = as.factor(rep(c(1:npops),NH_summary$N[valid_pops]))
    
    # convert to hierfstat format 
    if(!exists("NH_hierfstat"))
    {
      NH_hierfstat = genind2hierfstat(NH_genind)
    }
    
    # calculate Roger's r
    Dr = genet.dist(NH_hierfstat,method="Dr")
    
    # write distance matrix to disk
    write.table(full(Dr),paste0(file_name,"_Dr.csv"),sep=",",col.names = F, row.names = F)
  }
  
  if(genout_ans==FALSE)
  {
    file.remove(genout_file)
  } else {
    cat("Writing neighborhood genepop file...\n")
  }
  
  return(NH_summary)
  cat("Processing complete.\n")
  cat("\n")
}