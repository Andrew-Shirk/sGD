---
title: "sGD-tutorial"
author: "Andrew Shirk and Samuel Cushman"
date: "2018-03-06"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{sGD-tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

# sGD description

sGD provides spatially explicit estimates of genetic diversity indices and Wright's neighborhood size in continuous populations isolated by distance or resistance as described in Shirk & Cushman 2014. Inputs include a set of spatially referenced codominant marker genotypes, a matrix of pairwise distances between all sampled individuals (either in Euclidean units or as effective distances quantified by methods such as least-cost-path or circuit theory), and several parameters, including the minimum neighborhood size and the neighborhood radius (specified in Euclidean or effective distances). sGD defines a local neighborhood extent (based on a user-specified radius) around each sample location, and then calculates genetic diversity indices and/or Wright's neighborhood size (a local measure of effective population size that accounts for the continuous structure of populations isolated by distance). 

# Installing sGD

To install sGD from github, you will need the devtools package installed first.


```r
install.packages("devtools")
```

Then you can install sGD from the github respository and load the library. Note that sGD will also load any packages it depends on during the install, if you don't already have them in your library. These dependencies allow sGD to calculate genetic diversity indices, but if you want to also calculate Wright's Neighborhood Size (see tutorial below), you'll need to download the free program NeEstimator 2.0 from http://molecularfisherieslaboratory.com.au/neestimator-software.


```r
library("devtools")
install_github("Andrew-Shirk/sGD")
library("sGD")
```

Now you're ready to use the sGD functions and demo data. The tutorial that follows is designed to provide a basic introduction and guide you in the use of sGD.

# sGD Introduction and Tutorial
Maintenance of genetic diversity in a population is critical to long-term survival and adaptation in a changing environment (Frankham et al. 2002). Populations lose allelic diversity and heterozygosity through the processes of genetic drift and inbreeding, respectively (Amos and Harwood, 1998). If this diversity is not replaced by the processes of mutation and immigration, the population loses alleles and heterozygosity at a rate proportional to the number of breeding
individuals (Wright, 1931). Diminished genetic diversity may reduce a population’s capacity for adaptation in dynamic environments (though dominance, epistasis, and pleiotropy may strongly affect this relationship; Reed and Frankham, 2001), lower the frequency of heterozygotes (and therefore the prevalence of heterozygote vigor or heterosis), favor inbreeding depression, and ultimately increase the risk of extinction (Frankham et al., 2002).

Wright (1931) first defined the concept of effective population size (*N~e~*) as the number of breeding individuals in an ideal population experiencing the same rate of drift as the real population of interest. Factors like unequal sex ratios, a fluctuating population size, overlapping generations, high variability in the distribution of offspring, and non-random mating make the effective size of a population lower than the census size (*N~c~*) because they influence the number of individuals that contribute their genetic variation to the next generation (Frankham et al., 2002). Therefore, the risk of extinction posed by reduced genetic diversity is more closely related to *N~e~* than *N~c~*, making *N~e~* is an important parameter in
conservation genetics (Lande and Barrowclough, 1987). 

To measure indices of genetic diversity and effective size for a groups of individuals meeting the classical definition of a population (i.e. each subpopulation is internally panmictic and having little or no immigration from other subpopulations), one simply groups individuals by subpopulation and calculates genetic indices for each. There are a number of excellent R packages for this, including adegenet, hierfstat, and diveRsity. Standalone software such as Genepop, FSTAT, Arlequin, and Genalex also perform these calculations. An example of a population meeting the classical definition might be a species of fish inhabiting a series of lakes isolated from one another. Genetic analyses on populations meeting the classical definition of a population can be done without bias, so long as each subpopulation is appropriately delineated (e.g. not grouped together with other subpopulations). Furthermore, there is no expectation of spatial heterogeneity in these genetic indices with a discrete subpopulation, because they are panmictic.

Sample groupings are obvious for discrete populations, but what about continuous populations? Wright (1946) first posed an answer to this question. Under his concept, mating probabilities and relatedness among individuals across a landscape are a function of Euclidean distance. So individuals farther apart are more differentiated than individuals nearby, and there are no sharp discontinuities in this pattern of differentiation. Relatedness is continuously varying across space as a smooth function of distance. He called this dynamic 'Isolation by Distance' (Wright, 1931). He thought of these continuous populations as being comprised of many small subpopulations that overlap with each other in space, and termed these local breeding populations 'genetic neighborhoods'. Within each neighborhood, mating approaches (but never reaches) random. Under Wright's definition of a genetic neighborhood, the number of effective individuals (i.e. those contributing their genetic variation to the next generation) within the local extent of breeding (*NS*) where most matings occur is defined by:

>$NS = 4\pi \sigma ^2D$

where $\sigma$ is the mean squared parent-offspring dispersal distance along one axis in a two-dimensional habitat (assuming a Gaussian dispersal function), and *D* is the ideal population density (i.e. the number of ideal individuals per unit area). A genetic neighborhood under isolation by distance forms a circle with a radius of 2$\sigma$ (Figure 1A below), reflecting the local area within which gene flow is high relative to drift (i.e. the neighborhood, rather than the global population, approaches panmixia). An example of a population meeting the assumptions of IBD might be a common species of tuna in the Pacific. Isolation by distance is generally evident when a species that can disperse efficiently is nearly continuously distributed in habitats that are essentially uniform, such that dispersal and mating probabilities are equal in all directions, and only increasing distance reduces the probability. 

We can explore the concept of IBD in R. In the sGD package, you'll find simulated microsatellite genotypes for a population of 576 simulated individuals arrayed in a regular grid. We can think of IBD as a raster grid where all cells have a value of 1. Lets plot our simulated individuals on this simple IBD landscape:


```r
library(maptools)
library(raster)

# read in location data and landscape
xy_data <- read.csv(system.file("extdata","sGD_demo_xy.csv",package="sGD"))
proj <- "+proj=utm +zone=10 +datum=NAD83"   
xy_points <- SpatialPoints(xy_data[,c(2,3)],proj4string=CRS(proj))
IBD_landscape <- raster(nrows=256, ncols=256, crs=proj, xmn = 0, ymn = 0, xmx = 25600, ymx = 25600, vals=1)
plot(IBD_landscape)
plot(xy_points,add=T,pch=20,col="black")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

The black dots represent individuals and the x and y axes are distances in meters. Now let's assume $\sigma$ for this population is 3000 meters, and that dispersal follows a Gaussian function. Under Wright's model, a genetic neighborhood in this population is a circle with a radius of 2$\sigma$, or 6000 meters. Let's pick two individuals (1 and 2) in this population and visualize the neighborhoods surrounding them:


```r
pt1 <- xy_points[255]
pt2 <- xy_points[324]
pt1_NH <- buffer(pt1,6000)
pt2_NH <- buffer(pt2,6000)
pt1_NH <- spChFIDs(pt1_NH,"pt1") 
pt2_NH <- spChFIDs(pt2_NH,"pt2") 
pts_NH <- spRbind(pt1_NH,pt2_NH)
plot(IBD_landscape)
plot(pts_NH,add=T)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

Notice that the genetic neighborhoods of these two individuals overlap. That means individuals located where the overlap occurs have some probability of mating with either individual 1 or 2. By having overlapping neighborhoods, each of which defines a local breeding pool, we can capture the continuous nature of the population. Let's see how many individuals are in these two neighborhoods:


```r
NH1_N <- nrow(crop(xy_points,pt1_NH)@coords)
NH2_N <- nrow(crop(xy_points,pt2_NH)@coords)
NH1_N
```

```
## [1] 101
```

```r
NH2_N
```

```
## [1] 101
```

Notice how both neighborhoods have the same number of individuals. That's exactly what you'd expect when the neighborhoods are circles of the same size and the population is uniformly distributed in a regular grid. But notice what happens in neighborhoods located on the edge of a population:


```r
pt3 <- xy_points[12]
pt3_NH <- buffer(pt3,6000)
plot(IBD_landscape)
plot(xy_points,add=T,pch=20,col="black")
plot(pt3_NH,add=T)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

```r
NH3_N <- nrow(crop(xy_points,pt3_NH)@coords)
NH3_N
```

```
## [1] 56
```

This demonstrates how genetic neighborhoods near the periphery of an IBD population have fewer individuals. We can use sGD to calculate *N*, number of individuals (and here we're talking about the actual number *N~c~*, not the effective number *N~e~*) in each neighborhood, but first we need to calculate a pairwise Euclidean distance matrix between all individuals:


```r
library(sGD)
IBD_dist <- distmat(xy_points,method="ed")
IBD_dist [1:10,1:10]
```

```
##          [,1]    [,2]    [,3]    [,4]    [,5]    [,6]    [,7]    [,8]
##  [1,]    0.00 1043.48 2086.96 3130.43 4173.91 5217.39 6260.87 7304.35
##  [2,] 1043.48    0.00 1043.48 2086.96 3130.43 4173.91 5217.39 6260.87
##  [3,] 2086.96 1043.48    0.00 1043.48 2086.96 3130.43 4173.91 5217.39
##  [4,] 3130.43 2086.96 1043.48    0.00 1043.48 2086.96 3130.43 4173.91
##  [5,] 4173.91 3130.43 2086.96 1043.48    0.00 1043.48 2086.96 3130.43
##  [6,] 5217.39 4173.91 3130.43 2086.96 1043.48    0.00 1043.48 2086.96
##  [7,] 6260.87 5217.39 4173.91 3130.43 2086.96 1043.48    0.00 1043.48
##  [8,] 7304.35 6260.87 5217.39 4173.91 3130.43 2086.96 1043.48    0.00
##  [9,] 8347.83 7304.35 6260.87 5217.39 4173.91 3130.43 2086.96 1043.48
## [10,] 9391.30 8347.83 7304.35 6260.87 5217.39 4173.91 3130.43 2086.96
##          [,9]   [,10]
##  [1,] 8347.83 9391.30
##  [2,] 7304.35 8347.83
##  [3,] 6260.87 7304.35
##  [4,] 5217.39 6260.87
##  [5,] 4173.91 5217.39
##  [6,] 3130.43 4173.91
##  [7,] 2086.96 3130.43
##  [8,] 1043.48 2086.96
##  [9,]    0.00 1043.48
## [10,] 1043.48    0.00
```

This shows you just the pairwise distances between the first 10 individuals in the population, but you get the idea. sGD uses this information to determine the membership of neighborhoods surrounding every individual, based on a threshold distance we'll call 'radius'. We also need to provide sGD with the genotypes from each location, imported as a genind object by the adegenet package. Here we'll use genotypes simulated by the program CDPOP based on the same population of 576 individuals we've already described. The genotypes are in genepop format, with a 3 digit code for alleles (note the number of digits needs to be specified in the read.genepop function). Finally, we need to specify a minimum sample (min_N) size per neighborhood required to produce an estimate. This avoids calculating genetic indices from too few samples and producing unreliable estimates.


```r
IBD_genepop_file <- system.file("extdata","sGD_demo_IBD.gen",package="sGD") 
library(adegenet)
IBD_genind_obj = read.genepop(IBD_genepop_file,ncode = 3L,quiet=T)
radius <- 6000 # for this demo, units are in meters 6000
IBD_GD <- sGD(IBD_genind_obj,xy_data,IBD_dist,radius,min_N=20,metrics="GD")
```

```
## Reading input files...
## Input summary:
## 	 individuals: 576 
## 	 loci: 20 
## 	 minimum sample size: 20 
## 	 maximum sample size: NA 
## 	 neighborhood radius: 6000 
## Determining neighborhood membership from dist.mat and NH_radius...
## Calculating genetic diversity indices...
```

```r
library(ggplot2)

ggplot()  +
geom_point(data=IBD_GD, aes(x=X, y=Y,color=N),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

This clearly shows the gradients in local *N~c~*, with fewer individuals near the periphery. However, genetic diversity indices and local effective population size (i.e. Wright's *NS*) are usually of more interest than *N~c~*. Let's look at the percentage of alleles from the total population that are present in each neighborhood (this is the Ap variable in the IBD_sGD data frame):


```r
ggplot()  +
geom_point(data=IBD_GD, aes(x=X, y=Y,color=Ap),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

In this simple landscape, the pattern of allele usage looks very similar to *N*. It's not exactly like *N* because the simulated genotypes used to produce this pattern were  We'll see soon how such patterns can vary spatially in complex landscapes. For now, you should explore genetic diversity indices calculated by sGD using the code above but changing the color in the ggplots to match the metric of interest. Options include the number of alleles (A), allelic richness (Ar), expected heterozygosity (Hs), observed heterozygosity (Ho), and inbreeding coefficient (FIS).

Now let's use sGD to calculate Wright's *NS* in these local breeding neighborhoods. To do this, we need to point sGD towards the NeEstimator 2.0 executable directory. On my system, the path is "C:/NeEstimator" but you can vary this depending on your installation. We also need to set *NS_ans* to TRUE. Note that in the code below, I have turned off the genetic diversity metric calculations by setting *GD_ans* to FALSE. If you don't need to calculate both types of output, turning one off will speed up the processing time. 



```r
NeEstimator_dir <- "C:/NeEstimator_2.01"
IBD_NS <- sGD(IBD_genind_obj,xy_data,IBD_dist,radius,min_N=20,metrics="NS",NeEstimator_dir=NeEstimator_dir)
```

```
## Reading input files...
## Input summary:
## 	 individuals: 576 
## 	 loci: 20 
## 	 minimum sample size: 20 
## 	 maximum sample size: NA 
## 	 neighborhood radius: 6000 
## Determining neighborhood membership from dist.mat and NH_radius...
## Calculating NS...
```

```r
ggplot()  +
geom_point(data=IBD_NS, aes(x=X, y=Y,color=NS_ex0.02),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

You can see how some neighborhoods in the center of the landscape have a much larger *NS* than the rest. This is an expected property of a population of uniform density. If we explored a population of variable density, the results could be very different. Note that in the *IBD_NS* output data frame, there are several variables that quantify *NS*. They differ in the amount they remove rare alleles from the analysis prior to inferring *NS*. See Waples & Do (2008) for a full explanation. The idea is that rare alleles bias the estimate, and NeEstimator allows the user to removing alleles below a certain frequency. The outputs of *NS_ex0*, *NS_ex0.02*, *NS_ex0.05*, and *NS_ex0.10* correspond to no removal of rare alleles, or removing alleles with a frequency below 0.02, 0.05, and 0.01 respectively.

Now let's consider a more complex landscape, where landscape heterogeneity creates spatial variation in how the landscape resists movement and gene flow. The American marten (*Martes americana*) is an example of this dynamic, because this small forest carnivore is vulnerable to predation in open habitats. Marten prefer to move through high canopy cover habitats and in places where there are snags and coarse woody debris that allow concealment and escape from predators. Resistance models are used to quantify landscape resistance to movement and gene flow. Let's look at the one included with the sGD package:


```r
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
```

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png)

In our example, the dark patches of low resistance in this landscape might be areas of high canopy cover, and the white areas of high resistance might be clearcuts or recent burns. These landscapes reflect the concept of isolation by resistance, where mating probabilities are driven not by Euclidean distance, but by the effective or ecological distances between individuals, given the resistance in the landscape. We can still apply the neighborhood concepts discussed already to these complex landscapes for the purposes of calculating genetic diversity indices in continuous populations. We just have to make one adjustment. Rather than define neighborhoods by Euclidean distances, we'll do so using cost-weighted distance given the resistance model. We describe this concept more fully in Shirk et al. (2011) and Shirk & Cushman (2014). Let's calculate a cost-weighted distance matrix with sGD:


```r
IBR_dist <- distmat(xy_points,method="cd",landscape=IBR_landscape)
IBR_dist [1:10,1:10]
```

```
##           [,1]     [,2]     [,3]     [,4]     [,5]     [,6]     [,7]
##  [1,]     0.00  1182.84  2757.11  7200.61 10095.58 16050.36 23291.78
##  [2,]  1182.84     0.00  1574.26  6076.35  9016.30 14971.07 22212.49
##  [3,]  2757.11  1574.26     0.00  5026.35  8197.67 14152.44 21393.86
##  [4,]  7200.61  6076.35  5026.35     0.00  4550.61 10505.38 17746.80
##  [5,] 10095.58  9016.30  8197.67  4550.61     0.00  7750.00 15571.93
##  [6,] 16050.36 14971.07 14152.44 10505.38  7750.00     0.00  8018.38
##  [7,] 23291.78 22212.49 21393.86 17746.80 15571.93  8018.38     0.00
##  [8,] 32429.65 31350.36 30531.73 26884.67 24709.80 17156.24  9800.00
##  [9,] 39438.23 38358.94 37540.31 33893.25 31718.38 24164.82 17292.03
## [10,] 48330.87 47251.58 46432.95 42785.89 40611.02 33057.46 26184.67
##           [,8]     [,9]    [,10]
##  [1,] 32429.65 39438.23 48330.87
##  [2,] 31350.36 38358.94 47251.58
##  [3,] 30531.73 37540.31 46432.95
##  [4,] 26884.67 33893.25 42785.89
##  [5,] 24709.80 31718.38 40611.02
##  [6,] 17156.24 24164.82 33057.46
##  [7,]  9800.00 17292.03 26184.67
##  [8,]     0.00  7629.90 16522.54
##  [9,]  7629.90     0.00  9333.45
## [10,] 16522.54  9333.45     0.00
```

Notice how these distances are larger than the Euclidean distances in the IBD_dist matrix we looked at earlier. That's because these numbers reflect *effective* distances given the resistance model. A genetic neighborhood based on effective distances in a complex landscape isn't a perfect circle. Rather, it's a complex shape. Let's assume the marten population in our example has a mean parent-offspring dispersal distance ($\sigma$) of 8000m in cost-weighted distance (when we talk about distance under IBR, they need to be in cost-distance units). So 2$\sigma$ is 16000m. Let's look at the shape and size of a single neighborhood, given this resistance surface and the same uniform population of 576 individuals. First, we'll calculate the cost-weighted distance to pt1 which we used above when considering IBD neighborhoods:


```r
# calculate transition surface and do geocorrection
library(gdistance)
tr <- transition(IBR_landscape,transitionFunction = function(x) {1/mean(x)},directions=8) 
trCorrC<-geoCorrection(tr,type="c",multpl=FALSE,scl=FALSE)
cd_pt1 <- accCost(trCorrC, pt1)
plot(cd_pt1)
plot(pt1,col="red",pch=20,add=T)
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

New lets threshold that cost-weighted distance at 16000m to define the limits of the neighborhood surrounding that point:


```r
reclass_mat <- matrix(c(0,16000,1,16000,100000,NA),byrow=T,nrow=2)
cd_pt1_bin <- reclassify(cd_pt1,reclass_mat)
plot(cd_pt1_bin,col="black")
plot(pt1,col="red",pch=20,add=T)
```

![plot of chunk unnamed-chunk-14](figure/unnamed-chunk-14-1.png)

Notice the irregular shape of the neighborhood. In ecological space, this neighborhood is still a circle (i.e., the outer edge is still 2$\sigma$). By defining the neighborhood and the landscape distance in cost-weighted units, we can now use sGD to calculate genetic diversity indices and Wright's *NS*, based on a CDPOP simulation of our 576 individuals on this resistance landscape. First let's look at the pattern of Ap on the landscape:


```r
IBR_genepop_file <- system.file("extdata","sGD_demo_IBR.gen",package="sGD") 
IBR_genind_obj <- read.genepop(IBR_genepop_file,ncode = 3L,quiet=T)
radius <- 16000 # for this demo, units are in meters 6000
IBR_sGD = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,metrics=c("GD","NS"),NeEstimator_dir=NeEstimator_dir)
```

```
## Reading input files...
## Input summary:
## 	 individuals: 576 
## 	 loci: 20 
## 	 minimum sample size: 20 
## 	 maximum sample size: NA 
## 	 neighborhood radius: 16000 
## Determining neighborhood membership from dist.mat and NH_radius...
## Calculating genetic diversity indices...
## Calculating NS...
```

```r
ggplot()  +
geom_point(data=IBR_sGD, aes(x=X, y=Y,color=Ap),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

Notice the complex pattern where the most central areas with high allelic diversity is driven by the resistance model. Also, note that the black circles represent neighborhoods where the minimum sample size (here set to N=20) was not met, so no estimates were made. Now let's look at the distribution of *NS* estimates:


```r
boxplot(IBR_sGD$NS_ex0.02)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)

Notice how there is a major outlier with an estimated *NS* > 1000. Occasionally, the algorithm NeEstimator uses to estimate *NS* produces very large estimates, and it's important to recognize these are unlikely, particularly when they are surrounded by neighborhoods with much lower estimates. Let's remove that outlier, and then look at the pattern of *NS* on the landscape:


```r
# get rid of an outlier
IBR_sGD$NS_ex0.02[which(IBR_sGD$NS_ex0.02>1000)] = NA

ggplot()  +
geom_point(data=IBR_sGD, aes(x=X, y=Y,color=NS_ex0.02),size=5) + 
scale_color_gradient(low="blue", high="green",na.value = "black") + 
theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank())
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

Again, we see a complex pattern that corresponds to the pattern of resistance in the landscape. At this point, you should explore other genetic diversity indices in the IBD_sGD data frame and see how their patterns play out on this landscape. Notice that, with the added resistance in the landscape, the genetic patterns tend to be stronger than what we observed earlier for the IBD population, particularly for genetic indices of heterozygosity.

One thing we haven't considered yet is what to do when you don't know 2$\sigma$. In simulations, dispersal distances can be recorded and used to calculate $\sigma$. In real populations, however, it is more difficult to estimate this parameter. Ideally, one could infer it from the genetic data directly. One idea suggested in Waples et al(2013) was to set the neighborhood radius to the size that minimizes FIS departures from zero (i.e. minimize both positive and negative inbreeding coefficients). A non-zero FIS value is an indicator that the local neighborhood is not in Hardy-Weinberg equilibrium, so if FIS in neighborhoods approaches zero, one might assume the neighborhood is sized appropriately to meet the assumptions of a population required before calculating genetic diversity indices and *NS*. In Shirk & Cushman (2014), we explored this idea and found it to be a useful concept but not terribly precise in estimating the true value of 2$\sigma$. We are currently exploring new ideas for this, and hope to update sGD soon with alternative approaches. For now, however, we recommend iteratively modifying the radius of the neighborhood until the distribution of FIS values is centered on zero. We already know the IBR population has a radius 2$\sigma$ = 16000m from our simulations. This produced the following distribution of FIS values:


```r
hist(IBR_sGD$FIS)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)

```r
median(na.omit(IBR_sGD$FIS))
```

```
## [1] -0.006
```

Notice the distribution is centered on zero and the median FIS value is very close to zero. Now let's add 50% more to the radius and look at the distribution (this will take several minutes, depending on your computer, as the larger the radius, the more computations are made):


```r
radius <- 24000
IBR_sGD_24k = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,metrics="GD")
```

```
## Reading input files...
## Input summary:
## 	 individuals: 576 
## 	 loci: 20 
## 	 minimum sample size: 20 
## 	 maximum sample size: NA 
## 	 neighborhood radius: 24000 
## Determining neighborhood membership from dist.mat and NH_radius...
## Calculating genetic diversity indices...
```

```r
hist(IBR_sGD_24k$FIS)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

```r
median(na.omit(IBR_sGD_24k$FIS))
```

```
## [1] 0.0029
```

See how the distribution of FIS values shifted to the right (positive values) and the median FIS moved almost 5x farther from zero? This is exactly what you expect due to a Wahlund effect that occurs when the neighborhood extent is larger than the local extent of breeding. Now let's reduce the radius to 10000m and recalculate FIS:


```r
radius <- 10000
IBR_sGD_10k = sGD(IBR_genind_obj,xy_data,IBR_dist,radius,min_N=20,metrics="GD")
```

```
## Reading input files...
## Input summary:
## 	 individuals: 576 
## 	 loci: 20 
## 	 minimum sample size: 20 
## 	 maximum sample size: NA 
## 	 neighborhood radius: 10000 
## Determining neighborhood membership from dist.mat and NH_radius...
## Calculating genetic diversity indices...
```

```r
hist(IBR_sGD_10k$FIS)
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

```r
median(na.omit(IBR_sGD_10k$FIS))
```

```
## [1] -0.01035
```

Now the distribution has shifted to the left (negative FIS values), in what amounts to an inverse Wahlund effect. Repeating this evaluation over a range of radius values can provide a crude means of inferring $\sigma$ when it is not known empirically or from simulations.

By now, you should have an understanding of what sGD is designed to do and the basics of using it in R. Make sure you also explore the help sections for sGD commands (e.g. type ?sGD or ?distmat) to learn more about the inputs, parameters, and outputs.

# Literature Cited and Recommended Reading

Amos, and J. Harwood (1998). Factors affecting levels of genetic diversity in natural populations. Philosophical Transactions of the Royal Society of London Series B-Biological Sciences 353, 177-186.

G. Caughley (1994). Directions in conservation biology. Journal of animal ecology, 215-244.

C.C. Cockerham, and B.S. Weir (1977). Digenic descent measures for finite populations. Genet. Res 30, 121-147.

S.A. Cushman, and E.L. Landguth (2012). Multi-taxa population connectivity in the Northern Rocky Mountains. Ecological Modelling 231, 101-112.

S.A. Cushman, J.S. Lewis, and E.L. Landguth (2013a). Evaluating the intersection of a regional wildlife connectivity network with highways. Movement Ecology 1, 12.

S.A. Cushman, K.S. McKelvey, J. Hayden, and M.K. Schwartz (2006). Gene flow in complex landscapes: Testing multiple hypotheses with causal modeling. American Naturalist 168, 486-499.

S.A. Cushman, T.N. Wasserman, E.L. Landguth, and A.J. Shirk (2013b). Re-Evaluating Causal Modeling with Mantel Tests in Landscape Genetics. Diversity 5, 51-72.

C.G. Eckert, K.E. Samis, and S.C. Lougheed (2008). Genetic variation across species’ geographical ranges: the central–marginal hypothesis and beyond. Molecular Ecology 17, 1170-1188.

R. Frankham, D.A. Briscoe, and J.D. Ballou (2002). Introduction to conservation genetics. Cambridge University Press.

G. García-Ramos, and M. Kirkpatrick (1997). Genetic models of adaptation and gene flow in peripheral populations. Evolution, 21-28.

R.H. Gardner (1999). Rule: map generation and a spatial analysis program.

S.C. Goslee, and D.L. Urban (2007). The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22, 1-19.

W.G. Hill (1981). Estimation of effective population size from data on linkage disequilibrium. Genetical Research 38, 209-216.

T. Jombart (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics 24, 1403-1405.

R. Lande, and G.F. Barrowclough (1987). Effective population size, genetic variation, and their use in population management. Viable populations for conservation, 87-123.

E.L. Landguth, and S.A. Cushman (2010). CDPOP: A spatially explicit cost distance population genetics program. Molecular Ecology Resources 10, 156-161.

R. Leblois, F. Rousset, and A. Estoup (2004). Influence of spatial and temporal heterogeneities on the estimation of demographic parameters in a continuous population using individual microsatellite data. Genetics 166, 1081-1092.

B.H. McRae (2006). Isolation by resistance. Evolution 60, 1551-1561.

M.P. Miller (2005). Alleles In Space (AIS): computer software for the joint analysis of interindividual spatial and genetic information. Journal of Heredity 96, 722-724.

M.C. Neel, K.S. McKelvey, N. Ryman, M.W. Lloyd, R. Short Bull, F.W. Allendorf, M.K. Schwartz, and R.S. Waples (2013). "Estimation of effective population size in continuously distributed populations: There goes the neighborhood.," in,  (Heredity), (in press).

R Development Core Team (2013). "R: A language and environment for statistical computing. R Foundation for Statistical Computing," in,  (Vienna, Austria: R Development Core Team).

F.J. Rohlf, and G.D. Schnell (1971). An investigation of the isolation-by-distance model. American Naturalist, 295-324.

A.J. Shirk, and S.A. Cushman (2011). sGD: software for estimating spatially explicit indices of genetic diversity. Molecular Ecology Resources 11, 922-934.

A.J. Shirk, D.O. Wallin, S.A. Cushman, C.G. Rice, and K.I. Warheit (2010). Inferring landscape effects on gene flow: a new model selection framework. Molecular Ecology 19, 3603-3619.

A.J. Shirk, S.A. Cushman (2014). Spatially-explicit estimation of Wright's neighborhood size in continuous populations. Frontiers in Ecology and Evolution, 2, 62.

P. Sinnock (1975). The Wahlund effect for the two-locus model. American Naturalist, 565-570.

J. van Etten (2011). "gdistance: distances and routes on geographical grids. R package version 1.1–2," in.

S. Wahlund (1928). The combination of populations and the appearance of correlation examined from the standpoint of the study of heredity. Hereditas 11, 65-106.

R.S. Waples, and C. Do (2008). LDNE: a program for estimating effective population size from data on linkage disequilibrium. Molecular Ecology Resources 8, 753-756.

R.S. Waples, and P.R. England (2011). Estimating Contemporary Effective Population Size on the Basis of Linkage Disequilibrium in the Face of Migration. Genetics 189, 633-644.

T. Wasserman, S. Cushman, A. Shirk, E. Landguth, and J. Littell (2012a). Simulating the effects of climate change on population connectivity of American marten (Martes americana) in the northern Rocky Mountains, USA. Landscape ecology 27, 211-225.

T.N. Wasserman, S.A. Cushman, J.S. Littell, A.J. Shirk, and E.L. Landguth (2012b). Population connectivity and genetic diversity of American marten (Martes americana) in the United States northern Rocky Mountains in a climate change context. Conservation Genetics, 1-13.

B.S. Weir (1996). Genetic data analysis II: methods for discrete population genetic data. Sunderland, MA: Sinauer Associates, 91-132.

S. Wright (1931). Evolution in Mendelian populations. Genetics 16, 97.

S. Wright (1943). Isolation by distance. Genetics 28, 114-138.

S. Wright (1946). Isolation by distance under diverse systems of mating. Genetics 31, 39.

