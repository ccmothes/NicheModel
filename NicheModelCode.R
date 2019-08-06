
# Niche modeling in Maxent over many species ------------------------------



library(dismo)
library(rgdal)
library(maptools)
library(maps)
library(rJava)
library(cluster)
library(raster)
library(geosphere)
library(dplyr)
library(rworldmap)


# Create unique geographic extent for each species --------

##presence data for all species
data<- read.csv("Data.csv")

data.occ.sp <- SpatialPointsDataFrame(data[,c('long','lat'),], data = data)


##create list of unique species names
specs <-unique(data.occ.sp$Species)

save(specs, file = "species_names_vector.RData")


##calculate buffer distance for each species by measuring the euclidean distance
##between the two most spatially segregated clusters of points

dist <- list()
for (i in 1:length(unique(data.occ.sp$Species))){
  s <- subset(data.occ.sp, data.occ.sp$Species == unique(data.occ.sp$Species)[i])
  d <- distm(s, fun = distHaversine)
  ag <-agnes(d, diss = TRUE, metric = "euclidean", method = "single")
  dist[[i]] <- max(ag$height)/1000
  print(i)
}
save(dist, file = "native_range_buffers.RData")

buffer_dist <- data.frame(species = specs[1:31], distance_km = as.numeric(paste(dist[1:31])))
##only use half the distance to make geo extent so make new column with distance/2
buffer_dist <- mutate(buffer_dist, distance_km_half = distance_km/2)
write.csv(buffer_dist, file = "buffer_distances.csv")

##create buffer around native range localities for each species to make geographic extent polygons

range.output = list()
for (i in 1:length(unique(data.occ.sp$Species))){
  s <- subset(data.occ.sp, data.occ.sp$Species == unique(data.occ.sp$Species)[i])
  x <- circles(s@coords, d = buffer_dist$distance_km_half[i], lonlat = T)
  s.range <- gUnionCascaded(x@polygons)
  range.output[[i]]<- s.range
  print(i)
}


##crop any polygons with ranges extending into Florida

bhs<-countriesCoarse[which(countriesCoarse$ISO3 == 'BHS'),]
cub<-countriesCoarse[which(countriesCoarse$ISO3 == 'CUB'),]
##Anolis porcatus
range.output[[7]]<-crop(range.output[[7]], cub)
##Anolis distichus
extent<-c(-79.5, -67.5, 17.5, 28)
range.output[[10]]<-crop(range.output[[10]], extent)
##Anolis carinatus
map1<-gUnion(bhs, cub)
range.output[[14]]<-crop(range.output[[14]], map1)
###Tarentola annularis
extent2<- c(-20, 62, -36, 37)
range.output[[25]]<-crop(range.output[[25]], extent2)


save(range.output, file='bufferedRanges.RData')


##Load in background points (squamate occurences)
backgroundPoints <-read.csv("backgroundPoints.csv")
bg <- SpatialPointsDataFrame(backgroundPoints[,c('Long','Lat'),], data = backgroundPoints)



# Clip background points to geographic extent for each species ---------------

bg.output<-list()
for (i in 1:length(unique(data.occ.sp$Species))){
bg.output[[i]]<-crop(bg, y=range.output[[i]], snap="out")
print(i)
}

save(bg.output, file='backgroundPoints.RData')


# Download predictor variables and clip to geographic extent for each species --------


##Load predictor variables from WorldClim using the raster package
predictors <- getData('worldclim', var = 'bio', res=10)

###NOTE: for this demonstration, 2.5 minute resolution was used for reproducibility
###when using 30 second resoltion (as in Mothes et al.) global files were downloaded directly 
###from the website to an external hardrive due to large file size

##Clip predictors to range of each species

predictors.ranges<-list()
for (i in 1:length(range.output)){
	ext <- extent(range.output[[i]])
	predictors.crop <- crop(predictors, y = ext, snap="out") 
	out <- raster(predictors.crop)
	cropped <- setValues(out, NA) 
	range.r <- rasterize(range.output[[i]], cropped) 
	predictors.ranges[[i]] <- mask(x=predictors.crop, mask=range.r)
	writeRaster(predictors.ranges[[i]], filename = paste(specs[i], "pred",sep = "_"), overwrite=T)
	print(i)
}

save(predictors.ranges, file = "predictors.ranges.RData")

##clip predictors to Florida

fl <- readOGR("FL_Counties_Project.shp")
#set projection of shapefile to match that of the predictors
new.prj <- proj4string(predictors)
fl.prj <- spTransform(fl, new.prj)

c <- crop(predictors, extent(fl.prj), snap = 'out')
f <- rasterize(fl.prj, c)
predictors.fl <- mask(x = c, mask = f, filename = "predictors_fl.tif")
names(predictors.fl) <- names(predictors)

# Run maxent for all species ----------------------------------------------


#create output directories for each species

dir.create("maxent_AIC")
setwd("maxent_AIC")
for (i in 1:length(specs)){
  dir.create(paste(specs[i]))
}

folders <- list.dirs(path = getwd())


##make sequence of betas to use
betas <- array(dim = c(24, 2))
betas[,1] <- "betamultiplier="
betas[,2] <- c(seq(0.2,1, 0.2), seq(2, 20, 1))
betas <- data.frame(betas)
betas[,3] <- paste(betas[,1], betas[,2],sep = '')
names(betas)[3] <- 'setbeta'
betas[,2] <- as.character(betas[,2])


##build maxent models with all 19 layers at each value of Beta, calculate AICc

names(bg.output) <- specs
names(range.output) <- specs


for (i in 1:length(specs)){
  setwd(folders[grep(specs[i],folders)][1])  
  

  occ <- subset(data.occ.sp, data.occ.sp$Species == specs[i])
  occ <- coordinates(occ)
  
  background = bg.output[[i]]
  background <- coordinates(background)  
  
  allbetas <- data.frame(array(dim = c(24, 7)))
  names(allbetas) <- c('species','beta','no.params','no.pres','loglik','AIC','AICc')
  allbetas[,1] <- specs[i]
  
  for (k in 1:nrow(betas)){
    
 
	xm <- maxent(x = predictors.ranges[[i]], p = occ, a = background, args = c(betas[k,"setbeta"],"responsecurves=TRUE","writebackgroundpredictions=TRUE"), path = folders[grep(specs[i],folders)][1])
    
    save(xm, file = paste(specs[i], betas[k,"X2"], "xm.RData", sep = '_'))
    
    #get lambdas info out of xm
    lambdas <- data.frame(strsplit(xm@lambdas, ','))
    lambdas <- t(lambdas)
    lambdas <- data.frame(lambdas)
    rownames(lambdas) <- 1:nrow(lambdas)
    #vitalconstants <- lambdas[(nrow(lambdas)-3):nrow(lambdas),1:2]
    lambdas <- lambdas[-((nrow(lambdas)-3):nrow(lambdas)),]
    colnames(lambdas) <- c('variable','lambda_estimate','min','max')
    lambdas$lambda_estimate <- as.numeric(as.character(lambdas$lambda_estimate))
    
    #read and write predictions at presences
    pres.raw <- read.csv('species_samplePredictions.csv')
    write.csv(pres.raw,paste(specs[i], betas[k,"X2"], "pres.predictions.csv", sep = '_') )
    
    #calculate AICc
    no.params <- length(which(lambdas$lambda_estimate != 0))
    no.pres <- xm@results[1]
    loglik <- sum(log(pres.raw$Raw.prediction))
    AIC <- (2*no.params) -(2*loglik)
    AICc <- AIC + (((2*no.params)*(no.params+1))/(no.pres-no.params-1))
    
    
    allbetas[k,2:7] <-c(betas[k,"X2"], no.params, no.pres, loglik, AIC, AICc)

    
  
    print(k)
  }
  
  write.csv(allbetas, file = paste(specs[i],"allbetas.csv", sep="_"))
  print(i)
}


#make output data frame for beta results
bestbetas <- data.frame(array(dim = c(31, 5)))
names(bestbetas) <- c('species','no.pres','no.params default','no.params best','bestbeta')


for (i in 1:length(specs)){
  #load beta file
  inpath <- folders[grep(specs[i],folders)][1]
  f <- list.files(path = inpath, pattern = 'allbetas')
  setwd(inpath)
  data <- read.csv(f)
  
  
  bestbetas[i,1] <- as.character(data$species[1])
  bestbetas[i,2] <- data$no.pres[1]
  bestbetas[i,3] <- data[which(data$beta == 1),'no.params']
  
  if(max(data$no.pres) > min(data$no.params))
    bestbetas[i,5] <- min(data[which(data$AICc == min(data[which(data$no.params < data$no.pres&data$no.params >0),'AICc'])),'beta'])
  
  if((max(data$no.pres) > min(data$no.params)))
    bestbetas[i,4] <- data[which(data$beta == bestbetas[i,5]),'no.params']
  print(i)
}

setwd("..")
save(bestbetas, file = 'bestbetas.RData')
write.csv(bestbetas, file = "bestbetas.csv")

#make plots to check behavior of likelihoods

for (i in 1:length(specs)){
  #load beta file
  inpath <- folders[grep(specs[i],folders)][1]
  f <- list.files(path = inpath, pattern = 'allbetas')
  setwd(inpath)
  data <- read.csv(f)
  

  
  pdf(file = paste(specs[i],"betaplots.pdf",sep = '_'))
  par(mfrow = c(1,2))
  plot(data$no.params, data$loglik, main = specs[i], ylab = "LogLik", xlab = 'No. Parameters', sub = paste("Best Beta =",bestbetas$bestbeta[i], sep = ' '), type = 'o')
  plot(data$no.params, data$AICc, main = specs[i], ylab = "AICc", xlab = 'No. Parameters', sub = paste("Best Beta =",bestbetas$bestbeta[i], sep = ' '), type = 'o')
  dev.off()
  
  
}


##Now rerun maxent model for each species using the best beta value to get best model results
###create new directories for best beta results
dir.create("bestbeta")
setwd("bestbeta")
for (i in 1:length(specs)){
  dir.create(paste(specs[i]))
}

output_folders3 <- list.dirs(path = getwd())
output_folders3 <- output_folders3[-1]


for (i in 1:length(specs)){
  setwd(output_folders3[grep(specs[i],output_folders3)][1])  
  
  occ <- subset(data.occ.sp, data.occ.sp$Species == specs[i])
  occ <- coordinates(occ)
  
  background = bg.output[[i]]
  background <- coordinates(background)  
  
  xm <- maxent(x = predictors.ranges[[i]], p = occ, a = background, args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"outputformat=logistic", "responsecurves=true", "replicates=10"), path = output_folders3[grep(specs[i],output_folders3)])
  
  
  save(xm, file = paste(specs[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  #predict model onto native range with logistic output and save raster file
  native.range <- predict(xm, predictors.ranges[[i]], filename = paste(specs[i],'native.tif',sep = '_'), overwrite =T)
  
  ##average rasters from each replicate to get single suitability map
  f <- list.files(path = getwd(), pattern = 'native')
  ras <- lapply(f, raster)
  STACK1 <- stack(ras) 
  mean <- stackApply(STACK1, indices = rep(1,nlayers(STACK1)), fun = "mean", na.rm = T)
  writeRaster( x = mean, filename = paste(specs[i], 'native_range_suitability_map.tif', sep = '_'), overwrite = TRUE)
  
  #predict model onto florida with logistic output
  florida <- predict(xm, predictors.fl, filename = paste(specs[i],'florida.tif',sep = '_'), overwrite = TRUE, progress = 'text')
  
  ##average raster for florida for each replicate
  f <- list.files(path = getwd(), pattern = 'florida')
  ras <- lapply(f, raster)
  STACK1 <- stack(ras) 
  mean <- stackApply(STACK1, indices = rep(1,nlayers(STACK1)), fun = "mean", na.rm = T)
  writeRaster( x = mean, filename = paste(specs[i], 'florida_suitability_map.tif', sep = '_'), overwrite = TRUE)
  
  
  
  print(i)
}


