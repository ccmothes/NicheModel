###Script for Niche Modeling in Maxent Over Many Species##
###Original script created by Amber N. Wright
###Modified by Caitlin C. Mothes


library(dismo)
library(geovec)
library(rgdal)
library(maptools)
library(maps)
library(rJava)
library(cluster)
library(raster)

setwd("C:/Users/lab/Desktop/Niche_Project")

##presence data for all species
data<- read.csv("C:/Users/Caitlin/Desktop/NicheProject/NicheModel/Data.csv")

data.occ.sp <- SpatialPointsDataFrame(data[,c('long','lat'),], data = data)


##list of species
Species<- read.csv ("C:/Users/Caitlin/Downloads/Niche Project/Niche Project/NicheModel/Species.csv", head=F)

specs <-unique(data$Species)

save(specs, file = "species_names_vector.RData")


##calculate buffer distance for each species
dist<-agnes(data, metric="euclidean", method="single")

####OR just load file with buffer distances for each species
distances<-read.csv("C:/Users/Caitlin/Downloads/Niche Project/Niche Project/NicheModel/buffer.distances.csv")


##create buffer around native range localities for each species
dist.output=list()
	for(i in 1:length (unique(distances$Species))){
	dist<-distances$distance_m
	dist.output[[i]]<-dist
}

dist<-as.vector(dist.output[[i]], mode="numeric")

range.output = list()
for (i in 1:length(unique(data.occ.sp$Species))){
  s <- subset(data.occ.sp, data.occ.sp$Species == unique(data.occ.sp$Species)[i])
  x <- circles(s@coords, d = dist[i], lonlat = T)
  s.range <- gUnionCascaded(x@polygons)
  range.output[[i]]<- s.range
  print(i)
}


##crop any polygons with ranges extending into Florida

save(range.output, file='bufferedRanges.RData')


##Load in background points (all lizard localities)
all.lizards<-read.csv("C:/Users/Caitlin/Downloads/Niche Project/Niche Project/NicheModel/LizardLocalities.csv")
lizards <- SpatialPointsDataFrame(all.lizards[,c('Long','Lat'),], data = all.lizards)
bg<-lizards


##clip bg points to range buffers for each species
bg.output<-list()
for (i in 1:length(unique(data.occ.sp$Species))){
bg.output[[i]]<-crop(bg, y=range.output[[i]], snap="out")
print(i)
}

save(bg.output, file='backgroundPoints.RData')


##Load predictor files (Bioclim layers)
outpath <- 'C:/Users/Caitlin/Desktop/NicheProject/NicheModel/bio'

f <- list.files(path=outpath, pattern='tif$')
setwd(outpath)
predictors <- stack(f)


##Clip predictors to range of each species

predictors.ranges<-list()
for (i in 1:31){
	ext <- extent(range.output[[i]])
	predictors.crop <- crop(predictors, y = ext, snap="out") 
	out <- raster(predictors.crop)
	cropped <- setValues(out, NA) 
	range.r <- rasterize(range.output[[1]], cropped) 
	predictors.masked <- mask(x=predictors.crop, mask=range.r)
	predictors.ranges_1<- predictors.masked
	writeRaster(predictors.ranges[[i]], filename = paste(specs[i], "pred",sep = "_"), overwrite=T)
	print(i)
}



###Start Modeling Process###

##Load in all data if haven't already done so

#presence data

data<- read.csv("C:/Users/lab/Desktop/Niche_Project/Data.csv")
data.occ.sp <- SpatialPointsDataFrame(data[,c('long','lat'),], data = data)


#background data
load("C:/Users/lab/Desktop/Niche_Project/backgroundPoints.RData")

#list of species
load("C:/Users/lab/Desktop/Niche_Project/species_names_vector.RData")

#load range polygons
load("C:/Users/lab/Desktop/Niche_Project/bufferedRanges.RData")


#climate predictors
outpath <- "C:/Users/lab/Desktop/Niche_Project"
predictors.ranges <- list()
for (i in 1:31){
  f <- list.files(path=outpath, pattern='pred.grd')
  setwd(outpath)
  predictors.ranges[[i]] <- stack(f[i])
  print(i)
}


setwd("C:/Users/lab/Desktop/Niche_Project/maxent_AIC")
outpath = "C:/Users/lab/Desktop/Niche_Project/maxent_AIC"

folders <- list.dirs(path = outpath)


#make sequence of betas to use
betas <- array(dim = c(25, 2))
betas[,1] <- "betamultiplier="
betas[,2] <- c(seq(0,1, 0.2), seq(2, 20, 1))
betas <- data.frame(betas)
betas[,3] <- paste(betas[,1], betas[,2],sep = '')
names(betas)[3] <- 'setbeta'
betas[,2] <- as.character(betas[,2])


###build maxent models with all 19 layers at each value of Beta, calculate AICc

names(bg.output) <- specs
names(range.output) <- specs


for (i in 1:31){
  setwd(folders[grep(specs[i],folders)][1])  
  

  occ <- subset(data.occ.sp, data.occ.sp$Species == specs[i])
  occ <- coordinates(occ)
  
  background = bg.output[[i]]
  background <- coordinates(background)  
  
  allbetas <- data.frame(array(dim = c(25, 7)))
  names(allbetas) <- c('species','beta','no.params','no.pres','loglik','AIC','AICc')
  allbetas[,1] <- specs[i]
  
  for (k in 1:nrow(betas)){
    
 
	xm <- maxent(x = predictors.ranges[[i]], p = occ, a = background, args = c(betas[k,"setbeta"],"outputformat=raw","writebackgroundpredictions=TRUE"), path = folders[grep(specs[i],folders)][1])
    
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
bestbetas <- data.frame(array(dim = c(153, 5)))
names(bestbetas) <- c('species','no.pres','no.params default','no.params best','bestbeta')

folders <- list.dirs(path = outpath)

for (i in 1:31){
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

setwd("C:/Users/lab/Desktop/Niche_Project/maxent_AIC")
save(bestbetas, file = 'bestbetas.RData')

#make plots to check behavior of likelihoods

for (i in 1:31){
  #load beta file
  inpath <- folders[grep(specs[i],folders)][1]
  f <- list.files(path = inpath, pattern = 'allbetas')
  setwd(inpath)
  data <- read.csv(f)
  
  setwd("C:/Users/lab/Desktop/Niche_Project/maxent_AIC")



  
  pdf(file = paste(specs[i],"betaplots.pdf",sep = '_'))
  par(mfrow = c(1,2))
  plot(data$no.params, data$loglik, main = specs[i], ylab = "LogLik", xlab = 'No. Parameters', sub = paste("Best Beta =",bestbetas$bestbeta[i], sep = ' '), type = 'o')
  plot(data$no.params, data$AICc, main = specs[i], ylab = "AICc", xlab = 'No. Parameters', sub = paste("Best Beta =",bestbetas$bestbeta[i], sep = ' '), type = 'o')
  dev.off()
  
  
}



load("C:/Users/lab/Desktop/Niche_Project/maxent_AIC/bestbetas.RData")

##okay, running with the best betas now
outpath = "C:/Users/lab/Desktop/Niche_Project/maxent_AIC/BestBetas"
output_folders3 <- list.dirs(path = outpath)
output_folders3 <- output_folders3[-1]

#re-run each species with all 19 variables and the best beta value. seems redundant with above (could just re-load the best xm from above) but this starts a new set of folders for each species with just this model. this way will have the correct html file. 

for (i in 1:31){
  setwd(folders[grep(specs[i],folders)][1])  
  
  occ <- subset(data.occ.sp, data.occ.sp$Species == specs[i])
  occ <- coordinates(occ)
  
  background = bg.output[[i]]
  background <- coordinates(background)  
  
  xm <- maxent(x = predictors.ranges[[i]], p = occ, a = background, args = c(betas[which(betas$X2 == bestbetas[i,"bestbeta"]),"setbeta"],"outputformat=logistic", "responsecurves=true"), path = output_folders3[grep(specs[i],output_folders3)])
  
  setwd(output_folders2[grep(specs[i],output_folders2)])
  
  save(xm, file = paste(specs[i], bestbetas[i,"bestbeta"], "xm.RData", sep = '_'))
  
  write.csv(data.frame(xm@results[7:44,]), file = paste(specs[i],bestbetas[i,"bestbeta"],'variable_scores.csv', sep = '_'))
  
  print(i)
}


setwd("C:/Users/lab/Desktop/Niche_Project/maxent_AIC")
load("C:/Users/lab/Desktop/Niche_Project/maxent_AIC/bestbetas.RData")
write.csv(bestbetas, 'bestbetas.csv')
allbetas_best <- read.csv('bestbetas.csv')

####now predict to natvie range and Florida from model with best beta value for current climate

load("C:/Users/lab/Desktop/Niche_Project/predictors.fl.RData")

outpath <- ("C:/Users/lab/Desktop/Niche_Project/maxent_AIC/BestBetas")

output_folders <- list.dirs(path = outpath)

for (i in 1:31){
  setwd(output_folders[grep(specs[i],output_folders)][1])
  
  #read best model from Best Beta Maxent Output with best layers
  f <- list.files(path = getwd(), pattern = 'xm')    
  bestmodel <- f[grep(paste("_",allbetas_best[i,'bestbeta'],"_", sep = '') , f, fixed = TRUE)]
  load(paste(getwd(),bestmodel, sep = '/'))
  
  #load occurrence and background data for plotting
  
  occ <- subset(data.occ.sp, data.occ.sp$Species == specs[i])
  occ <- coordinates(occ)
  
  background = bg.output[[i]]
  background <- coordinates(background)
  
  #save predictions to Best layers and Best beta output
  setwd(output_folders[grep(specs[i],output_folders)][1])
  px.rangeL <- predict(xm, predictors.ranges[[i]], filename = paste(specs[i],'px.range.current.logistic',sep = '_'), overwrite =T)
  
  px.rangeR <- predict(xm, predictors.ranges[[i]], args = "outputformat=raw", filename = paste(specs[i],'px.range.current.raw',sep = '_'), overwrite = T)
  
  px.flL <- predict(xm, predictors.fl, filename = paste(specs[i],'px.FL.current.logistic',sep = '_'), overwrite = TRUE, progress = 'text')
  
  px.flR <- predict(xm, predictors.fl, args = "outputformat=raw", filename = paste(specs[i],'px.FL.current.raw',sep = '_'), overwrite = TRUE)
  
  #save values at presence points to Best layers Best beta
  
  px.flL.values <- getValues(px.flL)
  px.flL.pres <- px.flL.values[unique(cellFromXY(px.flL, occ))]
  write.csv(px.flL.pres, 'px.flL.pres.csv')
  
  px.flR.values <- getValues(px.flR)
  px.flR.pres <- px.flR.values[unique(cellFromXY(px.flR, occ))]
  write.csv(px.flR.pres, 'px.flR.pres.csv')
  
  px.rangeL.values <- getValues(px.rangeL)
  px.rangeL.pres <- px.rangeL.values[unique(cellFromXY(px.rangeL, occ))]
  write.csv(px.rangeL.pres, 'px.rangeL.pres.csv')
  
  px.rangeR.values <- getValues(px.rangeR)
  px.rangeR.pres <- px.rangeR.values[unique(cellFromXY(px.rangeR, occ))]
  write.csv(px.rangeR.pres, 'px.rangeR.pres.csv')
  
  print(i)
  
  
}


###Change GRI files into .asc to load into GIS
for( i in 1: length(specs)){
 setwd(output_folders[grep(specs[i],output_folders)][1])
 predict<-raster(paste(specs[i], 'px.FL.current.logistic.gri', sep = '_'))


writeRaster(predict, filename= specs[i],'Suitability_Map.asc',sep='_', overwrite=TRUE)
print(i)
} 



