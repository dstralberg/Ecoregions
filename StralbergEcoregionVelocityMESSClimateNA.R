library(raster)
#library(SDMTools)
library(yaImpute) # for ann function
#library(gdata) 
#library(rasterVis)
#library(dplyr) 
#library(stringr) #for str_sub function


LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
setwd(eco <- "D:/NorthAmerica/Ecoregions/")
ecocurr <- raster("currentlcc.tif")
# datlcc <- read.csv("CECEcoregionSampleLCC.csv")
#l3eco <- shapefile("NA_CEC_Eco_Level3_LCC.shp")
ecolu <- read.csv("ecoregion_lookup.csv")
ecolev3 <- raster("F:/GIS/ecoregions/CEC/ceclev3idlcc.asc")

future <- "E:/CMIP5/"
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

# sampleclim<-cbind(datlcc,extract(curclim,as.matrix(cbind(datlcc[,3],datlcc[,4]))))
# sc <- na.omit(sampleclim)
# sc$NA_L3CODE <- as.factor(as.character(sc$NA_L3CODE))

fattail <- function(x, alpha, c) {
	left <- (c/(2*alpha*gamma(1/c)))
	right <- exp(-1*((abs(x/alpha))^c))
	result <- left*right
	return(right)
	}

ftmean <- function(alpha, c) {
		result <-(alpha*gamma(2/c))/gamma(1/c)
		return(result)
		}

#Separate out individual ecoregion predictions
pres <- ecocurr
pres[pres!=1] <- NA
pres <- raster::trim(pres)
writeRaster(pres,filename=paste(eco,1,"pred_curr",sep=""), format="GTiff", overwrite=TRUE)
for (i in 2:181){
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- NA
  pres <- raster::trim(pres)
  writeRaster(pres,filename=paste(eco,i,"pred_curr",sep=""), format="GTiff", overwrite=TRUE)	
}

#Separate out individual mapped ecoregions
curr <- ecolev3
curr[curr!=1] <- NA
curr <- raster::trim(curr)
writeRaster(curr,filename=paste(eco,1,"eco_curr",sep=""), format="GTiff", overwrite=TRUE)
for (i in 2:181){
  curr <- ecolev3
  curr[curr==1] <- 0
  curr[curr==i] <- 1
  curr[curr!=1] <- NA
  curr <- raster::trim(curr)
  writeRaster(curr,filename=paste(eco,i,"eco_curr",sep=""), format="GTiff", overwrite=TRUE)	
}

#2080s RCP 8.5
setwd(eco)
eco2080 <- list.files(eco, pattern="2085.tif$")
eco2080 <- grep(pattern="rcp85",eco2080,value=TRUE) 
emptystack<-stack(raster(eco2080[1]))
emptystack<-dropLayer(emptystack,1)
for (i in 1:181){
	velstack <-emptystack
	pres <- ecocurr
	pres[pres==1] <- 0
	pres[pres==i] <- 1
	pres[pres!=1] <- 0	
	present<-as.data.frame(rasterToPoints(pres))
	names(present)[3] <- "current"
	for (j in 1:length(eco2080)) {
		fut <- raster(eco2080[j])
		fut[fut==1] <- 0
		fut[fut==i] <- 1
		fut[fut!=1] <- 0
		future<-as.data.frame(rasterToPoints(fut))
		names(future)[3] <- "EcoLev3_2080"
		p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$current)
		f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$EcoLev3_2080)
		p.xy2<-p.xy[p.xy[,4]==1,1:3,drop=FALSE]
		f.xy2<-f.xy[f.xy[,4]==1,1:3,drop=FALSE]
		if(nrow(f.xy2)>0){
			d.ann <- as.data.frame(ann(as.matrix(p.xy2[,-1,drop=FALSE]),as.matrix(f.xy2[,-1,drop=FALSE]),k=1, verbose=F)$knnIndexDist)
			d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
			names(d1b) <- c("ID","X","Y","bvel")
			f.xy <- as.data.frame(f.xy)
			names(f.xy) <- c("ID","X","Y","Pres")
			d1b<-left_join(f.xy,d1b,by=c("ID","X","Y"))
			d1b$fat <- fattail(d1b$bvel, 8333.3335, 0.5) # alpha value that results in mean of 50 km / century or 500 m/year
			sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
			sppref <- extend(sppref,ecocurr)
			sppref[is.na(sppref[])] <- 0
			} else {
			print(i)
			sppref <- ecocurr*0}
		velstack <- addLayer(velstack,sppref)
	}
	velmean <- calc(velstack,fun=mean,na.rm=TRUE)
	projection(velmean) <- LCC
	writeRaster(velmean,filename=paste(i,"refmean2080_rcp85",sep=""), format="GTiff", overwrite=TRUE)
}

#Combined ecoregion refugia
ref <- list.files(eco, pattern="2080_rcp85.tif$")	
refs <- stack(ref)
rmax <- refs[[1]]
for (i in 2:length(ref)) {
  s <- stack(rmax,refs[[i]])
  rmax <- calc(s,max)
}
rmax <- mask(rmax,ecocurr)
writeRaster(rmax,filename="_rcp85_2080_maxref1", format="GTiff", overwrite=TRUE)

#2050s RCP 8.5
eco2050 <- list.files(eco, pattern="2055.tif$")
eco2050 <- grep(pattern="rcp85",eco2050,value=TRUE) 
for (i in 1:181){
  velstack <-emptystack
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- 0	
  present<-as.data.frame(rasterToPoints(pres))
  names(present)[3] <- "current"
  for (j in 1:length(eco2050)) {
    fut <- raster(eco2050[j])
    fut[fut==1] <- 0
    fut[fut==i] <- 1
    fut[fut!=1] <- 0	
    future<-as.data.frame(rasterToPoints(fut))	
    names(future)[3] <- "EcoLev3_2050"
    p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$current)
    f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$EcoLev3_2050)
    p.xy2<-p.xy[p.xy[,4]==1,1:3,drop=FALSE]
    f.xy2<-f.xy[f.xy[,4]==1,1:3,drop=FALSE]
    if(nrow(f.xy2)>0){
      d.ann <- as.data.frame(ann(as.matrix(p.xy2[,-1,drop=FALSE]),as.matrix(f.xy2[,-1,drop=FALSE]),k=1, verbose=F)$knnIndexDist)
      d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
      names(d1b) <- c("ID","X","Y","bvel")
      f.xy <- as.data.frame(f.xy)
      names(f.xy) <- c("ID","X","Y","Pres")
      d1b<-left_join(f.xy,d1b,by=c("ID","X","Y"))
      d1b$fat <- fattail(d1b$bvel, 4166.667, 0.5) # alpha value that results in mean of 25 km / half-century or 500 m/year
      sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
      sppref <- extend(sppref,ecocurr)
      sppref[is.na(sppref[])] <- 0
    } else {
      print(i)
      sppref <- ecocurr*0}
    velstack <- addLayer(velstack,sppref)
  }
  velmean <- calc(velstack,fun=mean,na.rm=TRUE)
  projection(velmean) <- LCC
  writeRaster(velmean,filename=paste(i,"refmean2050_rcp85",sep=""), format="GTiff", overwrite=TRUE)
}

#Combined ecoregion refugia
ref <- list.files(eco, pattern="2050_rcp85.tif$")	
refs <- stack(ref)
rmax <- refs[[1]]
for (i in 2:length(ref)) {
  s <- stack(rmax,refs[[i]])
  rmax <- calc(s,max)
}
writeRaster(rmax,filename="_rcp85_2050_maxref1", format="GTiff", overwrite=TRUE)

#2080s RCP 4.5
eco2080 <- list.files(eco, pattern="2085.tif$")
eco2080 <- grep(pattern="rcp45",eco2080,value=TRUE) 
for (i in 1:181){
  velstack <-emptystack
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- 0	
  present<-as.data.frame(rasterToPoints(pres))
  names(present)[3] <- "current"
  for (j in 1:length(eco2080)) {
    fut <- raster(eco2080[j])
    fut[fut==1] <- 0
    fut[fut==i] <- 1
    fut[fut!=1] <- 0	
    future<-as.data.frame(rasterToPoints(fut))	
    names(future)[3] <- "EcoLev3_2080"
    p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$current)
    f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$EcoLev3_2080)
    p.xy2<-p.xy[p.xy[,4]==1,1:3,drop=FALSE]
    f.xy2<-f.xy[f.xy[,4]==1,1:3,drop=FALSE]
    if(nrow(f.xy2)>0){
      d.ann <- as.data.frame(ann(as.matrix(p.xy2[,-1,drop=FALSE]),as.matrix(f.xy2[,-1,drop=FALSE]),k=1, verbose=F)$knnIndexDist)
      d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
      names(d1b) <- c("ID","X","Y","bvel")
      f.xy <- as.data.frame(f.xy)
      names(f.xy) <- c("ID","X","Y","Pres")
      d1b<-left_join(f.xy,d1b,by=c("ID","X","Y"))
      d1b$fat <- fattail(d1b$bvel, 8333.3335, 0.5) # alpha value that results in mean of 50 km / century or 500 m/year
      sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
      sppref <- extend(sppref,ecocurr)
      sppref[is.na(sppref[])] <- 0
    } else {
      print(i)
      sppref <- ecocurr*0}
    velstack <- addLayer(velstack,sppref)
  }
  velmean <- calc(velstack,fun=mean,na.rm=TRUE)
  projection(velmean) <- LCC
  writeRaster(velmean,filename=paste(i,"refmean2080_rcp45",sep=""), format="GTiff", overwrite=TRUE)
}

#Combined ecoregion refugia
ref <- list.files(e, pattern="2080_rcp45.tif$")	
refs <- stack(ref)
rmax <- refs[[1]]
for (i in 2:length(ref)) {
  s <- stack(rmax,refs[[i]])
  rmax <- calc(s,max)
}
writeRaster(rmax,filename="_rcp45_2080_maxref1", format="GTiff", overwrite=TRUE)

#2050s RCP 4.5
eco2050 <- list.files(eco, pattern="2055.tif$")
eco2050 <- grep(pattern="rcp45",eco2050,value=TRUE)
for (i in 1:181){
  velstack <-emptystack
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- 0	
  present<-as.data.frame(rasterToPoints(pres))
  names(present)[3] <- "current"
  for (j in 1:length(eco2050)) {
    fut <- raster(eco2050[j])
    fut[fut==1] <- 0
    fut[fut==i] <- 1
    fut[fut!=1] <- 0	
    future<-as.data.frame(rasterToPoints(fut))	
    names(future)[3] <- "EcoLev3_2050"
    p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$current)
    f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$EcoLev3_2050)
    p.xy2<-p.xy[p.xy[,4]==1,1:3,drop=FALSE]
    f.xy2<-f.xy[f.xy[,4]==1,1:3,drop=FALSE]
    if(nrow(f.xy2)>0){
      d.ann <- as.data.frame(ann(as.matrix(p.xy2[,-1,drop=FALSE]),as.matrix(f.xy2[,-1,drop=FALSE]),k=1, verbose=F)$knnIndexDist)
      d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
      names(d1b) <- c("ID","X","Y","bvel")
      f.xy <- as.data.frame(f.xy)
      names(f.xy) <- c("ID","X","Y","Pres")
      d1b<-left_join(f.xy,d1b,by=c("ID","X","Y"))
      d1b$fat <- fattail(d1b$bvel, 4166.667, 0.5) # alpha value that results in mean of 25 km / half-century or 500 m/year
      sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
      sppref <- extend(sppref,ecocurr)
      sppref[is.na(sppref[])] <- 0
    } else {
      print(i)
      sppref <- ecocurr*0}
    velstack <- addLayer(velstack,sppref)
  }
  velmean <- calc(velstack,fun=mean,na.rm=TRUE)
  projection(velmean) <- LCC
  writeRaster(velmean,filename=paste(i,"refmean2050_rcp45",sep=""), format="GTiff", overwrite=TRUE)
}

#Combined ecoregion refugia
ref <- list.files(e, pattern="2050_rcp45.tif$")	
refs <- stack(ref)
rmax <- refs[[1]]
for (i in 2:length(ref)) {
  s <- stack(rmax,refs[[i]])
  rmax <- calc(s,max)
}
writeRaster(rmax,filename="_rcp45_2050_maxref1", format="GTiff", overwrite=TRUE)


##Mapping
library(raster)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
setwd(eco <- "D:/NorthAmerica/Ecoregions/")
ecolu <- read.csv("ecoregion_lookup.csv")
ecolev3 <- raster("F:/GIS/ecoregions/CEC/ceclev3idlcc.asc")

buff <- function(e,km) {
  e[1] <- e[1] - 1000*km
  e[2] <- e[2] + 1000*km
  e[3] <- e[3] - 1000*km
  e[4] <- e[4] + 1000*km
  return(e)
}

bluegreen.colors <- colorRampPalette(c("#FFFACD", "lemonchiffon","#FFF68F", "khaki1","#ADFF2F", "greenyellow", "#00CD00", "green3", "#48D1CC", "mediumturquoise", "#007FFF", "blue"), space="Lab", bias=2)
provstate <- shapefile("F:/GIS/basemaps/province_state_lcc.shp")
mexstate <- shapefile("F:/GIS/basemaps/mexstatesLCC.shp")
poly <- bind(provstate,mexstate)
reds <- colorRampPalette(c("red"))

pdf(file=paste(eco,"ecoregionrefugia_rcp85_2080s_MESS.pdf",sep=""), width=8, height=10.5)
i<-1
while (i<181) {
  par(mfcol=c(3,2), oma=c(0,0,0,0))
  par(mar=c(2,2,2,2), bty="n")
  par(xpd=TRUE)	
  map1 <- raster(paste(eco,i,"eco_curr.tif",sep=""))
  map1 <- extend(map1,buff(extent(map1),100))
  map2 <- raster(paste(eco,i+1,"eco_curr.tif",sep=""))
  map2 <- extend(map2,buff(extent(map2),100))
  map3 <- raster(paste(eco,i+2,"eco_curr.tif",sep=""))
  map3 <- extend(map3,buff(extent(map3),100))
  map4 <- raster(paste(eco,i+3,"eco_curr.tif",sep=""))
  map4 <- extend(map4,buff(extent(map4),100))
  map5 <- raster(paste(eco,i+4,"eco_curr.tif",sep=""))
  map5 <- extend(map5,buff(extent(map5),100))
  map6 <- raster(paste(eco,i+5,"eco_curr.tif",sep=""))
  
  map1a <- raster(paste(eco,i,"refmean2080_rcp85.tif",sep=""))
  map1a[map1a == 0] <- NA
  try(map1a <- trim(map1a))
  map2a <- raster(paste(eco,i+1,"refmean2080_rcp85.tif",sep=""))	
  map2a[map2a == 0] <- NA
  try(map2a <- trim(map2a))
  map3a <- raster(paste(eco,i+2,"refmean2080_rcp85.tif",sep=""))
  map3a[map3a == 0] <- NA
  try(map3a <- trim(map3a))
  map4a <- raster(paste(eco,i+3,"refmean2080_rcp85.tif",sep=""))
  map4a[map4a == 0] <- NA
  try(map4a <- trim(map4a))
  map5a <- raster(paste(eco,i+4,"refmean2080_rcp85.tif",sep=""))
  map5a[map5a == 0] <- NA
  try(map5a <- trim(map5a))
  map6a <- raster(paste(eco,i+5,"refmean2080_rcp85.tif",sep=""))
  map6a[map6a == 0] <- NA
  try(map6a <- trim(map6a))
  
  try(map1b <- raster(paste(eco,i,"novel_rcp85_2080s.tif",sep="")))
  try(map2b <- raster(paste(eco,i+1,"novel_rcp85_2080s.tif",sep="")))	
  try(map3b <- raster(paste(eco,i+2,"novel_rcp85_2080s.tif",sep="")))
  try(map4b <- raster(paste(eco,i+3,"novel_rcp85_2080s.tif",sep="")))
  try(map5b <- raster(paste(eco,i+4,"novel_rcp85_2080s.tif",sep="")))
  try(map6b <- raster(paste(eco,i+5,"novel_rcp85_2080s.tif",sep="")))
  
  map1c <- raster(paste(eco,i,"novel_curr.tif",sep=""))
  map2c <- raster(paste(eco,i+1,"novel_curr.tif",sep=""))	
  map3c <- raster(paste(eco,i+2,"novel_curr.tif",sep=""))
  map4c <- raster(paste(eco,i+3,"novel_curr.tif",sep=""))
  map5c <- raster(paste(eco,i+4,"novel_curr.tif",sep=""))
  map6c <- raster(paste(eco,i+5,"novel_curr.tif",sep=""))
  
  try(map1d <- map1b<quantile(map1c,0.01))
  try(map1d[map1d == 0] <- NA)
  try(map2d <- map2b<quantile(map2c,0.01))
  try(map2d[map2d == 0] <- NA)
  try(map3d <- map3b<quantile(map3c,0.01))
  try(map3d[map3d == 0] <- NA)
  try(map4d <- map4b<quantile(map4c,0.01))
  try(map4d[map4d == 0] <- NA)
  try(map5d <- map3b<quantile(map5c,0.01))
  try(map5d[map5d == 0] <- NA)
  try(map6d <- map4b<quantile(map6c,0.01))
  try(map6d[map6d == 0] <- NA)
  
  plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map1a)))
  plot(map1d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map1d)))
  scalebar((xmax(map1)-xmin(map1))/2)
  
  plot(map2,legend=FALSE,axes=FALSE,col="gray",main=paste(i+1, ". ", ecolu[i+1,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map2a*map2d, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  plot(map2d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map2d)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map3d)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10))
  scalebar((xmax(map3)-xmin(map3))/2)
  
  plot(map4,legend=FALSE,axes=FALSE,col="gray",main=paste(i+3, ". ", ecolu[i+3,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map4a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map4a)))
  plot(map4d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map4d)))
  scalebar((xmax(map4)-xmin(map4))/2)
  
  plot(map5,legend=FALSE,axes=FALSE,col="gray",main=paste(i+4, ". ", ecolu[i+4,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map5a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE,main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map5a)))
  plot(map5d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map5d)))
  scalebar((xmax(map5)-xmin(map5))/2)
  
  plot(map6,legend=FALSE,axes=FALSE,col="gray",main=paste(i+5, ". ", ecolu[i+5,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map6a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map6a)))
  plot(map6d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map6d)))
  scalebar((xmax(map6)-xmin(map6))/2)
  
  i<-i+6
}
map1 <- raster::trim(raster(paste(eco,i,"eco_curr.tif",sep="")))
map1a <- raster::trim(raster(paste(eco,i,"refmean2080_rcp85.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
plot(map1d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map1d)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()

pdf(file=paste(eco,"ecoregionrefugia_rcp85_2050s.pdf",sep=""), width=8, height=10.5)
i<-1
while (i<181) {
  par(mfcol=c(3,2), oma=c(0,0,0,0))
  par(mar=c(2,2,2,2), bty="n")
  par(xpd=TRUE)	
  map1 <- raster(paste(eco,i,"eco_curr.tif",sep=""))
  map1 <- extend(map1,buff(extent(map1),100))
  map2 <- raster(paste(eco,i+1,"eco_curr.tif",sep=""))
  map2 <- extend(map2,buff(extent(map2),100))
  map3 <- raster(paste(eco,i+2,"eco_curr.tif",sep=""))
  map3 <- extend(map3,buff(extent(map3),100))
  map4 <- raster(paste(eco,i+3,"eco_curr.tif",sep=""))
  map4 <- extend(map4,buff(extent(map4),100))
  map5 <- raster(paste(eco,i+4,"eco_curr.tif",sep=""))
  map5 <- extend(map5,buff(extent(map5),100))
  map6 <- raster(paste(eco,i+5,"eco_curr.tif",sep=""))
  map1a <- raster(paste(eco,i,"refmean2050_rcp85.tif",sep=""))
  map1a[map1a == 0] <- NA
  try(map1a <- trim(map1a))
  map2a <- raster(paste(eco,i+1,"refmean2050_rcp85.tif",sep=""))	
  map2a[map2a == 0] <- NA
  try(map2a <- trim(map2a))
  map3a <- raster(paste(eco,i+2,"refmean2050_rcp85.tif",sep=""))
  map3a[map3a == 0] <- NA
  try(map3a <- trim(map3a))
  map4a <- raster(paste(eco,i+3,"refmean2050_rcp85.tif",sep=""))
  map4a[map4a == 0] <- NA
  try(map4a <- trim(map4a))
  map5a <- raster(paste(eco,i+4,"refmean2050_rcp85.tif",sep=""))
  map5a[map5a == 0] <- NA
  try(map5a <- trim(map5a))
  map6a <- raster(paste(eco,i+5,"refmean2050_rcp85.tif",sep=""))
  map6a[map6a == 0] <- NA
  try(map6a <- trim(map6a))
  
  plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map1a)))
  scalebar((xmax(map1)-xmin(map1))/2)
  
  plot(map2,legend=FALSE,axes=FALSE,col="gray",main=paste(i+1, ". ", ecolu[i+1,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10))
  scalebar((xmax(map3)-xmin(map3))/2)
  
  plot(map4,legend=FALSE,axes=FALSE,col="gray",main=paste(i+3, ". ", ecolu[i+3,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map4a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map4a)))
  scalebar((xmax(map4)-xmin(map4))/2)
  
  plot(map5,legend=FALSE,axes=FALSE,col="gray",main=paste(i+4, ". ", ecolu[i+4,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map5a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE,main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map5a)))
  scalebar((xmax(map5)-xmin(map5))/2)
  
  plot(map6,legend=FALSE,axes=FALSE,col="gray",main=paste(i+5, ". ", ecolu[i+5,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map6a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map6a)))
  scalebar((xmax(map6)-xmin(map6))/2)
  
  i<-i+6
}
map1 <- raster::trim(raster(paste(eco,i,"eco_curr.tif",sep="")))
map1a <- raster::trim(raster(paste(eco,i,"refmean2050_rcp85.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map1a)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()

pdf(file=paste(eco,"ecoregionrefugia_rcp45_2080s.pdf",sep=""), width=8, height=10.5)
i<-1
while (i<181) {
  par(mfcol=c(3,2), oma=c(0,0,0,0))
  par(mar=c(2,2,2,2), bty="n")
  par(xpd=TRUE)	
  map1 <- raster(paste(eco,i,"eco_curr.tif",sep=""))
  map1 <- extend(map1,buff(extent(map1),100))
  map2 <- raster(paste(eco,i+1,"eco_curr.tif",sep=""))
  map2 <- extend(map2,buff(extent(map2),100))
  map3 <- raster(paste(eco,i+2,"eco_curr.tif",sep=""))
  map3 <- extend(map3,buff(extent(map3),100))
  map4 <- raster(paste(eco,i+3,"eco_curr.tif",sep=""))
  map4 <- extend(map4,buff(extent(map4),100))
  map5 <- raster(paste(eco,i+4,"eco_curr.tif",sep=""))
  map5 <- extend(map5,buff(extent(map5),100))
  map6 <- raster(paste(eco,i+5,"eco_curr.tif",sep=""))
  map1a <- raster(paste(eco,i,"refmean2080_rcp45.tif",sep=""))
  map1a[map1a == 0] <- NA
  try(map1a <- trim(map1a))
  map2a <- raster(paste(eco,i+1,"refmean2080_rcp45.tif",sep=""))	
  map2a[map2a == 0] <- NA
  try(map2a <- trim(map2a))
  map3a <- raster(paste(eco,i+2,"refmean2080_rcp45.tif",sep=""))
  map3a[map3a == 0] <- NA
  try(map3a <- trim(map3a))
  map4a <- raster(paste(eco,i+3,"refmean2080_rcp45.tif",sep=""))
  map4a[map4a == 0] <- NA
  try(map4a <- trim(map4a))
  map5a <- raster(paste(eco,i+4,"refmean2080_rcp45.tif",sep=""))
  map5a[map5a == 0] <- NA
  try(map5a <- trim(map5a))
  map6a <- raster(paste(eco,i+5,"refmean2080_rcp45.tif",sep=""))
  map6a[map6a == 0] <- NA
  try(map6a <- trim(map6a))
  
  plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map1a)))
  scalebar((xmax(map1)-xmin(map1))/2)
  
  plot(map2,legend=FALSE,axes=FALSE,col="gray",main=paste(i+1, ". ", ecolu[i+1,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10))
  scalebar((xmax(map3)-xmin(map3))/2)
  
  plot(map4,legend=FALSE,axes=FALSE,col="gray",main=paste(i+3, ". ", ecolu[i+3,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map4a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map4a)))
  scalebar((xmax(map4)-xmin(map4))/2)
  
  plot(map5,legend=FALSE,axes=FALSE,col="gray",main=paste(i+4, ". ", ecolu[i+4,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map5a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE,main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map5a)))
  scalebar((xmax(map5)-xmin(map5))/2)
  
  plot(map6,legend=FALSE,axes=FALSE,col="gray",main=paste(i+5, ". ", ecolu[i+5,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map6a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map6a)))
  scalebar((xmax(map6)-xmin(map6))/2)
  
  i<-i+6
}
map1 <- raster::trim(raster(paste(eco,i,"eco_curr.tif",sep="")))
map1a <- raster::trim(raster(paste(eco,i,"refmean2080_rcp45.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()

pdf(file=paste(eco,"ecoregionrefugia_rcp45_2050s.pdf",sep=""), width=8, height=10.5)
i<-1
while (i<181) {
  par(mfcol=c(3,2), oma=c(0,0,0,0))
  par(mar=c(2,2,2,2), bty="n")
  par(xpd=TRUE)	
  map1 <- raster(paste(eco,i,"eco_curr.tif",sep=""))
  map1 <- extend(map1,buff(extent(map1),100))
  map2 <- raster(paste(eco,i+1,"eco_curr.tif",sep=""))
  map2 <- extend(map2,buff(extent(map2),100))
  map3 <- raster(paste(eco,i+2,"eco_curr.tif",sep=""))
  map3 <- extend(map3,buff(extent(map3),100))
  map4 <- raster(paste(eco,i+3,"eco_curr.tif",sep=""))
  map4 <- extend(map4,buff(extent(map4),100))
  map5 <- raster(paste(eco,i+4,"eco_curr.tif",sep=""))
  map5 <- extend(map5,buff(extent(map5),100))
  map6 <- raster(paste(eco,i+5,"eco_curr.tif",sep=""))
  map1a <- raster(paste(eco,i,"refmean2050_rcp45.tif",sep=""))
  map1a[map1a == 0] <- NA
  try(map1a <- trim(map1a))
  map2a <- raster(paste(eco,i+1,"refmean2050_rcp45.tif",sep=""))	
  map2a[map2a == 0] <- NA
  try(map2a <- trim(map2a))
  map3a <- raster(paste(eco,i+2,"refmean2050_rcp45.tif",sep=""))
  map3a[map3a == 0] <- NA
  try(map3a <- trim(map3a))
  map4a <- raster(paste(eco,i+3,"refmean2050_rcp45.tif",sep=""))
  map4a[map4a == 0] <- NA
  try(map4a <- trim(map4a))
  map5a <- raster(paste(eco,i+4,"refmean2050_rcp45.tif",sep=""))
  map5a[map5a == 0] <- NA
  try(map5a <- trim(map5a))
  map6a <- raster(paste(eco,i+5,"refmean2050_rcp45.tif",sep=""))
  map6a[map6a == 0] <- NA
  try(map6a <- trim(map6a))
  
  plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map1a)))
  scalebar((xmax(map1)-xmin(map1))/2)
  
  plot(map2,legend=FALSE,axes=FALSE,col="gray",main=paste(i+1, ". ", ecolu[i+1,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10))
  scalebar((xmax(map3)-xmin(map3))/2)
  
  plot(map4,legend=FALSE,axes=FALSE,col="gray",main=paste(i+3, ". ", ecolu[i+3,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map4a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map4a)))
  scalebar((xmax(map4)-xmin(map4))/2)
  
  plot(map5,legend=FALSE,axes=FALSE,col="gray",main=paste(i+4, ". ", ecolu[i+4,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map5a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE,main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map5a)))
  scalebar((xmax(map5)-xmin(map5))/2)
  
  plot(map6,legend=FALSE,axes=FALSE,col="gray",main=paste(i+5, ". ", ecolu[i+5,3], sep=""),cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)  
  plot(map6a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map6a)))
  scalebar((xmax(map6)-xmin(map6))/2)
  
  i<-i+6
}
map1 <- raster::trim(raster(paste(eco,i,"eco_curr.tif",sep="")))
map1a <- raster::trim(raster(paste(eco,i,"refmean2050_rcp45.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()