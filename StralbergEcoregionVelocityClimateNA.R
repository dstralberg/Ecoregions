library(raster)
library(dismo) #for mess function
library(yaImpute) # for ann function

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
setwd(eco <- "D:/NorthAmerica/Ecoregions/")
ecocurr <- raster("currentlcc.tif")
datlcc <- read.csv("CECEcoregionSampleLCC.csv")
#l3eco <- shapefile("NA_CEC_Eco_Level3_LCC.shp")
ecolu <- read.csv("ecoregion_lookup.csv")
ecolev3 <- raster("F:/GIS/ecoregions/CEC/ceclev3idlcc.asc")

future <- "E:/CMIP5/"
cur <- "E:/CMIP5/baseline19812010/"
setwd(cur)
clim <- list.files(cur, pattern =".asc$")
curclim<-stack(clim)

sampleclim<-cbind(datlcc,extract(curclim,as.matrix(cbind(datlcc[,3],datlcc[,4]))))
sc <- na.omit(sampleclim)
sc$NA_L3CODE <- as.factor(as.character(sc$NA_L3CODE))

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

#Separate out individual ecoregion predictions and calculate baseline novel climates with MESS analysis
pres <- ecocurr
pres[pres!=1] <- NA
pres <- raster::trim(pres)
writeRaster(pres,filename=paste(eco,1,"pred_curr",sep=""), format="GTiff", overwrite=TRUE)
ptc <- crop(curclim,pres)
ptc <-mask(ptc,pres)
p <- cbind(sc,extract(pres,as.matrix(cbind(sc[,3],sc[,4]))))
names(p)[ncol(p)] <- "pres"
p1 <- p[p$pres==1,]
if (nrow(p1)>0) {
  x <- try(m1 <- mess(ptc, p1[,5:30], full=FALSE))
  writeRaster(m1,filename=paste(eco,1,"novel_curr.tif",sep=""),overwrite=TRUE)
}
for (i in 2:181){
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- NA
  pres <- raster::trim(pres)
  writeRaster(pres,filename=paste(eco,i,"pred_curr",sep=""), format="GTiff", overwrite=TRUE)	
  
  ptc <- crop(curclim,pres)
  ptc <-mask(ptc,pres)
  p <- cbind(sc,extract(pres,as.matrix(cbind(sc[,3],sc[,4]))))
  names(p)[ncol(p)] <- "pres"
  p1 <- p[p$pres==1,]
  if (nrow(p1)>0) {
    x <- try(m1 <- mess(ptc, p1[,5:30], full=FALSE))
    writeRaster(m1,filename=paste(eco,i,"novel_curr.tif",sep=""),overwrite=TRUE)
  }
}


#Calculate refugia index 
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

#Calculate climate dissimilarity for each future ecoregion projection
w  <- paste(future,"NA_ENSEMBLE_rcp85_2080s_Bioclim_ASCII/",sep="")
setwd(w)
futclim <- list.files(w,pattern=".asc$")
s <-stack(futclim)
m <- c(-1000, 0.999, 0,  0.999, 1.001, 1, 1.001, 1000,0)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
pres <- ecocurr
pres[pres!=1] <- 0	
p <- cbind(sc,extract(pres,as.matrix(cbind(sc[,3],sc[,4]))))
names(p)[ncol(p)] <- "pres"
p1 <- p[p$pres==1,]
if (nrow(p1)>0) {
  fut <- raster(paste(eco,"_pred_rcp85_2080s.tif",sep=""))
  fut[fut!=1] <- NA	
  x<- try(cellStats(fut,min))
  if (x!=Inf) {
    ft <- raster::trim(fut)
    f <- crop(s,ft)
    f1<-mask(f,ft)
    x <- try(m1 <- mess(f1, p1[,5:30], full=FALSE))
    writeRaster(m1,filename=paste(eco,1,"novel_rcp85_2080s.tif",sep=""),overwrite=TRUE)
  }
}
for (i in 2:181) {
  pres <- ecocurr
  pres[pres==1] <- 0
  pres[pres==i] <- 1
  pres[pres!=1] <- 0	
  p <- cbind(sc,extract(pres,as.matrix(cbind(sc[,3],sc[,4]))))
  names(p)[ncol(p)] <- "pres"
  p1 <- p[p$pres==1,]
  if (nrow(p1)>0) {
    fut <- raster(paste(eco,"_pred_rcp85_2080s.tif",sep=""))
    fut[fut==1] <- 0
    fut[fut==i] <- 1
    fut[fut!=1] <- NA	
    x<- try(cellStats(fut,min))
    if (x!=Inf) {
      ft <- raster::trim(fut)
      f <- crop(s,ft)
      f1<-mask(f,ft)
      x <- try(m1 <- mess(f1, p1[,5:30], full=FALSE))
      writeRaster(m1,filename=paste(eco,i,"novel_rcp85_2080s.tif",sep=""),overwrite=TRUE)
    }
  }
}
