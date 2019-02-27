library(raster)

LCC <- CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")
setwd(eco <- "D:/NorthAmerica/Ecoregions/")
ecolu <- read.csv("ecoregion_lookup.csv")
ecolev3 <- raster("F:/GIS/ecoregions/CEC/ceclev3idlcc.asc")
l3eco <- shapefile("NA_CEC_Eco_Level3_LCC.shp")

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
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  plot(map2d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map2d)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map3d)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10),smallplot=c(0.85,0.87,0.47,0.67),legend.args=list(text='refugia index',side=3,line=0.5,cex=0.8),axis.args=list(cex=0.7))
  plot(map3d, legend.only=TRUE, col=reds(1), alpha=0.3, smallplot=c(0.85,0.87,0.75,0.78), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='novel',side=3,line=0.3,cex=0.8))
  plot(map3, legend.only=TRUE,col="gray",smallplot=c(0.85,0.87,0.37,0.40), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='baseline',side=3,line=0.3,cex=0.8))
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

pdf(file=paste(eco,"ecoregionrefugia_rcp85_2050s_MESS.pdf",sep=""), width=8, height=10.5)
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
  
  try(map1b <- raster(paste(eco,i,"novel_rcp85_2050s.tif",sep="")))
  try(map2b <- raster(paste(eco,i+1,"novel_rcp85_2050s.tif",sep="")))	
  try(map3b <- raster(paste(eco,i+2,"novel_rcp85_2050s.tif",sep="")))
  try(map4b <- raster(paste(eco,i+3,"novel_rcp85_2050s.tif",sep="")))
  try(map5b <- raster(paste(eco,i+4,"novel_rcp85_2050s.tif",sep="")))
  try(map6b <- raster(paste(eco,i+5,"novel_rcp85_2050s.tif",sep="")))
  
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
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  plot(map2d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map2d)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map3d)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10),smallplot=c(0.85,0.87,0.47,0.67),legend.args=list(text='refugia index',side=3,line=0.5,cex=0.8),axis.args=list(cex=0.7))
  plot(map3d, legend.only=TRUE, col=reds(1), alpha=0.3, smallplot=c(0.85,0.87,0.75,0.78), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='novel',side=3,line=0.3,cex=0.8))
  plot(map3, legend.only=TRUE,col="gray",smallplot=c(0.85,0.87,0.37,0.40), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='baseline',side=3,line=0.3,cex=0.8))
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
map1a <- raster::trim(raster(paste(eco,i,"refmean2050_rcp85.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
plot(map1d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map1d)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()

pdf(file=paste(eco,"ecoregionrefugia_rcp45_2080s_MESS.pdf",sep=""), width=8, height=10.5)
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
  
  try(map1b <- raster(paste(eco,i,"novel_rcp45_2080s.tif",sep="")))
  try(map2b <- raster(paste(eco,i+1,"novel_rcp45_2080s.tif",sep="")))	
  try(map3b <- raster(paste(eco,i+2,"novel_rcp45_2080s.tif",sep="")))
  try(map4b <- raster(paste(eco,i+3,"novel_rcp45_2080s.tif",sep="")))
  try(map5b <- raster(paste(eco,i+4,"novel_rcp45_2080s.tif",sep="")))
  try(map6b <- raster(paste(eco,i+5,"novel_rcp45_2080s.tif",sep="")))
  
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
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  plot(map2d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map2d)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map3d)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10),smallplot=c(0.85,0.87,0.47,0.67),legend.args=list(text='refugia index',side=3,line=0.5,cex=0.8),axis.args=list(cex=0.7))
  plot(map3d, legend.only=TRUE, col=reds(1), alpha=0.3, smallplot=c(0.85,0.87,0.75,0.78), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='novel',side=3,line=0.3,cex=0.8))
  plot(map3, legend.only=TRUE,col="gray",smallplot=c(0.85,0.87,0.37,0.40), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='baseline',side=3,line=0.3,cex=0.8))
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
map1a <- raster::trim(raster(paste(eco,i,"refmean2080_rcp45.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
plot(map1d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map1d)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()

pdf(file=paste(eco,"ecoregionrefugia_rcp45_2050s_MESS.pdf",sep=""), width=8, height=10.5)
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
  
  try(map1b <- raster(paste(eco,i,"novel_rcp45_2050s.tif",sep="")))
  try(map2b <- raster(paste(eco,i+1,"novel_rcp45_2050s.tif",sep="")))	
  try(map3b <- raster(paste(eco,i+2,"novel_rcp45_2050s.tif",sep="")))
  try(map4b <- raster(paste(eco,i+3,"novel_rcp45_2050s.tif",sep="")))
  try(map5b <- raster(paste(eco,i+4,"novel_rcp45_2050s.tif",sep="")))
  try(map6b <- raster(paste(eco,i+5,"novel_rcp45_2050s.tif",sep="")))
  
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
  plot(map2a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map2a)))
  plot(map2d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map2d)))
  scalebar((xmax(map2)-xmin(map2))/2)
  
  plot(map3,legend=FALSE,axes=FALSE,col="gray",main=paste(i+2, ". ", ecolu[i+2,3], sep=""),legend.mar=10, cex.main=1.0)
  plot(poly,border="light grey",col="transparent",add=TRUE)
  plot(map3a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, maxpixels=min(1000000,ncell(map3a)))
  plot(map3d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map3d)))
  plot(map3a, legend.only=TRUE, zlim=c(0,1), col=bluegreen.colors(10),smallplot=c(0.85,0.87,0.47,0.67),legend.args=list(text='refugia index',side=3,line=0.5,cex=0.8),axis.args=list(cex=0.7))
  plot(map3d, legend.only=TRUE, col=reds(1), alpha=0.3, smallplot=c(0.85,0.87,0.75,0.78), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='novel',side=3,line=0.3,cex=0.8))
  plot(map3, legend.only=TRUE,col="gray",smallplot=c(0.85,0.87,0.37,0.40), axis.args = list(at=FALSE,labels=FALSE), legend.args=list(text='baseline',side=3,line=0.3,cex=0.8))
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
map1a <- raster::trim(raster(paste(eco,i,"refmean2050_rcp45.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
plot(map1d, add=TRUE, col=reds(1), alpha=0.3, zlim=c(1,1), axes=FALSE, legend=FALSE, maxpixels=min(1000000,ncell(map1d)))
scalebar((xmax(map1)-xmin(map1))/2)
dev.off()
