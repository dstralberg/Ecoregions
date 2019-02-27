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

pdf(file=paste(eco,"ecoregionrefugia_rcp85_2080s.pdf",sep=""), width=8, height=10.5)
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
map1a <- raster::trim(raster(paste(eco,i,"refmean2080_rcp85.tif",sep="")))
map1a[map1a == 0] <- NA
try(map1a <- trim(map1a))
plot(map1,legend=FALSE,axes=FALSE,col="gray",main=paste(i, ". ", ecolu[i,3], sep=""),legend.mar=0, cex.main=1.0)
plot(poly,border="light grey",col="transparent",add=TRUE)
plot(map1a, add=TRUE, col=bluegreen.colors(10), alpha=0.8, zlim=c(0,1), axes=FALSE, legend=FALSE, main=paste("Ecoregion ",i, sep=""),cex.main=1.0, min(1000000,ncell(map1a)))
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