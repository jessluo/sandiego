#
#      Interpolation of Environmental Data
#
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards


##{ Compute a distance coordinate -----------------------------------------

library("oce")

# ggplot(phy) + geom_point(aes(x=long, y=-depth, colour=temp), size=1) + facet_grid(transect~.)

# compute a distance from the Western-most point using geodDist
# for each transect, we project data on a plane defined by the mean latitude (that is not a problem here because the latitude varies so little)
lonRef <- min(phy$long)
phy <- ddply(phy, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
})

# ggplot(phy) + geom_point(aes(x=dist, y=-depth, colour=temp), size=1) + facet_grid(transect~.)
# -> OK

# }


##{ Evaluate anisotropy ---------------------------------------------------

library("sp")
library("gstat")

# subsample data because the computation of the variogram is too long otherwise
sub <- seq(1, nrow(phy), 15)
ggplot(phy[sub,]) + geom_point(aes(x=dist, y=-depth, colour=temp), size=1) + facet_grid(transect~.)


# regular variogram
points <- seq(0, 50, by=2)

# convert data to a spatial object
phySp <- phy[sub,]
coordinates (phySp) <- ~dist + depth

# compute the directional variogram
v <- variogram(temp~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(v, main="temp, directional, dist unchanged")
# -> as expected, big difference between the plots
#    0 is vertical, 90 is horizontal

# transform horizontal distance to test whether that captures the horizontal-vertical anisotropy
phy$distTr <- phy$dist / 1000

# redo the variogram
# convert the data into a spatial object
phySp <- phy[sub,]
coordinates(phySp) <- ~distTr+depth

v <- variogram(temp~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(v, main="temp, directional, dist x / 1000")
# -> we never reach a sill because there is a global trend in temperature accross the transect. The slopes are different though

v <- variogram(fluoro~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(v, main="fluoro, directional, dist x / 1000")
# -> here again there is a trend along the vertical (direction 0 here), but there is still an inflection point which seems to be around 20 in the vertical and 40 in the horizontal


phy$distTr <- phy$dist / 1800
phySp <- phy[sub,]
coordinates(phySp) <- ~distTr+depth

v <- variogram(temp~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(v, main="temp, directional, dist x / 1800")
# -> we never reach a sill because there is a global trend in temperature accross the transect. The slopes are different though

v <- variogram(fluoro~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(v, main="fluoro, directional, dist x / 1800")

# for salinity
V <- variogram(salinity~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(V, main="salinity, directional, distance / 1800")

# for oxygen
V <- variogram(oxygen~1, data=phySp, alpha=c(0,90), boundaries=points)
plot(V, main="oxygen, directional, distance / 1800")

# -> 1800 is approximately OK even though there is a strong trend in all variables that we cannot get rid of
#    need to experiment some more (detrend the data, test other factors)

# }


##{ Interpolate -----------------------------------------------------------

library("akima")
source("lib_manip.R")
# set plotting colors
jet.rainbow = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.rainbowrev = colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F", "cyan", "#007FFF", "blue", "#00007F"))

# Spectral colour map from ColorBrewer
spectral <- function(n=6) {
  library("RColorBrewer")
  rev(brewer.pal(name="Spectral", n=n))
}

# transform the distance to account for the anisotropy
f <- 1800
phy$distTr <- phy$dist / f

# new coordinate
n <- 50
xo <- seq(min(phy$distTr, na.rm=T), max(phy$distTr, na.rm=T), length.out=n*4)
yo <- seq(0, max(phy$depth, na.rm=T), length.out=n)

# subset phy for only the upcasts for the salinity / seawater density because of some problems in the data
phyup <- phy[phy$down.up=="up",]

# Interpolate temperature
ti <- ddply(phy, ~transect, function(x) {
  # have to subset otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 2),]

  # interpolate
  ti <- interp(x$distTr, x$depth, x$temp, xo=xo, yo=yo, duplicate="mean")

  # convert the result to data.frame form
  ti <- list2frame(ti)
  return(ti)
}, .progress="text")
names(ti) <- c("transect", "dist", "depth", "temp")
ti$dist <- ti$dist * f

ggplot(ti) + geom_tile(aes(x=dist, y=-depth, fill=temp)) + geom_contour(aes(x=dist, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbow(10), guide="colourbar")

ggplot(ti) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + geom_contour(aes(x=dist, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Temp (C)", na.value="grey80", low = "#2d669f", high = "#c8dcef")

write.csv(ti, file="data/interp_temp.csv", row.names=FALSE)


# Interpolate density
swi <- ddply(phyup, ~transect, function(x) {
  # have to subet otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 2),]
  
  # interpolate
  swi <- interp(x$distTr, x$depth, x$swRho, xo=xo, yo=yo, duplicate="mean")
  
  # convert the result to data.frame form
  swi <- list2frame(swi)
  return(swi)
}, .progress="text")
names(swi) <- c("transect", "dist", "depth", "swRho")
swi$dist <- swi$dist * f

ggplot(swi) + geom_tile(aes(x=dist, y=-depth, fill=swRho)) + geom_contour(aes(x=dist, y=-depth, z=swRho), colour="white", swize=0.5, alpha=0.5) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbowrev(10), guide="colourbar")

# greatest difference? 
ggplot(swi) + geom_tile(aes(x=dist, y=-depth, fill=c(NA, diff(swRho)))) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbowrev(10), guide="colourbar")

write.csv(swi, file="data/interp_swRho.csv", row.names=FALSE)


# Interpolate salinity
si <- ddply(phyup, ~transect, function(x) {
  # have to subet otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 2),]
  
  # interpolate
  si <- interp(x$distTr, x$depth, x$salinity, xo=xo, yo=yo, duplicate="mean")
  
  # convert the result to data.frame form
  si <- list2frame(si)
  return(si)
}, .progress="text")
names(si) <- c("transect", "dist", "depth", "salinity")
si$dist <- si$dist * f

ggplot(si) + geom_tile(aes(x=dist, y=-depth, fill=salinity)) + geom_contour(aes(x=dist, y=-depth, z=salinity), colour="white", size=0.5, alpha=0.5) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbow(10), guide="colourbar")

ggplot(si) + geom_tile(aes(x=dist, y=-depth, fill=c(NA, diff(salinity)))) + geom_contour(aes(x=dist, y=-depth, z=salinity), colour="white", size=0.5, alpha=0.5) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbow(10), guide="colourbar")


write.csv(si, file="data/interp_salinity.csv", row.names=FALSE)

# Interpolate fluorometry
fi <- ddply(phy, ~transect, function(x) {
  # have to subet otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 3),]
  
  # interpolate
  fi <- interp(x$distTr, x$depth, x$fluoro, xo=xo, yo=yo, duplicate="mean")
  
  # convert the result to data.frame form
  fi <- list2frame(fi)
  return(fi)
}, .progress="text")
names(fi) <- c("transect", "dist", "depth", "fluoro")
fi$dist <- fi$dist * f

ggplot(fi) + geom_tile(aes(x=dist, y=-depth, fill=fluoro)) + geom_contour(aes(x=dist, y=-depth, z=fluoro), colour="white", size=0.5, alpha=0.5) + facet_grid(transect~.) + scale_fill_gradientn(colours=jet.rainbow(10), guide="colourbar")

write.csv(fi, file="data/interp_fluoro.csv", row.names=FALSE)

# Interpolate oxygen

oi <- ddply(phy, ~transect, function(x) {
  # have to subet otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 2),]
  
  # interpolate
  oi <- interp(x$distTr, x$depth, x$oxygen, xo=xo, yo=yo, duplicate="mean")
  
  # convert the result to data.frame form
  oi <- list2frame(oi)
  return(oi)
}, .progress="text")
names(oi) <- c("transect", "dist", "depth", "oxygen")
oi$dist <- oi$dist * f

ggplot(oi) + geom_tile(aes(x=dist, y=-depth, fill=oxygen)) + facet_grid(transect~.) + scale_fill_gradientn("Oxygen", colours=spectral(), guide="colourbar", na.value=NA) + labs(x="Distance (km)", y="Depth (m)")

write.csv(oi, file="data/interp_oxygen.csv", row.names=FALSE)
# }
