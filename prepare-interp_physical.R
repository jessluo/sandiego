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
# phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")
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

# transform the distance to account for the anisotropy
f <- 1800
phy$distTr <- phy$dist / f

# new coordinate
n <- 50
xo <- seq(min(phy$distTr, na.rm=T), max(phy$distTr, na.rm=T), length.out=n*4)
yo <- seq(0, max(phy$depth, na.rm=T), length.out=n)

ti <- ddply(phy, ~transect, function(x) {
  # have to subet otherwise interp fails
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

ggplot(ti) + geom_tile(aes(x=dist, y=-depth, fill=temp)) + geom_contour(aes(x=dist, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.)

# tried to interpolate salinity but it is returning NAs


ti <- ddply(phy, ~transect, function(x) {
  # have to subet otherwise interp fails
  # no idea why...
  x <- x[seq(1, nrow(x), 2),]
  
  # interpolate
  ti <- interp(x$distTr, x$depth, x$swRho, xo=xo, yo=yo, duplicate="mean")
  
  # convert the result to data.frame form
  ti <- list2frame(ti)
  return(ti)
}, .progress="text")
names(ti) <- c("transect", "dist", "depth", "swRho")
ti$dist <- ti$dist * f

ggplot(ti) + geom_tile(aes(x=dist, y=-depth, fill=temp)) + geom_contour(aes(x=dist, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5) + facet_grid(transect~.)
# }

