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
library("sp")
library("gstat")
library("oce")

# { Read data and calculate a distance matrix -----------------------------

# read data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")


# longitude should be negative, we are in the Western hemisphere
d$long <- -d$long

# compute a distance matrix using geodDist
# TODO: this data frame is transect 2 only. put this in ddply and let it compute per transect
# geodDist computes in km, so need to multiply by 1000 to get distances in meters
d$dist <- geodDist(d$lat, d$long, d$lat[1], d$long[1]) * 1000

# convert data to a spatial object
dSp <- d
coordinates (dSp) <- ~dist + depth

# computing the directional variogram
V <- variogram(temp~1, data=dSp, alpha=c(0,90), cutoff=50)
plot(V, main="temp, directional, distance unchanged")
# as expected, big difference between the plots

# compute a new coordinate system
d$distnew <- d$dist / 1000
dSp <- d
coordinates (dSp) <- ~distnew + depth

# computing the directional variogram for temperature
V <- variogram(temp~1, data=dSp, alpha=c(0,90), cutoff=50)
plot(V, main="temp, directional, distance / 1000")

# for salinity
V <- variogram(salinity~1, data=dSp, alpha=c(0,90), cutoff=50)
plot(V, main="salinity, directional, distance / 1000")

# for fluorometry
V <- variogram(fluoro~1, data=dSp, alpha=c(0,90), cutoff=50)
plot(V, main="fluorometry, directional, distance / 1000")

# for oxygen
V <- variogram(oxygen~1, data=dSp, alpha=c(0,90), cutoff=50)
plot(V, main="oxygen, directional, distance / 1000")

# tried different divisors for distance. overall, 1000 is ok... but there is a significant trend in the horizontal in all variables and I can't seem to find an inflection point.
# why does it take so long to calculate and plot? 

# TODO: interpolate with interp() using this anisotrophy ratio
# need to first cast into wide format
