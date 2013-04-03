#
#      Plotting the biological data
#      
#
#--------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("foreach")
library("doParallel")
library("oce")
registerDoParallel(cores=detectCores())
parallel <- TRUE


##{ Read data ------------------------------------------------
# read binned data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")

# read physical data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")

bio <- read.csv("data/bio.csv", stringsAsFactors=FALSE)
bio$dateTime <- as.POSIXct(bio$dateTime, tz="GMT")

ti <- read.csv("data/interp_temp.csv", stringsAsFactors=FALSE)
si <- read.csv("data/interp_salinity.csv", stringsAsFactors=FALSE)
swi <- read.csv("data/interp_swRho.csv", stringsAsFactors=FALSE)
fi <- read.csv("data/interp_fluoro.csv", stringsAsFactors=FALSE)
oi <- read.csv("data/interp_oxygen.csv", stringsAsFactors=FALSE)

# }

##{ Delineate the front -----------------------------------------------------
##### how to define the frontal water mass?
ggplot(data=phy) + geom_point(aes(x=salinity, y=temp, colour=long)) + facet_grid(transect~.) + scale_colour_gradientn(colours=rainbow(10))

# initialize
d$front <- NA

# delineate the frontal region
d[d$transect==1 & d$cast <=11,]$front <- "east"
d[d$transect==1 & d$cast >=12 & d$cast <= 15,]$front <- "front"
d[d$transect==1 & d$cast >=16,]$front <- "west"
d[d$transect==2 & d$cast <=22,]$front <- "east"
d[d$transect==2 & d$cast >=23 & d$cast <= 28,]$front <- "front"
d[d$transect==2 & d$cast >=29,]$front <- "west"
d[d$transect==3 & d$cast <=8,]$front <- "west"
d[d$transect==3 & d$cast >=9 & d$cast <= 13,]$front <- "front"
d[d$transect==3 & d$cast >=14,]$front <- "east"
d$front <- factor(d$front, levels=c("west", "front", "east"))


# }
##{ Plot concentrations on top of interpolated phy vars----------------------
                                 
# calculate distances for the binned data matrix                                 
lonRef <- min(phy$long)
d <- ddply(d, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
}, .progress="text")

# plotting just the appendicularians only

Tplot <- ggplot(ti[which(ti$transect!=3),]) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Temp (C)", na.value="grey80", low = "#2d669f", high = "#c8dcef")

Tplot + geom_point(aes(x=dist/1000, y=-depth, size=concentration, colour=concentration>0), alpha=0.7, data=d[d$taxon=="appendicularians",]) + facet_grid(transect~.,) + scale_colour_manual("Presence / Absence", values=c("grey60", "black")) + scale_area("Density", range=c(1,10)) + labs(title="Appendicularians")


# just a simple app concentrations by depth plot
df <- d[which(d$taxon=="appendicularians"),]
ggplot(data=df[which(df$concentration > 0),]) + geom_point(mapping=aes(x=-depth, y=log(concentration)), stat="identity") + labs(title="Appendicularian densities depth profile", x="Depth", y="Log-Density") + xlim(-135, 0) + geom_smooth(aes(x=-depth, y=log(concentration)), method="loess") + coord_flip()

# ggplot(data=df[which(df$concentration > 0),]) + geom_point(mapping=aes(x=concentration, y=-depth), stat="identity") + labs(title="Appendicularian densities depth profile", x="Density", y="Depth") + ylim(-135, 0) + scale_x_continuous(guide=list(guide_axis(), guide_axis(position="top", trans = function(x) x * 2)))

# }

