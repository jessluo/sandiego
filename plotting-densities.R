#
#      Plotting the biological data
#      
#
#--------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")



##{ Read data ------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")

ti <- read.csv("data/interp_temp.csv", stringsAsFactors=FALSE)

# }

##{ Plot concentrations on top of interpolated phy vars----------------------
                                 
# calculate distances for the binned data matrix                                 
lonRef <- min(phy$long)
d <- ddply(d, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
})

# plotting just the appendicularians only

Tplot <- ggplot(ti[which(ti$transect!=3),]) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + geom_contour(aes(x=dist, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Temp (C)", na.value="grey80", low = "#2d669f", high = "#c8dcef")

Tplot + geom_point(aes(x=dist/1000, y=-depth, size=concentration, colour=concentration>0), alpha=0.7, data=d[d$taxon=="appendicularians",]) + facet_grid(transect~.,) + scale_colour_manual("Presence / Absence", values=c("grey60", "black")) + scale_area("Density", range=c(1,10)) + labs(title="Appendicularians")


# just a simple app concentrations by depth plot
df <- d[which(d$taxon=="appendicularians"),]
ggplot(data=df[which(df$concentration > 0),]) + geom_point(mapping=aes(x=log(concentration), y=-depth), stat="identity") + labs(title="Appendicularian densities depth profile", x="Log-Density", y="Depth") + ylim(-135, 0) + geom_smooth(aes(x=log(concentration), y=-depth), method="glm")

ggplot(data=df[which(df$concentration > 0),]) + geom_point(mapping=aes(x=concentration, y=-depth), stat="identity") + labs(title="Appendicularian densities depth profile", x="Density", y="Depth") + ylim(-135, 0)

# }

