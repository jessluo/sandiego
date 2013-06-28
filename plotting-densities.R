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
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="America/Los_Angeles")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards


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

# everything together

# identify explanatory variables of interest
locVars <- c("depth", "long")
hydroVars <- c("temp", "salinity", "fluoro", "oxygen", "swRho")
vars <- c(locVars, hydroVars)
# define new groups for the analysis
d$group2 <- d$group
# d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + front + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)


# casting data into wide format with group 2
dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + front + depth + long + temp + salinity + swRho + fluoro + oxygen ~ group2, value.var="concentration")


# violin plots
ggplot(d[d$concentration>0,]) + geom_violin(aes(x=group2, y=-depth, weight=concentration, colour=group2), alpha=0.4) + facet_grid(transect~front)

ggplot(d[d$concentration>0,]) + geom_violin(aes(x=front, y=-depth, weight=concentration, colour=group2), alpha=0.4) + facet_grid(transect~group2) + labs(colour="Taxon")


ggplot(d) + geom_path(aes(x=-depth, y=concentration, colour=group2), alpha=0.5) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

ggplot(d[d$taxon != "appendicularians",]) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + facet_grid(transect~.) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

dd <- data.frame(x=1, y=sort(rnorm(100)), z=1:100)
ggplot(dd) + geom_violin(aes(x=x, y=y))
ggplot(dd) + geom_violin(aes(x=x, y=y, weight=z))



ggplot(d[d$concentration>0 & d$group != "Tunicates",]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=group), alpha=0.4) + facet_grid(transect~.)


ggplot(d[d$concentration>0,]) + geom_bar(aes(x=-depth, fill=group2), stat="bin", position="stack") + coord_flip() + facet_grid(transect~.)

# why doesn't this work?
# need to sum the abundances
ggplot(d[d$concentration>0,]) + geom_bar(aes(x=-depth, y=concentration, fill=group), stat="identity", position="stack") + coord_flip() + facet_grid(transect~.)

# pred / prey
narco <- c("h10_Pegantha", "h3_Cunina", "h6_Solmundella", "h7_Pegantha", "r4_Aegina", "sol_large", "sol_small")
trachy <- c("h11_Haliscera", "h2_Haliscera", "h5_Liriope", "h5b", "h7_Rhopalonema", "h9_Aglaura", "h9_Arctapodema")
pred <- c(narco, "Beroida")

ggplot(d[d$taxon %in% narco,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

ggplot(d[d$taxon %in% trachy,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

ggplot() + geom_point(aes(x=long, y=-depth, size=concentration), colour="pink", alpha=0.7, data=d[d$group=="Ctenophores" & d$taxon !="Beroida" & d$concentration > 0,]) + facet_grid(transect~.) + scale_area(range=c(1,10))

ggplot() + geom_point(aes(x=long, y=-depth, size=concentration), colour="pink", alpha=0.7, data=d[d$group=="Ctenophores" & d$taxon !="Beroida" & d$concentration > 0,]) + geom_point(aes(x=long, y=-depth, size=concentration), colour="black", alpha=0.7, data=d[d$group=="Ctenophores" & d$taxon=="Beroida" & d$concentration >0,]) + facet_grid(transect~.) + scale_area(range=c(1,10))
# }


