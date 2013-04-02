#
#      Plots and figures for student seminar 2013
#      Jessica Luo
#
#--------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("oce")
`%ni%` <- Negate(`%in%`) 
library("pastecs")


##{ Read data ------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards


ti <- read.csv("data/interp_temp.csv", stringsAsFactors=FALSE)
si <- read.csv("data/interp_salinity.csv", stringsAsFactors=FALSE)
swi <- read.csv("data/interp_swRho.csv", stringsAsFactors=FALSE)
fi <- read.csv("data/interp_fluoro.csv", stringsAsFactors=FALSE)
oi <- read.csv("data/interp_oxygen.csv", stringsAsFactors=FALSE)

# }

##{ Define new groups & delineate the front ---------------------------------

# define new groups for the analysis
d$group2 <- d$group
# d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]

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


##{ Create interpolated physical variables plots ----------------------------

# calculate distances for the binned data matrix                                 
lonRef <- min(phy$long)
d <- ddply(d, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
}, .progress="text")

# temperature
tplot <- ggplot(ti) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="white", size=0.5, alpha=0.5, breaks=c(10, 15)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Temp (C)", na.value="grey80", low = "#2d669f", high = "#c8dcef")

# seawater density
swplot <- ggplot(swi) + geom_tile(aes(x=dist/1000, y=-depth, fill=swRho)) + geom_contour(aes(x=dist/1000, y=-depth, z=swRho), colour="white", size=0.5, alpha=0.5, breaks=c("1025", "1026", "1027")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Seawater Density", na.value="grey80", high = "#245280", low = "#a9c8e7")

# salinity
splot <- ggplot(si) + geom_tile(aes(x=dist/1000, y=-depth, fill=salinity)) + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="white", size=0.5, alpha=0.5, breaks=c("33", "33.25", "33.5", "33.75", "34")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Salinity", na.value="grey80", high = "#5b4581", low = "#d5cde4")

# fluorometry
fplot <- ggplot(fi) + geom_tile(aes(x=dist/1000, y=-depth, fill=fluoro)) + geom_contour(aes(x=dist/1000, y=-depth, z=fluoro), colour="white", size=0.5, alpha=0.5, breaks=c("0.2", "0.4", "0.6", "0.8")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Fluorometry", na.value="grey80", high = "#4b6b5d", low = "#dfe9e5")

# oxygen
oplot <- ggplot(oi) + geom_tile(aes(x=dist/1000, y=-depth, fill=oxygen)) + geom_contour(aes(x=dist/1000, y=-depth, z=oxygen), colour="white", size=0.5, alpha=0.5, breaks=c("2", "3", "4", "5")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth") + scale_fill_gradient("Oxygen", na.value="grey80", low = "#e9dfe3", high = "#825b6c")
# }


##{ Generate plots --------------------------------------------------------

# general density plots
ggplot(d) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

# density plots with three transects separated
ggplot(d) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + facet_grid(transect~.) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

# density plot without appendicularians, with and without transects separated
ggplot(d[d$taxon != "appendicularians",]) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

ggplot(d[d$taxon != "appendicularians",]) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + facet_grid(transect~.) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

# violin plot comparing depth distribution of all taxa with respect to transect and position relative to the front
ggplot(d[d$concentration>0,]) + geom_violin(aes(x=front, y=-depth, weight=concentration, colour=group2), alpha=0.4) + facet_grid(transect~group2) + labs(colour="Taxon")

# }

##{ Ctenophore assemblage analysis -------------------------------------------------
#  CA on ctenophores

# subsetting and casting data into wide format
ds <- d[d$group=="Ctenophores",]
dsC <- dcast (ds, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + swRho + fluoro + oxygen ~ taxon, value.var="concentration")
dCspp <- dsC[,names(dsC) %in% c("Beroida", "Bolinopsis", "Charistephane", "Dryodora glandiformis", "Haeckelia beehlri", "Hormiphora californiensis", "Juvenile Lobata", "Larval Lobata", "Mertensid", "Ocyropsis maculata", "Pleurobrachia", "Thalassocalycidae inconstans", "Velamen")]
colSums (dCspp)
# --> Pleurobrachia, Charistephane, Bolinopsis and Dryodora are the least common species

# select some taxa using Escoufier’s equivalent vectors -- escouf on dCspp 
ct.escouf <- escouf(dCspp)
plot(ct.escouf)
# --> these taxa are all different, just arbitrarily select some

# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# log transform
dCspp <- log(dCspp+1)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 29% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# axes 2 and 3
plot(allCA, scaling=2, main="CA biplot of species concentrations", choice=2:3)
text(allCA, dis="sp", col="red", choice=2:3)
#

# subset the ctenophores
rare <- c("Bolinopsis", "Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown")
dct <- d[d$group=="Ctenophores" & d$taxon %ni% rare,]

# initiate
dct$assemblage <- NA

# define the assemblages in the data frame
a1 <- c("Ocyropsis maculata", "Juvenile Lobata", "Thalassocalycidae inconstans")
a2 <- c("Beroida", "Larval Lobata", "Haeckelia beehlri", "Mertensid")
a3 <- "Hormiphora californiensis"
a4 <- "Velamen"
dct[dct$taxon %in% a1,]$assemblage=1
dct[dct$taxon %in% a2,]$assemblage=2
dct[dct$taxon %in% a3,]$assemblage=3
dct[dct$taxon %in% a4,]$assemblage=4

# order the taxa according to assemblage
dct$taxon <- factor(dct$taxon, levels=c("Ocyropsis maculata", "Juvenile Lobata", "Thalassocalycidae inconstans", "Beroida", "Larval Lobata", "Haeckelia beehlri", "Mertensid", "Hormiphora californiensis", "Velamen"))

# all data together, Ctenophore densities by taxon
ggplot(dct) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

# violin plots with constant widths (scaled to the same width) -- this is appropriate for presenting the CA analysis
ggplot(dct[dct$concentration>0,]) + geom_violin(aes(x=taxon, y=-depth, weight=concentration, colour=factor(assemblage)), alpha=0.7, scale="width") + labs(colour="Taxon")

# Ocyropsis and larval lobates only
ggplot(dct[dct$taxon %in% c("Ocyropsis maculata", "Larval Lobata"),]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))


# }

##{ Liriope tetraphylla DVM -----------------------------------------------

# plot simple Liriope plot 
ggplot(d[d$taxon=="h5_Liriope",]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

# plot Liriope concentrations on top of salinity profile
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h5_Liriope" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=13) + labs(title="Liriope tetraphylla")

# Calculation of Liriope rate of movement
# how fast are the liriope tetraphylla moving up in the water column?
dlt <- d[d$taxon=="h5_Liriope" & d$concentration >0,]
dlt$taxon <- "Liriope tetraphylla"

# how do you find the center of mass?
median(dlt[dlt$transect==1 & dlt$cast==1,]$depth)
# --> 55.47875
median(dlt[dlt$transect==1 & dlt$cast==2,]$depth)
# --> 42.01025

# use ddply to calculate median depth and time within each cast
medLiriope <- ddply(dlt, ~transect + cast, function(x){
  medDepth <- median(x$depth)
  medTime <- median(x$dateTime)
  med <- c(medDepth, medTime)
  return(med)
})

# rename column names
medLiriope <- rename(medLiriope, c("V1"="medDepth", "V2"="medTime"))

# convert to time format
medLiriope$medTime <- as.POSIXct(medLiriope$medTime, origin="1970-01-01 00:00:00.00 GMT", tz="GMT")

# calculate the depth and time differences
medLiriope$depthdiff <- c(NA, -diff(medLiriope$medDepth))
medLiriope$timediff <- c(NA, diff(medLiriope$medTime))

# assign NAs at the beginning of each transect
medLiriope[medLiriope$transect==2 & medLiriope$cast==5, which(names(medLiriope) %in% c("depthdiff", "timediff"))] <- NA
medLiriope[medLiriope$transect==3 & medLiriope$cast==1, which(names(medLiriope) %in% c("depthdiff", "timediff"))] <- NA

# calculate rate of change
medLiriope$rate <- medLiriope$depthdiff / medLiriope$timediff

# plot it up
ggplot(data=medLiriope) + geom_bar(aes(x=medTime, y=rate), stat="identity") + facet_grid(transect~., scales="free_x")
# --> why does the scales="free_x" not work?

# TODO: SCALE RATE OF CHANGE BY CONCENTRATION OF INDIVIDUALS. HOW TO DO THIS?


# }

