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
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
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
plottheme <- theme(title=element_text(size=rel(2)), plot.title=element_text(size=rel(1.8)), axis.text=element_text(size=rel(1.5)), strip.text=element_text(size=rel(2)), legend.text=element_text(size=rel(1.5)))
# plottheme <- theme(title=element_text(size=20), plot.title=element_text(size=24), axis.text=element_text(size=14), strip.text=element_text(size=14), legend.text=element_text(size=12))

# general density plots
ggplot(d) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

# density plots with three transects separated
ggplot(d) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + facet_grid(transect~.) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

# density plot without appendicularians, with and without transects separated
ggplot(d[d$taxon != "appendicularians",]) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")

ggplot(d[d$taxon != "appendicularians",]) + geom_density(aes(x=-depth, weight=concentration, colour=group2), alpha=0.5) + facet_grid(transect~.) + coord_flip() + labs(y="Density", x="Depth", title="Depth distribution", colour="Taxon")


# DISTRIBUTION PLOTS

# salinity
splot + labs(title="Salinity") + plottheme + theme(legend.position=c(.1, .35))

# solmaris, large size fraction
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="sol_large" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=13) + labs(title="Solmaris rhodoloma, large size fraction", y="Depth (m)") + plottheme + theme(legend.position=c(.1, .35))

# solmaris, small size fraction
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="sol_small" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=13) + labs(title="Solmaris rhodoloma, small size fraction", y="Depth (m)") + plottheme + theme(legend.position=c(.1, .35))

# appendicularians
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="appendicularians" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=15) + labs(title="Appendicularians", y="Depth (m)", x="Distance from start of transect (km)") + plottheme + theme(legend.position=c(.1, .35))

# doliolids
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="doliolids" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Doliolids", y="Depth (m)", x="Distance from start of transect (km)") + plottheme + theme(legend.position=c(.1, .35))

# hydromedusa
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration, colour=taxon), alpha=0.7, data=d[d$group=="Hydromedusae" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Hydromedusae", y="Depth (m)", x="Distance from start of transect (km)") + plottheme + theme(legend.position=c(.1, .35))


# violin plot comparing depth distribution of all taxa with respect to transect and position relative to the front
ggplot(d[d$concentration>0,]) + geom_violin(aes(x=front, y=-depth, weight=concentration, colour=group2), alpha=0.4) + facet_grid(transect~group2) + labs(colour="Taxon")


ggplot(d[d$taxon=="appendicularians",]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

# siphonophores
dS <- d[d$group=="Siphonophores",]

dS$taxon <- factor(dS$taxon, levels=c("Physonect", "Diphyidae", "Sphaeronectes", "Lilyopsis", "Prayidae"))

ggplot(dS) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

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
dct$taxon <- factor(dct$taxon, levels=c("Ocyropsis maculata", "Juvenile Lobata", "Thalassocalycidae inconstans", "Larval Lobata", "Haeckelia beehlri", "Beroida", "Mertensid", "Hormiphora californiensis", "Velamen"))

# all data together, Ctenophore densities by taxon
ggplot(dct) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))


# violin plots with constant widths (scaled to the same width) -- this is appropriate for presenting the CA analysis
ggplot(dct[dct$concentration>0,]) + geom_violin(aes(x=taxon, y=-depth, weight=concentration, colour=taxon), alpha=0.7, scale="width") + labs(colour="Taxa", y="Depth (m)", x="Taxa") + theme(title=element_text(size=18), plot.title=element_text(size=22), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_text(size=14), legend.text=element_text(size=11))

# Ocyropsis and larval lobates only
ggplot(dct[dct$taxon %in% c("Ocyropsis maculata", "Larval Lobata"),]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))


# }

##{ Liriope tetraphylla DVM -----------------------------------------------

# plot simple Liriope plot 
ggplot(d[d$taxon=="h5_Liriope",]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

# plot Liriope concentrations on top of salinity profile
splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h5_Liriope" & d$concentration>0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=13) + labs(title="Liriope tetraphylla") + plottheme + theme(legend.position=c(.1, .35))

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
# --> in transect 1, move from 55 m to 15 m in 2 hrs --> 33 cm / min migration rate
# --> max rate of change 


# }

##{ jelly biomass -----------------------------------------------
# calculating ballpark C biomass for different jelly taxa

# ballpark average lengths
# app length: 1 mm
# ctenophore length: 10 mm
# hydromedusae length: 2.5 mm
# solmaris large: 4 mm
# solmaris small: 0.5 mm
# siphonophore nectophore length: 7.5 mm
# doliolid length: 3.5 mm

Lapp <- 1.5 # 0

# low average
Lct <- 10
Lhy <- 3.5
Lsl <- 5
Lsm <- 1
Lsiph <- 8
Ldol <- 3.5

# high average length
Lct <- 25
Lhy <- 4.5
Lsl <- 8
Lsm <- 2
Lsiph <- 15
Ldol <- 4.5

# biomass calculations (in micrograms)
# from Lavaniegos and Ohman 2007
d$massSingle <- NA
d[d$taxon=="appendicularians",]$massSingle <- 29.49 * (Lapp ^ 2.88) # assume generalized Oikopleura (combining wet weight and C biomass together)
d[d$group=="Ctenophores",]$massSingle <- 4.8 * (Lct ^ 1.775) # this result in micrograms
d[d$group=="Hydromedusae",]$massSingle <- 1.8885 * (Lhy ^ 2.619)
d[d$taxon=="sol_large",]$massSingle <- 1.8885 * (Lsl ^ 2.619)
d[d$taxon=="sol_small",]$massSingle <- 1.8885 * (Lsm ^ 2.619)
d[d$group=="Siphonophores",]$massSingle <- 20.47 * (Lsiph ^ 0.834)
d[d$taxon=="doliolids",]$massSingle <- 0.51 * (Ldol ^ 2.28)

d$biomass <- d$abund * d$massSingle
massSum <- sum(d$biomass)

# scaled up to include potential full hydro, app, doliolid and ctenophore data
abund_appuncounted <- sum(d[d$group=="Tunicates" & d$transect==1,]$abund) +
  sum(d[d$group=="Tunicates" & d$transect==2,]$abund)
abund_uncounted <- sum(d[d$group %in% c("Tunicates", "Hydromedusae") & d$transect==1,]$abund) +
  sum(d[d$group %in% c("Tunicates", "Hydromedusae") & d$transect==2,]$abund) +
  sum(d[d$group=="Ctenophores",]$abund) * 1.615
biomass_uncounted <- sum(d[d$group %in% c("Tunicates", "Hydromedusae") & d$transect==1,]$biomass) +
  sum(d[d$group %in% c("Tunicates", "Hydromedusae") & d$transect==2,]$biomass) +
  sum(d[d$group=="Ctenophores",]$biomass) * 1.615

# convert to grams
(massSum + biomass_uncounted)* 10^-6

# }

