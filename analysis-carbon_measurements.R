
library("plyr")
library("reshape2")
library("ggplot2")


# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Import data --------------------------------------

d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

da <- read.csv("data/apps_binned_by_depth.csv", stringsAsFactors=FALSE)
da$dateTime <- as.POSIXct(da$dateTime, tz="America/Los_Angeles")

# append apps data to regular data
d <- rbind(d, da)

phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="America/Los_Angeles")

# import interpolated physical data
ti <- read.csv("data/interp_temp.csv", stringsAsFactors=FALSE)
si <- read.csv("data/interp_salinity.csv", stringsAsFactors=FALSE)
swi <- read.csv("data/interp_swRho.csv", stringsAsFactors=FALSE)
fi <- read.csv("data/interp_fluoro.csv", stringsAsFactors=FALSE)
oi <- read.csv("data/interp_oxygen.csv", stringsAsFactors=FALSE)

# calculate distances for the binned data matrix                                 
lonRef <- min(phy$long)
d <- ddply(d, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
}, .progress="text")

# Spectral colour map from ColorBrewer
spectral <- function(n=6) {
  library("RColorBrewer")
  rev(brewer.pal(name="Spectral", n=n))
}

# }

##{ analysis ----------------------------------------------------
# plotting

ggplot(d[d$taxon == "appendicularians" & d$concentration >0,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=taxon), alpha=0.6) + labs(x="Longitude", y="Depth (m)") + facet_grid(transect~.) + scale_colour_manual(values=c("black","red", "yellow")) + scale_area(range=c(1,13)) + guides(size=guide_legend(order=1), colour=guide_legend(order=2))

ggplot(d[d$taxon %in% c("appendicularians", "kowalevskiid", "fritillarid") & d$concentration >0,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=taxon), alpha=0.6) + labs(x="Longitude", y="Depth (m)") +  facet_grid(transect~.) + scale_colour_manual(values=c("black","red", "yellow")) + scale_area(range=c(1,13))

# import length measurements
l <- read.csv("raw_biological_data/length_measurements.csv", stringsAsFactors=F)
l$dateTime <- as.POSIXct(strptime(l$dateTime, format="%m-%d-%y%H%M%S", tz="America/Los_Angeles"))
# why does %OS not work? fix this later

l[l$group=="Siphonophore","group"] <- "Siphonophores"
l[l$group=="Kowaleskiidae","group"] <- "Kowalevskiidae"


# convert from pixels to mm
# 2048 px = 135 mm
# 1 px = 0.06592 mm

l$length <- l$length * 0.06592

ggplot(l[l$group %in% c("Siphonophores" ,"Ctenophores", "Hydromedusae", "Doliolids"),]) + geom_bar(aes(x=length), stat="bin") + facet_grid(~group, scales="free") + labs(x="Length (mm)")

Appendicularians <- c("Oikopleuridae", "Kowalevskiidae", "Fritillaridae")
ggplot(l[l$group %in% Appendicularians,]) + geom_bar(aes(x=length), stat="bin") + facet_grid(~group) + labs(x="Trunk Length (mm)")

# }


##{ jelly biomass -----------------------------------------------
# calculating ballpark C biomass for different jelly taxa

# ballpark average lengths

#weighted app mean
Lapp <- 
  mean(l[l$group == "Oikopleuridae","length"]) * .88 + mean(l[l$group == "Fritillaridae","length"]) * .08 + mean(l[l$group == "Kowalevskiidae","length"]) * .04 
# [1] 1.401407

Lct <- mean(l[l$group == "Ctenophores","length"])
# [1] 12.82307

Lhy <- mean(l[l$group =="Hydromedusae", "length"])
# [1] 2.664116

Lsiph <- mean(l[l$group =="Siphonophores", "length"])
# [1] 8.268967

Ldol <- mean(l[l$group =="Doliolids", "length"])
# [1] 3.540721


# biomass calculations (in micrograms)
# from Lavaniegos and Ohman 2007
d$massSingle <- NA
d[d$taxon %in% c("appendicularians", "kowalevskiid", "fritillarid"),]$massSingle <- 29.49 * (Lapp ^ 2.88) # assume generalized Oikopleura (combining wet weight and C biomass together)
d[d$group=="Ctenophores",]$massSingle <- 4.8 * (Lct ^ 1.775) # this result in micrograms
d[d$group=="Hydromedusae",]$massSingle <- 1.8885 * (Lhy ^ 2.619)
d[d$group=="Siphonophores",]$massSingle <- 20.47 * (Lsiph ^ 0.834)
d[d$taxon=="doliolids",]$massSingle <- 0.51 * (Ldol ^ 2.28)

d$biomass <- d$abund * d$massSingle
massSum <- sum(d$biomass)

# convert to grams
massSum * 10^-6

# }



