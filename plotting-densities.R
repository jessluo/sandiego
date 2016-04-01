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
library("gridExtra")

registerDoParallel(cores=detectCores())
parallel <- TRUE

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read data ------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="America/Los_Angeles")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards

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


phy <- ddply(phy, ~transect, function(x, lonR=lonRef) {
  latRef <- mean(x$lat)
  x$dist <- geodDist(x$lat, x$lon, latRef, lonR) * 1000
  # NB: geodDist computes in km
  return(x)
}, .progress="text")

# phy$dist <- phy$dist/1.157
# d$dist <- d$dist/1.157 # some weird correction - Jan 2016


# Spectral colour map from ColorBrewer
spectral <- function(n=6) {
  library("RColorBrewer")
  rev(brewer.pal(name="Spectral", n=n))
}

# }

##{ Create interpolated physical variables plots ----------------------------

# where's the pycnocline
p <- ggplot(phy) + geom_point(aes(x=long, y=-depth, colour=c(0,abs(diff(swRho))))) + facet_grid(transect~.) + scale_colour_gradient("Pycnocline strength", limits=c(0,0.12), low="white", high="black") + theme_bw()

# temperature
tplot <- ggplot(ti) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="black", size=0.5, alpha=0.5, breaks=c(10, 13, 16)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn(expression(paste("Temp. (",degree,"C)")), colours=spectral(), na.value=NA) + theme_bw()

# temperature no contours
tplot <- ggplot(ti) + geom_tile(aes(x=dist/1000, y=-depth, fill=temp)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Temp (C)", colours=spectral(), na.value=NA) + theme_bw()

# seawater density
swplot <- ggplot(swi) + geom_tile(aes(x=dist/1000, y=-depth, fill=swRho-1000)) + geom_contour(aes(x=dist/1000, y=-depth, z=swRho-1000), colour="black", size=1, alpha=0.7, breaks=c("24.5", "25.1")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Seawater Density", colours=spectral(), na.value=NA) + theme_bw() 

# swrho w/ salinity contours
swplot <- ggplot() + geom_tile(aes(x=dist/1000, y=-depth, fill=swRho-1000), data=swi) + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.5, breaks=c("33.3", "33.45"), data=si) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Seawater Density", colours=spectral(), na.value=NA) + theme_bw() 

# salinity
splot <- ggplot(si) + geom_tile(aes(x=dist/1000, y=-depth, fill=salinity)) + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.7, breaks=c("33.3", "33.45")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Salinity", na.value=NA, colours=spectral()) + theme_bw() + theme(legend.position=c(.1, .35))

ggplot(phy) + geom_point(aes(x=long, y=-depth, colour=salinity)) + facet_grid(transect~.) + scale_colour_gradientn(colours=spectral()) + labs(x="Longitude", y="Depth (m)") + theme_bw()

# salinity on top of sampling tracks (fig 2 modification)
# splot <- ggplot() + geom_tile(aes(x=dist/1000, y=-depth, fill=salinity), data=si) + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.7, breaks=c("33.3", "33.45"), data=si) + geom_point(aes(x=dist/1000, y=-depth), colour="grey80", alpha=0.35, size=0.75, data=phy) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Salinity", na.value=NA, colours=spectral()) + theme_bw() + theme(legend.position=c(.1, .35))

# fluorometry
fplot <- ggplot(fi) + geom_tile(aes(x=dist/1000, y=-depth, fill=fluoro)) + geom_contour(aes(x=dist/1000, y=-depth, z=fluoro), colour="black", size=0.5, alpha=0.5, breaks=c("0.2", "0.4", "0.6", "0.8")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Fluoro. (V)", na.value=NA, colours=spectral()) + theme_bw()

# w/o contour
fplot <- ggplot(fi) + geom_tile(aes(x=dist/1000, y=-depth, fill=fluoro)) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Fluoro. (V)", na.value=NA, colours=spectral()) + theme_bw()

# oxygen
oplot <- ggplot(oi) + geom_tile(aes(x=dist/1000, y=-depth, fill=oxygen)) + geom_contour(aes(x=dist/1000, y=-depth, z=oxygen), colour="black", size=0.5, alpha=0.5, breaks=c("2", "3", "4")) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn(expression(paste("Oxygen (ml ",L^-1,")")), na.value=NA, colours=spectral()) + theme_bw()

# oxygen w/o contour
oplot <- ggplot() + geom_tile(aes(x=dist/1000, y=-depth, fill=oxygen), data=oi) + facet_grid(transect~.) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Oxygen", na.value=NA, colours=spectral()) + theme_bw()
# }

##{ Summarize counts and abundances -----------------------------------------
res <- sort(daply(d[d$group=="Ctenophores",], ~taxon, function(x){sum(x$abund)}), decreasing=TRUE)
res <- sort(daply(d[d$group=="Hydromedusae",], ~taxon, function(x){sum(x$abund)}), decreasing=TRUE)
res <- sort(daply(d[d$group=="Siphonophores",], ~taxon, function(x){sum(x$abund)}), decreasing=TRUE)

res <- sort(daply(d[d$group %in% c("Hydromedusae", "Solmaris"),], ~taxon, function(x){sum(x$abund)}), decreasing=TRUE)
res/sum(res)
# }

##{ Plot concentrations on top of interpolated phy vars----------------------
# according to GBM results

# Plots for the paper

# solmundella - temperature
Solmundella <- tplot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="h6_Solmundella" & d$concentration > 0,]) + 
  facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + 
  labs(title=expression(paste(italic("Solmundella bitentaculata")))) + 
  theme(legend.position=c(.1, .35), plot.margin=unit(c(1,0,0.5,0), "lines")) + 
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

# Pegantha - temperature + fluorometry: fluoro base plot + isothermals
Pegantha <- fplot + 
  geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="black", size=1, alpha=0.4, breaks=c(10, 13, 16), data=ti) + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="h7_Pegantha" & d$concentration > 0,]) +
  scale_size_area("Density", max_size=7) + labs(title=expression(paste(italic("Pegantha")," sp."))) + 
  theme(legend.position=c(.1, .36), axis.title.y=element_blank()) + 
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

pdf("plots/bubble/Solmundella_Pegantha.pdf", height=8.5, width=11)
grid.arrange(Solmundella, Pegantha, nrow=1, widths=c(1,1.03))                                                                                                                                                                                                                                                                             
dev.off()

# solmaris - temperature
Solmaris <- tplot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="Solmaris" & d$concentration > 0,]) +
  facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + 
  labs(title=expression(paste(italic("Solmaris rhodoloma")))) + 
  theme(legend.position=c(.1, .35), plot.margin=unit(c(1,0,0.5,0), "lines")) +
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

# doliolids - fluoro
Doliolids <- fplot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="doliolids" & d$concentration > 0,]) +
  facet_grid(transect~.,) + scale_size_area("Density", max_size=7) + 
  labs(title="Doliolids") + theme(legend.position=c(.1, .38), axis.title.y=element_blank()) +
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

pdf("plots/bubble/Solmaris_doliolids.pdf", height=8.5, width=11)
grid.arrange(Solmaris, Doliolids, nrow=1, widths=c(1,1.03))                                                                                                                                                                                                                                                                             
dev.off()

# Sphaeronectes - temp + oxygen
Sphaeronectes <- oplot + 
  geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="black", size=1, alpha=0.4, breaks=c(10, 13, 16), data=ti) + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="Sphaeronectes" & d$concentration > 0,]) +
  facet_grid(transect~.,) + scale_size_area("Density", max_size=9) + 
  labs(title=expression(paste(italic("Sphaeronectes")," sp."))) + 
  theme(legend.position=c(.1, .39), plot.margin=unit(c(1,0,0.5,0), "lines")) +
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

# Ocyropsis - oxygen
Ocyropsis <- oplot +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="Ocyropsis maculata" & d$concentration > 0,]) + 
  facet_grid(transect~.,) + scale_size_area("Density", max_size=7) + 
  labs(title=expression(paste(italic("Ocyropsis maculata")))) + 
  theme(legend.position=c(.1, .35), axis.title.y=element_blank()) +
  guides(size=guide_legend(order=1), fill=guide_colourbar(order=2))

pdf("plots/bubble/Sphaeronectes_Ocyropsis.pdf", height=8.5, width=11)
grid.arrange(Sphaeronectes, Ocyropsis, nrow=1, widths=c(1,1.03))  
dev.off()

# h15 - temperature
h15 <- tplot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h15" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="h15")


# vsh - (depth) + temperature
vsh <- tplot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="vsh" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="Very Small Hydromedusae")

# apps - depth
apps <- ggplot() + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="appendicularians" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="Appendicularians") + theme_bw()
                                                                                                                                                                                                                                                                            
# Liriope - temperature + fluoro
Liriope <- tplot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.4, data=d[d$taxon=="h5_Liriope" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title=expression(paste(italic("Liriope tetraphylla"))))

Liriope2 <- tplot + 
  geom_contour(aes(x=dist/1000, y=-depth, z=fluoro), colour="black", size=1, alpha=0.5, breaks=c(0.2, 0.4, 0.6, 0.8), data=fi) + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h5_Liriope" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="Liriope tetraphylla")

Liriope2 <- fplot + 
  geom_contour(aes(x=dist/1000, y=-depth, z=temp), colour="black", size=1, alpha=0.5, breaks=c(10, 15), data=ti) + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h5_Liriope" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="Liriope tetraphylla")


# physonect - fluoro
Physonect <- fplot + 
  ggplot() + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="Physonect" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_area("Density", range=c(1,10)) + labs(title="Physonect")

# diphyidae - fluoro
Diphyidae <- fplot + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.7, breaks=c("33.3", "33.45"), data=si) +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="Diphyidae" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Diphyidae")

# NEW SIPHS
# muggiaea - fluoro
muat <- fplot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="muat" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Muggiaea")

# lensia - oxygen / salinity
lemu <- oplot +  geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.7, breaks=c("33.3", "33.45"), data=si) +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="lemu" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Lensia sp.")

# nanomia - oxygen / salinity
nabi <- oplot +  geom_contour(aes(x=dist/1000, y=-depth, z=fluoro), colour="black", size=1, alpha=0.5, breaks=c(0.2, 0.4, 0.6, 0.8), data=fi) + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="nabi" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Nanomia bijuga")

# agalma - fluoro + oxygen
agel <- fplot +  geom_contour(aes(x=dist/1000, y=-depth, z=oxygen), colour="black", size=1, alpha=0.5, breaks=c("2", "3", "4", "5"), data=oi) +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="agel" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Agalma elegans")


# Beroida - temp
Beroida <- tplot +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="Beroida" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density", max_size=10) + labs(title="Beroida")

# Cunina - temp
Cunina <- tplot +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.7, data=d[d$taxon=="h3_Cunina" & d$concentration > 0,]) + facet_grid(transect~.,) + scale_size_area("Density") + labs(title="Cunina")



# just a simple app concentrations by depth plot
df <- d[which(d$taxon=="appendicularians"),]
ggplot(data=df[which(df$concentration > 0),]) + geom_point(mapping=aes(x=-depth, y=log(concentration)), stat="identity") + labs(title="Appendicularian densities depth profile", x="Depth", y="Log-Density") + xlim(-135, 0) + geom_smooth(aes(x=-depth, y=log(concentration)), method="loess") + coord_flip()

# }

##{ Violin plots ------------------------------------------------------------

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

### VIOLIN PLOT
# calculate an average value for each bin by frontal region
# define which groups you want
d$group2 <- d$group
# d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]

d2 <- d[,names(d) %in% c("transect", "cast", "down.up", "dateTime", "dateTimeB", "depth", "front", "group", "group2", "taxon", "abund", "concentration")]

d2$depthBin <- floor(d$depth)

d2 <- ddply(d2, ~transect + front + depthBin + group2, function(x){
  avgConc <- mean(x$concentration)
  return(avgConc)
}, .progress="text")

d2 <- rename(d2, c("depthBin" = "depth", "V1" = "concentration"))

d2$group2 <- factor(d2$group2, levels=c("Solmaris", "Hydromedusae", "Siphonophores", "Ctenophores", "appendicularians", "doliolids"))

# }

## { TEMPORARY CODE to plot 4 taxa on top of salinity -----------------------------------
si <- si[si$transect==2,]
phy <- phy[phy$transect==2,]
d <- d[d$transect==2,]

mindist <- min(si[!is.na(si$salinity),]$dist)

d$dist <- d$dist - mindist
si$dist <- si$dist - mindist
phy$dist <- phy$dist - mindist

splot <- ggplot(si) + geom_tile(aes(x=dist/1000, y=-depth, fill=salinity)) + geom_contour(aes(x=dist/1000, y=-depth, z=salinity), colour="black", size=1, alpha=0.7, breaks=c("33.3", "33.45")) + labs(x="Distance (km)", y="Depth (m)") + scale_fill_gradientn("Salinity", na.value=NA, colours=spectral()) + theme_bw() + scale_x_continuous(limits=c(-5.5, 72)) + scale_y_continuous(limits=c(-138,0))

Solmaris <- splot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.4, data=d[d$taxon=="Solmaris" & d$concentration > 0,]) +
  scale_size_area("Density", max_size=15) + 
  labs(title=expression(paste(italic("Solmaris rhodoloma")))) + 
  theme(legend.position=c(.07, .35), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(size=guide_legend(order=1), fill="none")

Liriope <- splot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.4, data=d[d$taxon=="h5_Liriope" & d$concentration > 0,]) +
  scale_size_area("Density", max_size=9) + 
  labs(title=expression(paste(italic("Liriope tetraphylla")))) +
  theme(legend.position=c(.06, .35)) +  
  guides(size=guide_legend(order=1), fill="none")


Sphaeronectes <- splot + 
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.4, data=d[d$taxon=="Sphaeronectes" & d$concentration > 0,]) +
  scale_size_area("Density", max_size=7) + 
  labs(title=expression(paste(italic("Sphaeronectes")," sp."))) + 
  theme(legend.position=c(.06, .5), axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + 
  guides(size=guide_legend(order=1), fill="none")

Ocyropsis <- splot +
  geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.4, data=d[d$taxon=="Ocyropsis maculata" & d$concentration > 0,]) + 
  scale_size_area("Density", max_size=4) + 
  labs(title=expression(paste(italic("Ocyropsis maculata")))) + 
  theme(legend.position=c(.06, .45), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + guides(size=guide_legend(order=1), fill="none")

pdf("plots/bubble/Solmaris_Liriope_Sphaeronectes_Ocyropsis_Salinity.pdf", height=8.75, width=12.5)
grid.arrange(Solmaris, Sphaeronectes, Liriope, Ocyropsis, nrow=2, ncol=2)  
dev.off()

# one more - appendicularians
apps <- splot + geom_point(aes(x=dist/1000, y=-depth, size=concentration), alpha=0.5, data=d[d$taxon=="appendicularians" & d$concentration > 0,]) + scale_size_area("Density", max_size=9) + labs(title="Appendicularians") + theme_bw() + theme(legend.position=c(0.06, 0.45))

apps <- apps + scale_x_continuous(limits=c(11, 89))

ggsave("plots/bubble/Appendicularians_salinity_tr2.pdf", apps, height=4, width=11)


# }

