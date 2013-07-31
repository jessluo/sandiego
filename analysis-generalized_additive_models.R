#
#      Analysis of the species-environment relationships through
#      generalized additive models
#
#  (c) Copyright 2013 Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library(ggplot2)
library(plyr)
library(reshape2)
library(mgcv)
library(car)
`%ni%` <- Negate(`%in%`) 

##{ Read data --------------------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

# }

##{ Define new groups & delineate the front ---------------------------------
# identify explanatory variables of interest
locVars <- c("depth", "long")
hydroVars <- c("temp", "salinity", "fluoro", "oxygen", "swRho")
vars <- c(locVars, hydroVars)

# define new groups for the analysis
d$group2 <- d$group
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]
d[d$group2=="appendicularians","group2"] <- "Appendicularians"
d[d$group2=="doliolids","group2"] <- "Doliolids"

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + front + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

dg$transect <- factor(dg$transect)
dg$front <- factor(dg$front, levels=c("west", "front", "east"))
# }

##{ Scatterplot and GAM for solmaris--------------------------------------------
dsol <- d[d$taxon=="sol_large",]

dsol <- dsol[,which(names(dsol) %in% c(hydroVars, "swRho", "concentration", "transect"))]

dsol$transect <- factor(dsol$transect)

scatterplotMatrix(dsol, pch=19, cex=0.5, reg.line=F, lwd.smooth=1.25, spread=F, ellipse=T, col=c("gray60", "#2957FF", "#FF8000"), col.axis="gray50")

ggplot(aes(x=salinity, y=concentration), data=dsol) + geom_point(color="#FF8000", alpha=0.75) + facet_wrap(transect~front)

ggplot(aes(x=salinity, y=concentration), data=dsol) + geom_point(color="#FF8000", alpha=0.75) + geom_smooth(se=F, method='gam', formula=y~s(x, bs="cr")) + facet_grid(transect~.) + labs(title="Solmaris large GAM")

# gams will be developed for most abundant species:
# solmaris - large and small, Liriope, Sphaeronectes, and Ocyropsis maculata

# solmaris
sol_gam <- gam(concentration ~ s(temp, bs="cr") + s(salinity, bs="cr") + s(fluoro, bs="cr") + s(oxygen, bs="cr"), data=d[d$taxon=="sol_large",])
summary(sol_gam)
plot(sol_gam, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col="#FF8000", shade=T, shade.col="gray90")

# solmaris - splitting it up by location with respect to the front
sol_gam <- gam(concentration ~ s(temp, bs="cr") + s(salinity, bs="cr") + s(fluoro, bs="cr") + s(oxygen, bs="cr"), data=d[d$taxon=="sol_large" & d$front=="front",])
summary(sol_gam)
sol_gam <- update(sol_gam, . ~ . - s(salinity, bs="cr"))
summary(sol_gam)
sol_gam <- update(sol_gam, . ~ . - s(oxygen, bs="cr"))
summary(sol_gam)
plot(sol_gam, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col="#FF8000", shade=T, shade.col="gray90")
vis.gam(sol_gam, type="response", plot.type="persp", phi=30, theta=30, border=NA)

sol_gam <- gam(concentration ~ s(temp, bs="cr") + s(salinity, bs="cr") + s(fluoro, bs="cr") + s(oxygen, bs="cr"), data=d[d$taxon=="sol_large" & d$front=="west",])
summary(sol_gam)
sol_gam <- update(sol_gam, . ~ . - s(oxygen, bs="cr"))
summary(sol_gam)
plot(sol_gam, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col="#FF8000", shade=T, shade.col="gray90")

sol_gam <- gam(concentration ~ s(temp, bs="cr") + s(salinity, bs="cr") + s(fluoro, bs="cr") + s(oxygen, bs="cr"), data=d[d$taxon=="sol_large" & d$front=="east",])
summary(sol_gam)
sol_gam <- update(sol_gam, . ~ . - s(salinity, bs="cr"))
summary(sol_gam)
sol_gam <- update(sol_gam, . ~ . - s(oxygen, bs="cr") + oxygen)
summary(sol_gam)
plot(sol_gam, pages=1, residuals=T, pch=19, cex=0.25, scheme=1, col="#FF8000", shade=T, shade.col="gray90")
vis.gam(sol_gam, type="response", plot.type="persp")

ggplot(aes(x=swRho, y=concentration), data=d[d$taxon=="sol_large",]) + geom_point(color="#FF8000", alpha=0.75) + geom_rug(sides="b") + geom_smooth(se=F, method='gam', formula=y~s(x, bs="cr")) + facet_grid(transect~front ,scales="free_y")

ggplot(aes(x=salinity, y=concentration), data=d[d$taxon=="sol_large",]) + geom_point(color="#FF8000", alpha=0.75) + geom_rug(sides="b") + geom_smooth(se=F, method='gam', formula=y~s(x, bs="cr")) + facet_grid(transect~front ,scales="free_y")

# }


## { Generalized additive mixed models ------------------------------------------


# }