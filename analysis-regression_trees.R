#
#      Analysis of the species-environment relationships through
#      regression trees
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("vegan")
library("mvpart")
library("plyr")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())

##{ Read and reformat data ------------------------------------------------

# read data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards


# identify explanatory variables of interest
locVars <- c("depth", "long")
hydroVars <- c("temp", "salinity", "fluoro", "oxygen")
vars <- c(locVars, hydroVars)

# define new groups for the analysis
d$group2 <- d$group
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

ggplot(dg[dg$concentration > 0,]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group2, scale="free")
ggplot(dg[dg$concentration > 0,]) + geom_histogram(aes(x=log(concentration))) + facet_wrap(~group2, scale="free")

# }

##{ Per group trees -------------------------------------------------------

m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("none"))
m <- mvpart(concentration ~ depth + long, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)

# apps are not very interesting
m <- mvpart(concentration ~ depth + long + temp + salinity + fluoro + oxygen, data=dg[dg$group2=="appendicularians",], xv=("min"), xval=5, xvmult=100)

# hydromedusae
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Hydromedusae",], xv=("1se"), xval=5, xvmult=100)

# siphonophores
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Siphonophores",], xv=("min"), xval=5, xvmult=100)

# ctenophores - cannot explain ctenophores with this method because they are not numerous enough. can only use logistic regression
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Ctenophores",], xv=("min"), xval=5, xvmult=100)


# }

##{ Community structure ---------------------------------------------------

# let's look at the influence of the different sizes of Solmaris on each other
dSol <- d[d$group=="Solmaris",]
dSol <- dcast(dSol, dateTimeB ~ taxon, value.var="concentration")
dSol <- join(dSol, unique(d[,c("dateTimeB", vars)]))

m <- mvpart(sol_small ~ temp + salinity + fluoro + oxygen + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)
m <- mvpart(sol_small ~ depth + long + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)

dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + fluoro + oxygen ~ group2, value.var="concentration")

# multivariate regression tree with all taxa
m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("appendicularians", "Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# basically, when you look at the multivariate regression tree for everything, it is 27% explained by temperature (which is probably just a proxy for depth) and then nothing else

# multivariate regression tree without the apps
m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# much better. though weird/repeating color scale. and difficult to see as it is overlaying one of the branches

# multivariate regression tree with just the jellies and the ctenophores
m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)

# multivariate regression tree with just the hydrozoans
m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# the latter three are about the same
# think about the interpretation for these plots ...


# multivariate regression tree with Principal components analysis of the group means
m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100, pca=T)
# difficult to interpret

# }


##{ Unconstrained ordination - CA  -----------------------------

# Correspondence Analysis
# is it going to be better to transform the concentration values? probably will be fine on either transformed values or raw values. can do a sqrt transformation to dampen the effect of dominant species.

# start with CA with no transformation to see how it performs
dCspp <- dgC[,names(dgC) %in% c("appendicularians", "Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# CA requires all row sums to be >0 in the community data matrix
# check to see which row sums = 0 and remove them
which(rowSums(dCspp) == 0)
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
# remove the NAs
dCspp <- na.omit(dCspp)
# perform the CA using cca() from the vegan library
allCA <- cca(dCspp)
# check the summary
head(summary(allCA)) # first two CA axes explain 56% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")

# removing the appendicularians from the CA
dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 59% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# removing the appendicularians improves the proportion explained of the first CA a little bit, but not by much

# just including the hydrozoans + the ctenophores in the CA
dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 67% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")

# I'm not sure which one I like better


# CA with sqrt transformation
dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# square root transform
dCspp <- sqrt(dCspp)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 64% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# so the proportion of variance explained does not change too much based on the transformation or no transformation, but the biplot changes quite a bit
# maybe because the difference between the large values of the large solmaris compared to the 0's and 1's of the ctenophores and hydromedusae made the first set of plots on the untransformed data have this "V" like shape -- which is less pronounced in this plot.

# CA with log transformation
dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# log transform
dCspp <- log(dCspp+1)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 65% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# again, the proportion of variance explained does not change too much, but the biplot changes a lot. and also changes the interpretation of which species are closer to each other. Prefer one over the other?

# }


# TODO:
# CCA or RDA with species ~ environment
# logistic regression on P/A data
