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
library("packfor")
library("plyr")
library("reshape2")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())

##{ Read and reformat data ------------------------------------------------

# read data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards

# rename env var names
d <- rename(d, replace=c("salinity" = "sal", "oxygen" = "oxy", "fluoro" = "flu", "swRho" = "rho"))

# identify explanatory variables of interest
locVars <- c("depth", "long", "front")
hydroVars <- c("temp", "sal", "flu", "oxy")
vars <- c(locVars, hydroVars)

# define broad groups for the analysis
d$group2 <- d$group
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]
d[d$group2=="appendicularians","group2"] <- "Appendicularians"
d[d$group2=="doliolids","group2"] <- "Doliolids"

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# visualize
# ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group2, scale="free")
# ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=log(concentration))) + facet_wrap(~group2, scale="free")

# casing data into wide format with group 2
dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + sal + flu + oxy ~ group2, value.var="concentration")

# }

## { Define secondary groups for analysis -------------------------------
# use the ordination results to define groups for analysis
# hydromedusae: excluding the species that are present in less than 1% of the population (and keeping broad taxonomic groups), these are the following groups:
# 1) deep trachymedusae: h2-haliscera, h7-rhopalonema, h9-aglaura, h11-haliscera, 2) liriope and h5b, 3) shallow narcos: h7-pegantha, h6-solmundella, 4) deep narco's -- h3-cunina, and 5) other hydros: h1, h15, vsh
# ctenophores: keep broad taxonomic groups: cydippids, lobates, beroe, thalasso, velamen
# siphs: keep groups: physonectae, diphyidae, sphaeronectes, prayidae

# select downcasts for inclusion
dcchoose <- unique(d[d$group=="Hydromedusae" & d$transect==1,"cast"])

# define functional/taxonomic groups
Deep_Trachy <- c("h11_Haliscera", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura")
Shallow_Trachy <- c("h5_Liriope", "h5b")
Shallow_Narco <- c("h7_Pegantha", "h6_Solmundella")
Deep_Narco <- c("h3_Cunina")
Other_Hydro <- c("h1", "h15", "vsh")
Cydippida <- c("Haeckelia beehlri", "Hormiphora californiensis", "Mertensid")
Lobata <- c("Bolinopsis", "Ocyropsis maculata", "Juvenile Lobata", "Larval Lobata")
Prayidae <- c("Lilyopsis", "Prayidae")

# assign these groups into a different column
d$group3 <- d$taxon
d$group3[d$taxon %in% Deep_Trachy] <- "Deep Trachy"
d$group3[d$taxon %in% Shallow_Trachy] <- "Shallow Trachy"
d$group3[d$taxon %in% Shallow_Narco] <- "Shallow Narco"
d$group3[d$taxon %in% Deep_Narco] <- "Deep Narco"
d$group3[d$taxon %in% Other_Hydro] <- "Other Hydro"
d$group3[d$taxon %in% Cydippida] <- "Cydippida"
d$group3[d$taxon %in% Lobata] <- "Lobata"
d$group3[d$taxon %in% Prayidae] <- "Prayidae"

# exclude rare taxa
`%ni%` <- Negate(`%in%`) 

exclude <- c("Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown", "Annatiara", "h10_Pegantha", "h13", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia")
d <- d[d$group3 %ni% exclude,]

# select which portions of the water column to sample
drt <- d[d$transect==1 & d$cast %in% dcchoose,]
drt <- rbind(drt, d[d$transect==2,])

# compute total concentration per group
drt <- ddply(drt, ~transect + cast + front + dateTimeB + group3, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# }

##{ Creating some species plots ------------------------------------------
dO=d[d$taxon=="Ocyropsis maculata",]
dO$lconc=log1p(dO$concentration)
dO$pres <- dO$concentration > 0

ggplot(dO) + geom_point(aes(x=sal, y=temp, size=pres, shape=long>119.6, colour=depth> 50), alpha=0.7) + scale_size_manual(values=c(1, 2.5))

dO=d[d$group=="Solmaris",]
dO$lconc=log1p(dO$concentration)
dO$pres <- dO$concentration > 0

ggplot(dO) + geom_point(aes(x=sal, y=temp, size=concentration, color=taxon, alpha=pres)) 
#}

##{ Per group trees -------------------------------------------------------
# SOLMARIS - small
m <- mvpart(concentration ~ temp + sal + flu + oxy, data=dg[dg$group2=="sol_small",], xv=("none"))
m <- mvpart(concentration ~ depth + long, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + sal + flu, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/solsm_tsf_min.txt", sep="\n", append=T)
# --> tree with xv=min with temp, sal, flu only explains 79% of the variation. adding in the front status adds <1% of the variation explained. so describe this tree with the 4 env variables only.

# SOLMARIS - large
m <- mvpart(concentration ~ temp + sal + flu + oxy, data=dg[dg$group2=="sol_large",], xv=("none"))
m <- mvpart(concentration ~ depth + long, data=dg[dg$group2=="sol_large",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + sal + flu + oxy + front, data=dg[dg$group2=="sol_large",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + sal + flu + front, data=dg[dg$group2=="sol_large",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + sal + flu + front, data=dg[dg$group2=="sol_large",], xv=("1se"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/sollg_tsffr_1se.txt", sep="\n", append=T)
# --> tree with xv=1se with temp, sal, flu and front explains 37% of the variation but using xv="min" the tree explains 71.5% of the variation. the major divisions are the same, however. Might be useful to show both. 
# --> tree with xv="min" using temp, flu, oxy and front explains 57% of the variation.

# APPENDICULARIANS
m <- mvpart(concentration ~ depth + long + temp + sal + flu + oxy, data=dg[dg$group2=="appendicularians",], xv=("min"), xval=5, xvmult=100)
# --> NO splits possible with apps. 

# HYDROMEDUSAE
# all as a group
m <- mvpart(concentration ~ temp + sal + flu + oxy, data=dg[dg$group2=="Hydromedusae",], xv=("min"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/hydro_tsfo_min.txt", sep="\n", append=T)
# --> all hydro's as a group with xv=min: best to use sal (it cross correlates with the front). explains 64% of the variation
# --> all hydro's with xv=1se, explains 46% of the variation
# --> also tried running the regression tree w/ all the hydros except the rare species. and it does not change the tree at all.

# try all hydros splitting up into functional group
dh <- d[d$group3 %ni% exclude & d$group=="Hydromedusae",]
dh <- ddply(dh, ~transect + cast + down.up + dateTimeB + group3, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

dh <- dcast(dh, dateTimeB ~ group3, value.var="concentration")
dh <- join(dh, unique(d[,c("dateTimeB", vars)]))

# remove the all zeros
spp <- c("Deep Narco", "Deep Trachy", "Other Hydro", "Shallow Narco", "Shallow Trachy")
dhspp <- dh[,names(dh) %in% spp]
dhenv <- dh[, c(hydroVars, "front")]
keep <- which(rowSums(dhspp) > 0)
dhspp <- dhspp[keep,]
dhenv <- dhenv[keep,]

dh <- cbind(dhspp, dhenv)

m <- mvpart(data.matrix(decostand(dh[,names(dh) %in% spp], method="chi.square")) ~  temp + sal + flu + oxy + front, data=dh, xv=("pick"), xval=5, xvmult=100)
# --> why do i get a lower variance explained using a MRT? I get only 23% variance explained here. 

# SIPHONOPHORES
m <- mvpart(concentration ~ temp + sal + flu + front, data=dg[dg$group2=="Siphonophores",], xv=("min"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/siph_tsffr_min.txt", sep="\n", append=T)
# --> model with temp, sal, flu and front with xv="min" explains 42.5% of the variation, xv="1se" explains 39.3% of the variation
# --> model with temp, sal, and flu with xv="min" explain 43% of the variation

m <- mvpart(concentration ~ depth + temp + sal + flu + oxy + front, data=d[d$taxon=="Physonect",], xv=("min"), xval=5, xvmult=100)
# --> a model with diphyidae only explains 19% of the variation. a model with physonects only explain 19.8% of the variation. a model with sphaeronectes explains 41.9% of the variation. a model with Lilyopsis explains <1 % of the variation. And no splits possible with Prayidae

# test to see how the mrt is for siphonophores
ds <- d[d$group=="Siphonophores",]

ds <- dcast(ds, dateTimeB ~ taxon, value.var="concentration")
ds <- join(ds, unique(d[,c("dateTimeB", vars)]))

# remove the all zeros
spp <- c("Diphyidae", "Lilyopsis", "Physonect", "Prayidae", "Sphaeronectes")
dsspp <- ds[,names(ds) %in% spp]
dsenv <- ds[, c(hydroVars, "front")]
keep <- which(rowSums(dsspp) > 0)
dsspp <- dsspp[keep,]
dsenv <- dsenv[keep,]

ds <- cbind(dsspp, dsenv)

m <- mvpart(data.matrix(decostand(ds[,names(ds) %in% spp], method="chi.square")) ~  temp + sal + flu + oxy + front, data=ds, xv=("min"), xval=5, xvmult=100)
# --> results from the mrt is a model that explains 33% of the variation

# CTENOPHORES
m <- mvpart(concentration ~ sal + oxy + front, data=dg[dg$group2=="Ctenophores",], xv=("min"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/cteno_sofr_min.txt", sep="\n", append=T)
# --> original model build suggested that ctenophores cannot be described with this method because they were not numerous enough. Using the regression trees model, the whole model can describe 15% of the variation for the ctenophores

# try with the different groups of ctenophores

dct <- d[d$group3 %ni% exclude & d$group == "Ctenophores",]
dct <- ddply(dct, ~transect + cast + down.up + dateTimeB + group3, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

dct[dct$group3 == "Thalassocalycidae inconstans","group3"] <- "Thalassocalycidae"

m <- mvpart(concentration ~ temp + sal + flu + oxy + front, data=dct[dct$group3=="Lobata",], xv=("min"), xval=5, xvmult=100)
# --> 20.7% of the variation
m <- mvpart(concentration ~ temp + sal + flu + oxy + front, data=dct[dct$group3=="Cydippida",], xv=("min"), xval=5, xvmult=100)
# --> 6% of the variation
m <- mvpart(concentration ~ temp + sal + flu + oxy + front, data=dct[dct$group3=="Beroida",], xv=("min"), xval=5, xvmult=100)
# --> a complicated model that explains 27% of the variation
m <- mvpart(concentration ~ temp + sal + flu + oxy + front, data=dct[dct$group3=="Velamen",], xv=("min"), xval=5, xvmult=100)
# --> poor model that explains < 5% of the variation
# --> Thalasso also terrible, 2% of variation

# ctenophores together
dct <- dcast(dct, dateTimeB ~ group3, value.var="concentration")
dct <- join(dct, unique(d[,c("dateTimeB", vars)]))

# remove the all zeros
spp <- c("Beroida", "Cydippida", "Lobata", "Thalassocalycidae", "Velamen")
dcspp <- dct[,names(dct) %in% spp]
dcenv <- dct[, c(hydroVars, "front")]
keep <- which(rowSums(dcspp) > 0)
dcspp <- dcspp[keep,]
dcenv <- dcenv[keep,]

dct <- cbind(dcspp, dcenv)

m <- mvpart(data.matrix(decostand(dct[,names(dct) %in% spp], method="chi.square")) ~  temp + sal + flu + oxy + front, data=dct, xv=("min"), xval=5, xvmult=100)
# --> results from the mrt is a model that explains 20% of the variation (actually better than the bulk model)

# DOLIOLIDS
m <- mvpart(concentration ~ front + temp + oxy + flu, data=dg[dg$group2=="Doliolids",], xv=("min"), xval=5, xvmult=100)
out <- c(capture.output(printcp(m)), capture.output(summary(m)))
cat(out, file="plots/regression_trees/dol_toffr_min.txt", sep="\n", append=T)
# --> 50% of the variation explained by temperature, sal and oxy alone.

# }

##{ Community structure ---------------------------------------------------

# let's look at the influence of the different sizes of Solmaris on each other
dSol <- d[d$group=="Solmaris",]
dSol <- dcast(dSol, dateTimeB ~ taxon, value.var="concentration")
dSol <- join(dSol, unique(d[,c("dateTimeB", vars)]))

m <- mvpart(sol_small ~ temp + sal + flu + oxy + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)
m <- mvpart(sol_small ~ depth + long + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)

# in De'ath (2002) he says to standardize the species and then sites to the same mean. And also to standardize the data to give chi-squared distances between sites. this then gives an opportunity to compare with cca
# decostand
# setting environmental vars
dCenv <- dgC[, hydroVars]
spp <- c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")
dCspp <- dgC[,names(dgC) %in% spp]
keep <- which(rowSums(dCspp) > 0)
dCspp <- dCspp[keep,]
dCenv <- dCenv[keep,]
# check again
nrow(dCenv) == nrow(dCspp)

dCbind <- cbind(dCspp, dCenv)

m <- mvpart(data.matrix(decostand(dCbind[,names(dCbind) %in% spp], method="chi.square")) ~  temp + sal + flu + oxy, data=dCbind, xv=("min"), xval=5, xvmult=100)

# creating multivariate species plots with new groups
# cast drt data frame (drt data frame is initialized above - data divided into groups defined by ordination analysis/taxonomic group and including only those portions of the data where all the species data is available)
drt <- dcast(drt, dateTimeB ~ group3, value.var="concentration")
drt <- join(drt, unique(d[,c("dateTimeB", vars)]))

rename <- c("appendicularians" = "Appendicularians", "doliolids" = "Doliolids", "Physonect" = "Physonectae", "sol_large" = "Solmaris_Lg", "sol_small" = "Solmaris_Sm", "Thalassocalycidae inconstans" = "Thalassocalycidae")
drt <- rename(drt, rename)

spp <- c("Appendicularians", "Beroida", "Cydippida", "Deep Narco", "Deep Trachy", "Diphyidae", "Doliolids", "Lobata", "Other Hydro", "Physonectae", "Prayidae", "Shallow Narco", "Shallow Trachy", "Solmaris_Lg", "Solmaris_Sm", "Sphaeronectes", "Thalassocalycidae", "Velamen")

# separate into environmental vars and species vars
drtenv <- drt[, vars]
drtspp <- drt[,names(drt) %in% spp]

keep <- which(rowSums(drtspp) > 0)
drtenv <- drtenv[keep,]
drtspp <- drtspp[keep,]

drtbind <- cbind(drtspp, drtenv)

m <- mvpart(data.matrix(decostand(drtbind[,names(drtbind) %in% spp], method="chi.square")) ~  temp + sal + flu + oxy, data=drtbind, xv=("pick"), xval=5, xvmult=100)
# --> still only explains 4% of the variability in the data.

# }

##{ Logistic Regression  ------------------------------------------------
# use logistic regression on presence / absence data

# code in Presence or Absence in original data frame
d$pa <- 0
d$pa[which(d$concentration != 0)] <- 1


# cast data frame
dC <- dcast (d, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + sal + flu + oxy ~ taxon, value.var="pa")
# make names you can call
names(dC) <- make.names(names(dC))

# try logistic regression on Ocyropsis maculata (Ctenophore)
nrow(dC)

ocmalogit <- glm(Ocyropsis.maculata ~ (temp + sal + flu + oxy)^4 + I(temp^2) + I(sal^2) + I(flu^2) + I(oxy^2) , data=dC, family=binomial("logit"))
summary(ocmalogit) # intercept, sal, and oxy are significant explanatory variables
# make a model only out of significant explanatory variables

m <- update(ocmalogit, . ~ . - temp:sal:flu:oxy)
summary(m)
m <- update(m, . ~ . - temp:sal:flu)
summary(m)
m <- update(m, . ~ . - sal:flu:oxy)
summary(m)

m <- step(ocmalogit, direction="backward")
summary(m)
m <- update(m, . ~ . - sal:flu:oxy)
summary(m)
m <- update(m, . ~ . - temp:sal:oxy)
summary(m)
m <- update(m, . ~ . - sal:flu)
summary(m)
m <- update(m, . ~ . - I(flu^2))
summary(m)
m <- update(m, . ~ . - temp:sal:oxy)
summary(m)
m <- update(m, . ~ . - temp:sal:oxy)
summary(m)
m <- update(m, . ~ . - temp:sal:oxy)
summary(m)




ocmalogit <- glm(Ocyropsis.maculata ~ sal + oxy, data=dC, family=binomial("logit"))
summary(ocmalogit)
plot(predict(ocmalogit, type="response") ~ sal + oxy, data=dC)

# TODO: try others

# }

