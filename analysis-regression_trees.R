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

ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group2, scale="free")
ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=log(concentration))) + facet_wrap(~group2, scale="free")

# casing data into wide format with group 2
dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + fluoro + oxygen ~ group2, value.var="concentration")

# }

##{ Creating some species plots
dO=d[d$taxon=="Ocyropsis maculata",]
dO$lconc=log1p(dO$concentration)
dO$pres <- dO$concentration > 0

ggplot(dO) + geom_point(aes(x=salinity, y=temp, size=pres, shape=long>119.6, colour=depth> 50), alpha=0.7) + scale_size_manual(values=c(1, 2.5))

#}

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

m <- mvpart(data.matrix(decostand(dCbind[,names(dCbind) %in% spp], method="chi.square")) ~  temp + salinity + fluoro + oxygen, data=dCbind, xv=("min"), xval=5, xvmult=100)


# # multivariate regression tree with all taxa
# m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("appendicularians", "Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# # --> basically, when you look at the multivariate regression tree for everything, it is 27% explained by temperature (which is probably just a proxy for depth) and then nothing else
# 
# # multivariate regression tree without the apps
# m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# # --> much better. percentage of variance explained is 57%, which is pretty good. but most of the changes driven by temperature and salinity. bar plots appear to represent absolute abundances rather than proportion of abundances
# 
# # mvpart with scaled species abundances
# m <- mvpart(scale(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")])) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# 
# # multivariate regression tree with just the jellies and the ctenophores
# m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# 
# # multivariate regression tree with just the hydrozoans
# m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100)
# # --> the latter three are about the same
# # think about the interpretation for these plots ...

# maybe perform a square root transformation on the species data for the multivariate regression tree



# # multivariate regression tree with Principal components analysis of the group means
# m <- mvpart(data.matrix(dgC[,names(dgC) %in% c("Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]) ~ temp + salinity + fluoro + oxygen, data=dgC, xv=("min"), xval=5, xvmult=100, pca=T)
# # difficult to interpret

# }

##{ Logistic Regression  ------------------------------------------------
# use logistic regression on presence / absence data

# code in Presence or Absence in original data frame
d$pa <- 0
d$pa[which(d$concentration != 0)] <- 1


# cast data frame
dC <- dcast (d, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + fluoro + oxygen ~ taxon, value.var="pa")
# make names you can call
names(dC) <- make.names(names(dC))

# try logistic regression on Ocyropsis maculata (Ctenophore)
nrow(dC)

ocmalogit <- glm(Ocyropsis.maculata ~ (temp + salinity + fluoro + oxygen)^4 + I(temp^2) + I(salinity^2) + I(fluoro^2) + I(oxygen^2) , data=dC, family=binomial("logit"))
summary(ocmalogit) # intercept, salinity, and oxygen are significant explanatory variables
# make a model only out of significant explanatory variables

m <- update(ocmalogit, . ~ . - temp:salinity:fluoro:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:fluoro)
summary(m)
m <- update(m, . ~ . - salinity:fluoro:oxygen)
summary(m)

m <- step(ocmalogit, direction="backward")
summary(m)
m <- update(m, . ~ . - salinity:fluoro:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - salinity:fluoro)
summary(m)
m <- update(m, . ~ . - I(fluoro^2))
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)




ocmalogit <- glm(Ocyropsis.maculata ~ salinity + oxygen, data=dC, family=binomial("logit"))
summary(ocmalogit)
plot(predict(ocmalogit, type="response") ~ salinity + oxygen, data=dC)

# TODO: try others

# }

