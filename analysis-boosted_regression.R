# boosted regression trees

library(ggplot2)
library(reshape2)
library(plyr)
library(gbm)

d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

d$front <- factor(d$front, levels=c("west", "front", "east"))

hist(log(d[d$taxon=="h5_Liriope","concentration"]))
hist(log(d[d$taxon=="sol_large","concentration"]))

?gbm
# n.trees - want something in the 1000's
# the more trees you have the more variance you can explain
# cv.folds - cross validation; 
# interaction depth - how big each tree is. the higher this number is, the more complex interactions you can fit. the higher the interaction depth and the greater the learning rate the steeper the drop off of the CV curve. BUT you want higher the number of trees, so you would want lower interaction rate and lower shrinkage. 
# do you care about autocorrelation/independence of samples? no because it is dealt w/ the CV
# bag.fraction - take the default (50%)
# train.fraction or nTrain - this is when you have a training set


ds <- d[d$taxon=="sol_small",]

sol_gbm <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=ds,
               distribution="gaussian",
               n.trees=8000,
               shrinkage=0.075,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=10)

best.iter <- gbm.perf(sol_gbm, method="OOB")
print(best.iter)

best.iter <- gbm.perf(sol_gbm, method="test") # need finite ylim values?
print(best.iter)

best.iter <- gbm.perf(sol_gbm, method="cv")
print(best.iter)

summary(sol_gbm, n.trees=1)
summary(sol_gbm, n.trees=best.iter)

print(pretty.gbm.tree(sol_gbm,1))
print(pretty.gbm.tree(sol_gbm, sol_gbm$n.trees))

par(mfrow=c(3,2))
plot(sol_gbm,1,best.iter,cont=20)
plot(sol_gbm,2,best.iter,cont=20)
plot(sol_gbm,3,best.iter,cont=20)
plot(sol_gbm,4,best.iter,cont=20)
plot(sol_gbm,5,best.iter,cont=20)
plot(sol_gbm,6,best.iter,cont=20)
par(mfrow=c(1,1))