#
#       Analysis of community structure via correlation
#
#       Jessica Luo, 2013
#
#--------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read data -------------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

# choose a portion of the data to analyze
# when the dataset is complete this will not be necessary
dcchoose <- unique(d[d$group=="Hydromedusae" & d$transect==1,"cast"])

# there is still one downcast in transect 1 that has not been analyzed
dc <- d[d$transect==1 & d$cast %in% dcchoose,]

# all of transect 2 is complete
dc <- rbind(dc, d[d$transect==2,])

# transect 3 is not yet complete and thus is not included in this analysis

# }

## { Define new groups for analysis ----------------------------------
Deep_Trachy <- c("h11_Haliscera", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura")
Shallow_Trachy <- c("h5_Liriope", "h5b")
Shallow_Narco <- c("h7_Pegantha", "h6_Solmundella")
Deep_Narco <- c("h3_Cunina")
Other_Hydro <- c("h1", "h15", "vsh")
Cydippida <- c("Haeckelia beehlri", "Hormiphora californiensis", "Mertensid")
Lobata <- c("Bolinopsis", "Ocyropsis maculata", "Juvenile Lobata", "Larval Lobata")
Prayidae <- c("Lilyopsis", "Prayidae")

# assign these groups into a different column
d$group2 <- d$taxon
d$group2[d$taxon %in% Deep_Trachy] <- "Deep Trachy"
d$group2[d$taxon %in% Shallow_Trachy] <- "Shallow Trachy"
d$group2[d$taxon %in% Shallow_Narco] <- "Shallow Narco"
d$group2[d$taxon %in% Deep_Narco] <- "Deep Narco"
d$group2[d$taxon %in% Other_Hydro] <- "Other Hydro"
d$group2[d$taxon %in% Cydippida] <- "Cydippida"
d$group2[d$taxon %in% Lobata] <- "Lobata"
d$group2[d$taxon %in% Prayidae] <- "Prayidae"

# exclude rare taxa
`%ni%` <- Negate(`%in%`) 
exclude <- c("Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown", "Annatiara", "h10_Pegantha", "h13", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia")
d <- d[d$group2 %ni% exclude,]

# compute total concentration per group
dc <- ddply(d, ~transect + cast + front + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# }

## { Construct a correlation heatmap for all groups together---------------------------

# convert to wide format
dcW <- dcast(dc, dateTimeB~group2, value.var="concentration")

# rename columns
rename <- c("appendicularians" = "Appendicularians", "doliolids" = "Doliolids", "Physonect" = "Physonectae", "sol_large" = "Solmaris_Lg", "sol_small" = "Solmaris_Sm", "Thalassocalycidae inconstans" = "Thalassocalyce")
dcW <- rename(dcW, rename)

# set the order
# use the results of the correspondence analysis
# --> CA showed that you can group the taxa into four different groups, 1) deep narco and deep trachy, 2) appendicularians, large solmaris, shallow trachy (Liriope), Velamen, doliolids, other Hydros, Shallow Narco, and small solmaris. 3) Lobata and Thalassocalycidae, and 4) Beroida, Cydippida, Diphyidae, Prayidae, Physonect and Sphaeronectes.
levels <- c("Deep Narco", "Deep Trachy", "Appendicularians", "Solmaris_Lg", "Doliolids", "Other Hydro", "Shallow Narco", "Solmaris_Sm", "Shallow Trachy", "Velamen", "Lobata", "Thalassocalyce", "Beroida", "Cydippida", "Diphyidae", "Prayidae", "Physonectae", "Sphaeronectes")

# reorder the columns
dcW <- dcW[,levels]

# calculate the spearman's correlation coefficient and melt into a dataframe
dcWm <- melt(cor(dcW, use="complete.obs", method="spearman"))

# remove the identity values
dcWm[which(dcWm$value==1),"value"] <- NA

# plot the heatmap
cortheme <- theme(axis.text.x=element_text(angle=45, vjust=0.5, size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title=element_text(size=26))

ggplot(data=dcWm) + geom_tile(aes(x=factor(X1, levels=levels), y=factor(X2, levels=levels), fill=value >0, alpha=abs(value))) + labs(x="", y="", title="", fill="Positive / Negative", alpha="Strength of Correlation") + scale_fill_discrete(labels=c("Negative", "Positive")) + theme_bw() + cortheme

# adding in significant values
# try using corrplot package
library(corrplot)
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(cor(dcW, use="complete.obs", method="spearman"), 0.9)
corr <- cor(dcW, use="complete.obs", method="spearman")
corr[which(corr==1)] <- 0
corrplot(corr, method="square")

corrplot(cor(dcW, use="complete.obs", method="spearman"), p.mat = res1[[1]], sig.level = 0.1)
corrplot(cor(dcW, use="complete.obs", method="spearman"), p.mat = res1[[1]], insig="blank")

# }


## { Correlation heatmap for groups separately -------------------------------------

# define the groups separately
group1 <- c("Deep Narco", "Deep Trachy")
group2 <- c("Appendicularians", "Solmaris_Lg", "Doliolids", "Other Hydro", "Shallow Narco", "Solmaris_Sm", "Shallow Trachy", "Velamen")
group3 <- c("Lobata", "Thalassocalyce")
group4 <- c("Beroida", "Cydippida", "Diphyidae", "Prayidae", "Physonectae", "Sphaeronectes")

# recast dc dataframe into wide format
dcW <- dcast(dc, dateTimeB~group2, value.var="concentration")

# subset into 4 groups
dc1 <- dcW[, group1]
dc2 <- dcW[, group2]
dc3 <- dcW[, group3]
dc4 <- dcW[, group4]

# calculate the spearman's correlation coefficient and melt into a dataframe
dc1m <- melt(cor(dc1, use="complete.obs", method="spearman"))
dc2m <- melt(cor(dc2, use="complete.obs", method="spearman"))
dc3m <- melt(cor(dc3, use="complete.obs", method="spearman"))
dc4m <- melt(cor(dc4, use="complete.obs", method="spearman"))


# plot the heatmap
ggplot(data=dc1m) + geom_tile(aes(x=factor(X1, levels=group1), y=factor(X2, levels=group1), fill=value >0, alpha=abs(value))) + labs(x="", y="", title="", fill="Positive / Negative", alpha="Strength of Correlation") + scale_fill_discrete(labels=c("Negative", "Positive")) + theme_bw() + cortheme
ggplot(data=dc2m) + geom_tile(aes(x=factor(X1, levels=group2), y=factor(X2, levels=group2), fill=value >0, alpha=abs(value))) + labs(x="", y="", title="", fill="Positive / Negative", alpha="Strength of Correlation") + scale_fill_discrete(labels=c("Negative", "Positive")) + theme_bw() + cortheme
ggplot(data=dc3m) + geom_tile(aes(x=factor(X1, levels=group3), y=factor(X2, levels=group3), fill=value >0, alpha=abs(value))) + labs(x="", y="", title="", fill="Positive / Negative", alpha="Strength of Correlation") + scale_fill_discrete(labels=c("Negative", "Positive")) + theme_bw() + cortheme
ggplot(data=dc4m) + geom_tile(aes(x=factor(X1, levels=group4), y=factor(X2, levels=group4), fill=value >0, alpha=abs(value))) + labs(x="", y="", title="", fill="Positive / Negative", alpha="Strength of Correlation") + scale_fill_discrete(labels=c("Negative", "Positive")) + theme_bw() + cortheme

# --> splitting up the groups do not necessarily help the visualization. The correlation coefficients do not change. There are low (and insignificant) correlation coefficients for many of the ctenophores and siphonophores because they are not abundant in the water column. I could bin them in larger depth ranges but I'm not sure that it will be ecologically meaningful to do so. 

# }

