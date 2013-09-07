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

# }

## { Define new groups for analysis ----------------------------------
Deep_Hydro <- c("h11_Haliscera", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura", "h3_Cunina")
Liriope <- c("h5_Liriope", "h5b")
Shallow_Narco <- c("h7_Pegantha", "h6_Solmundella")
Other_Hydro <- c("h1", "h15", "vsh")
Cydippida <- c("Haeckelia beehlri", "Hormiphora californiensis", "Mertensid")
Lobata_Thalasso <- c("Ocyropsis maculata", "Juvenile Lobata", "Larval Lobata", "Thalassocalycidae inconstans")
Prayidae <- c("Lilyopsis", "Prayidae")
Solmaris <- c("sol_large", "sol_small")

# assign these groups into a different column
d$group2 <- d$taxon
d$group2[d$taxon %in% Deep_Hydro] <- "Deep Hydro"
d$group2[d$taxon %in% Liriope] <- "Liriope"
d$group2[d$taxon %in% Shallow_Narco] <- "Shallow Narco"
d$group2[d$taxon %in% Solmaris] <- "Solmaris"
d$group2[d$taxon %in% Other_Hydro] <- "Other Hydro"
d$group2[d$taxon %in% Cydippida] <- "Cydippida"
d$group2[d$taxon %in% Lobata_Thalasso] <- "Lobata_Thalasso"
d$group2[d$taxon %in% Prayidae] <- "Prayidae"

# exclude rare taxa
`%ni%` <- Negate(`%in%`) 
exclude <- c("Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown", "Annatiara", "h10_Pegantha", "h13", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia", "Bolinopsis", "h1")
d <- d[d$group2 %ni% exclude,]

# compute total concentration per group
dg <- ddply(d, ~transect + cast + front + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# }

## { Construct a correlation heatmap for all taxa separately---------------------------

# convert to wide format
dW <- dcast(dg, dateTimeB~group2, value.var="concentration")

# rename columns
rename <- c("appendicularians" = "Apps", "doliolids" = "Doliolids", "Physonect" = "Physonectae")
dW <- rename(dW, rename)

# set the order
# use the results of the correspondence analysis
# --> CA showed that you can group the taxa into six different groups, 1) velamen, doliolids, liriope, appendicularians, shallow narco, other hydro and solmaris, 2) lobata_thalasso, 3) beroida, 4) cydippida, diphyidae, physonect and sphaeronectes, 5) deep hydro, 6) prayidae
levels <- c("Velamen", "Doliolids", "Liriope", "Apps", "Shallow Narco", "Other Hydro", "Solmaris", "Lobata_Thalasso", "Beroida", "Cydippida", "Diphyidae", "Physonectae", "Sphaeronectes", "Deep Hydro", "Prayidae")

# reorder the columns
dW <- dW[,levels]

# calculate the spearman's correlation coefficient
dcorr <- cor(dW, use="complete.obs", method="spearman")
dc <- melt(dcorr)

# remove the identity values
dc[which(dc$value==1),"value"] <- NA

# save only the lower triangle's rho values (for text later)
dcorr[upper.tri(dcorr,diag=T)] <- NA
dctext <- melt(dcorr)

# plot the heatmap
cortheme <- theme(axis.text.x=element_text(angle=90, vjust=0.5, size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title=element_text(size=26))

ggplot(mapping=aes(x=factor(Var1, levels=levels), y=factor(Var2, levels=levels))) + geom_tile(aes(fill=value), data=dc) + geom_text(aes(label=round(value,2)), data=dctext) + labs(x="", y="", title="", fill="") + scale_fill_gradient2(limits=c(-1,1), low="red", high="blue", na.value="grey90") + theme_bw() + cortheme

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

res1 <- cor.mtest(cor(dW, use="complete.obs", method="spearman"), 0.9)
corr <- cor(dW, use="complete.obs", method="spearman")
corrplot(corr, method="square")

corrplot.mixed(corr, upper="color", lower="number")

corrplot(cor(dW, use="complete.obs", method="spearman"), p.mat = res1[[1]], sig.level = 0.1)
corrplot(cor(dW, use="complete.obs", method="spearman"), p.mat = res1[[1]], insig="blank")

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

