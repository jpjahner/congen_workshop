###########################################
## 01_beta_simulate.R
###########################################

## Contact: Josh Jahner, jpjahner@gmail.com
## 19 v 25


## Using the metbrewer package to select plot colors
	## All color palettes based on art in the the Metropolitan Museum of Art
	## https://github.com/BlakeRMills/MetBrewer
install.packages("MetBrewer")
library(MetBrewer)
plot_colors <- met.brewer("Johnson", 5) ## choosing 5 colors from the Isfahan2 palette



## Exploring how different values of alpha (shape1) and beta (shape2) effect the shape of a beta distribution

xvals <- seq(0, 1, by=0.01)

par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
## Uniform distribution
	points(xvals, dbeta(xvals, shape1=1, shape2=1), type="l", col=plot_colors[3], lwd=6) 		## yellow
## Symmetrical U-shaped distributions
	points(xvals, dbeta(xvals, shape1=0.5, shape2=0.5), type="l", col=plot_colors[2], lwd=6)	## light red
	points(xvals, dbeta(xvals, shape1=0.1, shape2=0.1), type="l", col=plot_colors[1], lwd=6)	## dark red
## Left and right skewed U-shaped distribution
	points(xvals, dbeta(xvals, shape1=0.8, shape2=0.2), type="l", col=plot_colors[4], lwd=6)	## light blue
	points(xvals, dbeta(xvals, shape1=0.2, shape2=0.8), type="l", col=plot_colors[5], lwd=6)	## dark blue



## Using the beta distribution to model the site frequency spectrum from empirical data
	## Rosenthal et al. (2024) Influence of dams on sauger population structure and hybridization with introduced walleye
	## Try out a few different values for shape1 and shape2 below to see if you can match the plot on the figure
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
points(xvals, dbeta(xvals, shape1=NA, shape2=NA), type="l", col="black", lwd=6)
	## change the NAs above


## Function for simulating genetic data from the beta distribution

single_pop_sim <- function(nloci=NA, ninds=NA, alpha=NA, beta=NA){
	## create a matrix to hold the genetic data
	genotypes <- matrix(NA, ninds, nloci)
	for (i in 1:nloci){
		## sample the beta distribution once to determine the allele frequency for this locus
		afreq <- rbeta(1, alpha, beta)
		for (j in 1:ninds){
			## use two random draws from the binomial distribution (p=afreq) to determine the individual's genotype
				## potential genotypes: 0, 1, 2
			genotypes[j,i] <- sum(rbinom(2, 1, afreq))
		}
	}
	return(genotypes)
}



## Creating an example genetic data matrix

set.seed(17)
geno_mat <- single_pop_sim(nloci=100, ninds=20, alpha=0.48, beta=0.21)
geno_mat


