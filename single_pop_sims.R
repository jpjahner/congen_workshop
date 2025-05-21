###########################################
## single_pop_sims.R
###########################################

## Contact: Josh Jahner, jpjahner@gmail.com
## 19 v 25



## Using the beta distribution to model the site frequency spectrum

xvals <- seq(0, 1, by=0.01)

## Uniform distribution
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
points(xvals, dbeta(xvals, shape1=1, shape2=1), type="l", col="black", lwd=4)

## Symmetrical U-shaped distributions
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
points(xvals, dbeta(xvals, shape1=0.5, shape2=0.5), type="l", col="black", lwd=4)
points(xvals, dbeta(xvals, shape1=0.1, shape2=0.1), type="l", col="red", lwd=4)

## Left and right skewed U-shaped distribution
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
points(xvals, dbeta(xvals, shape1=0.8, shape2=0.2), type="l", col="black", lwd=4)
points(xvals, dbeta(xvals, shape1=0.2, shape2=0.8), type="l", col="red", lwd=4)

## beta distribution from empirical data
	## Rosenthal et al. (2024) Influence of dams on sauger population structure and hybridization with introduced walleye
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=2, cex.axis=1.5, las=1)
points(xvals, dbeta(xvals, shape1=NA, shape2=NA), type="l", col="black", lwd=4)



## Function for simulating genetic data from the beta distribution

single_pop_sim <- function(nloci=NA, ninds=NA, alpha=NA, beta=NA){
	## create a matrix to hold the genetic data
	genotypes <- matrix(NA, ninds, nloci)
	for (i in 1:nloci){
		## sample the beta distribution once to determine the allele frequency for this locus
		afreq <- rbeta(1, 0.48, 0.21)
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


