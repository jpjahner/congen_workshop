###########################################
## two_pop_sims.R
###########################################

## Contact: Josh Jahner, jpjahner@gmail.com
## 19 v 25


## Plot beta distributions for the two populations

xvals <- seq(0, 1, by=0.01)
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=1.5, las=1)
points(xvals, dbeta(xvals, shape1=0.5, shape2=0.05), type="l", col="blue", lwd=2)
points(xvals, dbeta(xvals, shape1=0.05, shape2=0.5), type="l", col="red", lwd=2)



## Function for simulating genetic data

single_pop_sim <- function(nloci=NA, ninds=NA, alpha=NA, beta=NA){
	## create a matrix to hold the genetic data
	genotypes <- matrix(NA, ninds, nloci)
	for (i in 1:nloci){
		## sample the beta distribution once to determine the allele frequency for this locus
		afreq <- rbeta(1, 0.48, 0.21)
		for (j in 1:ninds){
			## use two random draws from the binomial distribution (p=afrea) to determine the individual's genotype
				## potential genotypes: 0, 1, 2
			genotypes[j,i] <- sum(rbinom(2, 1, afreq))
		}
	}
	return(genotypes)
}



## Generate genotype matrix for both populations

ninds_per_pop <- 25
pop1_genos <- single_pop_sim(nloci=1000, ninds=ninds_per_pop, alpha=0.5, beta=0.05)
pop2_genos <- single_pop_sim(nloci=1000, ninds=ninds_per_pop, alpha=0.05, beta=0.5)



## principal component analysis (PCA) to look at population structure

all_genos <- rbind(pop1_genos, pop2_genos)
all_pca <- prcomp(all_genos, center=TRUE, scale=FALSE)
	summary(all_pca)



## plotting PCA results
	 ## first pop plotted in blue
	 ## second pop plotted in red

par(mar=c(5,5,1,1))
plot(all_pca$x[,1], all_pca$x[,2], type="n", xlab="PC 1", ylab="PC 2", cex.lab=1.5, las=1)
points(all_pca$x[1:ninds_per_pop,1], all_pca$x[1:ninds_per_pop,2], pch=21, bg="blue", cex=1.5)
points(all_pca$x[(ninds_per_pop+1):(ninds_per_pop*2),1], all_pca$x[(ninds_per_pop+1):(ninds_per_pop*2),2], pch=21, bg="red", cex=1.5)



