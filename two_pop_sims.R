###########################################
## two_pop_sims.R
###########################################

## Contact: Josh Jahner, jpjahner@gmail.com
## 19 v 25

## The goal of this section is to use your population genetic simulations
## to carry out two population genetic analyses. First, you will conduct
## a principal components analysis (PCA) to see if there is genetic
## structure between two simulated populations. In the PCA plot, if 
## populations are strongly divided on the left and right hand sides, that
## is evidence for strong genetic structure. Try out a few different
## combinations of beta distribution shapes and numbers of loci in your
## simulations. Which combinations do and do not have strong genetic
## structure for the PCA?



## Plot beta distributions for two populations

xvals <- seq(0, 1, by=0.01)
par(mar=c(5,5,1,1))
plot(0, type="n", xlim=c(0,1), ylim=c(0,5), xlab="Allele frequency", ylab="Beta density", cex.lab=1.5, las=1)
points(xvals, dbeta(xvals, shape1=0.5, shape2=0.05), type="l", col="blue", lwd=2)
points(xvals, dbeta(xvals, shape1=0.05, shape2=0.5), type="l", col="red", lwd=2)



## Function for simulating genetic data from a beta distribution

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

ninds_per_pop <- 20 ## NOTE: please keep this number even, or the parentage sims below will break
nloci_per_pop <- 500
pop1_genos <- single_pop_sim(nloci=nloci_per_pop, ninds=ninds_per_pop, alpha=0.5, beta=0.05)
pop2_genos <- single_pop_sim(nloci=nloci_per_pop, ninds=ninds_per_pop, alpha=0.05, beta=0.5)



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



## Next, you will take one population and cross the individuals to create
## a second offspring generation to assess your ability to recover parent
## offspring relationships. Again, try out a few combinations of beta
## distribution shapes and numbers of loci. For which combinations are 
## you able to recover accurate parent offspring relationships?


## create parent population genotypes

pedigree_ninds_per_pop <- 20 ## NOTE: please keep this number even, or the parentage sims will break
pedigree_nloci_per_pop <- 500
parent_genos <- single_pop_sim(nloci=pedigree_nloci_per_pop, ninds=pedigree_ninds_per_pop, alpha=0.5, beta=0.05)



## cross individuals from population 1 from above to create a second generation
	## when pop size is 20:
	## ind1 crossed with ind11
	## ind2 crossed with ind12
	## ind3 crossed with ind13
	## ind4 crossed with ind14
	## ind5 crossed with ind15
	## ind6 crossed with ind16
	## ind7 crossed with ind17
	## ind8 crossed with ind18
	## ind9 crossed with ind19
	## ind10 crossed with ind20

offspring_genos <- matrix(NA, (pedigree_ninds_per_pop/2), pedigree_nloci_per_pop)
for (i in 1:(pedigree_ninds_per_pop/2)){
	for (j in 1:pedigree_nloci_per_pop){
		if		(parent_genos[i,j]==0 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==0) { offspring_genos[i,j] <- 0 }							## possible genotypes = 0
		else if	(parent_genos[i,j]==0 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==1) { offspring_genos[i,j] <- rbinom(1, 1, 0.5) }			## possible genotypes = 0,1
		else if	(parent_genos[i,j]==0 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==2) { offspring_genos[i,j] <- 1 }							## possible genotypes = 1
		else if	(parent_genos[i,j]==1 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==0) { offspring_genos[i,j] <- rbinom(1, 1, 0.5) }			## possible genotypes = 0,1
		else if	(parent_genos[i,j]==1 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==1) { offspring_genos[i,j] <- sum(rbinom(2, 1, 0.5)) }	## possible genotypes = 0,1,2
		else if	(parent_genos[i,j]==1 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==2) { offspring_genos[i,j] <- rbinom(1, 1, 0.5) + 1 }		## possible genotypes = 1,2
		else if	(parent_genos[i,j]==2 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==0) { offspring_genos[i,j] <- 1 }							## possible genotypes = 1
		else if	(parent_genos[i,j]==2 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==1) { offspring_genos[i,j] <- rbinom(1, 1, 0.5) + 1 }		## possible genotypes = 1,2
		else if	(parent_genos[i,j]==2 && parent_genos[i+(pedigree_ninds_per_pop/2),j]==2) { offspring_genos[i,j] <- 2 }							## possible genotypes = 2
	}
}


## parent offspring relationships with the sequoia package
#install.packages("sequoia")
library(sequoia)

## sequoia requires row names with individual ids, so let's make those here
parent_names <- vector()
for (i in 1:pedigree_ninds_per_pop) { parent_names[i] <- paste0("parent_", i) }
offspring_names <- vector()
for (i in 1:(pedigree_ninds_per_pop/2)) { offspring_names[i] <- paste0("offspring_", i) }

## merge the pop1 and offspring genotype matrices, add the individual ids, and run sequoia
pedigree_genos <- rbind(parent_genos, offspring_genos)
rownames(pedigree_genos) <- c(parent_names, offspring_names)
pedigree_out <- GetMaybeRel(pedigree_genos)
	pedigree_out$MaybeTrio



