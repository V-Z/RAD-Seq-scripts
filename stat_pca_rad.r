## Libraries
# Install packages
# install.packages(pkgs=c("stringi", "ade4", "adegenet", "adegraphics", "vcfR", "ape", "pegas", "StAMPP", "memuse", "pinfsc50", "tidyr"), lib="rpkgs", repos="https://mirrors.nic.cz/R/", dependencies=TRUE)
# Libraries
library(package=stringi, lib.loc="rpkgs")
library(package=ade4, lib.loc="rpkgs")
library(package=adegenet, lib.loc="rpkgs")
library(package=adegraphics, lib.loc="rpkgs")
library(package=vcfR, lib.loc="rpkgs")
library(package=ape, lib.loc="rpkgs")
library(package=pegas, lib.loc="rpkgs")
library(package=StAMPP, lib.loc="rpkgs")
library(package=memuse, lib.loc="rpkgs")

## Input file
vcf.file <- Sys.getenv("VCFR")

## Functions
# Conversion from vcfR object to genlight in tetraploids
vcfR2genlight.tetra <- function(x, n.cores=2) {
	bi <- is.biallelic(x)
	if (sum(!bi) > 0) {
		msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
		msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
		msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
		warning(msg)
		x <- x[bi, ]
		}
	x <- addID(x)
	CHROM <- x@fix[, "CHROM"]
	POS <- x@fix[, "POS"]
	ID <- x@fix[, "ID"]
	x <- extract.gt(x)
	x[x == "0|0"] <- 0
	x[x == "0|1"] <- 1
	x[x == "1|0"] <- 1
	x[x == "1|1"] <- 2
	x[x == "0/0"] <- 0
	x[x == "0/1"] <- 1
	x[x == "1/0"] <- 1
	x[x == "1/1"] <- 2
	x[x == "1/1/1/1"] <- 4
	x[x == "0/1/1/1"] <- 3
	x[x == "0/0/1/1"] <- 2
	x[x == "0/0/0/1"] <- 1
	x[x == "0/0/0/0"] <- 0
	if (requireNamespace("adegenet")) {
		x <- new("genlight", t(x), n.cores=n.cores)
		}
		else {
			warning("adegenet not installed")
			}
	adegenet::chromosome(x) <- CHROM
	adegenet::position(x) <- POS
	adegenet::locNames(x) <- ID
	return(x)
	}

# Patch for faster PCA calculation on genlight objects, see https://github.com/thibautjombart/adegenet/pull/150
glPcaFast <- function(x, center=TRUE, scale=FALSE, nf=NULL, loadings=TRUE, alleleAsUnit=FALSE, returnDotProd=FALSE) {
	if(!inherits(x, "genlight")) stop("x is not a genlight object")
	# Keep the original mean / var code, as it's used further down and has some NA checks
	if(center) {
		vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
		if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
		}
	if(scale) {
		vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
		if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
		}
	# Convert to full data, try to keep the NA handling as similar to the original as possible - dividing by ploidy keeps the NAs
	mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
	# Handle NAs
	NAidx <- which(is.na(mx), arr.ind=T)
	if (center) {
		mx[NAidx] <- vecMeans[NAidx[,2]]
		} else {
			mx[NAidx] <- 0
			}
	# Center and scale
	mx <- scale(mx, center=if (center) vecMeans else F, scale=if (scale) vecVar else F)
	# All dot products at once using underlying BLAS to support thousands of samples, this could be replaced by 'Truncated SVD', but it would require more changes in the code around
	allProd <- tcrossprod(mx) / nInd(x) # Assume uniform weights
	# Perform the analysis
	# Eigenanalysis
	eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
	rank <- sum(eigRes$values > 1e-12)
	eigRes$values <- eigRes$values[1:rank]
	eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
	# Scan nb of axes retained
	if(is.null(nf)) {
		barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
		cat("Select the number of axes: ")
		nf <- as.integer(readLines(n=1))
		}
	# Rescale PCs
	res <- list()
	res$eig <- eigRes$values
	nf <- min(nf, sum(res$eig>1e-10))
	# res$matprod <- allProd # For debugging
	# use: li=XQU=V\Lambda^(1/2)
	eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
	res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
	# Get loadings
	# Need to decompose X^TDV into a sum of n matrices of dim p*r but only two such matrices are represented at a time
	if(loadings) {
		if(scale) {
			vecSd <- sqrt(vecVar)
			}
		res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # Create empty matrix
		# use: c1=X^TDV
		# and X^TV=A_1 + ... + A_n
		# with A_k=X_[k-]^T v[k-]
		myPloidy <- ploidy(x)
		for(k in 1:nInd(x)){
			temp <- as.integer(x@gen[[k]]) / myPloidy[k]
			if(center) {
				temp[is.na(temp)] <- vecMeans[is.na(temp)]
				temp <- temp - vecMeans
				} else {
					temp[is.na(temp)] <- 0
					}
			if(scale){
				temp <- temp/vecSd
				}
			res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
			}
		res$loadings <- res$loadings / nInd(x) # Don't forget the /n of X_tDV
		res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
		}
	# Format output
	colnames(res$scores) <- paste("PC", 1:nf, sep="")
	if(!is.null(indNames(x))) {
		rownames(res$scores) <- indNames(x)
		} else {
			rownames(res$scores) <- 1:nInd(x)
			}
	if(!is.null(res$loadings)) {
		colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
		if(!is.null(locNames(x)) & !is.null(alleles(x))) {
			rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
			} else {
				rownames(res$loadings) <- 1:nLoc(x)
				}
		}
	if(returnDotProd) {
		res$dotProd <- allProd
		rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
		}
	res$call <- match.call()
	class(res) <- "glPca"
	return(res)
	}

## Import SNP data from VCF
vcf.data <- read.vcfR(vcf.file)
# Checks
head(vcf.data)
vcf.data@fix[1:10,1:5]

## Various checks
# In VCF R
# DP - quickcheck
vcf.dp <- extract.gt(vcf.data, element='DP', as.numeric=TRUE)

# Boxplot
pdf("dp_per_sample.pdf", width=35, height=7)
boxplot(vcf.dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7)
abline(h=mean(x=as.numeric(vcf.dp), na.rm=TRUE), col="red")
dev.off()

# N missing SNPs per sample
write.table(x=summary(t(as.matrix(rad.genlight)))[7,], file="missing_per_sample.txt", sep="\t")

## Convert to genlight
rad.genlight <- vcfR2genlight.tetra(x=vcf.data, n.cores=2)
pop(rad.genlight) <- substr(indNames(rad.genlight), start=1, stop=4) # Add pop names: here beginning of individuals names
# See names
indNames(rad.genlight)
pop(rad.genlight)

## PCA
rad.pca <- glPcaFast(rad.genlight, nf=5)

# Save nice figs
pdf("PCA_all_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(rad.pca$scores, pop(rad.genlight),  xax=1, yax=2, col=transp("blue", 0.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=TRUE, pgrid.draw=FALSE, plot=FALSE)
g2 <- s.label(rad.pca$scores, xax=1, yax=2, ppoints.col="red", plabels=list(box=list(draw=FALSE), optim=TRUE), paxes.draw=TRUE, pgrid.draw=FALSE, plabels.cex=1, plot=FALSE)
ADEgS(c(g1, g2), layout=c(1, 2))
dev.off()

pdf("PCA_all_SNPs_ax13.pdf", width=14, height=7)
g1 <- s.class(rad.pca$scores, pop(rad.genlight),  xax=1, yax=3, col=transp("blue", 0.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=TRUE, pgrid.draw=FALSE, plot=FALSE)
g2 <- s.label(rad.pca$scores, xax=1, yax=3, ppoints.col="red", plabels=list(box=list(draw=FALSE), optim=TRUE), paxes.draw=TRUE, pgrid.draw=FALSE, plabels.cex=1, plot=FALSE)
ADEgS(c(g1, g2), layout=c(1, 2))
dev.off()

# Ploidy - differentiated plots
pdf ("PCA_all_ploidycol_SNPs_ax12.pdf", width=14, height=7)
g1 <- s.class(rad.pca$scores, as.factor(as.vector(ploidy(rad.genlight))), xax=1, yax=2, col=transp(c("#FF0000", "#0000FF")), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=TRUE, pgrid.draw=FALSE, plab.cex=0, plot=FALSE)
g2 <- s.label(rad.pca$scores, xax=1, yax=2, ppoints.col="red", plabels=list(box=list(draw=FALSE), optim=TRUE), paxes.draw=TRUE, pgrid.draw=FALSE, plabels.cex=1, plot=FALSE)
ADEgS(c(g1, g2), layout=c(1, 2))
dev.off()

## StAMPP and distance-based analyses
# Calculate Nei's distances between individuals/pops
rad.d.ind <- stamppNeisD(rad.genlight, pop=FALSE) # Nei's 1972 distance between individuals
stamppPhylip(rad.d.ind, file="neis_dist_inds.phy.dst") # Export matrix for SplitsTree
rad.d.pop <- stamppNeisD(rad.genlight, pop=TRUE) # Nei's 1972 distance between populations
stamppPhylip(rad.d.pop, file="neis_dist_pops.phy.dst") # Export matrix for SplitsTree
