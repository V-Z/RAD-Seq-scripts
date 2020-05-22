# Author: Vojtěch Zeisek, https://trapa.cz/
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html

# 

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Do not exit on error
options(error=expression(NULL))

## Libraries
# Install packages
# install.packages(pkgs=c("adegenet", "adegraphics", "BH", "pegas", "permute", "plogr", "sf", "spData", "StAMPP", "testthat", "vcfR"), lib="rpackages", repos="https://mirrors.nic.cz/R/", dependencies="Imports")
# Libraries
library(package=adegenet, lib.loc="rpackages")
library(package=adegraphics, lib.loc="rpackages")
library(package=StAMPP, lib.loc="rpackages")
library(package=vcfR, lib.loc="rpackages")

## Input file
vcf.file <- Sys.getenv("VCFR")

## Import SNP data from VCF
vcf.data <- read.vcfR(vcf.file)
# Checks
head(vcf.data)
vcf.data@fix[1:10,1:5]

## Various checks
# In VCF R
# DP - quickcheck
vcf.dp <- extract.gt(vcf.data, element='DP', as.numeric=TRUE)

# Boxplot of DP
pdf("dp_per_sample_boxplot.pdf", width=21, height=7)
	boxplot(vcf.dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)", las=2, cex=0.7)
	abline(h=mean(x=as.numeric(vcf.dp), na.rm=TRUE), col="red")
	dev.off()

# Barplot of DP
pdf("dp_per_sample_barplot.pdf", width=21, height=7)
	barplot(apply(X=arabidopsis.vcf.dp, MARGIN=2, FUN=mean, na.rm=TRUE), las=3)
	title("Mean DP per specimen")
	abline(h=seq(from=0, to=60, by=10), col="red")
	dev.off()

# Heatmpa of DP
pdf("dp_per_sample_heatmap.pdf", width=28, height=28)
	heatmap.bp(x=arabidopsis.vcf.dp, col.ramp=rainbow(n=100, start=0.1))
	title("DP per specimens and loci")
	dev.off()

## Convert to genlight
rad.genlight <- vcfR2genlight(x=vcf.data, n.cores=1)
pop(rad.genlight) <- substr(indNames(rad.genlight), start=1, stop=4) # Add pop names: here beginning of individuals names
# See names
indNames(rad.genlight)
pop(rad.genlight)

## Checking data

# N missing SNPs per sample
write.table(x=summary(t(as.matrix(rad.genlight)))[7,], file="missing_per_sample.tsv", sep="\t")

# Plot of missing data (white) and number of 2nd alleles
pdf("missing_data_2nd_allele.pdf", width=21, height=7)
	glPlot(x=rad.genlight, legend=TRUE, posi="topleft")
	dev.off()

# Sum of the number of second allele in each SNP
rad.freq <- glSum(rad.genlight)
# Plot distribution of (second) allele frequencies
pdf("distrib_2nd_allele.pdf", width=21, height=7)
	hist(x=rad.freq, proba=TRUE, col="gold", xlab="Allele frequencies", main="Distribution of (second) allele frequencies")
	lines(x=density(rad.freq)$x, y=density(rad.freq)$y*1.5, col="red", lwd=3 )
	dev.off()

# Mean number of second allele in each SNP
rad.mean <- glMean(rad.genlight)
rad.mean <- c(rad.mean, 1-rad.mean)
# Plot distribution of allele frequencies
pdf("mean_num_2nd_allele.pdf", width=21, height=7)
	hist(x=rad.mean, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies", main="Distribution of allele frequencies", nclass=20)
	lines(x=density(rad.mean, bw=0.05)$x, y=density(rad.mean, bw=0.05)$y*2, lwd=3)
	dev.off()

# Location of the missing values (NA)
rad.na.density <- density(glNA(rad.genlight), bw=10)
# Set range of xlim parameter from 0 to the length of original alignment
pdf("distrib_missing_data.pdf", width=21, height=7)
	plot(x=rad.na.density, type="n", xlab="Position in the alignment", main="Location of the missing values (NA)", xlim=c(0, 1701))
	polygon(c(rad.na.density$x, rev(rad.na.density$x)), c(rad.na.density$y, rep(0, length(rad.na.density$x))), col=transp("blue", alpha=0.3))
	points(glNA(rad.genlight), rep(0, nLoc(rad.genlight)), pch="|", cex=2, col="blue")
	dev.off()

## PCA
rad.pca <- glPca(rad.genlight, nf=5, center=TRUE, scale=FALSE, loadings=TRUE)

# Basic plot
pdf("pcoa.pdf", width=14, height=7)
	scatter.glPca(x=rad.pca, posi="bottomright")
	title("PCoA")
	dev.off()

# Loadings
pdf("pcoa_loadings.pdf", width=14, height=7)
	loadingplot.glPca(x=rad.pca)
	dev.off()

# Nicer plots
# Axis 1 and 2
pdf("pcoa_ax_12.pdf", width=14, height=7)
	g1 <- s.class(rad.pca$scores, pop(rad.genlight),  xax=1, yax=2, col=transp("blue", 0.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=TRUE, pgrid.draw=FALSE, plot=FALSE)
	g2 <- s.label(rad.pca$scores, xax=1, yax=2, ppoints.col="red", plabels=list(box=list(draw=FALSE), optim=TRUE), paxes.draw=TRUE, pgrid.draw=FALSE, plabels.cex=1, plot=FALSE)
	ADEgS(c(g1, g2), layout=c(1, 2))
	dev.off()

# Axis 1 and 3
pdf("pcoa_ax_13.pdf", width=14, height=7)
	g1 <- s.class(rad.pca$scores, pop(rad.genlight),  xax=1, yax=3, col=transp("blue", 0.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=TRUE, pgrid.draw=FALSE, plot=FALSE)
	g2 <- s.label(rad.pca$scores, xax=1, yax=3, ppoints.col="red", plabels=list(box=list(draw=FALSE), optim=TRUE), paxes.draw=TRUE, pgrid.draw=FALSE, plabels.cex=1, plot=FALSE)
	ADEgS(c(g1, g2), layout=c(1, 2))
	dev.off()

# Ploidy - differentiated plots
pdf ("pcoa_ploidycol_ax_12.pdf", width=14, height=7)
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

