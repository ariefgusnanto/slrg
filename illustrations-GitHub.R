# This file contains illustration to perform sparse logistic regression fitting
# using TCGA ESCA data. This is done in UNIX/Linux OS.

# Requirement: minimum 8 cores with 32GB RAM, and BLAS (tested on OpenBLAS) installed

# Packages imported (need to be available first):
# DNAcopy from Bioconductor
# wavethresh from CRAN
# MASS (with base R)
# RhpcBLASctl from CRAN
# pROC from CRAN
# CNAseg from https://github.com/ariefgusnanto/slrg
# (download CNAseg_1.6.tar.gz manually from the GitHub page, and install it from local drive)

# The TCGA ESCA data (GFC-TCGA-ESCA-data.RData) can be downloaded from the GitHub page manually
# and the preparation steps to get the file (from UCSC Xenabrowser) is given in preparation-GitHub.R

load("GDC-TCGA-ESCA-data.RData")
source("functions-GitHub.R")

## Creating the fixed predictors (clinical predictors)
temp <- cbind(ifelse(clin$gender=="male", 1, 0),
as.numeric(clin$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived))
colnames(temp) <- c("Gender", "Age")
X <- cbind(1, temp)
# X is the matrix of clinical predictors
# Rows: patients, Columns: predictors
# If you have other variables, they can be included as additional columns in X.
# Note that the first column of X is a vector of ones (intercept).

# Z is the matrix of CNA
# Rows: patients, columns: genomic location (for CNA, this must be ordered)
Z <- as.matrix(t(dat))
rm(dat)
gc()

######### Calling RhpcBLASctl for multi-thread computation
library(RhpcBLASctl)
blas_set_num_threads(8) # 8 is the number of thread


# Running the main function
# Note: the arguments are explicitly shown
res.HLIG <- log.sparse(y=y,X=X,Z=Z, penalty="HLIG", lambda.seq=exp(seq(-2,6,0.5)), 
                       tolerance=1e-3, write=NULL, plot.all=FALSE, opt.crit = "AIC", 
                       alpha.mix=0.5, alpha=1.00001, zero.threshold=1e-4, fixed.b0=FALSE,
                       epsilon=1e-6)

res.HLIGnet <- log.sparse(y=y,X=X,Z=Z, penalty="HLIGnet", lambda.seq=exp(seq(-2,6,0.5)), 
                       tolerance=1e-3, write=NULL, plot.all=FALSE, opt.crit = "AIC", 
                       alpha.mix=0.5, alpha=1.00001, zero.threshold=1e-4, fixed.b0=FALSE,
                       epsilon=1e-6)


## Plotting lambda
par(mfrow=c(1,2))
plot.lambda(res.HLIG, "AIC", lwd=2, xlab="Log lambda", ylab="AIC", main="HLIG")
plot.lambda(res.HLIGnet, "AIC", lwd=2, xlab="Log lambda", ylab="AIC", main="HLIGnet")

## plotting random effects
par(mfrow=c(2,1))
plot.re(res.HLIG, ann$chrom, xlab="Genomic Locations (Chromosome)", main="HLIG", las=1, cex.axis=0.8)
plot.re(res.HLIGnet, ann$chrom, xlab="Genomic Locations (Chromosome)", main="HLIGnet", las=1, cex.axis=0.8)

# Annotation, list of genes
res.HLIG.ann <- ann.out(res.HLIG, ann, "AIC", 1e-3)
res.HLIGnet.ann <- ann.out(res.HLIGnet, ann, "AIC", 1e-3)

# writing the list of genes to file
write.table(res.HLIG.ann, file="ESCA-res.HLIG.ann.txt",  quote=F, col.names=T, row.names=T)
write.table(res.HLIGnet.ann, file="ESCA-res.HLIGnet.ann.txt",  quote=F, col.names=T, row.names=T)


