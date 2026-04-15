###################################################################################
# this file contains R codes to prepare the TCGA ESCA data 

# Main source of the data are fomr the UCSC Xena browser
# The file TCGA-ESCA.gene-level_ascat3.tsv can be downloaded from
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-ESCA.gene-level_ascat3.tsv.gz
# The file gencode.v36.annotation.gtf.gene.probemap can be downloaded from
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap
# The file TCGA-ESCA.clinical.tsv can be downloaded from 
# https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-ESCA.clinical.tsv.gz

# Of course, gunzip those .gz files first.

d <- read.table("TCGA-ESCA.gene-level_ascat3.tsv", sep="\t", header=T, row.names=1)
ann <- read.table("gencode.v36.annotation.gtf.gene.probemap", sep="\t", header=T, row.names=1)[-c(60624:60660),]
cc <- read.table("TCGA-ESCA.clinical.tsv", sep="\t", header=T, row.names=1, fill=T)

rownames(cc) <- gsub("-",".", rownames(cc))
common_names <- intersect(rownames(cc), colnames(d))

clin <- cc[common_names, , drop = FALSE]
dat <- d[, common_names, drop = FALSE]
                                    
table(clin[,2])
# Outcome:
#         Adenomas and Adenocarcinomas Cystic, Mucinous and Serous Neoplasms 
#                                   58                                     1 
#              Squamous Cell Neoplasms 
#                                   73 


# Need to remove on case of cyctic ESCA
pmatch("Cystic", as.character(clin[,2]))
#Outcome:
# [1] 114

id.remove <- pmatch("Cystic", as.character(clin[,2]))

clin <- clin[-id.remove,]
dat <- dat[,-id.remove]

# remove missing values for age
id.remove <- which(is.na(clin$age_at_earliest_diagnosis_in_years.diagnoses.xena_derived))
clin <- clin[-id.remove,]
dat <- dat[,-id.remove]

# remove genes with missing values in CNA
row.var.dat <- apply(dat,1,var)
id.remove <- which(is.na(row.var.dat))
dat <- dat[-id.remove,]
ann <- ann[-id.remove,]

# removing chrM, chrX, chrY
id.remove <- which(ann$chrom %in% c("chrM", "chrX", "chrY"))
dat <- dat[-id.remove,]
ann <- ann[-id.remove,]

y <- as.numeric(clin[,2]==clin[2,2])

save(y, dat, ann, clin, file="GDC-TCGA-ESCA-data.RData")

