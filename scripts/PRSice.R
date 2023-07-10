# Read covariate file and extract relevant columns
covariate <- read.table("final.data.fam", header = F)
cov <- data.frame(covariate$V1, covariate$V2, covariate$V5)
colnames(cov) <- c("FID", "IID", "Sex")

# Read principal component scores file and select first 4 columns (PC1-PC4)
pcs <- read.table("chr.data.pca.eigenvec", header = F)
pcs <- pcs[, 1:4]
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:2))

# Merge covariate and principal component data based on FID and IID columns
covar <- merge(cov, pcs, by = c("FID", "IID"))

# Write merged data to a new file named "file.covariate"
write.table(covar[, 1:5], "file.covariate", quote = F, row.names = F, col.names = T)