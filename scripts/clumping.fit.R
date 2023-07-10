p.threshold <- c(0.001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5)

# Read the phenotype file 
phenotype <- read.table( "final.data.fam", header = F)
phen <- data.frame ( phenotype$V1, phenotype$V2, phenotype$V6)
colnames(phen) <- c( "FID",	"IID",	"Phenotype" )
write.table( phen, "phenotype.phen", row.names = F, col.names = T, quote = F, sep = "\t") 

# Read the PCs
pcs <- read.table( "data.pca.eigenvec", header = T)
pcs <- pca[,1:4]
colnames(pcs) <- c( "FID", "IID", paste0( "PC", 1:2)) 

# Read the covariate
covariate <- read.table( "final.data.fam", header=F)
cov<-data.frame( covariate$V1, covariate$V2, covariate$V5)
colnames(cov)<-c("FID",	"IID",	"Sex")

write.table( cov, "sex.sex", row.names = F, col.names = F, quote = F, sep = "\t") 
pheno <- merge( merge( phen, cov, by = c( "FID", "IID" )), pcs, by = c( "FID", "IID" ))
null.model <- lm( Phenotype~., data = pheno[ , !colnames( pheno ) %in% c( "FID", "IID" )])
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL

for(i in p.threshold){
  
    # Go through each p-value threshold
    prs <- read.table( paste0("clumping.",i,".profile"), header=T)
    pheno.prs <- merge( pheno, prs[ , c("FID","IID", "SCORE")], by = c("FID", "IID"))
    model <- lm( Phenotype~., data = pheno.prs[ , !colnames( pheno.prs ) %in% c("FID","IID")])
    model.r2 <- summary(model)$r.squared
    
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    
    # Obtain the coeffcient and p-value of association of PRS
    prs.coef <- summary(model)$coeff["SCORE", ]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    prs.result <- rbind(prs.result, data.frame( Threshold = i, R2 = prs.r2, P = prs.p, BETA = prs.beta, SE = prs.se))
}
# Best result is:
prs.result[ which.max( prs.result$R2 ), ]
write.table( prs.result[ which.max ( prs.result$R2 ), ], "result.txt",row.names = F, col.names = T, quote = F, sep = "\t") 
write.table( prs.result,  "Clumping.prs.result.txt", row.names = F, col.names = T, quote = F,sep ="\t") 