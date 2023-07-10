# Load the file 
het <- read.table( "clear.het.valid.sample", header=T)
sex <- read.table( "clear.sexcheck", header=T)
valid <- subset(sex, STATUS == "OK" & FID %in% het$FID)
# Create a result file
write.table(valid[,c("FID", "IID")], "sex.valid", row.names=F, col.names=T, sep="\t", quote=F) 