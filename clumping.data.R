# Read in file
dat <- read.table( "LifetimeMDD", header = T)
# Add BETA column
dat$BETA <- log(dat$OR)
# Save the new file
write.table( dat, "BaseDataset", quote = F, row.names = F, col.names = T)