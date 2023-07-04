# Load the file
dat <- read.table( "clear.het", header = T)
# Calculate the mean from column F (Method-of-moments F coefficient estimate)
m <- mean(dat$F) 
# Calculate the standard deviation  from column F (Method-of-moments F coefficient estimate)
s <- sd(dat$F)
valid <- subset(dat, F <= m+3*s & F >= m-3*s) 
# Create a result file
write.table( valid[,c(1,2)] , "clear.het.valid.sample", quote = F, row.names = F, col.names = T)