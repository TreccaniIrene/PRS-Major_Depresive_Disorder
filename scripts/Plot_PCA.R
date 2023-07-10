library(ggfortify)
library(plotly)
library(ggplot2)
library(ggrepel)

# Read in file
all_phase <- read.table( "all_phase3.psam", header = F)
pca <- read.table( "data.pca.eigenvec", header = F)
m <- merge( pca, all_phase, by.x = c("V2"), by.y = c("V1"), all = T)
df <- m[3:4]
colnames(m)["V5.y"]

# Rename some columns
colnames(m)[ colnames(m) == "V5.y" ] <- "Population"
colnames(m)[ colnames(m) == "V3.x" ] <- "PC1"
colnames(m)[ colnames(m) == "V4.x" ] <- "PC2"
row.names(m) <- m$V2

# Plot the result
mpca_res <- prcomp(df,scale. = T)
p <- ggplot( data = m, aes( x = PC1, y = PC2) ) + geom_point( aes( color = Population )) 
ggplotly(p)

# Remove the erroneous value
name <- row.names( m[ m$V2 == "127CO", ])
m <- m[ which( m$V2 != name[]), ]
i <- data.frame( m$V1, m$V2 )
write.table( i, "newdataset.txt", quote = F, row.names = F, col.names = F)