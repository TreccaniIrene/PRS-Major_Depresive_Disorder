library(ggfortify)
library(plotly)
library(ggplot2)
library(ggrepel)

path1 = "path_dataset1/"
path2 = "path_dataset2/"

# Read in file
all_phase <- read.table( "all_phase3.psam" , header = F)
pca <- read.table( paste( phat1,"data.pca.eigenvec", sep = "", collapse = NULL), header = F)
pca.two<-read.table( paste( phat2,"data.pca.eigenvec", sep = "", collapse = NULL), header = F)

# Merge all_phase and pca by column: chr and SNP
m<-merge(pca,all_phase, by.x = c("V2"),by.y=c("V1"),all=T)
df <- m[3:5]
colnames(m)["V5.y"]

# Rename some columns
colnames(m)[colnames(m)=="V3.x"] <- "PC1"
colnames(m)[colnames(m)=="V4.x"] <- "PC2"
colnames(m)[colnames(m)=="V5.y"] <- "Population"
row.names(m)<-m$V2

# Merge all_phase and pca.two by column: chr and SNP
t <-merge( pca.two, all_phase, by.x = c("V2"), by.y = c("V1"), all = T)
df <- m[3:4]
colnames(t)["V5.y"]

# Rename the column population
colnames(t)[colnames(t) == "V5.y"] <- "Population"
row.names(t) <- t$V2

# Save the name of row which Population column is NA
name <- row.names( m[ which( is.na( m$Population == "NA" ) == T),])
name2 <- row.names( t[ which( is.na( t$Population == "NA" ) == T),])

# Rename the Population column in the db complete in DB and OLD
m$Population[ which( m$V2 == name[] ) ] <- "TargetData1"
m$Population[ which( t$V2 == name2[] ) ] <- "TargetData2"

# Plot the result
mpca_res <- prcomp( df, scale. = T )
p <- ggplot( data = m, aes( x = PC1, y = PC2 ) ) + geom_point( aes( color=Population ) ) 
ggplotly(p)

# Remove the erroneous value
name <- row.names( m[ m$V2 == "127CO", ] )
m <- m[ which( m$V2 != name[] ), ]
i <- data.frame(m$V1,m$V2)
write.table( i,"newdataset.txt", quote=F,row.names=F,col.names = F)
