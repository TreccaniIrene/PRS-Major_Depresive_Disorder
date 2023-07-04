library(ggplot2)
library(factoextra)

# Read in file
eigenval <- read.table( "data.pca.eigenval", header = F)
eigenval <- eigenval[1:5,]

# First convert to percentage variance explained
pve <- data.frame(PC = 1:5, pve = eigenval/sum(eigenval)*100)

# Plot the result
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()