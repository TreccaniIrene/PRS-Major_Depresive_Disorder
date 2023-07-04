path = "/chr_"
# Read in file
db <- read.table( "Final_Target_QC.fam", header=F )
# Create files for each chromosoma
k=1
for (e in 1:22) {
    target <- read.table("Final_Target_QC.fam", sep = "", collapse = NULL), header=F)
    chr <- read.table( paste( paste( path,k, sep = "", collapse = NULL), "/Chr.fam", sep = "", collapse = NULL), header=F)
    # Create a new file with information of ids
    i <- c(chr[,c(1,2)], target[, c(1,2)])
    write.table( i, paste(paste( path,k, sep = "", collapse = NULL), "/recodedids.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = F)
    # Create a new file with information of sex
    s <- target [,c(1,2,5)]
    write.table( s, paste(paste( path,k, sep = "", collapse = NULL), "/recodedsex.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = F)
    # Create a new file with information of phenotype
    p <- target[,c(1,2,6)]
    write.table( p, paste(paste( path,k, sep = "", collapse = NULL), "/recodepheno.txt", sep = "", collapse = NULL), quote=F, row.names=F,col.names = F)
    k <- k+1
}