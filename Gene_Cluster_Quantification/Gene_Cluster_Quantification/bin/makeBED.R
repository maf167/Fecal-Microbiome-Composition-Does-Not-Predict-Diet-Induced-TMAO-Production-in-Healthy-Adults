##
## Finding operon genes and junctions in KEGG
##

# Usage: Rscript makeBED.R <ssdb_table.csv> <output path> [KIDs]

# ssdb_table is the output of a KEGG SSDB Gene Cluster search (you can convert the html to csv with Excel)
# annotation table contains semicolon-delimited descriptors in a table with cells matcthing the SSDB table
# output path - includes file prefix

### Parsing Arguments
args = commandArgs(trailingOnly = TRUE)
if(!(length(args) %in% 2:3)) {print("makeBED.R: Wrong Number of Arguments"); break}

clusters <- read.csv(args[1], header=F)

outPath = args[2]

if(length(args) == 2) KID = NULL else KID <- strsplit(args[3], ",")[[1]]

### Look up KEGG Gene Data
library(KEGGREST)
GeneTable <- data.frame(chr = character(),Start=numeric(), End=numeric(),
                        Gname = character())
BGCTable  <- data.frame(chr = character(),Start=numeric(), End=numeric(), 
                        BGCname = character())

# For each gene cluster homolog
for(i in 1:nrow(clusters)){
  strain = clusters[i,1]
  # For each gene
  for(k in 2:ncol(clusters)){
    if(!is.na(clusters[i,k]) & nchar(as.character(clusters[i,k])) > 0){
      
      # Z is an ID for the strain-gene pair
      Z = paste(strain, as.character(clusters[i,k]), sep=":")
      
      if(any(grepl("\\(",Z))) Z = strsplit(Z,"\\(")[[1]][1]

      # Query KEGG GENOME for strain-gene pair and extract position data to GeneTable
      g = keggGet(Z)
      if("POSITION" %in% names(g[[1]])){
        pos = g[[1]]$POSITION
        if(length(grep("complement", pos)) > 0){
          pos  = substr(pos, 12, nchar(pos)-1)
        }
        s = as.numeric(strsplit(pos, "\\.\\.")[[1]][1])
        e = as.numeric(strsplit(pos, "\\.\\.")[[1]][2])
        GeneTable <- rbind(GeneTable, 
                           data.frame(chr = strain,Start=s,End=e,
                                      Gname = as.character(clusters[i,k])
                                      )
                           )
      }
    }
  }
}



Kstrains <- as.character(unique(GeneTable$chr))

# For each strain, extract gene cluster position from gene data

for(i in Kstrains){
  tmpTab <- GeneTable[GeneTable$chr == i,]
  s <- min(c(tmpTab$Start, tmpTab$End)); e <- max(c(tmpTab$Start, tmpTab$End))
  BGCTable <- rbind(BGCTable, data.frame(chr = i,Start=s, End=e, BGCname = i))
}

## Write out

write.table(GeneTable, file = paste(outPath, "genes.bed", sep=""),  quote=F,
            row.names=F, col.names=F, sep="\t")
write.table(BGCTable, file = paste(outPath, "bgc.bed", sep=""),  quote=F, 
            row.names=F, col.names=F, sep="\t")

write.table(Kstrains, file = paste(outPath, "Kstrains.txt", sep=""),  
            quote=F, row.names=F, col.names=F, sep="\t")

for(K in KID){
  Ktable <- GeneTable[grep(K,GeneTable$Gname),]
  write.table(Ktable, file = paste(outPath,"_",K, ".bed", sep=""),  
              quote=F, row.names=F, col.names=F, sep="\t")
}