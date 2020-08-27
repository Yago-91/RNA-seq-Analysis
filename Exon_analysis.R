#Scrip for Exon quantification

#Parsing arguments
args <- commandArgs(trailingOnly = TRUE)

#Setting variables
table_dir <- paste(args[1], "/merge_exon_quant.table", sep = "")
table_dir
n_samples <- as.integer(args[2])
n_samples
Names <- unlist(strsplit(args[3], ","))
Names
exon_table <- read.csv(table_dir, sep = " ", header = TRUE)

#Analysis
for (col in 6:(n_samples+5)){
  colsum <- sum(exon_table[,col])
  exon_table[,col] <- exon_table[,col]/colsum*100
}

colnames(exon_table) <- c(colnames(exon_table)[1:5],Names)

#Saving data

write.csv(exon_table, paste(args[1], "/by_Exon_analysis.csv", sep = ""), row.names = FALSE)



