#https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

suppressMessages(library(tximport))
suppressMessages(library(DESeq2))

#Parsing command line arguments
args <- commandArgs(trailingOnly = T)

#Setting variables
loc_dir <- args[1]
quant_dirs <- unlist(strsplit(args[2], ","))
Species <- args[3]
Names <- unlist(strsplit(args[4], ","))
tx2gene_path <- args[5]
Groups_R <- unlist(strsplit(args[6], ","))
Group_A_name <- Groups_R[1]
Group_A_number <- as.integer(Groups_R[2])
Group_B_name <- Groups_R[3]
Group_B_number <- as.integer(Groups_R[4])

#Directories of all the files
files <- file.path(loc_dir, quant_dirs, "quant.sf")
print("Location directories of quant.sf files")
files
names(files) <- Names
print(paste("All passed files exist? ",all(file.exists(files))))

#Correspondance file
tx2gene <- read.csv(tx2gene_path, header = T, colClasses = c(rep("character",2)))

#Importing data
txi.salmon <- tximport(files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#Configuring differential analysis
sampleTable <- data.frame(condition = factor(rep(c(Group_B_name,Group_B_name), c(Group_A_number,Group_B_number))))
rownames(sampleTable) <- colnames(txi.salmon$counts)
dds <- DESeqDataSetFromTximport(txi.salmon, sampleTable, ~condition)
#Stablish reference group
dds$condition <- relevel(dds$condition, ref = Group_A_name)

#Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)

#Save results
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=paste(Group_A_name,"_vs_",Group_B_name,".csv",sep = ""))

#Save gene-wise analysis data in raw TPM
out_1 <- paste(loc_dir, "/Salmon/Salmon_merge_genes.csv", sep = "")
write.csv(txi.salmon$counts, out_1)