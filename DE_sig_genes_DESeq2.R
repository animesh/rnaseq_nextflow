# Differential Expression Analysis with DESeq2
# Faithful to the Griffith Lab RNA-seq tutorial (without external ID mapping file)
# https://rnabio.org/module-03-expression/0003/03/03/Differential_Expression-DESeq2/

# 1. Load required libraries
library(DESeq2)
suppressPackageStartupMessages(library(data.table))
library(dplyr)

# 2. Read htseq-count matrix from current working directory
htseqCounts <- fread("gene_read_counts_table_all_final.tsv")

# 3. Format htseq counts data for DESeq2
htseqCounts <- as.matrix(htseqCounts)
rownames(htseqCounts) <- htseqCounts[, "GeneID"]
htseqCounts <- htseqCounts[, colnames(htseqCounts) != "GeneID"]
class(htseqCounts) <- "integer"

# 4. Filter raw counts (require at least 1 of 6 samples with counts >= 10)
htseqCounts <- htseqCounts[which(rowSums(htseqCounts >= 10) >= 1), ]

# 5. Specify experimental design
metaData <- data.frame("Condition" = c("UHR", "UHR", "UHR", "HBR", "HBR", "HBR"))
metaData$Condition <- factor(metaData$Condition, levels = c("HBR", "UHR"))
rownames(metaData) <- colnames(htseqCounts)

# 6. Construct DESeqDataSet and run analysis
dds <- DESeqDataSetFromMatrix(countData = htseqCounts, colData = metaData, design = ~ Condition)
dds <- DESeq(dds)
res <- results(dds)

# 7. Format results and merge raw counts (using Ensembl GeneIDs directly)
deGeneResult <- as.data.frame(res)
deGeneResult$ensemblID <- rownames(deGeneResult)
deGeneResult <- as.data.table(deGeneResult)

original_counts <- as.data.frame(htseqCounts)
original_counts[, "ensemblID"] <- rownames(htseqCounts)
deGeneResult <- merge(deGeneResult, original_counts, by = "ensemblID", all.x = TRUE)

# 8. Sort and filter significant genes (FDR < 0.05)
deGeneResultSorted <- deGeneResult[order(deGeneResult$padj), ]
deGeneResultSignificant <- deGeneResultSorted[deGeneResultSorted$padj < 0.05, ]

# 9. Save output files to current working directory
fwrite(deGeneResultSorted, file = "DE_all_genes_DESeq2.tsv", sep = "\t")
fwrite(deGeneResultSignificant, file = "DE_sig_genes_DESeq2.tsv", sep = "\t")
saveRDS(dds, "dds.rds")
saveRDS(res, "res.rds")

cat("Analysis complete. Generated DE_sig_genes_DESeq2.tsv and DE_all_genes_DESeq2.tsv in CWD.\n")