####Environment set up####
getwd()
setwd("/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R")
#working dir is location of featureCounts table: "/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R"

#STEP 5
#install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#install packages that will be used in this analysis
#DESeq2, clusterProfiler, enrichplot, org.Hs.eg.db, 
BiocPackages <- c("DESeq2","pheatmap","clusterProfiler","org.Hs.eg.db","enrichplot","ggplot2", "ggrepel")
BiocManager::install(BiocPackages)

#load libraries
library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggrepel)


####Importing and Cleaning up of Data####
#import data
ftCountsData <- read.csv("ftCounts.txt", sep="\t")
#ftCountsData is a data.frame

#remove unwanted columns: Chr, Start, End, Strand, Length
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "Chr")]
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "Start")]
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "End")]
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "Strand")]
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "Length")]

#establish that the row names of the data.frame ftCountsData are the Gene IDs for downstream analysis
rownames(ftCountsData) <- ftCountsData$Geneid

#can now also get rid of the Geneid column
ftCountsData <- ftCountsData[, -which(names(ftCountsData) == "Geneid")]

#Clean up the column names from the file paths in the IBU cluster to the sample names
newColumnNames <- c("HER21","HER22","HER23","NonTNBC1","NonTNBC2","NonTNBC3","Normal1", "Normal2","Normal3","TNBC1","TNBC2","TNBC3")
names(ftCountsData) <- newColumnNames

####Creation of the Metadata Table####
#README with metadata info missing from /data/courses/rnaseq/breastcancer_de
samples <- c("HER21","HER22","HER23","NonTNBC1","NonTNBC2","NonTNBC3","Normal1", "Normal2","Normal3","TNBC1","TNBC2","TNBC3")
groups <- factor(c("HER2","HER2","HER2","NonTNBC","NonTNBC","NonTNBC","Normal","Normal","Normal","TNBC","TNBC","TNBC"))
replicates <- c(1,2,3)
metadata<-data.frame(sample=samples, group=groups)
metadata
rownames(metadata) <- metadata$sample

#evaluate if there is a match or mismatch in the order of the sample names 
missing_in_counts <- setdiff(rownames(metadata), colnames(ftCountsData))
missing_in_metadata <- setdiff(colnames(ftCountsData), rownames(metadata))

if (length(missing_in_counts) > 0 || length(missing_in_metadata) > 0) {
  stop("Mismatch between ftCountsData and metadata: Check sample names.")
}

#Check for NA values in ftCountsData
if (any(is.na(ftCountsData))) {
  print("NA values found in the ftCountsData.")
}

####DESeq2####
#using DESeqDataSetFromMatrix because the count data is from featureCounts
dds <- DESeqDataSetFromMatrix (countData=ftCountsData, colData=metadata, design=~group) 

#Run DESeq2
#overrides previous value stored in dds
dds <- DESeq(dds)

####PCA plot####
#Removing the dependence of the variance on the mean
#vsd is now dds without dependence of the variance on the mean
vsd <- vst(dds, blind = TRUE)

#In some way, visualise how the samples cluster based on their gene expression profiles
#Assess how the samples cluster based on their gene expression profiles
#intgroup - group
plotPCA(vsd, intgroup="group")
##prettify the plot
##maybe save it as a variable?

#Briefly comment on the observed pattern and what it means for downstream analysis
#The Normal, TNBC and NonTNBC samples cluster close to themselves.
#The HER2 sample on the other hand has two replicates that cluster together but one replicate is a distance away
#This could affect downstream analysis as by making it difficult to perform the contrast analysis and draw out the differences between one sample and another

#STEP 6
#####Pairwise contrast####
#TNBC vs NonTNBC - chosen because they are clustered together per sample. 
#HER2 does not show expected clustering patterns
#can't do results on vsd because vsd is class DESeqTransform, but results requires class DESeqDataSet = dds
TNBCvsNonTNBC_contrast <- results(dds, contrast=c("group", "TNBC", "NonTNBC"))
#TNBCvsNonTNBC_contrast is of class data.frame.
#columns = baseMean, log2FoldChange, lfcSE, stat, pvalue, padj
#rows = geneids

#HOW many genes are differentially expressed in the pairwise comp? look at padj < 0.05
#using padj < 0.05 bc these are the statistically significant genes (unlikely to be due to random chance)
#these are the differentially expressed genes from the pairwise contrast
DEgenes_TNBCvsNonTNBC_contrast <- TNBCvsNonTNBC_contrast[!is.na(TNBCvsNonTNBC_contrast$padj) & TNBCvsNonTNBC_contrast$padj < 0.05, ]

#sorted from smallest padj to the greatest
ordered_DEgenes_TNBCvsNonTNBC_contrast <- DEgenes_TNBCvsNonTNBC_contrast[order(DEgenes_TNBCvsNonTNBC_contrast$padj, na.last=NA), ]

#ordered_DEgenes_NonTNBCvsTNBC_contrast is of class DESqesResults. Convert that to data.frame for the csv
write.csv(as.data.frame(ordered_DEgenes_TNBCvsNonTNBC_contrast), "/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R/ordered_DEgenes_TNBCvsNonTNBC_contrast.csv", row.names = TRUE)

DEgenes_TNBCvsNonTNBC_number <- nrow(DEgenes_TNBCvsNonTNBC_contrast)
DEgenes_TNBCvsNonTNBC_number #number of genes that are differentially expressed in TNBC vs NonTNBC = 1728

#HOW many genes are upregulated
#upregulated - positive log2FoldChange, >0
#genes with + log2FoldChange are upregulated in TNBC compared to NonTNBC
upregulatedInTNBC <- ordered_DEgenes_TNBCvsNonTNBC_contrast[ordered_DEgenes_TNBCvsNonTNBC_contrast$log2FoldChange > 0, ]
numberOfGenes_upregulatedInTNBC <- nrow(upregulatedInTNBC) #783

#HOW many genes are downregulated
#downregulated - negative log2FoldChange, <0
#genes with - log2FoldChange are downregulated in TNBC compared to NonTNBC
downregulatedInTNBC <- ordered_DEgenes_TNBCvsNonTNBC_contrast[ordered_DEgenes_TNBCvsNonTNBC_contrast$log2FoldChange < 0,]
numberOfGenes_downregulatedInTNBC <- nrow(downregulatedInTNBC) #945

####Visualization of up and downregulation of DE genes####
#Volcano plot
volcanoPlotData <- as.data.frame(TNBCvsNonTNBC_contrast)
volcanoPlotData$gene <- mapIds(org.Hs.eg.db, 
                               keys = rownames(volcanoPlotData), 
                               keytype = "ENSEMBL", 
                               column = "SYMBOL",
                               multiVals = "first")  #will use the first match

#Remove rows with pvalue=NA or padj=NA or gene=NA
volcanoPlotData <- volcanoPlotData[!is.na(volcanoPlotData$pvalue) & !is.na(volcanoPlotData$padj) & !is.na(volcanoPlotData$gene), ]

#Initialize diffExpressed column
volcanoPlotData$diffExpressed <- "NOTSIG"

#Define significance thresholds
log2FC_threshold <- 0 #should this be 0.6???????????
padj_threshold <- 0.05 #why not just using the DE genes data set (the 1728 that are padj<0.05 and are not NA)

#Classify genes as "UP", "DOWN", or "NOTSIG"
#upregulated
volcanoPlotData$diffExpressed[volcanoPlotData$padj < padj_threshold & 
                                volcanoPlotData$log2FoldChange > log2FC_threshold] <- "UP"

#downregulated
volcanoPlotData$diffExpressed[volcanoPlotData$padj < padj_threshold & 
                                volcanoPlotData$log2FoldChange < -log2FC_threshold] <- "DOWN"

#Write the full volcano plot data to CSV
write.csv(volcanoPlotData, "/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R/volcanoPlotData.csv", row.names = TRUE)


#top 20 genes based on absolute log2FoldChange and significant padj
top_genes <- volcanoPlotData[volcanoPlotData$padj < padj_threshold, ] ##statisitically sig
top_genes <- top_genes[order(abs(top_genes$log2FoldChange), decreasing = TRUE), ][1:20, ]

write.csv(top_genes, "/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R/top_genes_volcanoPlotData.csv", row.names = TRUE)

#Add a label column for the top 20 genes
volcanoPlotData$label <- ifelse(rownames(volcanoPlotData) %in% rownames(top_genes), 
                                volcanoPlotData$gene, 
                                NA)

#Create the volcano plot
ggplot(data=volcanoPlotData, aes(x = log2FoldChange, y = -log10(pvalue), col = diffExpressed)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  ggtitle("Volcano Plot of TNBC vs NonTNBC") +
  scale_color_manual(values = c("UP" = "orange", "DOWN" = "#00AFBB", "NOTSIG" = "grey"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  geom_text_repel(data = subset(volcanoPlotData, !is.na(label)), 
                  aes(label = label), 
                  size = 3, 
                  box.padding = 0.3, 
                  point.padding = 0.2)

#####
#HEATMAP Numero???
# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Select the top 20 most highly expressed genes
ddstopGenes <- rownames(normalized_counts)[order(rowMeans(normalized_counts), decreasing = TRUE)[1:20]]

# Subset the normalized counts for these top genes
ddstopGenes_counts <- normalized_counts[ddstopGenes, ]

# Define colors for annotation
annotation_col <- data.frame(Group = metadata$group)
#rownames(annotation_col) <- rownames(metadata)

# Generate heatmap
pheatmap(ddstopGenes_counts,
         annotation=metadata,
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         scale = "row",  # Standardize values to z-scores per row (gene)
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         fontsize_row = 8, 
         fontsize_col = 10)

#Based on the original publication, select 2-3 genes that are of particular interest 
#and investigate their expression level. 
#You could use, for example, the normalized counts (see DESeq2::counts) 
#where the effect of between-sample differences in sequencing depth has been removed.
##HAVENT DONE YET


####Gene Ontology####
#STEP 7
#identify Gene Ontology terms that contain more differentially expressed genes 
#than expected by chance for the pairwise comparison(s) you examined in step 6

#run overrepresentation analysis using clusterProfiler::enrichGO
TNBCvsNonTNBC_contrast_genes_universe <- rownames(TNBCvsNonTNBC_contrast) #all genes; not DE
TNBCvsNonTNBC_contrast_genes_DE <- rownames(DEgenes_TNBCvsNonTNBC_contrast) #DE genes


ego <- enrichGO(gene = TNBCvsNonTNBC_contrast_genes_DE,
                universe = TNBCvsNonTNBC_contrast_genes_universe,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                keyType = "ENSEMBL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH")
write.csv(as.data.frame(ego), "/Users/amoriki/Desktop/Code/MSc/rnaSeq/rna_seq_R/enrichGOTerms.csv", row.names = TRUE)

head(ego)

barplot(ego)
dotplot(ego)
sessionInfo()
