#==========================================
#Load required Packages
#==========================================
install.packages("DESeq2")
library("DESeq2")
library(zoo)
#==========================================
#set paths and Sample info
#==========================================
setwd("D:/bioinformatic/bioinfo_prj/Bulk RNA seq")
#load the count datas obtained fron Htseq; check.names FALSE to remain the same name
count_data = read.csv(file = "TCGA-count-data.csv", header = TRUE, sep = "," , row.names = 1, check.names = FALSE )
dim(count_data)
View(count_data)
rownames(count_data)
#Load the data where the samples and respective condition are labelled
col_data = read.csv(file = "TCGA-column-data.csv", header = T,sep = ",",row.names = 1)
View(col_data)
rownames(col_data)
colnames(count_data)
all(rownames(col_data)==colnames(count_data)) #ensure all samples are labelled

#==========================================
#Data cleaning
#==========================================
class(count_data)
is.na(count_data) #check missing values
which(is.na(count_data),arr.ind = T)
#mean of gene values across sample are taken to fill in missing value
count_data[] = t(na.aggregate(t(count_data)))
count_data

#==========================================
# DESeq - for differential expressed genes
#==========================================
#convert categorical data to factor
col_data$condition <- factor(col_data$condition)cou
#Create DESeq Dataset
dds = DESeqDataSetFromMatrix(countData = round(count_data),
                              colData = col_data,
                              design = ~ condition)
#Set reference level will compare cancer to normal
dds$condition <- relevel(dds$condition, ref = "normal")
dds

#Run DESeq
dds <- DESeq(dds)
res1 <- results(dds)
summary(res1)

#==========================================
#DEGs
#==========================================
#highly up-regulated genes 
resSigup <- subset(res1, padj <0.05 & log2FoldChange>1)
write.csv(resSigup, "upregulated_genes.csv")

#highly down regulated genes 
resSigdwn <- subset(res1, padj <0.05 & log2FoldChange< -1)
write.csv(resSigdwn, "downregulated_genes.csv")

#for all differential expressed genes
resSig <- subset(res1, padj<0.05 & log2FoldChange>1 | padj<0.05 & log2FoldChange< -1)
write.csv(resSig, "DE_genes.csv")

#==========================================
#Enseble IDs to gene symbols
#==========================================

#==========================================
#Data visualisation
#==========================================
if(!require("BiocManager", quietly = T ))
  install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

EnhancedVolcano(res1,                        
  lab = "",                    # Labels for each point (empty = no labels shown)
  x = 'log2FoldChange',        
  y = 'pvalue',                
  selectLab = NULL,            # label gene of interest
  title = NULL,                # plot title
  cutoffLineType = 'twodash',  # Style of threshold lines
  cutoffLineWidth = 0.8,       # Thickness of threshold lines
  xlim = c(-8,8),              # x-axis limits
  xlab = bquote(~Log[2]~ 'fold change'),  
  ylim = c(0,12),              # y-axis limits
  ylab = bquote(~-Log[10]~italic(P)),     
  pCutoff = 0.05,              # p-value threshold for significance
  FCcutoff = 1.0,              # Fold change threshold for biological relevance
  #transcriptPointSize = 0.5,  # size of points
  #transcriptLabSize = 4.0,    # size of labels
  colAlpha = 1,                # Transparency of points (1 = opaque)
  shape = 19,                  # Shape of points (19 = solid circle)
  subtitle = NULL,             # Optional: subtitle
  legendPosition = 'top',      # Position of legend
  legendLabSize = 12,          # Font size of legend labels
  legendIconSize = 4.0,        # Size of legend icons
  gridlines.major = FALSE,     # Hide major gridlines
  gridlines.minor = FALSE,     # Hide minor gridlines
  drawConnectors = FALSE,      # Don't draw lines from labels to points
  widthConnectors = 0.2,       # Connector line thickness (if enabled)
  colConnectors = 'grey50',    # Connector line color
  border = 'full'              # Border style around plot
) 

png(file = "DE_volcano.png")
