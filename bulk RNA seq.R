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

#for all differentially expressed genes
resSig <- subset(res1, padj<0.05 & log2FoldChange>1 | padj<0.05 & log2FoldChange< -1)
write.csv(resSig, "DE_genes.csv")
