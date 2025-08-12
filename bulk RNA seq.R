#data:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272922
#=================================
#Load required Packages
#=================================
install.packages("DESeq2")
library("DESeq2")
#=================================
#set paths and Sample info
#==================================
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

#===========================================
#Data cleaning
#===========================================
class(count_data)
is.na(count_data) #check missing values
which(is.na(count_data),arr.ind = T)
install.packages("zoo")
library(zoo)
#mean of gene values across sample are taken to fill in missing value
count_data[] = t(na.aggregate(t(count_data)))
count_data
