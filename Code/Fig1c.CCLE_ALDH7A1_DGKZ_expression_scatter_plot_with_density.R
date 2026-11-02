# ALDH7A1-DGKZ scatter plot with density
# Set working directory
setwd("E:/KHH/OneDrive/HSK Lab/Synthetic_Lethality/Tumor_type/")

#========================================#
# Read gene pair data
#========================================#
gene_pairs<-read.delim("Results/Tables/Gene_pairs_and_metabolite.txt", sep = "\t", header = T, stringsAsFactors = F)
SL_genes<-read.delim("Results/Tables/Genes_in_gene_pairs.txt", sep = "\t", header = T, stringsAsFactors = F)

# CCLE cell line information
CCLE_cell_lines<-read.csv("E:/KHH/OneDrive/HSK Lab/Data/Cancer_Dependency_Map/sample_info.csv", stringsAsFactors = F, header = T)
sort(unique(CCLE_cell_lines$primary_disease))

# Read gene expression data
library(data.table)
library(preprocessCore)
library(ggplot2)

# Log2 scale
CCLE_exp<-fread("E:/KHH/OneDrive/HSK Lab/Data/Cancer_Dependency_Map/CCLE_expression.csv",
                data.table = F)
rownames(CCLE_exp)<-CCLE_exp[,1]
CCLE_exp<-CCLE_exp[,-1]

for(i in 1:nrow(CCLE_exp)){
  tmp<-subset(CCLE_cell_lines, CCLE_cell_lines$DepMap_ID==rownames(CCLE_exp)[i])
  if(nrow(tmp)!=0){
    rownames(CCLE_exp)[i]<-tmp$CCLE_Name
  }else{
    print(rownames(CCLE_exp)[i])
  }
}
#rownames(CCLE_exp)[which(rownames(CCLE_exp)=="ACH-001316")]<-"PEDS005TSUSP_KIDNEY"
CCLE_exp<-CCLE_exp[-which(rownames(CCLE_exp)=="ACH-001316"),]

cancer_type<-sort(unique(sapply(1:nrow(CCLE_exp),
                                function(x) paste0(strsplit(rownames(CCLE_exp)[x], split = "_", fixed = T)[[1]][-1], collapse = "_"))))
idx<-c(grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", colnames(CCLE_exp)),
       grep("ENGINEERED", rownames(CCLE_exp)),
       grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", rownames(CCLE_exp)))
CCLE_exp<-CCLE_exp[-idx,]


colnames(CCLE_exp)<-sapply(1:ncol(CCLE_exp), function(x) strsplit(colnames(CCLE_exp)[x], split = " (", fixed = T)[[1]][1])
CCLE_exp<-t(CCLE_exp)

exp<-CCLE_exp



# Quantile normalization
exp_quantile<-normalize.quantiles(exp)
rownames(exp_quantile)<-rownames(exp)
colnames(exp_quantile)<-colnames(exp)

exp<-exp_quantile
rm(exp_quantile)

# Z-score normalization
exp<-t(scale(t(exp)))

df<-data.frame(ALDH7A1 = exp["ALDH7A1",], DGKZ = exp["DGKZ",])


plot(x = exp["ALDH7A1",], y = exp["DGKZ",], pch = 19, cex = 1.2,
     xlab = "ALDH7A1", ylab = "DGKZ", main = "CCLE")
dev.off()
col1<-rgb(178, 34, 34, max = 255, alpha = 180)
col2<-rgb(65, 105, 225, max = 255, alpha = 85)
col3<-rgb(0, 0, 0, max = 255, alpha = 100)

p <- ggplot(df, aes(ALDH7A1, DGKZ)) + geom_point(size = 4, color = "black") + theme_classic() +
  labs(title = "CCLE", xlab = "ALDH7A1 gene expression", ylab = "DGKZ gene expression", size = 15)
ggExtra::ggMarginal(p, type = "density", color = col1, fill = col1)

df<-data.frame(AGPAT1 = exp["AGPAT1",], AGPAT2 = exp["AGPAT2",])


plot(x = exp["AGPAT1",], y = exp["AGPAT2",], pch = 19, cex = 1.2,
     xlab = "AGPAT1", ylab = "AGPAT2", main = "CCLE")
dev.off()
col1<-rgb(178, 34, 34, max = 255, alpha = 180)
col2<-rgb(65, 105, 225, max = 255, alpha = 85)
col3<-rgb(0, 0, 0, max = 255, alpha = 100)

p <- ggplot(df, aes(AGPAT1, AGPAT2)) + geom_point(size = 4, color = col3) + theme_classic() +
  labs(title = "CCLE", xlab = "AGPAT1 gene expression", ylab = "AGPAT2 gene expression", size = 15)
ggExtra::ggMarginal(p, type = "density", color = col1, fill = col1)


AGPAT1<-exp["AGPAT1",]
AGPAT2<-exp["AGPAT2",]

max_real<-colMaxs(exp[c("AGPAT1", "AGPAT2"),])
plot(x = 1:length(max_real), y = max_real, pch = 19, col = col3, ylab = "Max real", xlab = "Samples")

sorted.max_real<-sort(max_real,decreasing = F)
plot(x = 1:length(max_real), y = sorted.max_real, pch = 19, col = col1, ylab = "Sorted max real", xlab = "Sorted samples")




# 1. Get a random mixed vector of geneA and geneB, then get max(ran-geneA,ran-geneB) matrix (mean of 1000 iterations)
ran_geneA<-t(apply(matrix(nrow = 1000, ncol=length(AGPAT1)), 1, function(x) sample(AGPAT1)))
ran_geneB<-t(apply(matrix(nrow = 1000, ncol=length(AGPAT2)), 1, function(x) sample(AGPAT2)))

ran_max_geneA_geneB<-matrix(nrow=nrow(ran_geneA), ncol=ncol(ran_geneA))
for(i in 1:1000){
  tmp<-rbind(ran_geneA[i,], ran_geneB[i,])
  ran_max_geneA_geneB[i,]<-sort(colMaxs(tmp))
}
# 2. Sort max values in order and subject each values with the same rank
mean_ran_geneA_geneB<-colMeans(ran_max_geneA_geneB)
points(1:length(mean_ran_geneA_geneB), y = mean_ran_geneA_geneB, col = col2, pch = 19)

# 3. Get lowest 25% of each max vectors
Q1_max_geneA_geneB<-sorted.max_real[1:floor(length(sorted.max_real)/4)]
Q1_mean_ran_geneA_geneB<-mean_ran_geneA_geneB[1:floor(length(mean_ran_geneA_geneB)/4)]
plot(x = 1:length(Q1_max_geneA_geneB), y = Q1_max_geneA_geneB,
     pch = 19, col = col1, ylab = "Sorted max real", xlab = "Sorted samples")
points(1:length(Q1_mean_ran_geneA_geneB), y = Q1_mean_ran_geneA_geneB, col = col2, pch = 19)

# 4. Perform KS-test and draw plot
result<-ks.test(Q1_mean_ran_geneA_geneB, Q1_max_geneA_geneB, alternative = c("greater"))
distance<-KS.diss(Q1_mean_ran_geneA_geneB, Q1_max_geneA_geneB)



















