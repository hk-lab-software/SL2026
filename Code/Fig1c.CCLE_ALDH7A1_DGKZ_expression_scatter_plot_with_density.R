#=========================================
# Fig 1c
# ALDH7A1-DGKZ scatter plot with density
# Set working directory
#=========================================
setwd("Fig1c")

#========================================#
# Read data
#========================================#
# CCLE cell line information
CCLE_cell_lines<-read.csv("sample_info.csv", stringsAsFactors = F, header = T)
sort(unique(CCLE_cell_lines$primary_disease))

# Read gene expression data
library(data.table)
library(preprocessCore)
library(ggplot2)

# Log2 scale
CCLE_exp<-fread("CCLE_expression.csv", data.table = F)
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