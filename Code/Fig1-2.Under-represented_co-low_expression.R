#============================================================================#
# Fig1a
# 2.Under-represented_co-low_expression
#============================================================================#
# Set working directory
setwd("Fig1")

#========================================#
# FUNCTION
#========================================#
library(matrixStats)
#========================================#
# Calculate horizontal distance(From package "provenance")
#========================================#
KS.diss<-function(x, y){
  xy<-c(x, y)
  cad1<-stats::ecdf(x)
  cad2<-stats::ecdf(y)
  max(abs(cad2(xy) - cad1(xy)))
}
#========================================#
# Main Synthetic Lethality Test
# exp_of_geneA: vector of geneA expression
# exp_of_geneB: vector of geneB expression
# name_of_geneA: name of geneA
# name_of_geneB: name of geneB
#========================================#
SL_KS_test<-function(exp_of_geneA, exp_of_geneB, name_of_geneA, name_of_geneB){
  
  outfile=paste(name_of_geneA, name_of_geneB, sep = "-")
  pdf(file=paste(outfile,'.pdf',sep=''),width=15,height=15)
  plot(exp_of_geneA, exp_of_geneB, xlab = name_of_geneA, ylab = name_of_geneB,
       main = paste(name_of_geneA, name_of_geneB, sep = "-"), pch=19, cex = 1.5, cex.lab = 2, las = 1, cex.axis = 2, cex.main = 2)
  dev.off()
  
  max_geneA_geneB<-sapply(1:length(exp_of_geneA), function(x) max(exp_of_geneA[x], exp_of_geneB[x]))
  # 1. Get a random mixed vector of geneA and geneB, then get max(ran-geneA,ran-geneB) matrix (mean of 1000 iterations)
  ran_geneA<-t(apply(matrix(nrow = 1000, ncol=length(exp_of_geneA)), 1, function(x) sample(exp_of_geneA)))
  ran_geneB<-t(apply(matrix(nrow = 1000, ncol=length(exp_of_geneA)), 1, function(x) sample(exp_of_geneB)))
  
  ran_max_geneA_geneB<-matrix(nrow=nrow(ran_geneA), ncol=ncol(ran_geneA))
  for(i in 1:1000){
    tmp<-rbind(ran_geneA[i,], ran_geneB[i,])
    ran_max_geneA_geneB[i,]<-sort(colMaxs(tmp))
  }
  # 2. Sort max values in order and subject each values with the same rank
  max_geneA_geneB<-sort(max_geneA_geneB)
  mean_ran_geneA_geneB<-colMeans(ran_max_geneA_geneB)
  
  # 3. Get lowest 25% of each max vectors
  Q1_max_geneA_geneB<-max_geneA_geneB[1:floor(length(max_geneA_geneB)/4)]
  Q1_mean_ran_geneA_geneB<-mean_ran_geneA_geneB[1:floor(length(mean_ran_geneA_geneB)/4)]
  
  # 4. Perform KS-test and draw plot
  result<-ks.test(Q1_mean_ran_geneA_geneB, Q1_max_geneA_geneB, alternative = c("greater"))
  distance<-KS.diss(Q1_mean_ran_geneA_geneB, Q1_max_geneA_geneB)
  
  outfile=paste(name_of_geneA, name_of_geneB, sep = "-")
  pdf(paste("CDF_of_", outfile,'.pdf',sep=''), width=15, height=7)
  
  par(mfrow=c(1,2))
  plot(exp_of_geneA, exp_of_geneB,
       xlab= name_of_geneA, ylab= name_of_geneB, main= paste(name_of_geneA, name_of_geneB, sep = "-"),
       pch=19, cex = 1.5, cex.lab = 2, las = 1, cex.axis = 2, cex.main = 2)
  plot(ecdf(Q1_mean_ran_geneA_geneB), lty = 1, col="firebrick1", xlab = "Values", ylab = "CDF",
       main = paste("CDF of ", name_of_geneA, "-", name_of_geneB, sep = ""), cex.lab = 2, cex.axis = 2, cex.main = 2)
  lines(ecdf(Q1_max_geneA_geneB),lty = 1, col="darkblue")
  text(Q1_mean_ran_geneA_geneB[length(Q1_mean_ran_geneA_geneB)/4], y = 0.6,
       labels = paste0("P-value=", formatC(result$p.value,format = "e",digits=2)))
  text(Q1_mean_ran_geneA_geneB[length(Q1_mean_ran_geneA_geneB)/4], y = 0.5,
       labels = paste0("Distance=", round(distance, digits=4)))
  legend(x = Q1_mean_ran_geneA_geneB[1], y = 1, legend = c("Random", "Real"), col = c("firebrick1", "darkblue"), lty = 1)
  dev.off()
  
  # 5. Get horizontal max value of Real-Random and position
  real_random<-(Q1_max_geneA_geneB - Q1_mean_ran_geneA_geneB)
  
  return(c(result$p.value, distance, max(real_random), which(real_random==max(real_random))[1]))
}
#========================================#
# Read gene pair data
#========================================#
# Newly added gene pairs
gene_pairs<-read.delim("Gene_pairs_all_genes_paired_in_pathway.txt",
                       sep = "\t", header = T, stringsAsFactors = F)
genes<-sort(unique(c(gene_pairs$GeneA, gene_pairs$GeneB)))
SL_genes<-read.delim("Genes_in_gene_pairs.txt",
                     sep = "\t", header = T, stringsAsFactors = F)
SL_genes<-subset(SL_genes, SL_genes$hgnc_symbol%in%genes)

# CCLE cell line information
CCLE_cell_lines<-read.csv("sample_info.csv", stringsAsFactors = F, header = T)
sort(unique(CCLE_cell_lines$primary_disease))

# Read gene expression data
library(data.table)
CCLE_exp<-fread("CCLE_RNAseq_genes_rpkm_20180929.gct", data.table = F)
CCLE_exp[,-c(1:2)]<-log2(CCLE_exp[,-c(1:2)]+1)

# Convert metabolic SL genes to KEGG id
idx<-which(SL_genes$hgnc_symbol%in%CCLE_exp$Description)

# Genes not in CCLE data
tmp_genes<-SL_genes[-idx,]

# Convert genes that are not in the table to ENSEMBL ID
library(biomaRt)
mart<-useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

tmp_genes1<-getBM(attributes = c("hgnc_symbol", "entrezgene_id", "ensembl_gene_id"), filters = "hgnc_symbol",
                  values = tmp_genes$hgnc_symbol, mart = mart)
tmp_genes[which(!tmp_genes$hgnc_symbol%in%tmp_genes1$hgnc_symbol),]

for(i in 1:nrow(tmp_genes1)){
  tmp<-grep(tmp_genes1$ensembl_gene_id[i], CCLE_exp$Name)
  if(length(tmp)!=0){
    CCLE_exp$Description[tmp]<-tmp_genes1$hgnc_symbol[i]
  }
}
rm(tmp ,tmp_genes, tmp_genes1)

# Check
all(SL_genes$hgnc_symbol%in%CCLE_exp$Description)

tmp<-as.data.table(CCLE_exp)
tmp<-tmp[,-1]
tmp<-tmp[, lapply(.SD, mean), by = 'Description']
CCLE_exp<-as.data.frame(tmp)
rm(tmp)

rownames(CCLE_exp)<-CCLE_exp$Description
CCLE_exp<-CCLE_exp[,-1]
dim(CCLE_exp[SL_genes$hgnc_symbol,])

#=========#
# SL test
#=========#
set.seed(1111)
# Multicore
library(foreach)
library(doParallel)
library(preprocessCore)

exp<-CCLE_exp

idx<-grep("HAEMATOPOIETIC_AND_LYMPHOID_TISSUE", colnames(exp))

exp<-exp[,-idx]

# Quantile normalization
exp_quantile<-normalize.quantiles(as.matrix(exp))
rownames(exp_quantile)<-rownames(exp)
colnames(exp_quantile)<-colnames(exp)

exp<-exp_quantile
rm(exp_quantile)

# Z-score normalization
exp<-t(scale(t(exp)))

exp<-subset(exp, rownames(exp)%in%SL_genes$hgnc_symbol)

cl<-makeCluster(2)
registerDoParallel(cl)
results<-data.frame("GeneA"=NA, "GeneB"=NA, "Pathway"=NA,
                    "P_value"=NA, "Distance"=NA, "Horizontal_Max"=NA, "Position"=NA, stringsAsFactors = F)
start_time<-Sys.time()
results<-foreach(j = 1:nrow(gene_pairs), .packages = "matrixStats", .combine = rbind, .init = results)%dopar%{
  geneA<-as.numeric(subset(exp, rownames(exp)==gene_pairs$GeneA[j]))
  geneB<-as.numeric(subset(exp, rownames(exp)==gene_pairs$GeneB[j]))
  
  if(all(is.na(geneA))| all(is.na(geneB))){
    c(as.character(gene_pairs[j, 1:3]), rep(NA, 4))
  }else{
    c(as.character(gene_pairs[j, 1:3]), SL_KS_test(geneA, geneB, gene_pairs$GeneA[j], gene_pairs$GeneB[j]))
  }
  
}
stopCluster(cl)
end_time<-Sys.time()
print(end_time - start_time)
gc()

results<-results[-1,]
results$P_value<-as.numeric(results$P_value)
results$Distance<-as.numeric(results$Distance)
results$Horizontal_Max<-as.numeric(results$Horizontal_Max)
results$Position<-as.numeric(results$Position)

write.table(results, "All_cell_lines_all_gene_pair_results.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
significant_results<-subset(results, results$FDR<=0.25 &
                              results$Horizontal_Max>=0.8 &
                              results$Position<=10)
write.table(significant_results, "All_cell_lines_significant_gene_pair_results.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)
