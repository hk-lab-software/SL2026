#================================
# Fig1a
# 1.Bimodally_expressed_genes
#================================
# Set working directory
setwd("Fig1")

# Load data
gene_pairs<-read.delim("Gene_pairs_all_genes_paired_in_pathway", sep = '\t', header = T, stringsAsFactors = F)
genes<-unique(c(gene_pairs$GeneA, gene_pairs$GeneB))

#=============================================================
# Use Mclust package to fit a Gaussian finite mixture model
# Fit by EM algorithm
#=============================================================
library(mclust)
set.seed(1111)
#======================#
# FUNCTIONS
# Perform GMM fitting
# exp: expression vector
# gene: name of gene
#======================#
mclust<-function(exp, gene){
  gmm.1<-densityMclust(exp, G = 1)
  gmm.2<-densityMclust(exp, G = 2)
  
  gmm_res.1<-summary(gmm.1)
  gmm_res.2<-summary(gmm.2)
  
  density.1<-gmm.1$density
  density.2<-gmm.2$density
  
  tmp<-data.frame("x"=exp, "y.1"=density.1, "y.2"=density.2)
  tmp<-tmp[order(tmp$x),]
  
  x<-seq(min(exp)-diff(range(exp))/10, max(exp)+diff(range(exp))/10, length = 200)
  cdens<-predict(gmm.2, x, what = "cdens")
  cdnes<-t(apply(cdens, 1, function(d) d*gmm.2$parameters$pro))
  
  hist(exp, main = gene, freq = F, xlab = gene, col = 'grey', border = F, cex.axis = 1.6, ylim = c(0, max(cdens)))
  lines(tmp$x, tmp$y.1, lwd = 2, col = "seagreen4")
  matplot(x, cdens, type = "l", lwd = 3, add = T, lty = 1, col = c("firebrick", "royalblue"))
  loglikratio<-gmm_res.2$loglik - gmm_res.1$loglik
  pro<-gmm_res.2$pro
  legend("topright", legend = c("1 component", "Group 1", "Group 2", round(c(loglikratio, pro[1], pro[2]), digits = 4)),
        col = c("seagreen", "firebrick", "royalblue", "grey", "grey", "grey"), lty = c(1,1), lwd = c(2,2), cex = 0.7)

  return(c(loglikratio, pro, gmm_res.2$mean, gmm_res.2$variance))
}

#===========================================================#
# Perform GMM for each tumor type and each cancer groups
# 1. TCGA cancer data
# 2. GTEx normal data
#===========================================================#
library(data.table)
setDTthreads(3)
cancer_types<-list.files("Data/TCGA_gene_expression/Z-score_normalized/")
cancer_types<-sapply(1:length(cancer_types), function(x) strsplit(cancer_types[x], split = ".", fixed = T)[[1]][1])

# Clinically relevant cancer groups
cancer_group<-cancer_types[-grep("TCGA", cancer_types)]
cancer_types<-cancer_types[grep("TCGA", cancer_types)]

#===============#
# cancer types
#================#
for(i in 1:length(cancer_types)){
  print(cancer_types[i])
  
  # Make result table
  mclust_results<-matrix(NA, ncol = 16, nrow = length(genes))
  colnames(mclust_results)<-c("TCGA_Loglikelihoodratio", "TCGA_Group1", "TCGA_Group2", "TCGA_Group1_mean", "TCGA_Group2_mean",
                              "TCGA_Group1_variance", "TCGA_Group2_variance", "TCGA_mean",
                              "GTEx_Loglikelihoodratio", "GTEx_Group1", "GTEx_Group2", "GTEx_Group1_mean", "GTEx_Group2_mean",
                              "GTEx_Group1_variance", "GTEx_Group2_variance", "GTEx_mean")
  rownames(mclust_results)<-genes
  
  # Read TCGA and GTEx gene expression data
  TCGA_exp<-fread(paste0("Data/TCGA_gene_expression/Z-score_normalized/", cancer_types[i], ".txt"),
                  data.table = F, stringsAsFactors = F)
  GTEx_exp<-fread(paste0("Data/GTEx_gene_expression/Z-score_normalized/", cancer_types[i], "_normal.txt"),
                  data.table = F, stringsAsFactors = F)
  rownames(TCGA_exp)<-TCGA_exp$V1;rownames(GTEx_exp)<-GTEx_exp$V1
  TCGA_exp<-TCGA_exp[,-1];GTEx_exp<-GTEx_exp[,-1]
  
  # Create figure directory for each cancer type
  dir.create(paste0("Figures/GMM/", cancer_types[i]))
  # Perform GMM for each gene
  for(j in 1:length(genes)){
    TCGA<-as.numeric(TCGA_exp[genes[j],])
    GTEx<-as.numeric(GTEx_exp[genes[j],])
    
    if(all(is.na(TCGA)) | all(is.na(GTEx))){
      mclust_results[j,]<-NA
    }else{
      pdf(paste0("Figures/GMM/", cancer_types[i], "/", genes[j], ".pdf"))
      par(mfrow = c(1,2))
      
      mclust_results[j,8]<-mean(TCGA)
      mclust_results[j,16]<-mean(GTEx)
      
      tryCatch({
        mclust_results[j, 1:7]<-mclust(TCGA, paste0("TCGA_", genes[j]))
      }, error = function(e) print(e))
      tryCatch({
        mclust_results[j, 9:15]<-mclust(GTEx, paste0("GTEx_", genes[j]))
      }, error = function(e) print(e))
      
      par(mfrow = c(1,1))
      tmp_TCGA<-hist(TCGA, breaks = 25, plot = F)
      tmp_GTEx<-hist(GTEx, breaks = 25, plot = F)
      
      hist(TCGA, main = genes[j], freq = F, xlab = "Gene expression",
           col = rgb(178, 34, 34, max = 255, alpha = 125), border = F, cex.axis = 1.6, breaks = 25,
           xlim = range(c(TCGA, GTEx)), ylim = range(c(tmp_TCGA$density, tmp_GTEx$density)))
      hist(GTEx, add = T, freq = F, col = rgb(65, 105, 225, max = 255, alpha = 125), border = F, cex.axis = 1.6, breaks = 25)
      legend("topright", legend = c("Cancer(TCGA)", "Normal(GTEx)"),
             col = c("firebrick", "royalblue"), pch = 19, cex = 0.7)
      dev.off()
    }
  }
  write.table(mclust_results, paste0("Results/GMM/Cancer_type/", cancer_types[i], ".txt"), sep = "\t", quote = F, col.names = T, row.names = T)
}
#================#
# 5 cancer groups
#================#
for(i in 1:length(cancer_group)){
  print(cancer_group[i])
  
  # Make result table
  mclust_results<-matrix(NA, ncol = 16, nrow = length(genes))
  colnames(mclust_results)<-c("TCGA_Loglikelihoodratio", "TCGA_Group1", "TCGA_Group2", "TCGA_Group1_mean", "TCGA_Group2_mean",
                              "TCGA_Group1_variance", "TCGA_Group2_variance", "TCGA_mean",
                              "GTEx_Loglikelihoodratio", "GTEx_Group1", "GTEx_Group2", "GTEx_Group1_mean", "GTEx_Group2_mean",
                              "GTEx_Group1_variance", "GTEx_Group2_variance", "GTEx_mean")
  rownames(mclust_results)<-genes
  
  # Read TCGA and GTEx gene expression data
  TCGA_exp<-fread(paste0("Data/TCGA_gene_expression/Z-score_normalized/", cancer_group[i], ".txt"),
                  data.table = F, stringsAsFactors = F)
  GTEx_exp<-fread(paste0("Data/GTEx_gene_expression/Z-score_normalized/", cancer_group[i], "_normal.txt"),
                  data.table = F, stringsAsFactors = F)
  rownames(TCGA_exp)<-TCGA_exp$V1;rownames(GTEx_exp)<-GTEx_exp$V1
  TCGA_exp<-TCGA_exp[,-1];GTEx_exp<-GTEx_exp[,-1]
  
  # Create figure directory for each cancer type
  dir.create(paste0("Figures/GMM/", cancer_group[i]))
  # Perform GMM for each gene
  for(j in 1:length(genes)){
    TCGA<-as.numeric(TCGA_exp[genes[j],])
    GTEx<-as.numeric(GTEx_exp[genes[j],])
    
    if(all(is.na(TCGA)) | all(is.na(GTEx))){
      mclust_results[j,]<-NA
    }else{
      pdf(paste0("Figures/GMM/", cancer_group[i], "/", genes[j], ".pdf"))
      par(mfrow = c(1,2))
      
      mclust_results[j,8]<-mean(TCGA)
      mclust_results[j,16]<-mean(GTEx)
      
      tryCatch({
        mclust_results[j, 1:7]<-mclust(TCGA, paste0("TCGA_", genes[j]))
      }, error = function(e) print(e))
      tryCatch({
        mclust_results[j, 9:15]<-mclust(GTEx, paste0("GTEx_", genes[j]))
      }, error = function(e) print(e))
      
      par(mfrow = c(1,1))
      tmp_TCGA<-hist(TCGA, breaks = 25, plot = F)
      tmp_GTEx<-hist(GTEx, breaks = 25, plot = F)
      
      hist(TCGA, main = genes[j], freq = F, xlab = "Gene expression",
           col = rgb(178, 34, 34, max = 255, alpha = 125), border = F, cex.axis = 1.6, breaks = 25,
           xlim = range(c(TCGA, GTEx)), ylim = range(c(tmp_TCGA$density, tmp_GTEx$density)))
      hist(GTEx, add = T, freq = F, col = rgb(65, 105, 225, max = 255, alpha = 125), border = F, cex.axis = 1.6, breaks = 25)
      legend("topright", legend = c("Cancer(TCGA)", "Normal(GTEx)"),
             col = c("firebrick", "royalblue"), pch = 19, cex = 0.7)
      dev.off()
    }
  }
  write.table(mclust_results, paste0("Results/GMM/Cancer_type/", cancer_group[i], ".txt"), sep = "\t", quote = F, col.names = T, row.names = T)
}
#===============================================================#
# Interpret results
# Get two component GMM following genes with the below criteria
# Log-likelihood ratio>=3
# Mean expression difference between low and high group>=1.5
# Proportion of low expressing group is between 2%~50%
#===============================================================#
GMM_results<-c()
cancer_types<-c(cancer_types, cancer_group)
for(i in 1:length(cancer_types)){
  print(cancer_types[i])
  mclust_results<-read.delim(paste0("Results/GMM/Cancer_type/", cancer_types[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
  mclust_results$Cancer_type<-cancer_types[i]
  mclust_results$Genes<-rownames(mclust_results)
  
  GMM_results<-rbind(GMM_results, mclust_results)
}

GMM_following_genes<-c()
for(i in 1:length(cancer_types)){
  print(cancer_types[i])
  mclust_results<-read.delim(paste0("Results/GMM/Cancer_type/", cancer_types[i], ".txt"), sep = "\t", header = T, stringsAsFactors = F)
  tmp<-subset(mclust_results, mclust_results$TCGA_Loglikelihoodratio>3 &
                (mclust_results$TCGA_Group1>0.02 & mclust_results$TCGA_Group1<0.5) &
                (mclust_results$TCGA_Group2_mean - mclust_results$TCGA_Group1_mean)>=1.5)
  if(nrow(tmp)>=1){
    tmp$Cancer_type<-cancer_types[i]
    tmp$Genes<-rownames(tmp)
    GMM_following_genes<-rbind(GMM_following_genes, tmp)
  }
}
rownames(GMM_following_genes)<-NULL

#===============================================================#
# Compare two component GMM following genes with GTEx normal data
# Calculate difference of mean of TCGA high expressing group and GTEx expression
# GTEx normal tissue data should not follow a two component GMM
#===============================================================#
GMM_following_genes_normal<-subset(GMM_following_genes, !(GMM_following_genes$GTEx_Loglikelihoodratio>3) |
                                     !(GMM_following_genes$GTEx_Group1>0.02 | GMM_following_genes$GTEx_Group1<0.5) |
                                     !((GMM_following_genes$GTEx_Group2_mean - GMM_following_genes$GTEx_Group1_mean)>=1.5))
tmp<-c()
for(i in 1:length(cancer_types)){
  print(cancer_types[i])
  
  tmp1<-subset(GMM_following_genes_normal, GMM_following_genes_normal$Cancer_type==cancer_types[i])
  tmp1$t.test_pval<-NA
  tmp1$KS.test_pval<-NA
  tmp1$Mean_GTEx_TCGA_high<-NA
  
  # TCGA and GTEx gene expression data
  # Read TCGA and GTEx gene expression data
  TCGA_exp<-fread(paste0("Data/TCGA_gene_expression/Z-score_normalized/", cancer_types[i], ".txt"),
                  data.table = F, stringsAsFactors = F)
  GTEx_exp<-fread(paste0("Data/GTEx_gene_expression/Z-score_normalized/", cancer_types[i], "_normal.txt"),
                  data.table = F, stringsAsFactors = F)
  rownames(TCGA_exp)<-TCGA_exp$V1;rownames(GTEx_exp)<-GTEx_exp$V1
  TCGA_exp<-TCGA_exp[,-1];GTEx_exp<-GTEx_exp[,-1]
  
  for(j in 1:nrow(tmp1)){
    TCGA<-as.numeric(TCGA_exp[tmp1$Genes[j],])
    GTEx<-as.numeric(GTEx_exp[tmp1$Genes[j],])
    
    gmm.2<-densityMclust(TCGA, G = 2)
    gmm_res.2<-summary(gmm.2)
    TCGA_high_group<-TCGA[gmm_res.2$classification==2]
    tmp1$t.test_pval[j]<-t.test(TCGA_high_group, GTEx)$'p.value'
    tmp1$KS.test_pval[j]<-ks.test(TCGA_high_group, GTEx)$'p.value'
    tmp1$Mean_GTEx_TCGA_high[j]<-mean(GTEx)-mean(TCGA_high_group)
  }

  tmp<-rbind(tmp, tmp1)
}
GMM_following_genes_normal<-tmp;rm(tmp, tmp1)

# Remove genes that have high groups that differ from normal expression
GMM_following_genes_normal_filtered<-subset(GMM_following_genes_normal, GMM_following_genes_normal$t.test_pval>=0.05)
GMM_following_genes_normal_not_filtered<-subset(GMM_following_genes_normal, GMM_following_genes_normal$t.test_pval<0.05)

# Save table
write.table(GMM_results, "Results/GMM/All_GMM_results_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_following_genes, "Results/GMM/Two_component_GMM_following_genes_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_following_genes_in_each_cancer_type, "Results/GMM/Two_component_GMM_following_genes_in_each_cancer_type_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)

write.table(GMM_following_genes_normal, "Results/GMM/Two_component_GMM_following_genes_compared_with_normal_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_following_genes_normal_filtered, "Results/GMM/Two_component_GMM_following_genes_compared_with_normal_filtered_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_following_genes_normal_not_filtered, "Results/GMM/Two_component_GMM_following_genes_compared_with_normal_not_filtered_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_following_genes_in_each_cancer_type_normal, "Results/GMM/Two_component_GMM_following_genes_in_each_cancer_type_normal_filtered_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
write.table(GMM_genes_normal, "Results/GMM/Number_of_cancer_types_by_each_two_component_GMM_following_genes_normal_filtered_190617.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)
