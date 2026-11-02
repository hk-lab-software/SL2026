#====================================================================================================================================
# Fig1_Collecting_tumor_and_normal_gene_expression_data
#====================================================================================================================================
# Set working directory
setwd("Fig1")

# Package 'TCGAbiolinks' and 'SummarizedExperiment'
#library(devtools)
#devtools::install_github('BioinformaticsFMRP/TCGAbiolinks')
library(TCGAbiolinks)
library(SummarizedExperiment)

project<-TCGAbiolinks::getGDCprojects()
TCGA_project<-project[grep("TCGA", project$project_id),]
TCGA_project<-as.data.frame(TCGA_project[order(TCGA_project$project_id),], stringsAsFactor = F)
TCGA_project$primary_site<-as.character(TCGA_project$primary_site)
TCGA_project$disease_type<-as.character(TCGA_project$disease_type)
write.table(TCGA_project, "TCGA_projects.txt",
            sep = "\t", quote = F, col.names = T, row.names = F)

# Cancer types that have primary tissues in GTEx normal data.
# Reference file: TCGA_projects_compared_with_GTEx.xlsx
TCGA_project_id<-c("TCGA-ACC", "TCGA-BLCA", "TCGA-BRCA", "TCGA-CESC", "TCGA-COAD",
                   "TCGA-ESCA", "TCGA-GBM", "TCGA-LGG", "TCGA-LIHC", "TCGA-LUAD",
                   "TCGA-LUSC", "TCGA-OV", "TCGA-PAAD", "TCGA-PRAD", "TCGA-SKCM",
                   "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA", "TCGA-UCEC", "TCGA-UCS")
#===========================================#
# FUNCTION
# Prepare and download data from GDC portal
#===========================================#
download_TCGA_exp<-function(project_id){
  query<-GDCquery(project = project_id,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM",
                  experimental.strategy = "RNA-Seq")
  GDCdownload(query)
  exp<-GDCprepare(query, save = T, save.filename = paste0("RData/", project_id, ".RData"))
  exp_1<-assay(exp)
  rownames(exp_1)<-exp@rowRanges@elementMetadata$external_gene_name
  write.table(exp_1, paste0("Files/", project_id, ".txt"), sep = "\t", quote = F, row.names = T, col.names = T)
}
#===============================================#
# FUNCTION
# Error handling function for download_TCGA_exp
#===============================================#
download_TCGA_exp_error_handling<-function(project_id){
  tryCatch(download_TCGA_exp(project_id),
           error = function(project_id) download_TCGA_exp(project_id))
}

for(i in 1:length(TCGA_project_id)){
  print(TCGA_project_id[i])
  download_TCGA_exp_error_handling(TCGA_project_id[i])
}

#===============================================#
# Import GTEx normal tissue data
#===============================================#
TCGA_GTEx_tb<-read.delim("TCGA-GTEx.txt", sep = "\t", header = T, stringsAsFactors = F)
# Remove HNSC since it was prepared in another script
TCGA_GTEx_tb<-subset(TCGA_GTEx_tb, TCGA_GTEx_tb$TCGA!="TCGA-HNSC")
# Compare TCGA and GTEx gene symbol
# TCGA
tmp_query<-GDCquery(project = "TCGA-ACC",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "HTSeq - FPKM",
                experimental.strategy = "RNA-Seq")
exp<-GDCprepare(tmp_query, save = F)
TCGA_ENSEMBL_gene_ID<-data.frame("ENSEMBL_ID"=exp@rowRanges@elementMetadata$ensembl_gene_id,
                                "Gene_Name"=exp@rowRanges@elementMetadata$external_gene_name,
                                stringsAsFactors = F)
TCGA_ENSEMBL_gene_ID_unique<-unique(TCGA_ENSEMBL_gene_ID)
rm(tmp_query, exp)
# GTEx
GTEx_ENSEMBL_gene_ID<-read.delim("E:/KHH/HSK Lab/Data/GTEx/V7/GTEx_v7_gene_ensembl_tb.txt",
                            sep = "\t", header = T, stringsAsFactors = F)
GTEx_ENSEMBL_gene_ID$ENSEMBL_ID<-sapply(1:nrow(GTEx_ENSEMBL_gene_ID),
                                   function(x) strsplit(GTEx_ENSEMBL_gene_ID$ENSEMBL_ID[x], split = ".", fixed = T)[[1]][1])
GTEx_ENSEMBL_gene_ID<-GTEx_ENSEMBL_gene_ID[order(GTEx_ENSEMBL_gene_ID$ENSEMBL_ID),]
GTEx_ENSEMBL_gene_ID<-unique(GTEx_ENSEMBL_gene_ID)

ENSEMBL_gene_ID<-merge(TCGA_ENSEMBL_gene_ID_unique, GTEx_ENSEMBL_gene_ID, by = "ENSEMBL_ID")
colnames(ENSEMBL_gene_ID)<-c("ENSEMBL_ID", "TCGA_gene_name", "GTEx_gene_name")
write.table(ENSEMBL_gene_ID, "E:/KHH/OneDrive/HSK Lab/Synthetic_Lethality/Tumor_type/Data/TCGA_GTEx_ENSEMBL_ID.txt", sep = "\t", quote = F, col.names = T, row.names = F)
ENSEMBL_gene_ID<-read.delim("E:/KHH/OneDrive/HSK Lab/Synthetic_Lethality/Tumor_type/Data/TCGA_GTEx_ENSEMBL_ID.txt",
                            sep = "\t", header = T, stringsAsFactors = F)

# Log2(TPM + 1)
# Merge TCGA and GTEx data and perform quantile normalization
# Z-score
library(data.table)
#setDTthreads()
library(preprocessCore)
for(i in 1:nrow(TCGA_GTEx_tb)){
  TCGA<-TCGA_GTEx_tb$TCGA[i]
  GTEx<-strsplit(TCGA_GTEx_tb$GTEx[i], split = " | ", fixed = T)[[1]]
  
  # TCGA
  TCGA_exp<-fread(paste0("Files/", TCGA, ".txt"), data.table = T, stringsAsFactors = F)
  TCGA_exp$V1<-TCGA_ENSEMBL_gene_ID$ENSEMBL_ID
  TCGA_exp<-as.data.frame(TCGA_exp[,lapply(.SD, mean), by = "V1"], stringsAsFactors = F)
  rownames(TCGA_exp)<-TCGA_exp$V1
  TCGA_exp<-TCGA_exp[,-1]
  TCGA_exp<-TCGA_exp[ENSEMBL_gene_ID$ENSEMBL_ID,]
  
  # Remove normal derived samples
  # 11: Solid Tissue Normal
  # 12: Buccal Cell Normal
  # 13: EBV Immortalized Normal
  # 14: Bone Marrow Normal
  tmp<-c()
  for(k in 1:ncol(TCGA_exp)){
    tmp1<-strsplit(colnames(TCGA_exp)[k], split = "-", fixed = T)[[1]][4]
    tmp1<-c(grep("11", tmp1),
            grep("12", tmp1),
            grep("13", tmp1),
            grep("14", tmp1))
    
    if(length(tmp1)!=0){
      tmp<-c(tmp, k)
    }
  }
  
  if(length(tmp)!=0){
    TCGA_exp<-TCGA_exp[,-tmp]
  }
  rm(tmp, tmp1)
  TCGA_samples<-colnames(TCGA_exp)
  TCGA_exp<-apply(TCGA_exp, 2, function(x) x / sum(as.numeric(x)) * 10^6)
  
  # GTEx
  GTEx_exp<-matrix(NA, nrow = 56202)
  for(j in 1:length(GTEx)){
    tmp<-fread(paste0("E:/KHH/HSK Lab/Data/GTEx/V7/By_tissue/", GTEx[j], ".txt"), data.table = F, stringsAsFactors = F)
    tmp$V1<-sapply(1:nrow(tmp), function(x) strsplit(tmp$V1[x], split = ".", fixed = T)[[1]][1])
    rownames(tmp)<-tmp$V1
    tmp<-tmp[,-1]
    
    GTEx_exp<-cbind(GTEx_exp, tmp)
  }
  GTEx_exp<-GTEx_exp[,-1]
  GTEx_exp<-GTEx_exp[ENSEMBL_gene_ID$ENSEMBL_ID,]
  GTEx_samples<-colnames(GTEx_exp)
  
  # Merge TCGA and GTEx
  exp_tb<-cbind(TCGA_exp, GTEx_exp)
  exp_tb$Gene<-ENSEMBL_gene_ID$TCGA_gene_name
  
  # Get mean of duplicated genes
  exp_tb<-as.data.frame(as.data.table(exp_tb)[,lapply(.SD, mean), by = "Gene"], stringsAsFactors = F)
  rownames(exp_tb)<-exp_tb$Gene
  exp_tb<-exp_tb[,-1]
  exp_tb<-log2(exp_tb + 1)
  
  # Quantile normalization
  exp_quantile_normalized_tb<-normalize.quantiles(as.matrix(exp_tb))
  rownames(exp_quantile_normalized_tb)<-rownames(exp_tb)
  colnames(exp_quantile_normalized_tb)<-c(TCGA_samples, GTEx_samples)
  
  # Z-score
  TCGA_exp_quantile_normalized<-exp_quantile_normalized_tb[,TCGA_samples]
  TCGA_exp_z_score<-t(scale(t(TCGA_exp_quantile_normalized), center = T, scale = T))
  
  GTEx_exp_quantile_normalized<-exp_quantile_normalized_tb[,GTEx_samples]
  GTEx_exp_z_score<-t(scale(t(GTEx_exp_quantile_normalized), center = T, scale = T))
  
  # Save TCGA, GTEx, and merged gene expression data
  TCGA_exp_quantile_normalized<-as.data.frame(TCGA_exp_quantile_normalized)
  GTEx_exp_quantile_normalized<-as.data.frame(GTEx_exp_quantile_normalized)
  
  TCGA_exp_z_score<-as.data.frame(TCGA_exp_z_score)
  GTEx_exp_z_score<-as.data.frame(GTEx_exp_z_score)
  
  exp_quantile_normalized_tb<-as.data.frame(exp_quantile_normalized_tb)
  
  fwrite(TCGA_exp_quantile_normalized, paste0("Quantile_normalized/", TCGA, ".txt"),
              sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(GTEx_exp_quantile_normalized,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/GTEx_gene_expression/Quantile_normalized/", TCGA, "_normal.txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(TCGA_exp_z_score, paste0("Z-score_normalized/", TCGA, ".txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(GTEx_exp_z_score,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/GTEx_gene_expression/Z-score_normalized/", TCGA, "_normal.txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(exp_quantile_normalized_tb,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/TCGA_and_GTEx_merged_gene_expression/", TCGA, ".txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  
  print(round(i/nrow(TCGA_GTEx_tb)*100, digits = 2))
}

#=========================================================#
# Make another data with clinically relevant tumor types
# Gland origin: ACC, BRCA, PRAD, THCA
# Gynecological: CESC, UCEC, UCS
# Brain: GBM, LGG
# GI tract: COAD, ESCA, PAAD, STAD
# Germ cell: OV, TGCT
#=========================================================#
cancer_groups<-list("Gland_origin"=c("TCGA-ACC", "TCGA-BRCA", "TCGA-PRAD", "TCGA-THCA"),
                    "Gynecological"=c("TCGA-CESC", "TCGA-UCEC", "TCGA-UCS"),
                    "Brain"=c("TCGA-GBM", "TCGA-LGG"),
                    "GI_tract"=c("TCGA-COAD", "TCGA-ESCA", "TCGA-PAAD", "TCGA-STAD"),
                    "Germ_cell"=c("TCGA-OV", "TCGA-TGCT"))
for(i in 1:length(cancer_groups)){
  TCGA<-cancer_groups[[i]]
  GTEx<-c()
  for(j in 1:length(TCGA)){
    GTEx<-c(GTEx, strsplit(subset(TCGA_GTEx_tb, TCGA_GTEx_tb$TCGA==TCGA[j])$GTEx, split = " | ", fixed = T)[[1]])
  }
  
  # TCGA
  TCGA_exp<-c()
  for(j in 1:length(TCGA)){
    tmp<-fread(paste0("Files/", TCGA[j], ".txt"), data.table = T, stringsAsFactors = F)
    tmp$V1<-TCGA_ENSEMBL_gene_ID$ENSEMBL_ID
    tmp<-as.data.frame(tmp[, lapply(.SD, mean), by = "V1"], stringsAsFactors = F)
    rownames(tmp)<-tmp$V1
    tmp<-tmp[,-1]
    tmp<-tmp[ENSEMBL_gene_ID$ENSEMBL_ID,]
    
    if(j==1){
      TCGA_exp<-tmp
    }else{
      TCGA_exp<-cbind(TCGA_exp, tmp)
    }
  }
  
  # Remove normal derived samples
  # 11: Solid Tissue Normal
  # 12: Buccal Cell Normal
  # 13: EBV Immortalized Normal
  # 14: Bone Marrow Normal
  tmp<-c()
  for(k in 1:ncol(TCGA_exp)){
    tmp1<-strsplit(colnames(TCGA_exp)[k], split = "-", fixed = T)[[1]][4]
    tmp1<-c(grep("11", tmp1),
            grep("12", tmp1),
            grep("13", tmp1),
            grep("14", tmp1))
    
    if(length(tmp1)!=0){
      tmp<-c(tmp, k)
    }
  }
  
  if(length(tmp)!=0){
    TCGA_exp<-TCGA_exp[,-tmp]
  }
  rm(tmp, tmp1)
  
  TCGA_samples<-colnames(TCGA_exp)
  
  # GTEx
  for(j in 1:length(GTEx)){
    tmp<-fread(paste0("E:/KHH/HSK Lab/Data/GTEx/V7/By_tissue/", GTEx[j], ".txt"), data.table = F, stringsAsFactors = F)
    tmp$V1<-sapply(1:nrow(tmp), function(x) strsplit(tmp$V1[x], split = ".", fixed = T)[[1]][1])
    rownames(tmp)<-tmp$V1
    tmp<-tmp[,-1]
    
    if(j==1){
      GTEx_exp<-tmp
    }else{
      GTEx_exp<-cbind(GTEx_exp, tmp)
    }
  }
  GTEx_exp<-GTEx_exp[,-1]
  GTEx_exp<-GTEx_exp[ENSEMBL_gene_ID$ENSEMBL_ID,]
  GTEx_samples<-colnames(GTEx_exp)
  
  # Merge TCGA and GTEx
  exp_tb<-cbind(TCGA_exp, GTEx_exp)
  exp_tb$Gene<-ENSEMBL_gene_ID$TCGA_gene_name
  
  # Get mean of duplicated genes
  exp_tb<-as.data.frame(as.data.table(exp_tb)[,lapply(.SD, mean), by = "Gene"], stringsAsFactors = F)
  rownames(exp_tb)<-exp_tb$Gene
  exp_tb<-exp_tb[,-1]
  exp_tb<-log2(exp_tb + 1)
  
  # Quantile normalization
  exp_quantile_normalized_tb<-normalize.quantiles(as.matrix(exp_tb))
  rownames(exp_quantile_normalized_tb)<-rownames(exp_tb)
  colnames(exp_quantile_normalized_tb)<-c(TCGA_samples, GTEx_samples)
  
  # Z-score
  TCGA_exp_quantile_normalized<-exp_quantile_normalized_tb[,TCGA_samples]
  TCGA_exp_z_score<-t(scale(t(TCGA_exp_quantile_normalized), center = T, scale = T))
  
  GTEx_exp_quantile_normalized<-exp_quantile_normalized_tb[,GTEx_samples]
  GTEx_exp_z_score<-t(scale(t(GTEx_exp_quantile_normalized), center = T, scale = T))
  
  # Save TCGA, GTEx, and merged gene expression data
  TCGA_exp_quantile_normalized<-as.data.frame(TCGA_exp_quantile_normalized)
  GTEx_exp_quantile_normalized<-as.data.frame(GTEx_exp_quantile_normalized)
  
  TCGA_exp_z_score<-as.data.frame(TCGA_exp_z_score)
  GTEx_exp_z_score<-as.data.frame(GTEx_exp_z_score)
  exp_quantile_normalized_tb<-as.data.frame(exp_quantile_normalized_tb)
  
  fwrite(TCGA_exp_quantile_normalized, paste0("Quantile_normalized/", names(cancer_groups)[i], ".txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(GTEx_exp_quantile_normalized,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/GTEx_gene_expression/Quantile_normalized/", names(cancer_groups)[i], "_normal.txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(TCGA_exp_z_score, paste0("Z-score_normalized/", names(cancer_groups)[i], ".txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(GTEx_exp_z_score,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/GTEx_gene_expression/Z-score_normalized/", names(cancer_groups)[i], "_normal.txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
  fwrite(exp_quantile_normalized_tb,
         paste0("E:/KHH/HSK Lab/Synthetic_Lethality/Tumor_type/Data/TCGA_and_GTEx_merged_gene_expression/", names(cancer_groups)[i], ".txt"),
         sep = "\t", quote = F, row.names = T, col.names = T)
}
