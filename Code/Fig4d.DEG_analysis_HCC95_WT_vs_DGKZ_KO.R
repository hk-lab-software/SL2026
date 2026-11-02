#============================================
# Fig 4d
# DEG_analysis_HCC95_WT_vs_DGKZ_KO_241002
#============================================
# Set working directory
setwd("Fig4d")

# Read count table
library(data.table)
library(matrixStats)
count_tb<-fread("HCC95_WT_vs_DGKZ_KO_count.txt",
                sep = "\t", header = T, data.table = F)

# Convert ENSG ID to gene symbol
library(biomaRt)
ensembl<-useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
tmp<-listAttributes(ensembl)
genes_tb<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                   filters="ensembl_gene_id",
                   values=count_tb$ENSG_ID, mart=ensembl)
count_tb$Gene_symbol<-NA
idx<-c()
for(i in 1:nrow(count_tb)){
  if(i%%1000==0){
    print(i)
  }
  if(count_tb$ENSG_ID[i]%in%genes_tb$ensembl_gene_id){
    tmp<-subset(genes_tb, genes_tb$ensembl_gene_id==count_tb$ENSG_ID[i])
    tmp<-tmp$hgnc_symbol
    if(""%in%tmp){
      count_tb$Gene_symbol[i]<-count_tb$ENSG_ID[i]
    }else if(length(tmp)>1){
      idx<-c(idx, i)
    }else{
      count_tb$Gene_symbol[i]<-tmp
    }
  }else{
    count_tb$Gene_symbol[i]<-count_tb$ENSG_ID[i]
  }
  
}

# Manually curate genes
genes_tb[which(genes_tb$ensembl_gene_id%in%count_tb$ENSG_ID[idx]),]


count_tb<-count_tb[order(count_tb$Gene_symbol),]
count_tb<-cbind(count_tb[,8], count_tb[,1:7])
colnames(count_tb)[1:2]<-c("Gene_symbol", "ENSG_ID")

# Duplicate gene symbols
tmp<-count_tb$Gene_symbol[which(duplicated(count_tb$Gene_symbol))]

tmp1<-subset(count_tb, !count_tb$Gene_symbol%in%tmp)
tmp2<-subset(count_tb, count_tb$Gene_symbol%in%tmp)

tmp3<-data.frame("Gene_symbol"=tmp, "ENSG_ID"=NA,
                 "WT1"=NA, "WT2"=NA, "WT3"=NA,
                 "KO1"=NA, "KO2"=NA, "KO3"=NA,
                 stringsAsFactors = F)
for(i in 1:nrow(tmp3)){
  tmp4<-subset(tmp2, tmp2$Gene_symbol==tmp3$Gene_symbol[i])
  tmp3$ENSG_ID[i]<-tmp4$ENSG_ID[1]
  tmp3[i,3:8]<-colMaxs(as.matrix(tmp4[,3:8]))
}
count_tb<-rbind(tmp1, tmp3)
count_tb<-count_tb[order(count_tb$Gene_symbol),]

write.table(count_tb, "HCC95_WT_vs_DGKZ_KO_count.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

# Read FPKM table
FPKM_tb<-fread("genes.fpkm_table",
                sep = "\t", header = T, data.table = F)
colnames(FPKM_tb)<-c("ENSG_ID", "WT1", "WT2", "WT3", "KO1", "KO2", "KO3")
genes_tb<-getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                filters="ensembl_gene_id",
                values=FPKM_tb$ENSG_ID, mart=ensembl)
FPKM_tb$Gene_symbol<-NA
idx<-c()
for(i in 1:nrow(FPKM_tb)){
  if(i%%1000==0){
    print(i)
  }
  if(FPKM_tb$ENSG_ID[i]%in%genes_tb$ensembl_gene_id){
    tmp<-subset(genes_tb, genes_tb$ensembl_gene_id==FPKM_tb$ENSG_ID[i])
    tmp<-tmp$hgnc_symbol
    if(""%in%tmp){
      FPKM_tb$Gene_symbol[i]<-FPKM_tb$ENSG_ID[i]
    }else if(length(tmp)>1){
      idx<-c(idx, i)
    }else{
      FPKM_tb$Gene_symbol[i]<-tmp
    }
  }else{
    FPKM_tb$Gene_symbol[i]<-FPKM_tb$ENSG_ID[i]
  }
  
}

# Manually curate genes
genes_tb[which(genes_tb$ensembl_gene_id%in%FPKM_tb$ENSG_ID[idx]),]

FPKM_tb$Gene_symbol[29444]<-"LINC00856"
FPKM_tb$Gene_symbol[56629]<-"EIF1B-AS1"

FPKM_tb<-FPKM_tb[order(FPKM_tb$Gene_symbol),]
FPKM_tb<-cbind(FPKM_tb[,8], FPKM_tb[,1:7])
colnames(FPKM_tb)[1:2]<-c("Gene_symbol", "ENSG_ID")
FPKM_tb[,3:8]<-log2(FPKM_tb[,3:8]+1)

write.table(FPKM_tb, "HCC95_WT_vs_DGKZ_KO_Log2_FPKM.txt",
            sep = "\t", quote = F, row.names = F, col.names = T)

boxplot(list("WT"=as.numeric(FPKM_tb[6356,3:5]), "KO"=as.numeric(FPKM_tb[6356,6:8])),
        ylim = c(2,5), ylab = "DGKZ expression (Log2(FPKM+1)",
        col = c("firebrick", "royalblue"))
tmp_test<-t.test(as.numeric(FPKM_tb[6356,3:5]), as.numeric(FPKM_tb[6356,6:8]))
text(x = 1.5, y = 4.8,
     substitute(paste(italic("P"), " = 0.0002")), cex = 2)

# DESEQ2
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# Perform DESeq2 analysis
coldata<-data.frame("Sample"=c("WT1", "WT2", "WT3", "KO1", "KO2", "KO3"), "Group"=NA)
rownames(coldata)<-coldata$Sample
coldata$Group<-c(rep("WT", 3), rep("KO", 3))
coldata$Group<-as.factor(coldata$Group)

# Check if names match
rownames(count_tb)<-count_tb$Gene_symbol
count_tb_1<-count_tb[,-c(1:2)]
all(rownames(coldata)==colnames(count_tb_1))

# Round values to integers
gene_exp_count_tb<-round(gene_exp_count_tb)

dds<-DESeqDataSetFromMatrix(countData = count_tb_1,
                            colData = coldata,
                            design = ~ Group)
dds

dds$Group<-factor(dds$Group, levels = c("WT", "KO"))

# Run Differential expression analysis
dds<-DESeq(dds)
results<-results(dds)
results
results<-results[order(results$padj),]

summary(results)

# Log fold change shrinkage for visualization and ranking
resultsNames(dds)

results_LFC<-lfcShrink(dds, coef = "Group_KO_vs_WT", type = "apeglm")
results_LFC<-results_LFC[order(results_LFC$padj),]

# MA-plot
plotMA(results_LFC)

plotCounts(dds, gene = rownames(results)[which.min(results$padj)], intgroup = "Group")

write.csv(as.data.frame(results),
          file = "DESeq2_results.csv")

# Volcano plot
results<-as.data.frame(results)

#
results<-read.csv("DESeq2_results.csv",
                  header =  T, stringsAsFactors = F)
results$FDR<-p.adjust(results$pvalue, method = "fdr")
results<-subset(results, !is.na(results$FDR))

col1<-rgb(178, 34, 34, max = 255, alpha = 255, names = "firebrick50")
col2<-rgb(65, 105, 225, max = 255, alpha = 255, names = "royalblue50")
col3<-rgb(0, 0, 0, max = 255, alpha = 255)

results$col<-sapply(1:nrow(results), function(x){
  if(results$log2FoldChange[x]>1 & results$FDR[x]<0.1){
    return(col1)
  }else if(results$log2FoldChange[x]<(-1) & results$FDR[x]<0.1){
    return(col2)
  }else{
    return(col3)
  }
})

col4<-rgb(13.3, 54.5, 13.3,max = 255, alpha = 255)
# SPRR2 genes
SPRR2_genes<-c("SPRR2A","SPRR2B", "SPRR2C", "SPRR2D", "SPRR2E", "SPRR2F", "SPRR2G")
results$col[which(results$X%in%SPRR2_genes)]<-col4

order(col1, col2, col3, col4)
tmp<-rbind(subset(results,results$col==col3), subset(results, results$col==col2))
tmp<-rbind(tmp, subset(results, results$col==col1))
tmp<-rbind(tmp ,subset(results, results$col==col4))
results<-tmp


#install.packages("plotrix")
library(plotrix)

plot(x = results$log2FoldChange, y = -log10(results$FDR), pch = 19,
     col = results$col, main = "HCC95 WT vs DGKZ KO",
     xlab = "Log2FoldChange(WT/KO)", ylab = "-log10(FDR)", cex = 2)
abline(h = -log10(0.1), lty = 2, lwd = 2, col = col1)
abline(v = 1, lty = 2, col = col1)
abline(v = -1, lty = 2, col = col1)

pdf("Results/HCC95_WT_vs_DGKZ_KO/HCC95_WT_vs_DGKZ_KO_volcano_plot.pdf", height = 6, width = 6)
par(mar = c(5.1, 5.1, 3.1, 2.1))
gap.plot(x = results$log2FoldChange,
         y = -log10(results$FDR),
         col = results$col,
         gap = c(62, 135),
         ylim = c(0, 150),
         gap.axis = "y",
         ytics = seq(10, from = 0, to = 150),
         yticlab = seq(10, from = 0, to = 150),
         las = 0,
         xlab = "Log2FoldChange(WT/KO)", ylab = "-log10(FDR)", cex = 2, pch = 19,
         main = substitute(paste("HCC95 WT vs ", italic("DGKZ"), " KO")), cex.main = 1.8, cex.lab = 1.5, cex.axis = 1.5)
abline(h = -log10(0.1), lty = 2, lwd = 2, col = col1)
abline(v = 1, lty = 2, col = col1)
abline(v = -1, lty = 2, col = col1)
dev.off()

significant_DEGs<-subset(results, abs(results$log2FoldChange)>1 &
                           results$FDR<0.1)



FPKM_tb<-read.delim("Data/RNA-seq/HCC95_sgNT_vs_sgDGKZ-1/HCC95_WT_vs_DGKZ_KO_Log2_FPKM.txt",
                    sep = "\t", header = T, stringsAsFactors = F)
Focal_Adhesion_genes<-read.delim("Data/RNA-seq/HCC95_sgNT_vs_sgDGKZ-1/KEGG_FOCAL_ADHESION.v2023.2.Hs.gmt.txt",
                                 sep = "\t", header = F, stringsAsFactors = F)
Focal_Adhesion_genes_FPKM_tb<-subset(FPKM_tb, FPKM_tb$Gene_symbol%in%Focal_Adhesion_genes$V1)
rownames(Focal_Adhesion_genes_FPKM_tb)<-Focal_Adhesion_genes_FPKM_tb$Gene_symbol
Focal_Adhesion_genes_FPKM_tb<-Focal_Adhesion_genes_FPKM_tb[,-c(1:2)]
Focal_Adhesion_genes_FPKM_tb<-t(Focal_Adhesion_genes_FPKM_tb)
library(pheatmap)
library(RColorBrewer)
pheatmap(Focal_Adhesion_genes_FPKM_tb, cluster_cols = T, cluster_rows = F, clustering_method = 'average',
         main = 'KEGG Focal Adhesion Genes',
         color = c(rep("#2166AC", 5), colorRampPalette(rev(brewer.pal(n=7, name = 'RdBu')))(50),
                   rep("#B2182B", 5)),
         fontsize_col = ,
         clustering_distance_rows = 'euclidean')

Anoikis_genes<-read.delim("Data/RNA-seq/HCC95_sgNT_vs_sgDGKZ-1/GOBP_ANOIKIS.v2023.2.Hs.gmt.txt",
                                 sep = "\t", header = F, stringsAsFactors = F)
Anoikis_genes_FPKM_tb<-subset(FPKM_tb, FPKM_tb$Gene_symbol%in%Anoikis_genes$V1)
rownames(Anoikis_genes_FPKM_tb)<-Anoikis_genes_FPKM_tb$Gene_symbol
Anoikis_genes_FPKM_tb<-Anoikis_genes_FPKM_tb[,-c(1:2)]
Anoikis_genes_FPKM_tb<-t(Anoikis_genes_FPKM_tb)
pheatmap(Anoikis_genes_FPKM_tb, cluster_cols = T, cluster_rows = F, clustering_method = 'average',
         main = 'Gene Ontology Anoikis Genes',
         color = c(rep("#2166AC", 5), colorRampPalette(rev(brewer.pal(n=7, name = 'RdBu')))(50),
                   rep("#B2182B", 5)),
         fontsize = 20,
         clustering_distance_rows = 'euclidean')





















