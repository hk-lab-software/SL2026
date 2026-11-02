#==================================================
# Fig 3c
# CCLE LUSC and HNSC cell lines
# ALDH7A1 expression and methylation correlation
#==================================================
# Set working directory
setwd("Fig3c")

# Data was downloaded from the Depmap portal.
# LUSC: Cell lines with Oncotreesubtype = Lung Squamous Cell Carcinoma
# HNSC: Cell liens with OncotreePrimaryDiseae = Head and Neck Squamous Carcinoma

# LUSC
LUSC_data<-read.csv("ALDH7A1_5_125931082_125932082 Methylation Fraction vs ALDH7A1 log2(TPM+1) filtered by LUSC.csv",
                    header = T, stringsAsFactors = F)
LUSC_data<-LUSC_data[-which(is.na(LUSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.)),]

pdf("LUSC.pdf", height = 5, width =5)
plot(x = LUSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.,
     y = LUSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "LUSC", cex = 1.5,
     xlim = c(0, 1), las = 1, ylim = c(0, 7))
abline(lm(LUSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3 ~ LUSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(LUSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.,
         LUSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3, method = "spearman")
dev.off()

# HNSC
HNSC_data<-read.csv("ALDH7A1_5_125931082_125932082 Methylation Fraction vs ALDH7A1 log2(TPM+1) filtered by HNSC.csv",
                    header = T, stringsAsFactors = F)
HNSC_data<-HNSC_data[-which(is.na(HNSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.)),]
HNSC_data<-HNSC_data[-which(is.na(HNSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3)),]
pdf("HNSC.pdf", height = 5, width =5)
plot(x = HNSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.,
     y = HNSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "HNSC", cex = 1.5,
     xlim = c(0, 1), las = 1, ylim = c(0, 8))
abline(lm(HNSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3 ~ HNSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(HNSC_data$ALDH7A1_5_125931082_125932082...Methylation.Fraction...Methylation..1kb.upstream.TSS.,
         HNSC_data$ALDH7A1...log2.TPM.1....Expression..Short.read..Public.26Q3, method = "spearman")
dev.off()

