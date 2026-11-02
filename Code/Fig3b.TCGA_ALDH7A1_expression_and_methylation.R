#=================================================
# Fig. 3b
# TCGA squamous carcinomas
# ALDH7A1 expression and methylation correlation
#=================================================
# Set working directory
setwd("Fig3b")

# Load package
library(data.table)

# LUSC
LUSC_data<-read.csv("LUSC.csv")
LUSC_data$Aggregate<-rowMeans(LUSC_data[,-c(1:2)])

pdf("LUSC.pdf", height = 5, width =5)
plot(x = LUSC_data$Aggregate, y = LUSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "LUSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(LUSC_data$ALDH7A1 ~ LUSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(LUSC_data$Aggregate, LUSC_data$ALDH7A1, method = "spearman")
dev.off()

HNSC_data<-read.csv("HNSC.csv")
HNSC_data$Aggregate<-rowMeans(HNSC_data[,-c(1:2)])

pdf("HNSC.pdf", height = 5, width =5)
plot(x = HNSC_data$Aggregate, y = HNSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "HNSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(HNSC_data$ALDH7A1 ~ HNSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(HNSC_data$Aggregate, HNSC_data$ALDH7A1, method = "spearman")
dev.off()

CESC_data<-read.csv("CESC.csv")
CESC_data$Aggregate<-rowMeans(CESC_data[,-c(1:2)])

pdf("CESC.pdf", height = 5, width =5)
plot(x = CESC_data$Aggregate, y = CESC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "CESC", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8), las = 1)
abline(lm(CESC_data$ALDH7A1 ~ CESC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(CESC_data$Aggregate, CESC_data$ALDH7A1, method = "spearman")
dev.off()

BLCA_data<-read.csv("BLCA.csv")
BLCA_data$Aggregate<-rowMeans(BLCA_data[,-c(1:2)])

pdf("BLCA.pdf", height = 5, width =5)
plot(x = BLCA_data$Aggregate, y = BLCA_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "BLCA", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8.5), las = 1)
abline(lm(BLCA_data$ALDH7A1 ~ BLCA_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(BLCA_data$Aggregate, BLCA_data$ALDH7A1, method = "spearman")
dev.off()

pdf("Merged.pdf", width = 20, height = 5)
par(mfrow = c(1,4))
plot(x = LUSC_data$Aggregate, y = LUSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "LUSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(LUSC_data$ALDH7A1 ~ LUSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = HNSC_data$Aggregate, y = HNSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "HNSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(HNSC_data$ALDH7A1 ~ HNSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = CESC_data$Aggregate, y = CESC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "CESC", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8), las = 1)
abline(lm(CESC_data$ALDH7A1 ~ CESC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = BLCA_data$Aggregate, y = BLCA_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "BLCA", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8.5), las = 1)
abline(lm(BLCA_data$ALDH7A1 ~ BLCA_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
dev.off()

# Add histology types to CESC and BLCA
# Read clinical data
CESC_clinical<-fread("CESC-clinical.tsv", data.table = F)

# CESC
unique(CESC_clinical$cases.disease_type)
CESC_clinical<-CESC_clinical[,c("cases.submitter_id","diagnoses.primary_diagnosis")]
CESC_clinical<-unique(CESC_clinical)


CESC_data$Type<-NA
for(i in 1:nrow(CESC_data)){
  tmp<-paste0(strsplit(CESC_data$Samples[i], split = "-")[[1]][1:3], collapse = "-")
  if(tmp%in%CESC_clinical$cases.submitter_id){
    tmp1<-unique(subset(CESC_clinical, CESC_clinical$cases.submitter_id==tmp)$diagnoses.primary_diagnosis)
    if(any(grepl("squamous", tmp1, ignore.case = T))){
      CESC_data$Col[i]<-"firebrick"
    }else{
      CESC_data$Col[i]<-"black"
    }
  }
}
CESC_data<-CESC_data[order(CESC_data$Col),]
pdf("CESC.pdf", height = 5, width =5)
plot(x = CESC_data$Aggregate, y = CESC_data$ALDH7A1,
     pch = 19, xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "CESC", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8), las = 1,
     col = CESC_data$Col)
abline(lm(CESC_data$ALDH7A1 ~ CESC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(CESC_data$Aggregate, CESC_data$ALDH7A1, method = "spearman")
legend("topright", col = c("firebrick", "black"), pch = 19, legend = c("Squamous", "Others"))
dev.off()

# BLCA
BLCA_clinical<-fread("BLCA-clinical.txt", sep = "\t", data.table = F)

BLCA_data$Type<-NA
BLCA_data$Col<-NA
for(i in 1:nrow(BLCA_data)){
  tmp<-paste0(strsplit(BLCA_data$Samples[i], split = "-")[[1]][1:3], collapse = "-")
  if(tmp%in%BLCA_clinical$`Case ID`){
    tmp1<-unique(subset(BLCA_clinical, BLCA_clinical$`Case ID`==tmp)$`mRNA cluster`)
    if(tmp1=="Basal_squamous"){
      BLCA_data$Col[i]<-"firebrick"
      BLCA_data$Type[i]<-tmp1
    }else{
      BLCA_data$Col[i]<-"black"
    }
  }
}
BLCA_data<-BLCA_data[order(BLCA_data$Col),]

pdf("BLCA.pdf", height = 5, width =5)
plot(x = BLCA_data$Aggregate, y = BLCA_data$ALDH7A1,
     pch = 19, xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "BLCA", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8), las = 1,
     col = BLCA_data$Col)
abline(lm(BLCA_data$ALDH7A1 ~ BLCA_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(BLCA_data$Aggregate, BLCA_data$ALDH7A1, method = "spearman")
legend("topright", col = c("firebrick", "black"), pch = 19, legend = c("Squamous", "Others"))
dev.off()


#=========================================#
# TCGA-ESCA
#=========================================#
# Gene expression
ESCA_exp<-fread("ESCA-expression.cct",
                data.table = F)
ESCA_exp<-subset(ESCA_exp, ESCA_exp$attrib_name=="ALDH7A1")
rownames(ESCA_exp)<-"ALDH7A1"
ESCA_exp<-ESCA_exp[,-1]

# Gene methylation
ESCA_methylation<-fread("ESCA-methylation.cct",
                        data.table = F)
ESCA_methylation<-subset(ESCA_methylation, ESCA_methylation$attrib_name=="ALDH7A1")
rownames(ESCA_methylation)<-"ALDH7A1"
ESCA_methylation<-ESCA_methylation[,-1]

ESCA_samples<-intersect(colnames(ESCA_exp), colnames(ESCA_methylation))


ESCA_data<-data.frame("Samples"=ESCA_samples, "ALDH7A1"=NA, "Aggregate"=NA)
ESCA_data$Aggregate<-as.numeric(ESCA_methylation[,ESCA_data$Samples])
ESCA_data$ALDH7A1<-as.numeric(ESCA_exp[,ESCA_data$Samples])
ESCA_data$Samples<-sapply(1:length(ESCA_data$Samples), function(x) gsub(".", "-", ESCA_data$Samples[x], fixed = T))
ESCA_data$Aggregate<-ESCA_data$Aggregate+0.5

pdf("ESCA.pdf", height = 5, width =5)
plot(x = ESCA_data$Aggregate, y = ESCA_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "ESCA", cex = 1.5,
     xlim = c(0, 1), ylim = c(4, 12), las = 1)
abline(lm(ESCA_data$ALDH7A1 ~ ESCA_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(ESCA_data$Aggregate, ESCA_data$ALDH7A1, method = "spearman")
dev.off()

# Add squamous
ESCA_clinical<-fread("ESCA-clinical.txt", data.table = F)
ESCA_data$Type<-NA
for(i in 1:nrow(ESCA_data)){
  ESCA_data$Type[i]<-ESCA_clinical[6,gsub("-", ".", ESCA_data$Samples[i], fixed = T)]
}
ESCA_data$Col<-NA
for(i in 1:nrow(ESCA_data)){
  if(ESCA_data$Type[i]=="esophagussquamouscellcarcinoma"){
    ESCA_data$Col[i]<-"firebrick"
  }else{
    ESCA_data$Col[i]<-"black"
  }
}
ESCA_data<-ESCA_data[order(ESCA_data$Col),]

plot(x = ESCA_data$Aggregate, y = ESCA_data$ALDH7A1,
     pch = 19, xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "ESCA", cex = 1.5,
     xlim = c(0, 1), ylim = c(4, 12), las = 1,
     col = ESCA_data$Col)
abline(lm(ESCA_data$ALDH7A1 ~ ESCA_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
cor.test(ESCA_data$Aggregate, ESCA_data$ALDH7A1, method = "spearman")
legend("topright", col = c("firebrick", "black"), pch = 19, legend = c("Squamous", "Others"))



CESC_squamous_data<-subset(CESC_data, CESC_data$Col=="firebrick")
BLCA_squamous_data<-subset(BLCA_data, BLCA_data$Col=="firebrick")
ESCA_squamous_data<-subset(ESCA_data, ESCA_data$Col=="firebrick")

pdf("Merged-squamous.pdf", width = 25, height = 5)
par(mfrow = c(1,5))
plot(x = LUSC_data$Aggregate, y = LUSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "LUSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(LUSC_data$ALDH7A1 ~ LUSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = HNSC_data$Aggregate, y = HNSC_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "HNSC", cex = 1.5,
     xlim = c(0.1, 0.8), las = 1)
abline(lm(HNSC_data$ALDH7A1 ~ HNSC_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = CESC_squamous_data$Aggregate, y = CESC_squamous_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "CESC (squamous)", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8), las = 1)
abline(lm(CESC_squamous_data$ALDH7A1 ~ CESC_squamous_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = BLCA_squamous_data$Aggregate, y = BLCA_squamous_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "BLCA (squamous)", cex = 1.5,
     xlim = c(0.1, 0.8), ylim = c(0, 8.5), las = 1)
abline(lm(BLCA_squamous_data$ALDH7A1 ~ BLCA_squamous_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)

plot(x = ESCA_squamous_data$Aggregate, y = ESCA_squamous_data$ALDH7A1,
     pch = 19, col = "black", xlab = "Aggregated beta value",
     ylab = "ALDH7A1", main = "ESCA (squamous)", cex = 1.5,
     xlim = c(0, 1), ylim = c(4, 12), las = 1)
abline(lm(ESCA_squamous_data$ALDH7A1 ~ ESCA_squamous_data$Aggregate),
       lwd = 5, col = "firebrick", lty = 2)
dev.off()

LUSC_res<-cor.test(LUSC_data$Aggregate, LUSC_data$ALDH7A1, method = "spearman")
HNSC_res<-cor.test(HNSC_data$Aggregate, HNSC_data$ALDH7A1, method = "spearman")
CESC_sq_res<-cor.test(CESC_squamous_data$Aggregate, CESC_squamous_data$ALDH7A1, method = "spearman")
BLCA_sq_res<-cor.test(BLCA_squamous_data$Aggregate, BLCA_squamous_data$ALDH7A1, method = "spearman")
ESCA_sq_res<-cor.test(ESCA_squamous_data$Aggregate, ESCA_squamous_data$ALDH7A1, method = "spearman")

results_tb<-matrix(NA, nrow = 2, ncol = 5)
rownames(results_tb)<-c("rho", "p-value")
colnames(results_tb)<-c("LUSC", "HNSC", "CESC (sq)", "BLCA (sq)", "ESCA (sq)")
results_tb[1,1]<-LUSC_res$estimate
results_tb[2,1]<-LUSC_res$p.value
results_tb[1,2]<-HNSC_res$estimate
results_tb[2,2]<-HNSC_res$p.value
results_tb[1,3]<-CESC_sq_res$estimate
results_tb[2,3]<-CESC_sq_res$p.value
results_tb[1,4]<-BLCA_sq_res$estimate
results_tb[2,4]<-BLCA_sq_res$p.value
results_tb[1,5]<-ESCA_sq_res$estimate
results_tb[2,5]<-ESCA_sq_res$p.value
write.table(results_tb, "Results_tb.txt", sep = "\t", row.names = T, col.names = T)



