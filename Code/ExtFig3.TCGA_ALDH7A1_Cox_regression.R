#=============================================
# Extended Fig 3
# TCGA squamous carcinoma
# ALDH7A1 Cox regression
#=============================================
# Set working directory
setwd("ExtendedFig3")

# Load packages
library(data.table)
library(survival)

# LUSC
LUSC_expression<-fread("LUSC-expression.txt",
                       data.table = F)
LUSC_exposure<-fread("LUSC-exposure.tsv",
                     data.table = F)
LUSC_clinical<-fread("LUSC-clinical.tsv",
                     data.table = F)
tmp_clinical<-data.frame("Sample"=unique(LUSC_clinical$cases.submitter_id),
                         "Vital_status"=NA, "Time"=NA,
                         "Stage"=NA, "Sex"=NA, "Age"=NA,
                         "Smoking"=NA,
                         "ALDH7A1_exp"=NA,
                         stringsAsFactors = F)
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(LUSC_clinical, LUSC_clinical$cases.submitter_id==tmp_clinical$Sample[i])
  tmp_days_to_death<-as.numeric(unique(tmp$demographic.days_to_death))
  if(!is.na(tmp_days_to_death)){
    tmp_clinical$Vital_status[i]<-1
    tmp_clinical$Time[i]<-tmp_days_to_death
  }else{
    tmp_clinical$Vital_status[i]<-0
    tmp_clinical$Time[i]<-unique(as.numeric(tmp$diagnoses.days_to_last_follow_up))
  }
  tmp<-subset(tmp, tmp$diagnoses.diagnosis_is_primary_disease)
  
  tmp_clinical$Stage[i]<-unique(tmp$diagnoses.ajcc_pathologic_stage)
  tmp_clinical$Sex[i]<-ifelse(unique(tmp$demographic.sex_at_birth)=="male", 1, 2)
  tmp_clinical$Age[i]<-round(as.numeric(unique(tmp$diagnoses.age_at_diagnosis))/365, 0)
  
  tmp<-subset(LUSC_exposure, LUSC_exposure$cases.submitter_id==tmp_clinical$Sample[i])
  tmp_clinical$Smoking[i]<-tmp$exposures.tobacco_smoking_status
}
tmp_clinical$Stage[which(tmp_clinical$Stage=="'--")]<-NA
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IA")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IB")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIA")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIB")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIIA")]<-"Stage III"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIIB")]<-"Stage III"

tmp_clinical$Stage[which(tmp_clinical$Smoking=="Unknown")]<-NA
tmp_clinical$Stage[which(tmp_clinical$Smoking=="Not Reported")]<-NA
tmp_clinical$Smoking <- factor(tmp_clinical$Smoking, 
                              levels = c("Lifelong Non-Smoker",
                                         "Current Reformed Smoker, Duration Not Specified",
                                         "Current Reformed Smoker for < or = 15 yrs",
                                         "Current Reformed Smoker for > 15 yrs",
                                         "Current Smoker"))


ALDH7A1_LUSC_exp<-LUSC_expression[which(LUSC_expression$V1=="ALDH7A1"), -1]
for(i in 1:nrow(tmp_clinical)){
  if(any(grepl(tmp_clinical$Sample[i], colnames(ALDH7A1_LUSC_exp)))){
    tmp_clinical$ALDH7A1_exp[i]<-mean(as.numeric(ALDH7A1_LUSC_exp[,grep(tmp_clinical$Sample[i], colnames(ALDH7A1_LUSC_exp))]))
  }
}

LUSC_clinical<-tmp_clinical
LUSC_clinical<-LUSC_clinical[!is.na(LUSC_clinical$Vital_status),]
LUSC_clinical$Vital_status<-LUSC_clinical$Vital_status+1

LUSC_clinical<-LUSC_clinical[!is.na(LUSC_clinical$Time),]
LUSC_clinical<-LUSC_clinical[!is.na(LUSC_clinical$ALDH7A1_exp),]

res<-coxph(Surv(Time, Vital_status) ~ ALDH7A1_exp + Stage + Age + Sex + Smoking, data = LUSC_clinical)
LUSC_res<-summary(res)

# HNSC
HNSC_expression<-fread("HNSC-expression.cct",
                         data.table = F)
ALDH7A1_HNSC_exp<-HNSC_expression[which(HNSC_expression$attrib_name=="ALDH7A1"), -1]
colnames(ALDH7A1_HNSC_exp)<-gsub(".", "-", colnames(ALDH7A1_HNSC_exp), fixed = T)
ALDH7A1_HNSC_exp[1,]<-scale(as.numeric(ALDH7A1_HNSC_exp))[,1]

HNSC_clinical<-fread("HNSC-clinical.tsv",
                     data.table = F)
tmp_clinical<-data.frame("Sample"=unique(HNSC_clinical$cases.submitter_id),
                         "Vital_status"=NA, "Time"=NA,
                         "Stage"=NA, "Sex"=NA, "Age"=NA,
                         "Subsite_oropharnyx"=NA,
                         "ALDH7A1_exp"=NA,
                         stringsAsFactors = F)
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(HNSC_clinical, HNSC_clinical$cases.submitter_id==tmp_clinical$Sample[i])
  tmp_days_to_death<-as.numeric(unique(tmp$demographic.days_to_death))
  if(!is.na(tmp_days_to_death)){
    tmp_clinical$Vital_status[i]<-2
    tmp_clinical$Time[i]<-tmp_days_to_death
  }else{
    tmp_clinical$Vital_status[i]<-1
    tmp_clinical$Time[i]<-unique(as.numeric(tmp$diagnoses.days_to_last_follow_up))
  }
  tmp<-subset(tmp, tmp$diagnoses.diagnosis_is_primary_disease)
  
  tmp_clinical$Stage[i]<-unique(tmp$diagnoses.ajcc_pathologic_stage)
  tmp_clinical$Sex[i]<-ifelse(unique(tmp$demographic.sex_at_birth)=="male", 1, 2)
  tmp_clinical$Age[i]<-round(as.numeric(unique(tmp$diagnoses.age_at_diagnosis))/365, 0)
  
  tmp_clinical$Subsite_oropharnyx[i]<-ifelse(unique(tmp$cases.primary_site)=="Oropharynx", 2, 1)
}
tmp_clinical$Stage[which(tmp_clinical$Stage=="'--")]<-NA
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IVA")]<-"Stage IV"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IVB")]<-"Stage IV"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IVC")]<-"Stage IV"

for(i in 1:nrow(tmp_clinical)){
  if(any(grepl(tmp_clinical$Sample[i], colnames(ALDH7A1_HNSC_exp)))){
    tmp_clinical$ALDH7A1_exp[i]<-mean(as.numeric(ALDH7A1_HNSC_exp[,grep(tmp_clinical$Sample[i], colnames(ALDH7A1_HNSC_exp))]))
  }
}

# HPV status
HNSC_HPV<-read.delim("HNSC_HPV_status_499.tsv",
                     sep = "\t", header = T, stringsAsFactors = F)
tmp_clinical$HPV<-NA
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(HNSC_HPV, HNSC_HPV$case_id==tmp_clinical$Sample[i])
  if(nrow(tmp)!=0){
    tmp_clinical$HPV[i]<-tmp$HPV_status
  }
}

HNSC_clinical<-tmp_clinical
HNSC_clinical<-HNSC_clinical[!is.na(HNSC_clinical$Vital_status),]
HNSC_clinical<-HNSC_clinical[!is.na(HNSC_clinical$Time),]
HNSC_clinical<-HNSC_clinical[!is.na(HNSC_clinical$ALDH7A1_exp),]
HNSC_clinical<-HNSC_clinical[!is.na(HNSC_clinical$HPV),]

res<-coxph(Surv(Time, Vital_status) ~ ALDH7A1_exp + Stage + Age + Sex + Subsite_oropharnyx + HPV, data = HNSC_clinical)
HNSC_res<-summary(res)


# CESC
CESC_expression<-fread("CESC-expression.txt", data.table = F)
ALDH7A1_CESC_exp<-CESC_expression[which(CESC_expression$V1=="ALDH7A1"), -1]
CESC_clinical<-fread("CESC-clinical.tsv", data.table = F)
CESC_clinical<-subset(CESC_clinical, CESC_clinical$cases.disease_type=="Squamous Cell Neoplasms")
tmp_clinical<-data.frame("Sample"=unique(CESC_clinical$cases.submitter_id),
                         "Vital_status"=NA, "Time"=NA,
                         "Stage"=NA, "Sex"=NA, "Age"=NA,
                         "ALDH7A1_exp"=NA,
                         stringsAsFactors = F)
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(CESC_clinical, CESC_clinical$cases.submitter_id==tmp_clinical$Sample[i])
  tmp_days_to_death<-as.numeric(unique(tmp$demographic.days_to_death))
  if(!is.na(tmp_days_to_death)){
    tmp_clinical$Vital_status[i]<-1
    tmp_clinical$Time[i]<-tmp_days_to_death
  }else{
    tmp_clinical$Vital_status[i]<-0
    tmp_clinical$Time[i]<-unique(as.numeric(tmp$diagnoses.days_to_last_follow_up))
  }
  tmp<-subset(tmp, tmp$diagnoses.diagnosis_is_primary_disease)
  
  tmp_clinical$Stage[i]<-unique(tmp$diagnoses.figo_stage)
  tmp_clinical$Sex[i]<-ifelse(unique(tmp$demographic.sex_at_birth)=="male", 1, 2)
  tmp_clinical$Age[i]<-round(as.numeric(unique(tmp$diagnoses.age_at_diagnosis))/365, 0)
}

CESC_new_clinical<-read.delim("CESC-clinical-1.txt", sep = "\t", header = T)
for(i in 1:nrow(tmp_clinical)){
  if(is.na(tmp_clinical$Time[i])){
    tmp_clinical$Time[i]<-subset(CESC_new_clinical, CESC_new_clinical$Sample==tmp_clinical$Sample[i])$Time
  }
}

tmp_clinical$Stage[which(tmp_clinical$Stage=="'--")]<-NA
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IA")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IA1")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IA2")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IB")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IB1")]<-"Stage I"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IB2")]<-"Stage I"

tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIA")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIA1")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIA2")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIB")]<-"Stage II"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIIA")]<-"Stage III"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IIIB")]<-"Stage III"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IVA")]<-"Stage IV"
tmp_clinical$Stage[which(tmp_clinical$Stage=="Stage IVB")]<-"Stage IV"

for(i in 1:nrow(tmp_clinical)){
  if(any(grepl(tmp_clinical$Sample[i], colnames(ALDH7A1_CESC_exp)))){
    tmp_clinical$ALDH7A1_exp[i]<-mean(as.numeric(ALDH7A1_CESC_exp[,grep(tmp_clinical$Sample[i], colnames(ALDH7A1_CESC_exp))]))
  }
}

CESC_clinical<-tmp_clinical
CESC_clinical<-CESC_clinical[!is.na(CESC_clinical$Vital_status),]
CESC_clinical$Vital_status<-CESC_clinical$Vital_status+1

CESC_clinical<-CESC_clinical[!is.na(CESC_clinical$Time),]
CESC_clinical<-CESC_clinical[!is.na(CESC_clinical$ALDH7A1_exp),]

# squamous
# Read clinical data
CESC_tmp_clinical<-fread("CESC-clinical.tsv", data.table = F)

# CESC
unique(CESC_tmp_clinical$cases.disease_type)
CESC_tmp_clinical<-CESC_tmp_clinical[,c("cases.submitter_id","diagnoses.primary_diagnosis")]
CESC_tmp_clinical<-unique(CESC_tmp_clinical)

CESC_clinical$Type<-NA
for(i in 1:nrow(CESC_clinical)){
  if(CESC_clinical$Sample[i]%in%CESC_tmp_clinical$cases.submitter_id){
    tmp1<-unique(subset(CESC_tmp_clinical, CESC_tmp_clinical$cases.submitter_id==CESC_clinical$Sample[i])$diagnoses.primary_diagnosis)
    if(any(grepl("squamous", tmp1, ignore.case = T))){
      CESC_clinical$Type[i]<-"squamous"
    }else{
      CESC_clinical$Type[i]<-"Not squamous"
    }
  }
}

CESC_sq_clinical<-subset(CESC_clinical, CESC_clinical$Type=="squamous")


res<-coxph(Surv(Time, Vital_status) ~ ALDH7A1_exp + Stage + Age, data = CESC_sq_clinical)
CESC_res<-summary(res)

# BLCA
BLCA_expression<-fread("BLCA-expression.txt", data.table = F)
ALDH7A1_BLCA_exp<-BLCA_expression[which(BLCA_expression$V1=="ALDH7A1"), -1]
BLCA_clinical<-fread("BLCA-clinical.tsv", data.table = F)
BLCA_tmp<-fread("BLCA-clinical-1.txt", sep = "\t", data.table = F)
BLCA_tmp<-subset(BLCA_tmp, BLCA_tmp$`mRNA cluster`=="Basal_squamous")
BLCA_clinical<-subset(BLCA_clinical, BLCA_clinical$cases.submitter_id%in%BLCA_tmp$`Case ID`)
tmp_clinical<-data.frame("Sample"=unique(BLCA_clinical$cases.submitter_id),
                         "Vital_status"=NA, "Time"=NA,
                         "Stage"=NA, "Sex"=NA, "Age"=NA,
                         "ALDH7A1_exp"=NA,
                         stringsAsFactors = F)
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(BLCA_clinical, BLCA_clinical$cases.submitter_id==tmp_clinical$Sample[i])
  tmp_days_to_death<-as.numeric(unique(tmp$demographic.days_to_death))
  if(!is.na(tmp_days_to_death)){
    tmp_clinical$Vital_status[i]<-1
    tmp_clinical$Time[i]<-tmp_days_to_death
  }else{
    tmp_clinical$Vital_status[i]<-0
    tmp_clinical$Time[i]<-unique(as.numeric(tmp$diagnoses.days_to_last_follow_up))
  }
  tmp<-subset(tmp, as.logical(tmp$diagnoses.diagnosis_is_primary_disease))
  
  tmp_clinical$Stage[i]<-unique(tmp$diagnoses.ajcc_pathologic_stage)
  tmp_clinical$Sex[i]<-ifelse(unique(tmp$demographic.sex_at_birth)=="male", 1, 2)
  tmp_clinical$Age[i]<-round(as.numeric(unique(tmp$diagnoses.age_at_diagnosis))/365, 0)
}

BLCA_new_clinical<-read.delim("BLCA-clinical-2.txt", sep = "\t", header = T)
for(i in 1:nrow(tmp_clinical)){
  if(is.na(tmp_clinical$Time[i])){
    tmp_clinical$Time[i]<-subset(BLCA_new_clinical, BLCA_new_clinical$Sample==tmp_clinical$Sample[i])$Time
  }
}
for(i in 1:nrow(tmp_clinical)){
  if(any(grepl(tmp_clinical$Sample[i], colnames(ALDH7A1_BLCA_exp)))){
    tmp_clinical$ALDH7A1_exp[i]<-mean(as.numeric(ALDH7A1_BLCA_exp[,grep(tmp_clinical$Sample[i], colnames(ALDH7A1_BLCA_exp))]))
  }
}

BLCA_clinical<-tmp_clinical
BLCA_clinical<-BLCA_clinical[!is.na(BLCA_clinical$Vital_status),]
BLCA_clinical$Vital_status<-BLCA_clinical$Vital_status+1

BLCA_clinical<-BLCA_clinical[!is.na(BLCA_clinical$Time),]
BLCA_clinical<-BLCA_clinical[!is.na(BLCA_clinical$ALDH7A1_exp),]

BLCA_tmp_clinical<-fread("BLCA-clinical-1.txt", sep = "\t", data.table = F)

BLCA_clinical$Type<-NA
for(i in 1:nrow(BLCA_clinical)){
  if(BLCA_clinical$Sample[i]%in%BLCA_tmp_clinical$`Case ID`){
    tmp1<-unique(subset(BLCA_tmp_clinical, BLCA_tmp_clinical$`Case ID`==BLCA_clinical$Sample[i])$`mRNA cluster`)
    if(tmp1=="Basal_squamous"){
      BLCA_clinical$Type[i]<-"squamous"
    }else{
      BLCA_clinical$Type[i]<-"Not squamous"
    }
  }
}

res<-coxph(Surv(Time, Vital_status) ~ ALDH7A1_exp + Stage + Age, data = BLCA_clinical)
BLCA_res<-summary(res)

#===========#
# ESCA
#===========#
ESCA_expression<-fread("ESCA-expression.cct", data.table = F)
ALDH7A1_ESCA_exp<-ESCA_expression[which(ESCA_expression$attrib_name=="ALDH7A1"), -1]
colnames(ALDH7A1_ESCA_exp)<-sapply(1:ncol(ALDH7A1_ESCA_exp), function(x) gsub(".", "-", colnames(ALDH7A1_ESCA_exp)[x], fixed = T))

# scale
ALDH7A1_ESCA_exp<-t(scale(t(ALDH7A1_ESCA_exp[1,])))

ESCA_clinical<-fread("ESCA-clinical.txt", data.table = F)
rownames(ESCA_clinical)<-ESCA_clinical$attrib_name
ESCA_clinical<-ESCA_clinical[,-1]
ESCA_clinical<-t(ESCA_clinical)
ESCA_clinical<-as.data.frame(ESCA_clinical)
ESCA_clinical<-subset(ESCA_clinical, ESCA_clinical$histological_type=="esophagussquamouscellcarcinoma")
rownames(ESCA_clinical)<-sapply(1:nrow(ESCA_clinical), function(x) gsub(".", "-", rownames(ESCA_clinical)[x], fixed = T))
tmp_clinical<-data.frame("Sample"=unique(rownames(ESCA_clinical)),
                         "Vital_status"=NA, "Time"=NA,
                         "Stage"=NA, "Sex"=NA, "Age"=NA,
                         "ALDH7A1_exp"=NA,
                         stringsAsFactors = F)
for(i in 1:nrow(tmp_clinical)){
  tmp<-subset(ESCA_clinical, rownames(ESCA_clinical)==tmp_clinical$Sample[i])
  if(nrow(tmp)!=0){
    tmp_clinical$Vital_status[i]<-as.numeric(tmp$status)
    tmp_clinical$Time[i]<-as.numeric(tmp$overall_survival)
    tmp_clinical$Stage[i]<-unique(tmp$pathologic_stage)
    tmp_clinical$Sex[i]<-ifelse(unique(tmp$gender)=="male", 1, 2)
    tmp_clinical$Age[i]<-as.numeric(unique(tmp$years_to_birth))
  }
  if(tmp_clinical$Sample[i]%in%colnames(ALDH7A1_ESCA_exp)){
    tmp_clinical$ALDH7A1_exp[i]<-ALDH7A1_ESCA_exp[1,tmp_clinical$Sample[i]]
  }
}

ESCA_clinical<-tmp_clinical
ESCA_clinical<-ESCA_clinical[!is.na(ESCA_clinical$Vital_status),]
ESCA_clinical$Vital_status<-ESCA_clinical$Vital_status+1

ESCA_clinical<-ESCA_clinical[!is.na(ESCA_clinical$Time),]
ESCA_clinical<-ESCA_clinical[!is.na(ESCA_clinical$ALDH7A1_exp),]

res<-coxph(Surv(Time, Vital_status) ~ ALDH7A1_exp + Stage + Age, data = ESCA_clinical)
ESCA_res<-summary(res)

#==================================#
LUSC_coe<-LUSC_res$coefficients
LUSC_conf<-LUSC_res$conf.int

HNSC_coe<-HNSC_res$coefficients
HNSC_conf<-HNSC_res$conf.int

CESC_coe<-CESC_res$coefficients
CESC_conf<-CESC_res$conf.int

BLCA_coe<-BLCA_res$coefficients
BLCA_conf<-BLCA_res$conf.int

ESCA_coe<-ESCA_res$coefficients
ESCA_conf<-ESCA_res$conf.int

# Forest plot of multivariate Cox hazard ratios
ost_cox_plot<-data.frame("model"=c("LUSC", "HNSC", "CESC", "BLCA", "ESCA"),
                         "HR"=c(LUSC_conf[1,1], HNSC_conf[1,1], CESC_conf[1,1], BLCA_conf[1,1], ESCA_conf[1,1]),
                         "CI_low"=c(LUSC_conf[1,3], HNSC_conf[1,3], CESC_conf[1,3], BLCA_conf[1,3], ESCA_conf[1,3]),
                         "CI_high"=c(LUSC_conf[1,4], HNSC_conf[1,4], CESC_conf[1,4], BLCA_conf[1,4], ESCA_conf[1,4]),
                         "significance"=c(LUSC_coe[1,5], HNSC_coe[1,5], CESC_coe[1,5], BLCA_coe[1,5], ESCA_coe[1,5]))
write.table(ost_cox_plot, "Cox_results.txt", sep = "\t", quote = F, row.names = F, col.names = T)

# ---- 1. USER SETTINGS --------------------------------------
infile    <- "Cox_results.txt"
outstem   <- "forest_plot"
left_lab  <- "ALDH7A1 protective"   # text left of the null line
right_lab <- "ALDH7A1 harmful"      # text right of the null line
row_lab   <- "Cohort"               # header for the first column
x_title   <- "Hazard Ratio"

hr_limits <- c(0.6, 1.4)                 # x axis range, in HR units
hr_ticks  <- seq(0.6, 1.4, by = 0.1)     # tick positions, in HR units

# "log"    = log-spaced positions (standard for forest plots)
# "linear" = evenly spaced positions
x_spacing <- "linear"

fig_w  <- 11.33   # figure width  (inches)
row_h  <- 0.41    # height per row (inches)
base_h <- 1.00    # height for header + axis (inches)

# ---- 2. READ AND PREPARE ----------------------------------
d <- ost_cox_plot

d$hr_txt <- sprintf("%.2f (%.2f-%.2f)", d$HR, d$CI_low, d$CI_high)
d$p_txt  <- ifelse(d$significance < 0.01, "<0.01",
                   sprintf("%.2f", d$significance))
n <- nrow(d)

clipped <- d$CI_low < hr_limits[1] | d$CI_high > hr_limits[2]
if (any(clipped))
  warning("CI extends beyond hr_limits for: ",
          paste(d$model[clipped], collapse = ", "),
          " - widen hr_limits or the whiskers will be cut off.")

# ---- 3. AXIS ----------------------------------------------
tf     <- if (x_spacing == "log") log else identity
ax_lo  <- tf(hr_limits[1])
ax_hi  <- tf(hr_limits[2])
fmt_hr <- function(v) sub("\\.?0+$", "", sprintf("%.2f", v))

# ---- 4. LAYOUT GEOMETRY -----------------------------------
panel_L <- 0.42
panel_R <- 0.85
x_model <- 0.005
x_hr    <- 0.075
x_p     <- 0.895

map_hr <- function(v)
  panel_L + (tf(v) - ax_lo) / (ax_hi - ax_lo) * (panel_R - panel_L)

x_null <- map_hr(1)
y      <- n:1
y_head <- n + 1.1
y_axis <- 0.25

# ---- 5. DRAWING ROUTINE -----------------------------------
draw_forest <- function() {
  par(mar = c(0.4, 0.4, 0.4, 0.4), xpd = NA)
  plot(NA, xlim = c(0, 1), ylim = c(-1.05, y_head + 0.45),
       axes = FALSE, xlab = "", ylab = "")
  
  text(x_model, y_head, row_lab,                 adj = 0, font = 2, cex = 0.95)
  text(x_hr,    y_head, "Hazard Ratio (95% CI)", adj = 0, font = 2, cex = 0.95)
  text(x_p,     y_head, "p-value",               adj = 0, font = 2, cex = 0.95)
  text(x_null - 0.012, y_head, left_lab,  adj = 1, cex = 0.95)
  text(x_null + 0.012, y_head, right_lab, adj = 0, cex = 0.95)
  
  segments(x_null, y_axis, x_null, y_head - 0.45, lty = 2, lwd = 1.4)
  
  segments(map_hr(d$CI_low), y, map_hr(d$CI_high), y, lwd = 1.6)
  points(map_hr(d$HR), y, pch = 15, cex = 1.7)
  
  text(x_model, y, d$model,  adj = 0, font = 2, cex = 0.95)
  text(x_hr,    y, d$hr_txt, adj = 0, cex = 0.95)
  text(x_p,     y, d$p_txt,  adj = 0, cex = 0.95)
  
  segments(panel_L, y_axis, panel_R, y_axis, lwd = 1.2)
  tick_x <- map_hr(hr_ticks)
  segments(tick_x, y_axis, tick_x, y_axis - 0.12, lwd = 1.2)
  text(tick_x, y_axis - 0.34, fmt_hr(hr_ticks), cex = 0.9)
  text(mean(c(panel_L, panel_R)), y_axis - 0.85, x_title, cex = 1.0)
}

# ---- 6. OUTPUT --------------------------------------------
fig_h <- base_h + row_h * n

pdf(paste0(outstem, ".pdf"), width = fig_w, height = fig_h, useDingbats = FALSE)
draw_forest()
dev.off()

png(paste0(outstem, ".png"), width = fig_w * 150, height = fig_h * 150, res = 150)
draw_forest()
dev.off()

cat("axis:", paste(hr_limits, collapse = "-"), "| spacing:", x_spacing,
    "| ticks:", paste(fmt_hr(hr_ticks), collapse = " "), "\n")
cat("wrote", paste0(outstem, ".pdf"), "and", paste0(outstem, ".png"), "\n")
