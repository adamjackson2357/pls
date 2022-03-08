rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
setwd("../")

library(lme4)
library(stringr)
library(RColorBrewer)
library(tableone)

#### Load the data

covars=readRDS("Data/Covariates.rds")
proteins=readRDS("Data/Proteins.rds")
dim(proteins)
dim(covars)
print(all(rownames(covars)==rownames(proteins)))
str(covars)

# Challenge 1

myVars <- c("cohort", "phase", "sex", "age", "Imputed_alcohol", "Imputed_bms", "Imputed_smoking_status",
            "Imputed_education", "Imputed_activity", "LY_subtype")
## Vector of categorical variables that need transformation
catVars <- c("status", "trt", "ascites", "hepato",
             "spiders", "edema", "stage")
## Create a TableOne object
tab2 <- CreateTableOne(vars = myVars, data = covars,strata = "Status")

# Challenge 2
dir.create("Figures", showWarnings = FALSE)

pdf("Figures/Distributions.pdf")
par(mfrow=c(3,3), mar=c(3,4.5,1,1))
for (k in 1:ncol(proteins)){
  xfull=density(proteins[,k])
  x0=density(proteins[as.character(covars$type)=="0",k])
  x1=density(proteins[as.character(covars$type)=="1",k])
  plot(density(proteins[,k]), col="skyblue", xlab="", main="", 
       ylab=colnames(proteins)[k], las=1, ylim=range(c(xfull$y, x0$y, x1$y)))
  lines(x0, col="darkgreen")
  lines(x1, col="tomato")
  if (k==1){
    legend("topleft", lwd=2, col=c("skyblue", "darkgreen", "tomato"), 
           legend = c("Full sample", "Controls", "Cases"), cex=0.7)
  }
}
dev.off()

### Linear Mixed Models

#### Non adjusted on WBC

### Pooled LY cases (Challenge 3)

denoised=NULL
Beta_pooled=NULL
pvalue_pooled=NULL

f0='proteins[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + (1 | plate)'
f1=paste(f0, '+ type')

for (k in seq(1:ncol(proteins))){
  print(k)
  model=lmer(as.formula(f1), data=covars, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  model0=lmer(as.formula(f0), data=covars, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  pvalue_pooled=c(pvalue_pooled, anova(model, model0)$'Pr(>Chisq)'[2])
  
  beta=fixef(model)[c('(Intercept)', 'type1')]
  Beta_pooled=c(Beta_pooled, fixef(model)['type1'])
  
  X=cbind(rep(1, length(covars$type)), as.numeric(covars$type))
  denoised=cbind(denoised, (X%*%beta + resid(model)))
}

colnames(denoised)=colnames(proteins)
rownames(denoised)=rownames(proteins)
saveRDS(denoised, "Data/Proteins_denoised.rds")

Table_pooled=cbind(Beta_pooled, pvalue_pooled)
rownames(Table_pooled)=colnames(proteins)


### By subtype ( pooled Challenge 5)

f0='proteins_subtype[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + (1 | plate)'
f1=paste(f0, '+ type')

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  print(subtype)
  
  ids=c(covars$egm_id[covars$LY_subtype==""], covars$egm_id[covars$LY_subtype==subtype])
  covars_subtype=covars[covars$egm_id%in%ids,]
  proteins_subtype=proteins[covars_subtype$egm_id, ]
  Beta_subtype=NULL
  pvalue_subtype=NULL
  
  for (k in seq(1:ncol(proteins_subtype))){
    model=lmer(as.formula(f1), data=covars_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    model0=lmer(as.formula(f0), data=covars_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    pvalue_subtype=c(pvalue_subtype, anova(model, model0)$'Pr(>Chisq)'[2])
    
    Beta_subtype=c(Beta_subtype, fixef(model)['type1'])
  }
  
  Table_subtype=cbind(Beta_subtype, pvalue_subtype)
  rownames(Table_subtype)=colnames(proteins_subtype)
  
  assign(paste0("Table_", subtype), Table_subtype)
  
  Table_subtype=as.data.frame(Table_subtype)
  Table_subtype=cbind(round(Table_subtype$Beta_subtype, digits = 2), sprintf("%.2e", Table_subtype$pvalue_subtype))
  rownames(Table_subtype)=colnames(proteins)
  colnames(Table_subtype)=c("Beta", "pvalue")
  
  assign(paste0("Table_ready_", subtype), Table_subtype)
}

Table=cbind(Table_pooled, Table_BCLL, Table_DLBL, Table_FL, Table_MM)
colnames(Table)=paste(rep(c("Pooled", "CLL", "DLBCL", "FL", "MM"), each=2), rep(c("Beta", "Pvalue"), 5))
dir.create("Tables", showWarnings = FALSE)
write.table(Table, "Tables/Univariate_full_non_adjusted_wbc.txt", sep="\t")
# > supplementary table 5


#### Adjusted on WBC

### Pooled cases (Challenge 4)

wbc=readRDS("Data/wbc.rds")
covars_wbc=merge(covars, wbc, by = "row.names", all = FALSE)
rownames(covars_wbc)=covars_wbc$Row.names
proteins_wbc=proteins[rownames(covars_wbc),]

Beta_pooled_wbc=NULL
pvalue_pooled_wbc=NULL

f0='proteins_wbc[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + CD4 + NK + Bcells + Mono + CD8 + (1 | plate)'
f1=paste(f0, '+ type')

for (k in seq(1:ncol(proteins_wbc))){
  print(k)
  model=lmer(as.formula(f1), data=covars_wbc, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  model0=lmer(as.formula(f0), data=covars_wbc, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  pvalue_pooled_wbc=c(pvalue_pooled_wbc, anova(model, model0)$'Pr(>Chisq)'[2])
  
  Beta_pooled_wbc=c(Beta_pooled_wbc, fixef(model)['type1'])
}

Table_pooled_wbc=cbind(Beta_pooled_wbc, pvalue_pooled_wbc)
rownames(Table_pooled_wbc)=colnames(proteins_wbc)

### By subtype (subtype Challenge 5)

f0='proteins_wbc_subtype[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + CD4 + NK + Bcells + Mono + CD8 + (1 | plate)'
f1=paste(f0, '+ type')

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  print(subtype)
  
  ids=c(covars_wbc$egm_id[covars_wbc$LY_subtype==""], covars_wbc$egm_id[covars_wbc$LY_subtype==subtype])
  covars_wbc_subtype=covars_wbc[covars_wbc$egm_id%in%ids,]
  proteins_wbc_subtype=proteins_wbc[rownames(covars_wbc_subtype),]
  
  Beta_subtype=NULL
  pvalue_subtype=NULL
  
  for (k in seq(1:ncol(proteins_wbc_subtype))){
    # print(k)
    model=lmer(as.formula(f1), data=covars_wbc_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    model0=lmer(as.formula(f0), data=covars_wbc_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    pvalue_subtype=c(pvalue_subtype, anova(model, model0)$'Pr(>Chisq)'[2])
    
    Beta_subtype=c(Beta_subtype, fixef(model)['type1'])
  }
  
  Table_subtype=cbind(Beta_subtype, pvalue_subtype)
  rownames(Table_subtype)=colnames(proteins_wbc_subtype)
  
  assign(paste0("Table_", subtype, "_wbc"), Table_subtype)
}

Table_wbc=cbind(Table_pooled_wbc, Table_BCLL_wbc, Table_DLBL_wbc, Table_FL_wbc, Table_MM_wbc)
colnames(Table_wbc)=paste(rep(c("Pooled", "CLL", "DLBCL", "FL", "MM"), each=2), rep(c("Beta", "Pvalue"), 5))
write.table(Table_wbc, "Tables/Univariate_subset_adjusted_wbc.txt")
# > supplementary table 6


#### Non adjusted on 224 pairs

### Pooled 

Beta_pooled_no_wbc=NULL
pvalue_pooled_no_wbc=NULL

f0='proteins_wbc[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + (1 | plate)'
f1=paste(f0, '+ type')

for (k in seq(1:ncol(proteins_wbc))){
  print(k)
  model=lmer(as.formula(f1), data=covars_wbc, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  model0=lmer(as.formula(f0), data=covars_wbc, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
  pvalue_pooled_no_wbc=c(pvalue_pooled_no_wbc, anova(model, model0)$'Pr(>Chisq)'[2])
  
  Beta_pooled_no_wbc=c(Beta_pooled_no_wbc, fixef(model)['type1'])
}

Table_pooled_no_wbc=cbind(Beta_pooled_no_wbc, pvalue_pooled_no_wbc)
rownames(Table_pooled_no_wbc)=colnames(proteins_wbc)


### By subtype

f0='proteins_wbc_subtype[,k] ~ age + sex + cohort + phase + Imputed_bmi + Imputed_education + Imputed_activity + Imputed_smoking_status + Imputed_alcohol + (1 | plate)'
f1=paste(f0, '+ type')

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  print(subtype)
  
  ids=c(covars_wbc$egm_id[covars_wbc$LY_subtype==""], covars_wbc$egm_id[covars_wbc$LY_subtype==subtype])
  covars_wbc_subtype=covars_wbc[covars_wbc$egm_id%in%ids,]
  proteins_wbc_subtype=proteins_wbc[rownames(covars_wbc_subtype),]
  Beta_subtype=NULL
  pvalue_subtype=NULL
  
  for (k in seq(1:ncol(proteins_wbc_subtype))){
    model=lmer(as.formula(f1), data=covars_wbc_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    model0=lmer(as.formula(f0), data=covars_wbc_subtype, REML=FALSE, control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-04)))
    pvalue_subtype=c(pvalue_subtype, anova(model, model0)$'Pr(>Chisq)'[2])
    
    Beta_subtype=c(Beta_subtype, fixef(model)['type1'])
  }
  
  Table_subtype=cbind(Beta_subtype, pvalue_subtype)
  rownames(Table_subtype)=colnames(proteins_wbc_subtype)
  
  assign(paste0("Table_", subtype, "_no_wbc"), Table_subtype)
}

Table_no_wbc=cbind(Table_pooled_no_wbc, Table_BCLL_no_wbc, Table_DLBL_no_wbc, Table_FL_no_wbc, Table_MM_no_wbc)
colnames(Table_no_wbc)=paste(rep(c("Pooled", "CLL", "DLBCL", "FL", "MM"), each=2), rep(c("Beta", "Pvalue"), 5))
write.table(Table_no_wbc, "Tables/Univariate_subset_non_adjusted_wbc.txt")


#### Manhattan plots

Table=read.table("Tables/Univariate_full_non_adjusted_wbc.txt")
Table_wbc=read.table("Tables/Univariate_subset_adjusted_wbc.txt")
Table_no_wbc=read.table("Tables/Univariate_subset_non_adjusted_wbc.txt")

MyPal=brewer.pal("Paired", n = 12)
MyPalPairsExt <- colorRampPalette(MyPal)(28)

Groups=c("EGF", "FGF2", "GCSF", "VEGF", "GMSCF", "TGFa",
         "Eotaxin", "Fractalkine", "GRO", "MCP1", "MCP3", "MDC", "MIP1a", "MIP1b", "IP10", "IL8", 
         "IL1b", "IL2", "IL4", "IL5", "IL6", "IL7", "IL10", "IL13", "INFa", "INFg", "TNFa", "sCD40L")

Table=Table[Groups, ]
Table_wbc=Table_wbc[Groups, ]
Table_no_wbc=Table_no_wbc[Groups, ]

dir.create("Figures", showWarnings = FALSE)
pdf("Figures/Fig1_A.pdf", height = 7, width = 10)
par(mar=c(6,5,3,3))
plot(-log10(Table$Pooled.Pvalue), pch=17, col=ifelse(Table$Pooled.Beta<0, yes=MyPal[2], no=MyPal[8]), 
     xaxt="n", ylab=expression(-log[10](italic(p))), xlab="", las=1,
     ylim=c(min(c(-log10(Table$Pooled.Pvalue), -log10(Table_wbc$Pooled.Pvalue))), max(c(-log10(Table$Pooled.Pvalue), -log10(Table_wbc$Pooled.Pvalue)))))
points(-log10(Table_wbc$Pooled.Pvalue), pch=19, col=ifelse(Table_wbc$Pooled.Beta<0, yes=MyPal[2], no=MyPal[8]))
points(-log10(Table_no_wbc$Pooled.Pvalue), pch=18, col=ifelse(Table_no_wbc$Pooled.Beta<0, yes=MyPal[2], no=MyPal[8]))
abline(h=-log10(0.05/28), col="red")
abline(v=seq(1, 28), lty=3, col="grey")
axis(1, at=1:28, ifelse(Table_wbc$Pooled.Beta<0, yes=rownames(Table), no=''), las=2, col.axis=MyPal[2])
axis(1, at=1:28, ifelse(Table_wbc$Pooled.Beta<0, yes='', no=rownames(Table)), las=2, col.axis=MyPal[8])
legend("topright", pch=c(19, 17, 18, 19, 19), col=c("black", "black", "black", MyPal[2], MyPal[8]),
       legend = c("Subset adjusted on WBC", "Full not adjusted on WBC", "Subset not adjusted on WBC", paste("Negative", expression(beta)), paste("Positive", expression(beta))))
dev.off()


### By subtype

for(k in 1:4){
  subtype=c("BCLL", "DLBL", "FL", "MM")[k]
  subtype_name=c("CLL", "DLBCL", "FL", "MM")[k]
  
  Table1=eval(parse(text=paste0("Table_", subtype)))
  Table2=eval(parse(text=paste0("Table_", subtype, "_wbc")))
  Table3=eval(parse(text=paste0("Table_", subtype, "_no_wbc")))
  
  Table1=Table1[Groups, ]
  Table2=Table2[Groups, ]
  Table3=Table3[Groups, ]
  
  pdf(paste0("Figures/Fig1_", subtype_name, ".pdf"), height = 7, width = 10)
  par(mar=c(6,5,3,3))
  plot(-log10(Table1[,2]), pch=17, col=ifelse(Table1[,1]<0, yes=MyPal[2], no=MyPal[8]), 
       xaxt="n", ylab=expression(-log[10](italic(p))), xlab="", las=1,
       ylim=c(min(c(-log10(Table$MM.Pvalue), -log10(Table_wbc$MM.Pvalue))), max(c(-log10(Table$MM.Pvalue), -log10(Table_wbc$MM.Pvalue)))),
       main=paste0(subtype_name, " (N=", sum(covars$type==0)+sum(covars$LY_subtype==subtype), ")"))
  points(-log10(Table2[,2]), pch=19, col=ifelse(Table2[,1]<0, yes=MyPal[2], no=MyPal[8]))
  points(-log10(Table3[,2]), pch=18, col=ifelse(Table3[,1]<0, yes=MyPal[2], no=MyPal[8]))
  abline(h=-log10(0.05/28), col="red")
  abline(v=seq(1, 28), lty=3, col="grey")
  axis(1, at=1:28, ifelse(Table2[,1]<0, yes=rownames(Table), no=''), las=2, col.axis=MyPal[2])
  axis(1, at=1:28, ifelse(Table2[,1]<0, yes='', no=rownames(Table)), las=2, col.axis=MyPal[8])
  legend("topright", pch=c(19, 17, 18, 19, 19), col=c("black", "black", "black", MyPal[2], MyPal[8]),
         legend = c("Subset adjusted on WBC", "Full not adjusted on WBC", "Subset not adjusted on WBC", paste("Negative", expression(beta)), paste("Positive", expression(beta))))
  dev.off()
}

