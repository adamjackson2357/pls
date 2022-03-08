rm(list=ls())
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
setwd("../")

library(lme4)
library(sgPLS)
library(utils)
library(pheatmap)
library(RColorBrewer)
source("Scripts/pls_functions.R")


### Loading and preparing the data

Groups=c("EGF", "FGF2", "GCSF", "VEGF", "GMSCF", "TGFa",
         "Eotaxin", "Fractalkine", "GRO", "MCP1", "MCP3", "MDC", "MIP1a", "MIP1b", "IP10", "IL8", 
         "IL1b", "IL2", "IL4", "IL5", "IL6", "IL7", "IL10", "IL13", "INFa", "INFg", "TNFa", "sCD40L")

X_pooled=readRDS("Data/Proteins_denoised.rds")
covars=readRDS("Data/Covariates.rds")
X_pooled=X_pooled[,Groups]
Y_pooled=covars$type
print(all(rownames(X_pooled)==rownames(Y_pooled)))

Xgroups=c(6, 16)

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  sets=covars$egm_set[covars$LY_subtype==subtype]
  X=X_pooled[covars$egm_set%in%sets,]
  Y=Y_pooled[covars$egm_set%in%sets]
  assign(paste0("X_", subtype), X)
  assign(paste0("Y_", subtype), Y)
}

Y_diagnostic=covars$Time_to_diag[as.character(covars$type)=="1"] ## time to diagnostic in days
X_diagnostic=X_pooled[as.character(covars$type)=="1", ]
names(Y_diagnostic)=rownames(X_diagnostic)

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  ids=covars$egm_id[covars$LY_subtype==subtype]
  X=X_diagnostic[ids,]
  Y=Y_diagnostic[ids]
  assign(paste0("X_diag_", subtype), X)
  assign(paste0("Y_diag_", subtype), Y)
}

MyParameters=readRDS("Data/MyParameters.rds")

MyPal=brewer.pal("Paired", n = 12)


### Pooled cases (Fig 3)

MyPLSDA_pooled <- plsda(X_pooled, Y_pooled, ncomp=1, )
MysPLSDA_pooled=splsda(X_pooled, Y_pooled, keepX=MyParameters$Pooled$sPLSDA$NVar, ncomp=1, )
MygPLSDA_pooled <- gPLSda(X_pooled, Y_pooled, ncomp = 1, ind.block.x = Xgroups, keepX = MyParameters$Pooled$gPLSDA$NGroups)
MysgPLSDA_pooled <- sgPLSda(X_pooled, Y_pooled, ncomp = 1, ind.block.x = Xgroups, 
                            keepX = MyParameters$Pooled$sgPLSDA$NGroups, alpha.x = MyParameters$Pooled$sgPLSDA$alpha)


### Loadings plot

Loadings=cbind(MysPLSDA_pooled$loadings$X, MysgPLSDA_pooled$loadings$X, rep(NA, 28), rep(NA, 28))
Loadings=as.vector(t(Loadings))
Loadings=Loadings[-c(length(Loadings)-1, length(Loadings))]

pdf(paste0("Figures/Fig2_A.pdf"), height=6, width=10)
par(mar=c(10,5,3,3))
plot(Loadings, col=c(MyPal[6], MyPal[10], NA, NA), xaxt="n", ylab="Loadings Coefficients",  type="h", lwd=3, xlab="")
axis(1, at=seq(1.5, 28*4, by=4), labels = colnames(MyPLSDA_pooled$X), las=2)
axis(1, at=c(0, Xgroups, 28)*4, line = 6, labels = NA)
axis(1, at=c(3, 10.5, 21.5)*4, labels = c("Growth Factors", "Chemokines", "Cytokines"), line=6, tick = FALSE)
abline(v=c(0, Xgroups, 28)*4, lty=3, col="black")
legend("bottomright", legend = c("sPLS-DA", "sgPLS-DA"), lty=1, lwd=3,
       col=c(MyPal[6], MyPal[10]), cex=0.75)
dev.off()


### Compute scores from loadings:

S_X_pooled=X_pooled%*%MyPLSDA_pooled$loadings$X

### Misclassification rates plot

MyPredict=predict(MysPLSDA_pooled, newdata = X_pooled)
fitted=MyPredict$class$max.dist
table(fitted)

## Challenge 6

MSEP_sPLSDA=NULL
for(subtype in c("", "BCLL", "DLBL", "FL", "MM")){
  idx=which(covars$LY_subtype==subtype)
  MSEP_sPLSDA[[subtype]]=1-sum(diag(table(Y_pooled[idx], MyPredict$class$max.dist[idx])))/length(idx)
}

MyPredict=predict(MysgPLSDA_pooled, newdata = X_pooled)

MSEP_sgPLSDA=NULL
for(subtype in c("", "BCLL", "DLBL", "FL", "MM")){
  idx=which(covars$LY_subtype==subtype)
  MSEP_sgPLSDA[[subtype]]=1-sum(diag(table(Y_pooled[idx], MyPredict$class$max.dist[idx])))/length(idx)
}

MSEP=cbind(MSEP_sPLSDA, MSEP_sgPLSDA)
rownames(MSEP)=c("Controls", "CLL", "DLBCL", "FL", "MM")

MSEP=cbind(MSEP, rep(NA, 5), rep(NA, 5))
MSEP=unlist(t(MSEP)[1:(nrow(t(MSEP))*ncol(t(MSEP)))])
#MSEP=as.vector(t(MSEP))
MSEP=MSEP[-c(length(MSEP)-1, length(MSEP))]

pdf("Figures/Fig2_C.pdf", height = 5, width=9)
plot(MSEP, type = "h", lwd=3, xaxt="n", xlab="", ylab="Misclassification Rates", col=c(MyPal[6], MyPal[10], NA, NA),
     ylim=c(0, max(MSEP[!is.na(MSEP)])), las=1)
axis(1, at=seq(1.5, 5*4, by=4), labels = c("Controls", "CLL", "DLBCL", "FL", "MM"))
dev.off()


### Stability throughout cut-offs: pooled cases

Stability_results=StabilityPlot(X = X_pooled, Y = Y_pooled, NIter = 100)
saveRDS(Stability_results, "Results/stability_splsda.rds")

pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, filename = "Figures/Fig2_B.pdf", height=5, width=10)


### Run PLS by subtype # Challenge 7

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  X=eval(parse(text=paste0("X_", subtype)))
  Y=eval(parse(text=paste0("Y_", subtype)))
  
  MysPLSDA <- splsda(X, Y, keepX=eval(parse(text=paste0("MyParameters$", subtype, "$sPLSDA$NVar"))), ncomp=1, )
  assign(paste0("MysPLSDA_", subtype), MysPLSDA)
  
  MysgPLSDA <- sgPLSda(X, Y, ncomp = 1, ind.block.x = Xgroups, 
                       keepX = eval(parse(text=paste0("MyParameters$", subtype, "$sgPLSDA$NGroups"))), 
                       alpha.x = eval(parse(text=paste0("MyParameters$", subtype, "$sgPLSDA$alpha"))))
  assign(paste0("MysgPLSDA_", subtype), MysgPLSDA)
}

Loadings=cbind(MysPLSDA_BCLL$loadings$X, MysPLSDA_DLBL$loadings$X, MysPLSDA_MM$loadings$X, rep(NA, 28), rep(NA, 28))
Loadings=as.vector(t(Loadings))

pdf(paste0("Figures/Fig3_A.pdf"), height=6, width=10)
par(mar=c(10,5,3,3))
plot(Loadings, col=c(MyPal[1], MyPal[2], MyPal[3], NA, NA), xaxt="n", ylab="Loadings Coefficients",  type="h", lwd=3, xlab="", 
     ylim=c(min(Loadings[!is.na(Loadings)]), 0.9))
axis(1, at=seq(1.5, 28*5, by=5), labels = colnames(X_pooled), las=2)
axis(1, at=c(0, Xgroups, 28)*5, line = 6, labels = NA)
axis(1, at=c(3, 10.5, 21.5)*5, labels = c("Growth Factors", "Chemokines", "Cytokines"), line=6, tick = FALSE)
abline(v=c(0, Xgroups, 28)*5, lty=3, col="black")
legend("topright", lty=1, lwd=3, col=c(MyPal[1], MyPal[2], MyPal[3]), legend = c("CLL", "DLBCL", "MM"))
dev.off()

### Challenge 9
Loadings=cbind(MysgPLSDA_BCLL$loadings$X, MysgPLSDA_DLBL$loadings$X, MysgPLSDA_MM$loadings$X, rep(NA, 28), rep(NA, 28))
Loadings=as.vector(t(Loadings))

pdf(paste0("Figures/Fig3_B.pdf"), height=6, width=10)
par(mar=c(10,5,3,3))
plot(Loadings, col=c(MyPal[1], MyPal[2], MyPal[3], NA, NA), xaxt="n", ylab="Loadings Coefficients",  type="h", lwd=3, xlab="", 
     ylim=c(min(Loadings[!is.na(Loadings)]), 0.9))
axis(1, at=seq(1.5, 28*5, by=5), labels = colnames(X_pooled), las=2)
axis(1, at=c(0, Xgroups, 28)*5, line = 6, labels = NA)
axis(1, at=c(3, 10.5, 21.5)*5, labels = c("Growth Factors", "Chemokines", "Cytokines"), line=6, tick = FALSE)
abline(v=c(0, Xgroups, 28)*5, lty=3, col="black")
legend("topright", lty=1, lwd=3, col=c(MyPal[1], MyPal[2], MyPal[3]), legend = c("CLL", "DLBCL", "MM"))
dev.off()


### Stability throughout cut-offs: by subtype # Challenge 8

for (subtype in c("BCLL", "DLBL", "FL", "MM")){
  print(subtype)
  
  MyStab=NULL
  
  X=eval(parse(text=paste0("X_", subtype)))
  Y=eval(parse(text=paste0("Y_", subtype)))
  
  Stability_results=StabilityPlot(X = X, Y = Y, NIter = 100)
  
  pheatmap(Stability_results, cluster_rows = FALSE, cluster_cols = FALSE, 
           display_numbers = TRUE, filename = paste0("Figures/Supp_fig2_", subtype, ".pdf"), height=5, width=10)
}


### Time to diagnosis # Challenge 10
MyParameters$Time$sPLS$NVar

MysPLS_diag=spls(X = X_diagnostic, Y=Y_diagnostic, ncomp = 1, keepX = MyParameters$Time$sPLS$NVar)
MysPLS_diag_BCLL=spls(X = X_diag_BCLL, Y=Y_diag_BCLL, ncomp = 1,  keepX = 1)
MysPLS_diag_DLBL=spls(X = X_diag_DLBL, Y=Y_diag_DLBL, ncomp = 1,  keepX = 1)
MysPLS_diag_MM=spls(X = X_diag_MM, Y=Y_diag_MM, ncomp = 1,  keepX = 1)


### Loadings plot

Loadings=cbind(MysPLS_diag$loadings$X, MysPLS_diag_BCLL$loadings$X, 
               MysPLS_diag_DLBL$loadings$X, MysPLS_diag_MM$loadings$X, rep(NA, 28), rep(NA, 28))
Loadings=as.vector(t(Loadings))
Loadings=Loadings[-c(length(Loadings)-1, length(Loadings))]

pdf(paste0("Figures/Fig4_A.pdf"), height=5, width=10)
par(mar=c(6,5,3,3))
plot(Loadings, col=c(MyPal[4], MyPal[1], MyPal[2], MyPal[3], NA, NA), xaxt="n", ylab="Loadings Coefficients",  type="h", lwd=3, xlab="")
axis(1, at=seq(1.5, 28*6, by=6), labels = rownames(MysPLS_diag$loadings$X), las=2)
legend("topleft", lty=1, lwd=3, col=c(MyPal[4], MyPal[1], MyPal[2], MyPal[3]), 
       legend = c("Pooled cases", "CLL", "DLBCL", "MM"))
dev.off()


### Stability # Challenge 11

NVar=c(1,1,1,1)
X_diag_pooled=X_diagnostic
Y_diag_pooled=Y_diagnostic

NIter=100

for(i in 1:4){
  print(i)
  subtype=c("pooled", "DLBL", "BCLL", "MM")[i]
  X=eval(parse(text=paste0("X_diag_", subtype)))
  Y=eval(parse(text=paste0("Y_diag_", subtype)))
  
  TmpStab=NULL
  
  for (k in 1:NIter){
    s=sample(seq(1, nrow(X)), size=nrow(X), replace = TRUE) # BOOTSTRAP
    X_boot=X[s,]
    rownames(X_boot)=seq(1:nrow(X))
    Y_boot=Y[s]
    
    TmpsPLS=spls(X_boot, Y_boot, keepX=NVar[i], ncomp=1 )
    TmpStab=rbind(TmpStab, TmpsPLS$loadings$X[,1])
  }
  
  Stab=apply(TmpStab, 2, FUN=function(x){sum(x!=0)})/NIter
  assign(paste0("Stab_", subtype), Stab)
}

Stab=cbind(Stab_pooled, Stab_BCLL, Stab_DLBL, Stab_MM, rep(NA, 28), rep(NA, 28))
Stab=as.vector(t(Stab))

saveRDS(Stab, "Results/stability_ttd.rds")

pdf("Figures/Fig4_B.pdf", height=5, width=10)
par(mar=c(6,5,3,3))
plot(Stab, col=c(MyPal[4], MyPal[1], MyPal[2], MyPal[3], NA, NA), xaxt="n", ylab="Selection Proportion",  
     type="h", lwd=3, xlab="")
axis(1, at=seq(1.5, 28*6, by=6), labels = rownames(MysPLS_diag$loadings$X), las=2)
dev.off()


