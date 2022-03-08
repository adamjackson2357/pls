rm(list=ls())
dev.off()
path=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)
setwd("../")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tableone))
suppressPackageStartupMessages(library(lme4))

#### Descriptive analyses and tables ####

covars = readRDS("Data/Covariates.rds")
proteins = readRDS("Data/Proteins.rds")
dim(proteins)
dim(covars)
print(all(rownames(covars) == rownames(proteins)))

# convert to cases and controls
covars$type <- factor(ifelse(covars$type == 0, "control", "case"), levels = c("control", "case"))

# vector of the variable names
var_names <- c("cohort", "phase", "sex", "age", "Imputed_alcohol", "Imputed_bmi",
               "Imputed_smoking_status", "Imputed_education", "Imputed_activity",
               "LY_subtype")

# vector of the categorical variable names
cat_vars <- c("cohort", "phase", "sex",
              "Imputed_smoking_status", "Imputed_education", "Imputed_activity",
              "LY_subtype")

# use table one function
table_one <- CreateTableOne(vars = var_names,
                            strata = "type",
                            data = covars,
                            factorVars = cat_vars)
print(table_one, showAllLevels = TRUE, formatOptions = list(big.mark = ","))

# get the case and control proteins
case_proteins <- proteins[covars$type == "case",]
control_proteins <- proteins[covars$type == "control",]

# plot all the protein distributions
par(mfrow=c(3, 3))
for (i in names(proteins)){
  plot(density(case_proteins[, i]), main = paste("Cases", i))
  plot(density(control_proteins[, i]), main = paste("Controls", i))
  plot(density(proteins[, i]), main = paste("Full Sample", i))
}

# lots of instances of binomial distributions

#### Linear Mixed Models ####


# create a linear mixed models function
non_adjusted_lrt <- function(X, covars) {
  model0 <- lmer(X ~  age + sex + cohort + phase + Imputed_bmi +
                   Imputed_education + Imputed_activity + Imputed_smoking_status +
                   Imputed_alcohol + (1 | plate), data = covars, REML = FALSE)
  
  model1 <- lmer(X ~  age + sex + cohort + phase + Imputed_bmi +
                   Imputed_education + Imputed_activity + Imputed_smoking_status +
                   Imputed_alcohol + (1 | plate) + type, data = covars, REML = FALSE)
  
  summary(model1)$coefficients
  
  res = c(summary(model1)$coefficients["typecase", 1:2],
          anova(model0, model1)$'Pr(>Chisq)'[2],
          silent = TRUE)
  names(res) = c("coef", "coef.se", "pval")
  res = round(res, 2)
  
  return(res)
}

# apply across all proteins
all_lmm = t(apply(proteins, 2, FUN = non_adjusted_lrt, covars))
all_lmm

#### Adjust by WBC counts ####

# get wbc data
wbc <- readRDS("Data/wbc.rds")

# join the covars and wbc
covars_wbc <- merge(covars, wbc, by = 0)
rownames(covars_wbc) <- covars_wbc$Row.names
covars_wbc <- covars_wbc[,-1]

# subset the proteins as well
proteins_wbc <- subset(proteins, rownames(proteins) %in% rownames(wbc))

head(covars_wbc)
head(wbc)
dim(proteins_wbc)

# create a linear mixed models function with 
wbc_adjusted_lrt <- function(X, covars) {
  model0 <- lmer(X ~  age + sex + cohort + phase + Imputed_bmi +
                   Imputed_education + Imputed_activity + Imputed_smoking_status +
                   Imputed_alcohol + (1 | plate) +
                   CD8 + CD4 + NK + Bcells + Mono, data = covars, REML = FALSE)
  
  model1 <- lmer(X ~  age + sex + cohort + phase + Imputed_bmi +
                   Imputed_education + Imputed_activity + Imputed_smoking_status +
                   Imputed_alcohol + (1 | plate) + 
                   CD8 + CD4 + NK + Bcells + Mono + type, data = covars, REML = FALSE)
  
  summary(model1)$coefficients
  
  res = c(summary(model1)$coefficients["typecase", 1:2],
          anova(model0, model1)$'Pr(>Chisq)'[2],
          silent = TRUE)
  names(res) = c("coef", "coef.se", "pval")
  res = round(res, 2)
  
  return(res)
}

# run when including wbc
wbc_lmm = t(apply(proteins_wbc, 2, FUN = wbc_adjusted_lrt, covars_wbc))
wbc_lmm

# create plot
dim(all_lmm)[1]
not_adjusted = all_lmm[,"pval"]
not_adjusted <- cbind(p_val = -log10(all_lmm[,"pval"]), model = rep("not_adjusted", dim(all_lmm)[1]))
wbc_adjusted <- cbind(p_val = -log10(wbc_lmm[,"pval"]), model = rep("wbc_adjusted", dim(all_lmm)[1]))

p_vals <- data.frame((rbind(not_adjusted, wbc_adjusted)))
p_vals$p_val <- as.numeric(p_vals$p_val)
p_vals

ggplot(data = p_vals, aes(x = rownames(p_vals), y = p_val)) +
  geom_point()


#### Stratification Analysis ####

table(covars$LY_subtype)
table_one

# function to run the stratified log likelihood analysis
strat_lrt <- function(covars, proteins, subtype, lrt){
  subtype_egm <- subset(covars, LY_subtype == subtype)$egm_set
  subtype_covars <- data.frame(subset(covars, egm_set %in% subtype_egm))
  subtype_proteins <- data.frame(subset(proteins, rownames(proteins) %in% rownames(subtype_covars)))
  
  subtype_lmm <- t(apply(subtype_proteins, 2, FUN = lrt, subtype_covars))
  return(subtype_lmm)
}

# run for each subtype for the non-wbc-adjusted data
strat_lrt(covars, proteins, "BCLL", non_adjusted_lrt)
strat_lrt(covars, proteins, "DLBL", non_adjusted_lrt)
strat_lrt(covars, proteins, "FL", non_adjusted_lrt)
strat_lrt(covars, proteins, "MM", non_adjusted_lrt)

# run for each subtype for the wbc adjusted data
strat_lrt(covars_wbc, proteins_wbc, "BCLL", wbc_adjusted_lrt)
strat_lrt(covars_wbc, proteins_wbc, "DLBL", wbc_adjusted_lrt)
strat_lrt(covars_wbc, proteins_wbc, "FL", wbc_adjusted_lrt)
strat_lrt(covars_wbc, proteins_wbc, "MM", wbc_adjusted_lrt)