library(metafor)
library(tidyr)
library(dplyr)
library(readr)

setwd("../data/Mstat_P5.8e-14")

load("BF_allchr_10qtls_noDHW_Mstat.RData")
BF_ds1m_21c <- ds1m
BF_ds1 <- ds1
BF_rem <- BF_ds1 %>% group_by(SNP) %>% tally() %>% filter(n>=20)
BF_ds2 <- BF_ds1 %>% filter((SNP %in% BF_rem$SNP) == F)
ds2m_BF <- getmstatistic(BF_ds2$b, BF_ds2$SE, BF_ds2$SNP_Probe, BF_ds2$file, save_dir=paste0("BF_noRepeatSNP_10qtls_20cohorts"))
save(ds2m_BF, BF_ds2, file = paste0("BF_noRepeatSNP_10qtls_20cohorts/BF_noRepeatSNP_10qtls_20cohorts.RData"))

load("DRM_allchr_10qtls_noDHW_Mstat.RData")
DRM_ds1m_21c <- ds1m
DRM_ds1 <- ds1
DRM_rem <- DRM_ds1 %>% group_by(SNP) %>% tally() %>% filter(n>=20)
DRM_ds2 <- DRM_ds1 %>% filter((SNP %in% DRM_rem$SNP) == F)
ds2m_DRM <- getmstatistic(DRM_ds2$b, DRM_ds2$SE, DRM_ds2$SNP_Probe, DRM_ds2$file, save_dir=paste0("DRM_noRepeatSNP_10qtls_20cohorts"))
save(ds2m_DRM, DRM_ds2, file = paste0("DRM_noRepeatSNP_10qtls_20cohorts/DRM_noRepeatSNP_10qtls_20cohorts.RData"))

load("SVLM_allchr_10qtls_noDHW_Mstat.RData")
SVLM_ds1m_21c <- ds1m
SVLM_ds1 <- ds1
SVLM_rem <- SVLM_ds1 %>% group_by(SNP) %>% tally() %>% filter(n>=20)
SVLM_ds2 <- SVLM_ds1 %>% filter((SNP %in% SVLM_rem$SNP) == F)
ds2m_SVLM <- getmstatistic(SVLM_ds2$b, SVLM_ds2$SE, SVLM_ds2$SNP_Probe, SVLM_ds2$file, save_dir=paste0("SVLM_noRepeatSNP_10qtls_20cohorts"))
save(ds2m_SVLM, SVLM_ds2, file = paste0("SVLM_noRepeatSNP_10qtls_20cohorts/SVLM_noRepeatSNP_10qtls_20cohorts.RData"))

c <- read_delim("../covariates_mregression.csv") %>% select(-DNAm_Array)

M_regression <- function(mstat){
  tmp <- mstat %>% select(study_names_in, M, M_se) %>% unique()
  ds_df1 <- merge(tmp, c, by.x="study_names_in")
  vs <- colnames(c)[2:ncol(c)]
  results <- lapply(vs, function(x) {
    print(x)
    mod <- rma(yi = M, vi = M_se^2, mod = as.formula(paste("~", x)), data = ds_df1)
    data.frame(variable = x, pval = mod$pval[2])
  })
  out <- do.call(rbind, results)
  return(out)
}

o_21c <- rbind(M_regression(ds2m_BF$M_dataset) %>% mutate(method="BF"),
           M_regression(ds2m_DRM$M_dataset) %>% mutate(method="DRM"),
           M_regression(ds2m_SVLM$M_dataset) %>% mutate(method="SVLM"))
write.table(o_21c, "Mregression_20cohorts_10qtls_noRepeatSNPs_pval5.8e-14.csv", col=T, row=F, sep=",", quote=F)

#o_20c <- rbind(M_regression(ds2m_BF$M_dataset %>% filter(study_names_in != "Dutch_Hunger_Winter_Families_Study")) %>% mutate(method="BF"),
#           M_regression(ds2m_DRM$M_dataset %>% filter(study_names_in != "Dutch_Hunger_Winter_Families_Study")) %>% mutate(method="DRM"),
#           M_regression(ds2m_SVLM$M_dataset %>% filter(study_names_in != "Dutch_Hunger_Winter_Families_Study")) %>% mutate(method="SVLM"))
#write.table(o_20c, "Mregression_20cohorts_nobib_15qtls_pval5.8e-14.csv", col=T, row=F, sep=",", quote=F)

#DRM <- DRM_ds1m_21c$M_dataset %>% select(study_names_in, M, M_se) %>% unique()
#DRM_df1 <- merge(DRM, c, by.x="study_names_in")
#rma(yi = M, vi = M_se^2, mod = M ~ DNAm_Array + Relatedness + is_Netherlands + Imputation_reference_panel, data = DRM_df1) # pval - DNAmArray: 0.03; imputational panel: 0.005

#SVLM <- SVLM_ds1m_21c$M_dataset %>% select(study_names_in, M, M_se) %>% unique()
#SVLM_df1 <- merge(SVLM, c, by.x="study_names_in")
#rma(yi = M, vi = M_se^2, mod = M ~ DNAm_Array + Relatedness + is_Netherlands + Imputation_reference_panel, data = SVLM_df1) # pval - DNAmArray: 0.045; imputational panel: 0.0015

#M_regression_noDHW <- function(mstat){
#  tmp <- mstat %>% select(study_names_in, M, M_se) %>% unique()
#  ds_df1 <- merge(tmp, c, by.x="study_names_in") %>% select(-DNAm_Array)
#  vs <- colnames(ds_df1)[4:28]
#  results <- lapply(vs, function(x) {
#    print(x)
#    mod <- rma(yi = M, vi = M_se^2, mod = as.formula(paste("~", x)), data = ds_df1)
#    data.frame(variable = x, pval = mod$pval[2])
#  })
#  out <- do.call(rbind, results)
#  return(out)
#}

#o_20c_1 <- rbind(M_regression_noDHW(BF_ds1m_21c$M_dataset %>% filter((study_names_in %in% c("Dutch_Hunger_Winter_Families_Study")==F))) %>% mutate(method="BF"),
#           M_regression_noDHW(DRM_ds1m_21c$M_dataset %>% filter((study_names_in %in% c("Dutch_Hunger_Winter_Families_Study")==F))) %>% mutate(method="DRM"),
#           M_regression_noDHW(SVLM_ds1m_21c$M_dataset %>% filter((study_names_in %in% c("Dutch_Hunger_Winter_Families_Study")==F))) %>% mutate(method="SVLM"))
#write.table(o_20c_1, "Mregression_20cohorts_noDHW_15qtls_pval5.8e-14.csv", col=T, row=F, sep=",", quote=F)

## rerun M stat
library(getmstatistic)
load("BF_chr1_1kg_Mstat.RData")
ds1_BF <- ds1 %>% filter(file != "GS_unrelated")
asso <- ds1_BF %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_BF <- ds1_BF %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_BF <- getmstatistic(ds1_BF$b, ds1_BF$SE, ds1_BF$SNP_Probe, ds1_BF$file, save_dir=paste0("/BF", "_noGS"))
save(ds1_BF, ds1m_BF, file = "BF_chr1_1kg_noGS_Mstat.RData")

load("DRM_chr1_1kg_Mstat.RData")
ds1_DRM <- ds1 %>% filter(file != "Dutch_Hunger_Winter_Families_Study")
asso <- ds1_DRM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_DRM <- ds1_DRM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_DRM <- getmstatistic(ds1_DRM$b, ds1_DRM$SE, ds1_DRM$SNP_Probe, ds1_DRM$file, save_dir=paste0("/", "DRM", "_noDHW"))
save(ds1_DRM, ds1m_DRM, file = "DRM_chr1_1kg_noDHW_Mstat.RData")

load("SVLM_chr1_1kg_Mstat.RData")
ds1_SVLM <- ds1 %>% filter(file != "Dutch_Hunger_Winter_Families_Study")
asso <- ds1_SVLM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_SVLM <- ds1_SVLM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_SVLM <- getmstatistic(ds1_SVLM$b, ds1_SVLM$SE, ds1_SVLM$SNP_Probe, ds1_SVLM$file, save_dir=paste0("/", "SVLM", "_noDHW"))

M_regression_BF <- function(mstat){
  tmp <- mstat %>% select(study_names_in, M, M_se) %>% unique()
  ds_df1 <- merge(tmp, c, by.x="study_names_in") %>% select(-is_Scotland)
  vs <- colnames(ds_df1)[4:28]
  results <- lapply(vs, function(x) {
    print(x)
    mod <- rma(yi = M, vi = M_se^2, mod = as.formula(paste("~", x)), data = ds_df1)
    data.frame(variable = x, pval = mod$pval[2])
  })
  out <- do.call(rbind, results)
  return(out)
}

M_regression_DRM_SVLM <- function(mstat){
  tmp <- mstat %>% select(study_names_in, M, M_se) %>% unique()
  ds_df1 <- merge(tmp, c, by.x="study_names_in") %>% select(-DNAm_Array)
  vs <- colnames(ds_df1)[4:28]
  results <- lapply(vs, function(x) {
    print(x)
    mod <- rma(yi = M, vi = M_se^2, mod = as.formula(paste("~", x)), data = ds_df1)
    data.frame(variable = x, pval = mod$pval[2])
  })
  out <- do.call(rbind, results)
  return(out)
}

o_redo <- rbind(M_regression_BF(ds1m_BF$M_dataset) %>% mutate(method="BF"),
           M_regression_DRM_SVLM(ds1m_DRM$M_dataset) %>% mutate(method="DRM"),
           M_regression_DRM_SVLM(ds1m_SVLM$M_dataset) %>% mutate(method="SVLM"))
write.table(o_redo, "Mregression_20cohorts_BFnoGS_DRMnoDHW_SVLMnoDHW.csv", col=T, row=F, sep=",", quote=F)

## rerun M stat without GS and DHW
library(getmstatistic)
load("BF_chr1_1kg_Mstat.RData")
ds1_BF <- ds1 %>% filter(file %in% c("GS_unrelated","Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_BF %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_BF <- ds1_BF %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_BF <- getmstatistic(ds1_BF$b, ds1_BF$SE, ds1_BF$SNP_Probe, ds1_BF$file)
save(ds1_BF, ds1m_BF, file = "BF_chr1_1kg_noGS_noDHW_Mstat.RData")

load("DRM_chr1_1kg_Mstat.RData")
ds1_DRM <- ds1 %>% filter(file %in% c("GS_unrelated","Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_DRM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_DRM <- ds1_DRM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_DRM <- getmstatistic(ds1_DRM$b, ds1_DRM$SE, ds1_DRM$SNP_Probe, ds1_DRM$file, save_dir=paste0("/", "DRM", "_noDHW"))
save(ds1_DRM, ds1m_DRM, file = "DRM_chr1_1kg_noGS_noDHW_Mstat.RData")

load("SVLM_chr1_1kg_Mstat.RData")
ds1_SVLM <- ds1 %>% filter(file %in% c("GS_unrelated","Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_SVLM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_SVLM <- ds1_SVLM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_SVLM <- getmstatistic(ds1_SVLM$b, ds1_SVLM$SE, ds1_SVLM$SNP_Probe, ds1_SVLM$file, save_dir=paste0("/", "SVLM", "_noDHW"))
save(ds1_SVLM, ds1m_SVLM, file = "SVLM_chr1_1kg_noGS_noDHW_Mstat.RData")

M_regression_1 <- function(mstat){
  tmp <- mstat %>% select(study_names_in, M, M_se) %>% unique()
  ds_df1 <- merge(tmp, c, by.x="study_names_in") %>% select(-DNAm_Array, -is_Scotland)
  vs <- colnames(ds_df1)[4:27]
  results <- lapply(vs, function(x) {
    print(x)
    mod <- rma(yi = M, vi = M_se^2, mod = as.formula(paste("~", x)), data = ds_df1)
    data.frame(variable = x, pval = mod$pval[2])
  })
  out <- do.call(rbind, results)
  return(out)
}

o_redo1 <- rbind(M_regression_1(ds1m_BF$M_dataset) %>% mutate(method="BF"),
           M_regression_1(ds1m_DRM$M_dataset) %>% mutate(method="DRM"),
           M_regression_1(ds1m_SVLM$M_dataset) %>% mutate(method="SVLM"))
write.table(o_redo1, "Mregression_19cohorts_noGS_noDHW.csv", col=T, row=F, sep=",", quote=F)

