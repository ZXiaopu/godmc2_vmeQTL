library(metafor)
library(tidyr)
library(dplyr)
library(readr)
library(getmstatistic)

setwd("../data/Mstatistics/Mstat_P5.8e-14")
load("BF_noRepeatSNP_10qtls_22cohorts/BF_noRepeatSNP_10qtls_22cohorts.RData")
ds1_BF <- BF_ds2 %>% filter(file %in% c("Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_BF %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_BF <- ds1_BF %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_BF <- getmstatistic(ds1_BF$b, ds1_BF$SE, ds1_BF$SNP_Probe, ds1_BF$file)
save(ds1_BF, ds1m_BF, file = "BF_noRepeatSNP_10qtls_22cohorts/BF_noRepeatSNP_10qtls_21cohorts_noDHW.RData")

load("DRM_noRepeatSNP_10qtls_22cohorts/DRM_noRepeatSNP_10qtls_22cohorts.RData")
ds1_DRM <- DRM_ds2 %>% filter(file %in% c("Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_DRM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_DRM <- ds1_DRM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_DRM <- getmstatistic(ds1_DRM$b, ds1_DRM$SE, ds1_DRM$SNP_Probe, ds1_DRM$file)
save(ds1_DRM, ds1m_DRM, file = "DRM_noRepeatSNP_10qtls_22cohorts/DRM_noRepeatSNP_10qtls_21cohorts_noDHW.RData")

load("SVLM_noRepeatSNP_10qtls_22cohorts/SVLM_noRepeatSNP_10qtls_22cohorts.RData")
ds1_SVLM <- SVLM_ds2 %>% filter(file %in% c("Dutch_Hunger_Winter_Families_Study") == F)
asso <- ds1_SVLM %>% group_by(SNP_Probe) %>% tally() %>% filter(n>=10)
ds1_SVLM <- ds1_SVLM %>% filter(SNP_Probe %in% asso$SNP_Probe)
ds1m_SVLM <- getmstatistic(ds1_SVLM$b, ds1_SVLM$SE, ds1_SVLM$SNP_Probe, ds1_SVLM$file)
save(ds1_SVLM, ds1m_SVLM, file = "SVLM_noRepeatSNP_10qtls_22cohorts/SVLM_noRepeatSNP_10qtls_21cohorts_noDHW.RData")

c <- read_delim("../covariates_mregression.csv") %>% select(-DNAm_Array)

M_regression_1 <- function(mstat){
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

o_redo1 <- rbind(M_regression_1(ds1m_BF$M_dataset) %>% mutate(method="BF"),
           M_regression_1(ds1m_DRM$M_dataset) %>% mutate(method="DRM"),
           M_regression_1(ds1m_SVLM$M_dataset) %>% mutate(method="SVLM"))
write.table(o_redo1, "Mregression_21cohorts_noDHW.csv", col=T, row=F, sep=",", quote=F)

