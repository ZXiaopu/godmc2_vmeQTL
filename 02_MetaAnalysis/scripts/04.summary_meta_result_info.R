library(tidyr)
library(readr)
library(dplyr)

info <- read_delim("../process_reports/meta_res_check_report/meta_res_check_all_chunk.out",col_names=F, delim=" - ")
info$X2 <- as.numeric(gsub("; cpg","",info$X2))
info$Method <- rep(c("BF","DRM","SVLM"),119)

BF_sum <- info %>% filter(Method=="BF")
DRM_sum <- info %>% filter(Method=="DRM")
SVLM_sum <- info %>% filter(Method=="SVLM")

BF_sum_asso <- sum(BF_sum$X2)
BF_sum_cpg <- sum(BF_sum$X3)
DRM_sum_asso <- sum(DRM_sum$X2)
DRM_sum_cpg <- sum(DRM_sum$X3)
SVLM_sum_asso <- sum(SVLM_sum$X2)
SVLM_sum_cpg <- sum(SVLM_sum$X3)

message(paste0("association included in BF method is:", BF_sum_asso, ", DRM method is:", DRM_sum_asso, ", SVLM method is:", SVLM_sum_asso))
message(paste0("cpg included in BF method is:", BF_sum_cpg, ", DRM method is:", DRM_sum_cpg, ", SVLM method is:", SVLM_sum_cpg))

'''
association included in BF method is:12146723082, DRM method is:11730880971, SVLM method is:11730877264
cpg included in BF method is:850309, DRM method is:819472, SVLM method is:819472
'''
