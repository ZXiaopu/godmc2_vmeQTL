library(tidyr)
library(dplyr)
library(readr)

arguments <- commandArgs(T)
input_path <- arguments[1]
cohort <- arguments[2]
g_chunk <- arguments[3]
cpg_chunk <- arguments[4]
home_path <- arguments[5]
samplesize <- arguments[6]

setwd(input_path)

main <- function(){
    file_BF <- paste0("vQTL_BF_cis_genetic",g_chunk,"_cpg",cpg_chunk,"_pval1_output")
    file_DRM <- paste0("vQTL_drm_cis_genetic",g_chunk,"_cpg",cpg_chunk,"_1_1_pval1_output")
    file_SVLM <- paste0("vQTL_svlm_cis_genetic",g_chunk,"_cpg",cpg_chunk,"_1_1_pval1_output")

    BFres <- read_delim(file_BF)
    DRMres <- read_delim(file_DRM)
    SVLMres <- read_delim(file_SVLM)
    
    cpgref <- read_delim(paste0(home_path,"/data/EPIC_probe_ref_file_index.txt")) %>% filter(ID %in% BFres$Probe)
    
    for (i in unique(cpgref$file_i)){
        cpgs <- cpgref %>% filter(file_i == i)
        BFtmp <- BFres %>% filter(Probe %in% cpgs$ID) %>% mutate(SNP_Probe = paste(SNP,Probe,sep="_"), N = samplesize)
        DRMtmp <- DRMres %>% filter(Probe %in% cpgs$ID) %>% mutate(SNP_Probe = paste(SNP,Probe,sep="_"), N = samplesize)
        SVLMtmp <- SVLMres %>% filter(Probe %in% cpgs$ID) %>% mutate(SNP_Probe = paste(SNP,Probe,sep="_"), N = samplesize)

        write.table(BFtmp, paste0(home_path,"/data/BF_",cohort,"_gchunk_",g_chunk,"_cpg_chunk_",cpg_chunk,"_cpg_index",i,"_meta_input.txt"), col=F, row=F, sep="\t", quote=F)
        write.table(DRMtmp, paste0(home_path,"/data/DRM_",cohort,"_gchunk_",g_chunk,"_cpg_chunk_",cpg_chunk,"_cpg_index",i,"_meta_input.txt"), col=F, row=F, sep="\t", quote=F)
        write.table(SVLMtmp, paste0(home_path,"/data/SVLM_",cohort,"_gchunk_",g_chunk,"_cpg_chunk_",cpg_chunk,"_cpg_index",i,"_meta_input.txt"), col=F, row=F, sep="\t", quote=F)
}
}
main()

#df_BF <- input_files_BF %>% lapply(read_delim) %>% bind_rows() %>%
#  select(SNP, A1, A2, Freq, Probe, b, SE, p)
#combined_df_drm <- input_files_drm %>% lapply(read_delim) %>% bind_rows() %>%
#  select(SNP, A1, A2, Freq, Probe, b, SE, p)
#combined_df_svlm <- input_files_svlm %>% lapply(read_delim) %>% bind_rows() %>%
#  select(SNP, A1, A2, Freq, Probe, b, SE, p)

#write.table(combined_df_BF, file = paste0(input_path, "/vQTL_cis_BF_results_summary.txt"), sep="\t", quote=F, col=T, row=F)
#write.table(combined_df_drm, file = paste0(input_path, "/vQTL_cis_drm_results_summary.txt"), sep="\t", quote=F, col=T, row=F)
#write.table(combined_df_svlm, file = paste0(input_path, "/vQTL_cis_svlm_results_summary.txt"), sep="\t", quote=F, col=T, row=F)

#count_results <- function(cohort, threshold, BFres, DRMres, SVLMres){
#  BFres_thres <- BFres %>% filter(p<threshold) %>% mutate(method = "BF")
#  DRMres_thres <- DRMres %>% filter(p<threshold) %>% mutate(method = "DRM")
#  SVLMres_thres <- SVLMres %>% filter(p<threshold) %>% mutate(method = "SVLM")
#  temp <- rbind(BFres_thres, DRMres_thres, SVLMres_thres) %>%
#    group_by(SNP, Probe,method) %>% tally() %>% group_by(SNP, Probe) %>% tally()
#  
#  BFres_thres1 <- BFres %>%
#    group_by(Probe) %>% tally() %>% mutate(method = "BF")
#  DRMres_thres1 <- DRMres %>% filter(p<threshold) %>%
#    group_by(Probe) %>% tally() %>% mutate(method = "DRM")
#  SVLMres_thres1 <- SVLMres %>% filter(p<threshold) %>%
#    group_by(Probe) %>% tally() %>% mutate(method = "SVLM")
#
#  temp1 <- rbind(BFres_thres1, DRMres_thres1, SVLMres_thres1) %>%
#    group_by(Probe,method) %>% tally() %>% group_by(Probe) %>% tally()
#
#  sum <- data.frame(Cohort = cohort,
#                    threshold = threshold,
#                    asso_BF = nrow(BFres_thres),
#                    asso_DRM = nrow(DRMres_thres),
#                    asso_SVLM = nrow(SVLMres_thres),
#                    asso_three_methods = temp %>% filter(n==3) %>% nrow(),
#                    asso_two_methods = temp %>% filter(n==2) %>% nrow(),
#                    asso_one_method = temp %>% filter(n==1) %>% nrow(),
#                    vCpG_BF = nrow(BFres_thres1),
#                    vCpG_DRM = nrow(DRMres_thres1),
#                    vCpG_SVLM = nrow(SVLMres_thres1),
#                    vCpG_three_methods = temp1 %>% filter(n==3) 

#out <- rbind(count_results(cohort, 0.05, BFres, DRMres, SVLMres),
#                            count_results(cohort, 1e-3, BFres, DRMres, SVLMres),
#                            count_results(cohort, 1e-4, BFres, DRMres, SVLMres),
#                            count_results(cohort, 1e-5, BFres, DRMres, SVLMres),
#                            count_results(cohort, 1e-6, BFres, DRMres, SVLMres),
