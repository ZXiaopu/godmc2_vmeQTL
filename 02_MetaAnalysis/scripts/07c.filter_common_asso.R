library(tidyr)
library(readr)
library(dplyr)

args <- commandArgs(T)
chr <- args[1]
index <- args[2]

setwd(paste0("../data/Meta_Results/chr",chr))

BF_index <- read_delim(paste0("BF_random_meta_results_cpg_index",index,"1_fixed_meta_pval5e-8.txt"))
DRM_index <- read_delim(paste0("DRM_random_meta_results_cpg_index",index,"1_fixed_meta_pval5e-8.txt"))
SVLM_index <- read_delim(paste0("SVLM_random_meta_results_cpg_index",index,"1_fixed_meta_pval5e-8.txt"))

cohort_info <- read_delim("../../../GoDMC_vmeQTL_phase1_data/vQTL_results/Cohorts_with_all_results/Cohort_Information_short.csv") %>% mutate(cohort=Study) %>% select(cohort, Nsamples_3a)

calSampleSize <- function(df, method){
    lines <- readLines(paste0(method,"_random_meta_results_cpg_index",index,"1.txt.info"))
    input_lines <- lines[grepl("Input File", lines)]

    indices <- as.numeric(sub(paste0(".*Input File ([0-9]+) : ",method,"_.*"), "\\1", input_lines))
    cohorts <- sub(paste0(".*Input File [0-9]+ : ",method,"_(.*)_cpg_index.*"), "\\1", input_lines)

    df_cohorts <- data.frame(
                index = indices,
                cohort = cohorts,
                stringsAsFactors = FALSE
    )

    df_cohorts1 <- merge(df_cohorts, cohort_info, by.x="cohort")
    df_cohorts1 <- df_cohorts1[order(as.numeric(df_cohorts1$index)),]
    
    df$valid_indices <- lapply(df$Direction, function(x) {
        chars <- strsplit(x, "")[[1]]
        which(chars != "?")
    })

    df$valid_counts <- lapply(df$Direction, function(x) {
        chars <- strsplit(x, "")[[1]]
        length(which(chars != "?"))
    })
    
    df$valid_counts <- as.numeric(df$valid_counts)

    df$total_sample_size <- sapply(df$valid_indices, function(idx) {
        if (length(idx) == 0) return(0)
        sum(df_cohorts1$Nsamples_3a[idx], na.rm = TRUE)
    })

    return(df)
}

BF_index1 <- calSampleSize(BF_index, "BF") %>% filter(valid_counts>=3 & total_sample_size>=10000)
DRM_index1 <- calSampleSize(DRM_index, "DRM") %>% filter(valid_counts>=3 & total_sample_size>=10000)
SVLM_index1 <- calSampleSize(SVLM_index, "SVLM") %>% filter(valid_counts>=3 & total_sample_size>=10000)

BF_bonf <- BF_index1 %>% filter(Pvalue<5.3e-10)
DRM_bonf <- DRM_index1 %>% filter(Pvalue<5.3e-10)
SVLM_bonf <- SVLM_index1 %>% filter(Pvalue<5.3e-10)

BF_lead <- BF_bonf %>% filter(MarkerName %in% DRM_index1$MarkerName & MarkerName %in% SVLM_index1$MarkerName)
DRM_lead <- DRM_bonf %>% filter(MarkerName %in% BF_index1$MarkerName & MarkerName %in% SVLM_index1$MarkerName)
SVLM_lead <- SVLM_bonf %>% filter(MarkerName %in% DRM_index1$MarkerName & MarkerName %in% BF_index1$MarkerName)

asso <- c(BF_lead$MarkerName, DRM_lead$MarkerName, SVLM_lead$MarkerName) %>% unique()
BF_o <- BF_index1 %>% filter(MarkerName %in% asso) %>% 
            mutate(SNP_Probe = MarkerName, method="BF") %>% select(-valid_indices) %>% 
            separate(MarkerName, sep="_", into=c("SNP","A1","A2","CpG"))
DRM_o <- DRM_index1 %>% filter(MarkerName %in% asso) %>% 
            mutate(SNP_Probe = MarkerName, method="DRM") %>% select(-valid_indices) %>%
            separate(MarkerName, sep="_", into=c("SNP","A1","A2","CpG"))
SVLM_o <- SVLM_index1 %>% filter(MarkerName %in% asso) %>% 
            mutate(SNP_Probe = MarkerName, method="SVLM") %>% select(-valid_indices) %>%
            separate(MarkerName, sep="_", into=c("SNP","A1","A2","CpG"))

combined_out <- rbind(BF_o, DRM_o, SVLM_o)
message(paste0("For results bonferroni significance in at least one method: n(asso) = ",length(unique(combined_out$SNP_Probe))," n(cpg) = ", length(unique(combined_out$CpG))))
write.table(combined_out, paste0("One_method_5.3e-10_two_methods_5e-8_chr",chr,"_index_",index,"_cohort3_sample10k_new_KORA_output.txt_sep21"), col=T, row=F, sep="\t", quote=F)
