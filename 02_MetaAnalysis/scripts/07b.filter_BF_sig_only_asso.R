library(tidyr)
library(readr)
library(dplyr)

args <- commandArgs(T)
chr <- args[1]
index <- args[2]

setwd(paste0("../data/Meta_Results/chr",chr))

BF_index <- read_delim(paste0("BF_random_meta_results_cpg_index",index,"1_fixed_meta_pval5e-8.txt"))
run_by_BF_only <- read_delim(paste0("Missing_asso_run_by_BF_only_index",index,".txt"))

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
BF_only <- BF_index1 %>% filter(MarkerName %in% run_by_BF_only$SNP_Probe) %>% 
            mutate(SNP_Probe = MarkerName) %>% select(-valid_indices) %>% 
            separate(MarkerName, sep="_", into=c("SNP","A1","A2","CpG")) 

message(paste0("Results detected by BF but not DRM/SVLM are: n(asso) = ",length(unique(BF_only$SNP_Probe))," n(cpg) = ", length(unique(BF_only$CpG))))
write.table(BF_only, file=paste0("Missing_asso_BF_5e-8_index",index,"_new_KORA.txt"), col=T, row=F, sep="\t", quote=F)

