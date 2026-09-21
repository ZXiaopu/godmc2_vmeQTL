library(getmstatistic)
library(gridExtra)       
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)

input_path <- "../data/Mstatistics/M_sta_input"

BFfiles <- list.files(path = input_path, pattern = "^BF.*\\_all_pval5e-8_output")
DRMfiles <- list.files(path = input_path, pattern = "^DRM.*\\_all_pval5e-8_output")
SVLMfiles <- list.files(path = input_path, pattern = "^SVLM.*\\_all_pval5e-8_output")

run_m <- function(files, method, pval, include_two_cohorts){
    df <- lapply(files, function(x) read_delim(paste0(input_path,"/",x), col_names=c("SNP","Chr","BP","A1","A2","Freq","Probe","Probe_chr","Probe_bp","Gene","Orientation","b","SE","p")) %>% 
            mutate(file = gsub(paste0(method,"_"), "", x)) %>% 
            mutate(file = gsub("_all_pval5e-8_output","", file))) %>% 
            bind_rows()
    
    o <- df %>% filter(p<pval)    # filter(SNP %in% ref_1kg$snp) %>% 
    o$SNP_Probe <- paste0(o$SNP, "_", o$Probe) 

    setDT(o)
    n_probes <- uniqueN(o[file == "bib_eur_mother", Probe])

    # Pre-filter
    ref <- o[file == "bib_eur_mother"]

    counts_all <- o[, .N, by = .(SNP, Probe, SNP_Probe)]
    counts_ref <- ref[, .N, by = .(SNP, Probe, SNP_Probe)]

    keep <- unique(c(counts_all[N >= 10, SNP_Probe],counts_ref[N >= 10, SNP_Probe]))
    o3 <- o[J(keep), on = "SNP_Probe"]

    setDT(o3)
    asso <- o3[o3[order(Probe, file != "bib_eur_mother"),.SD[1],by = Probe],on = "SNP_Probe"]

    ## bib_eur_mother cohort does not have as many association available as others
#    ds1 <- asso
#    ds1m <- getmstatistic(ds1$b, ds1$SE, ds1$SNP_Probe, ds1$file, save_dir=paste0("../data/Mstatistics/Mstat_P",pval))
#    save(ds1m, ds1, file = paste0(input_path, "/Mstat_P", pval, "/", method, "_allchr_10qtls_Mtat.RData"))
}

run_m(BFfiles, "BF", 5.8e-14, TRUE)
run_m(DRMfiles, "DRM", 5.8e-14, TRUE)
run_m(SVLMfiles, "SVLM", 5.8e-14, TRUE)

# filter the distances between CpGs and SNPs to avoid the SNPs/CpGs are correlated
setwd("../data/Mstatistics/Mstat_P5.8e-14")
library(tidyr)
library(dplyr)
library(readr)
library(meffil)

selectIndep <- function(ds1m){
    M_df <- ds1m$M_dataset
    asso_df <- M_df %>% mutate(pair = variant_names_in) %>% 
        dplyr::select(pair, variant_names_in) %>% unique() %>% 
        separate(variant_names_in, sep="_", into=c("SNP","A1","A2","CpG")) %>% 
        separate(SNP, sep=":", into=c("chr","pos"))
    
    annots <- meffil.get.features("epic") %>% filter(name %in% asso_df$CpG) %>%
        mutate(CpG=name) %>% dplyr::select(CpG, chromosome, position)
    
    asso_df1 <- merge(asso_df, annots, by.x="CpG")
    colnames(asso_df1) <- c("CpG","Pair","SNP_chr","SNP_pos","A1","A2","CpG_chr","CpG_pos")
    
    asso_df1 <- asso_df1 %>% mutate(SNP_chr = as.numeric(SNP_chr), SNP_pos = as.numeric(SNP_pos),
                                    CpG_chr = as.numeric(gsub("chr","",CpG_chr)), CpG_pos = as.numeric(CpG_pos))
    
    o <- data.frame()
    for (chr in c(1:22)){
        asso_df2 <- asso_df1 %>% filter(SNP_chr==chr)
        cpg_dist_df <- asso_df2 %>% distinct(CpG, CpG_chr, CpG_pos) %>%
            mutate(min_dist_other_CpG = min_dist_to_others(CpG_pos)) %>% ungroup() %>%
            dplyr::select(CpG, min_dist_other_CpG)
    
        snp_dist_df <- asso_df2 %>% distinct(SNP_chr, SNP_pos) %>%
            mutate(min_dist_other_SNP = min_dist_to_others(SNP_pos)) %>% ungroup()

        asso_df2 <- asso_df2 %>% left_join(cpg_dist_df, by = "CpG") %>%
            left_join(snp_dist_df, by = c("SNP_chr", "SNP_pos")) %>%
            arrange(CpG_chr, CpG_pos)
        
        o <- rbind(o, asso_df2)
    }

    return(o)
}

min_dist_to_others <- function(pos) {
  n <- length(pos)
  ord <- order(pos)
  s_pos <- pos[ord]

  left_dist <- c(Inf, diff(s_pos))
  right_dist <- c(diff(s_pos), Inf)
  min_d <- pmin(left_dist, right_dist)

  res <- numeric(n)
  res[ord] <- min_d
  return(res)
}

load("BF_allchr_10qtls_Mstat.RData")
BF_ds1 <- ds1
BF_ds1m_dist <- selectIndep(ds1m)

load("DRM_allchr_10qtls_Mstat.RData")
DRM_ds1 <- ds1
DRM_ds1m_dist <- selectIndep(ds1m)

load("SVLM_allchr_10qtls_Mstat.RData")
SVLM_ds1 <- ds1
SVLM_ds1m_dist <- selectIndep(ds1m)

BF_indep <- BF_ds1m_dist %>% filter(min_dist_other_CpG>50000 & min_dist_other_SNP>50000)
DRM_indep <- DRM_ds1m_dist %>% filter(min_dist_other_CpG>50000 & min_dist_other_SNP>50000)
SVLM_indep <- SVLM_ds1m_dist %>% filter(min_dist_other_CpG>50000 & min_dist_other_SNP>50000)

BF_ds1 <- BF_ds1 %>% filter(SNP_Probe %in% BF_indep$Pair)
DRM_ds1 <- DRM_ds1 %>% filter(SNP_Probe %in% DRM_indep$Pair)
SVLM_ds1 <- SVLM_ds1 %>% filter(SNP_Probe %in% SVLM_indep$Pair)

library(getmstatistic)
runMsta <- function(ds1, method){
    ds1m <- getmstatistic(ds1$b, ds1$SE, ds1$SNP_Probe, ds1$file, save_dir=paste0(method, "_noRepeatSNP_10qtls_22cohorts"))
    save(ds1m, ds1, file = paste0(method, "_allchr_10qtls_indep_Mtat.RData"))
}

runMsta(BF_ds1, "BF")
runMsta(DRM_ds1, "DRM")
runMsta(SVLM_ds1, "SVLM")

