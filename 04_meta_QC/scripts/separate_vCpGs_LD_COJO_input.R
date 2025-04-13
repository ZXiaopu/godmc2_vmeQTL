library(tidyr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)
chr = args[1]
pair_list = args[2]
BFres = args[3]
DRMres = args[4]
SVLMres = args[5]
outdir = args[6]

#chr="22"
#godmc_path="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC"
#pair_list=paste0(godmc_path, "/data/vmeQTL_vCpG_pairs_1e-5_three_methods_chr22.txt")
#BFres="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/03_MetaAnalysis/meta_input/chr22/vQTL_cis_BF_results_summary_birthcohort1946_chr22.txt"
#DRMres="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/03_MetaAnalysis/meta_input/chr22/vQTL_cis_drm_results_summary_birthcohort1946_chr22.txt"
#SVLMres="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/03_MetaAnalysis/meta_input/chr22/vQTL_cis_svlm_results_summary_birthcohort1946_chr22.txt"
#outdir="vmeQTL_meta_QC/LDinput"

#pairs = read_delim(pair_list) %>% mutate(pair=paste(SNP, CpG, sep="_"))
#BF_res <- read_delim(BFres) %>% filter(SNP_CpG_pair %in% pairs$pair)
#DRM_res <- read_delim(DRMres) %>% filter(SNP_CpG_pair %in% pairs$pair)
#SVLM_res <- read_delim(SVLMres) %>% filter(SNP_CpG_pair %in% pairs$pair)

separate_cpgs <- function(res, chr, method){
    cpg_list <- unique(res$Probe)
    write.table(cpg_list, file = paste0(outdir, "/chr", chr, "/cpglist_", method))
    for (cpg in cpg_list){
        tmp <- res %>% filter(Probe == cpg) %>%
                    mutate(SNP1 = SNP, P = p) %>% 
                    separate(SNP1, sep="_", into=c("chr:pos","a1","a2")) %>%
                    separate(`chr:pos`, sep=":", into=c("CHR","BP")) %>%
                    select(CHR, SNP, BP, P)
        write.table(tmp, file = paste0(outdir, "/chr", chr, "/", cpg, "_", method, ".LDinput"), 
                    col=T, row=F, sep=" ", quote=F)

        tmp1 <- res %>% filter(Probe == cpg) %>% mutate(freq=Freq, se=SE) %>%
                    select(SNP, A1, A2, freq, b, se, p, N)
        write.table(tmp1, file = paste0(outdir, "/chr", chr, "/", cpg, "_", method, ".cojo.ma"),
                    col=T, row=F, sep=" ", quote=F)    
    }
}


pairs = read_delim(pair_list) %>% mutate(pair=paste(SNP, CpG, sep="_"))
BF_res <- read_delim(BFres) %>% filter(SNP_CpG_pair %in% pairs$pair)
DRM_res <- read_delim(DRMres) %>% filter(SNP_CpG_pair %in% pairs$pair)
SVLM_res <- read_delim(SVLMres) %>% filter(SNP_CpG_pair %in% pairs$pair)

separate_cpgs(BF_res, chr, "BF")
separate_cpgs(DRM_res, chr, "DRM")
separate_cpgs(SVLM_res, chr, "SVLM")

