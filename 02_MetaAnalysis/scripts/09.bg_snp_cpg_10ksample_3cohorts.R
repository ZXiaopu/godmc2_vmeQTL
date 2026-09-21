library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(data.table)
library(meffil)

files <- list.files(path="../../../GoDMC_vmeQTL_phase1_data/module1-3", recursive=TRUE)

setwd("../../../GoDMC_vmeQTL_phase1_data/module1-3")

genotype_files <- files[grep("results/02/data.allele_codes.gz",files)]
cpg_files <- files[grep("results/03/methylation_summary.RData",files)]

load_cpg <- function(x){
    load(x)
    df <- meth_summary %>% mutate(filename=x)
    return(df)
}

cpg_df <- map_df(cpg_files, function(x) load_cpg(x))
cpg_df$cohort <- gsub("/results/03/methylation_summary.RData","",cpg_df$filename)
cpg_df$n_sample <- cpg_df$`outlier.outliers.lower` + cpg_df$`outlier.outliers.upper` + cpg_df$`outlier.n`
n_cohort <- cpg_df %>% group_by(cpg) %>% tally()
cpg_remove1 <- n_cohort %>% filter(n<3)
cpg_sub <- n_cohort %>% filter(n>=3) %>% filter(n<15)

setDT(cpg_df)
cpg_remove2 <- cpg_df[cpg %in% cpg_sub$cpg,
                      .(total_n = sum(n_sample, na.rm = TRUE)),
                      by = cpg][total_n < 10000, cpg]

cpg_keep <- cpg_df %>% filter(!(cpg %in% c(cpg_remove1$cpg, cpg_remove2)))
annots <- meffil.get.features("epic") %>% filter(chromosome %in% paste0("chr",c(1:22)))
cpg_keep_autosomal <- cpg_keep %>% filter(cpg %in% annots$name)
df <- data.frame(unique(cpg_keep_autosomal$cpg))
write.table(df, file="../../../godmc2_vmeQTL/02_MetaAnalysis/data/Background_cpgs_cohort3_sample10k.txt", col=T, row=F, sep="\t", quote=F)

samplesize <- read_delim("../GoDMC_vmeQTL_phase1_data/vQTL_results/Cohorts_with_all_results/Cohort_Information_short.csv") %>% mutate(cohort=Study, n_sample=Nsamples_3a) %>% dplyr::select(cohort, n_sample)
genotype_df <- map_df(genotype_files, function(x) read_delim(x) %>% mutate(filename=x))
genotype_df$cohort <- gsub("/results/02/data.allele_codes.gz","",genotype_df$filename)
n_cohort1 <- genotype_df %>% group_by(SNP) %>% tally()
snp_remove1 <- n_cohort1 %>% filter(n<3)
snp_sub <- n_cohort1 %>% filter(n>=3 & n<15)

setDT(genotype_df)
setDT(samplesize)
genotype_df1 <- merge(genotype_df, samplesize, by = "cohort")
snp_remove2 <- genotype_df1[SNP %in% snp_sub$SNP,
                      .(total_n = sum(n_sample, na.rm = TRUE)),
                      by = SNP][total_n < 10000, SNP]

genotype_keep <- genotype_df1 %>% filter(!(SNP %in% c(snp_remove1$SNP, snp_remove2)))
genotype_keep_autosomal <- genotype_keep[!grep("23:",genotype_keep$SNP)]
df1 <- data.frame(unique(genotype_keep_autosomal$SNP))
write.table(df1, file="../../../godmc2_vmeQTL/02_MetaAnalysis/data/Background_snps_cohort3_sample10k.txt", col=T, row=F, sep="\t", quote=F)

