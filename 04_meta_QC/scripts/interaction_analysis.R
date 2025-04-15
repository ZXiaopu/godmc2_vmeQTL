library(tidyr)
library(dplyr)
library(readr)
args <- commandArgs(trailingOnly=TRUE)

#methylation_file <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/GoDMC_NSHD99/godmc_phase2/processed_data/methylation_data/untransformed_methylation_adjusted_pcs.RData"
#genotype_file <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/genetic_data/vmeQTL.chr21.filtered.traw"
#pair_file <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/adjust_mean_effect/vmeQTL_vCpG_adjust_mean_effect_chr21.txt"
#cov_file <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/GoDMC_NSHD99/godmc_phase2/processed_data/methylation_data/all_covariates.txt"
#out_dir <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/interaction_analysis"
#chr <- "21"

methylation_file <- args[1]
genotype_file <- args[2]
pair_file <- args[3]
cov_file <- args[4]
out_dir <- args[5]
chr <- args[6]

load(methylation_file)
geno <- read_delim(genotype_file)
pairs <- read_delim(pair_file)
pairs_uniq <- pairs %>% group_by(SNP,cpg) %>% tally()
cov <- read_delim(cov_file) %>% select(-Slide_factor)
cov <- cov[!grepl("pc",colnames(cov))]
cov$Sex_factor[cov$Sex_factor == "F"] <- 0
cov$Sex_factor[cov$Sex_factor == "M"] <- 1

index<-sapply(cov,function(.col){all(is.na(.col) | .col[1L] == .col)})
index[is.na(index)] <- FALSE
cov <- cov[,!index]

ge <- function(snp, beta, env){
    mol <- lm(as.numeric(beta) ~ as.numeric(snp)*as.numeric(env))
    sum <- summary(mol)$coef
    sum1 <- sum[rownames(sum)=="as.numeric(snp):as.numeric(env)",]
    
    if (length(sum1) == 0){
        sum1_tmp <- data.frame(col1 = NA,col2 = NA,col3 = NA,col4 = NA)
        colnames(sum1_tmp) <- colnames(sum)
        sum1 <- sum1_tmp
    } else {
        sum1 <- sum[rownames(sum)=="as.numeric(snp):as.numeric(env)",] %>% t() %>% data.frame()
        colnames(sum1) <- colnames(sum)
    }

    return(sum1)
}

runGE <- function(snp, cpg){
    geno1 <- geno %>% filter(SNP == snp) %>% select(-`(C)M`, -POS, -COUNTED, -ALT, -CHR) %>%
            t() %>% data.frame() %>%
            janitor::row_to_names(row_number=1)
    geno1$FID_IID <- rownames(geno1)
    geno1 <- geno1 %>% separate(FID_IID, sep="_", into=c("FID","IID")) %>% select(-FID)

    meth <- norm.beta[rownames(norm.beta) == cpg,] %>% data.frame()
    colnames(meth) <- "beta_value"
    meth$IID <- rownames(meth)

    d <- merge(meth, geno1, by.x="IID") %>% na.omit()
    colnames(d) <- c("IID","beta_value","genotype")

    d1 <- merge(d, cov, by.x="IID") %>% na.omit()

    out <- lapply(c(4:ncol(d1)), function(x) {
               ge(d1$genotype, d1$beta_value, d1[,x]) %>% mutate(statistics = paste0("Gx",colnames(d1)[x]))
               }) %>% bind_rows()
    return(out)
}

#gg <- function(snp1, beta, snp2){
#    mol1 <- lm(as.numeric(beta) ~ as.numeric(snp1)*as.numeric(snp2))
#    sum1 <- data.frame(summary(mol1)$coef)[-1,]
#    return(sum1)
#}

#runGG <- function(probe){
#    snplist <- pairs_uniq %>% filter(cpg == probe)
#    geno1 <- geno %>% filter(SNP %in% snplist$SNP) %>% select(-`(C)M`, -POS, -COUNTED, -ALT, -CHR) %>%
#            t() %>% data.frame() %>%
#            janitor::row_to_names(row_number=1)
#    geno1$FID_IID <- rownames(geno1)
#    geno1 <- geno1 %>% separate(FID_IID, sep="_", into=c("FID","IID")) %>% select(-FID)
     
#    geno2 <- geno1 %>% select(-IID)

#    out <- data.frame()
#    if(ncol(geno2)>1){
#        col_pairs <- combn(names(geno2), 2)

#        meth <- norm.beta[rownames(norm.beta) == probe,] %>% data.frame()
#        colnames(meth) <- "beta_value"
#        meth$IID <- rownames(meth)

#        d <- merge(meth, geno1, by.x="IID") %>% na.omit()

#        out <- lapply(c(1:ncol(col_pairs)), function(x) {
#                d1 <- d %>% select(col_pairs[1,x], col_pairs[2,x], beta_value) %>% na.omit()
#                gg(d1[,1], d1[,3], d1[,2]) %>% 
#                mutate(statistics = c(col_pairs[1,x], col_pairs[2,x], paste0(col_pairs[2,x],"*",col_pairs[2,x])))
#        }) %>% bind_rows()
#    }
#    return(out)
    
#}

GEout <- lapply(c(1:nrow(pairs_uniq)), function(x) {
                runGE(pairs_uniq$SNP[x], pairs_uniq$cpg[x]) %>%
                    mutate(vmeQTL = pairs_uniq$SNP[x], vCpG = pairs_uniq$cpg[x])
                }) %>% bind_rows()

#vCpG_uniq <- unique(pairs_uniq$cpg)
#GGout <- lapply(vCpG_uniq, function(x) {runGG(x) %>% mutate(vCpG = x)}) %>% bind_rows()


write.table(GEout, paste0(out_dir,"/GEresult_chr",chr,".txt"), col=T, row=F, sep="\t", quote=F)
#write.table(GGout, paste0(out_dir,"/GGresult_chr",chr,".txt"), col=T, row=F, sep="\t", quote=F)

