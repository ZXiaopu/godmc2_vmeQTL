library(readr)
library(dplyr)
library(tidyr)
library(janitor)
library(car)

args <- commandArgs(trailingOnly=TRUE)
norm.beta <- args[1]
pairs <- args[2]
genotype <- args[3]
outfile <- args[4]
#norm.beta <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/GoDMC_NSHD99/godmc_phase2/processed_data/methylation_data/untransformed_methylation_adjusted_pcs.RData"
#pairs <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/LD_COJO_input/chr21/vmeQTL_vCpG_pair_after_LD_COJO.txt"
#genotype <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/genetic_data/vmeQTL.chr21.filtered.traw"
#osca <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/GoDMC_NSHD99/godmc_phase2/resources/bin/osca"
#outfile <- "/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/04_meta_QC/processed_data/adjust_mean_effect/vmeQTL_vCpG_adjust_mean_effect_chr21.txt"

load(norm.beta)
pair <- read_delim(pairs)
pair1 <- pair %>% group_by(SNP, cpg) %>% tally()
geno <- read_delim(genotype)

adjust_mean <- function(snp, cpg){
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
    
    d$resi <- residuals(lm(beta_value ~ genotype, data = d))
    res <- run_vQTL(d$genotype, d$resi) %>% mutate(vmeQTL = snp, vCpG = cpg)
    return(res)
}

BF <- function(x, y){
  BF.P <- leveneTest(y ~ as.factor(x), center = median)$`Pr(>F)`[1]
  return(BF.P)
}

DRM <- function(x, y){
  Y.i <- tapply(y, x, median)
  DRM <- abs(y - Y.i[as.factor(x)])
  DRMres <- summary(lm(DRM ~ as.numeric(x)))
  return(DRMres)
}

SVLM <- function(x, y){
  resid.SVLM <- residuals(summary(lm(y ~ as.numeric(x))))
  SVLMres <- summary(lm((resid.SVLM^2) ~ as.numeric(x)))
  return(SVLMres)
}

## SNP is the genotype data; trait is the liver fat PDFF
run_vQTL <- function(snp, trait){
  BF_Pval <- BF(snp, trait)
  DRMres <- DRM(snp, trait)
  SVLMres <- SVLM(snp, trait)
  out <- data.frame(BF_P = BF_Pval,
                    DRM_P = DRMres$coef[2,4],
                    SVLM_P = SVLMres$coef[2,4])
  return(out)
}

out <- lapply(1:nrow(pair1), function(x) adjust_mean(pair1$SNP[x], pair1$cpg[x])) %>% bind_rows() 
out1 <- out %>%
    reshape2::melt(id.vars=c("vmeQTL", "vCpG")) 
out1 <- out1 %>% 
    mutate(variable = gsub("_P","",out1$variable), P_adjMean=value) %>%
    mutate(info = paste(vmeQTL, vCpG, variable, sep="_")) %>% select(info, P_adjMean)
pair <- pair %>% mutate(info = paste(SNP, cpg, method, sep="_"))

out2 <- merge(out1, pair, by.x="info")
write.table(out2, file = outfile, col=T, row=F, sep="\t", quote=F)
