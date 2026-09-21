library(tidyr)
library(readr)
library(dplyr)
library(meffil)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

cohortlist <- read_delim("../../../GoDMC_vmeQTL_phase1_data/vQTL_results/Cohorts_with_all_results/Cohort_list.csv", delim="\t", col_names="Cohort")

minfi_annots <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep_cpgs <- !(minfi_annots$chr %in% c("chrX", "chrY"))
minfi_keep <- minfi_annots[keep_cpgs, c("Name","chr", "pos", "Probe_rs", "Probe_maf", "CpG_rs", "CpG_maf", "SBE_rs", "SBE_maf")] %>% data.frame()
minfi_problem_cpgs <- minfi_keep %>% filter((!is.na(CpG_rs)) | (!is.na(SBE_rs)))
minfi_keep_cpgs <- minfi_keep %>% filter((is.na(CpG_rs)) & (is.na(SBE_rs)))

o <- data.frame()
for (cohort in cohortlist$Cohort){
    load(paste0(cohort, "/results/03/methylation_summary.RData"))
    df <- meth_summary[order(meth_summary$sd, decreasing=T),] %>% filter(cpg %in% minfi_keep_cpgs$Name)
    topCpGs <- head(df,10000) %>% mutate(Study=cohort)
    o <- rbind(o, topCpGs)
}

write.table(o, "../../03_FeatureCheck/data/VMP/vmeQTL_21cohort_top20k_CpGs_remove_snp_in_probe.txt", col=T, row=F, sep="\t", quote=F)

o <- read_delim("../../03_FeatureCheck/data/VMP/vmeQTL_21cohort_top20k_CpGs_remove_snp_in_probe.txt")

count10 <- o %>% group_by(cpg) %>% tally() %>% filter(n>=10)
write.table(count10, "../../03_FeatureCheck/data/VMP/vmeQTL_21cohort_top20k_CpGs_remove_snp_in_probe_repeated10cohorts.txt", col=T, row=F, sep="\t", quote=F)

o_count10 <- o %>% filter(cpg %in% count10$cpg)

count10$weighted_sd <- 0
count10$mean_sd <- 0

for (idx in c(1:nrow(count10))){
    print(idx)
    probe <- count10$cpg[idx]
    tmp <- o_count10 %>% filter(cpg == probe)
    sum_n <- sum(tmp$outlier.n)
    tmp$weight <- tmp$outlier.n/sum_n
    count10$weighted_sd[idx] <- sum(tmp$weight*tmp$sd)
    count10$mean_sd[idx] <- mean(tmp$sd)
}

count10_weighted_ordered <- count10[order(count10$weighted_sd, decreasing=T),]
count10_weighted_ordered_sd0.15 <- count10_weighted_ordered %>% filter(weighted_sd>0.15) #863
#count10_weighted_ordered_sd0.12 <- count10_weighted_ordered %>% filter(weighted_sd>0.12) #2,225 / sd 0.13 1599 / sd 0.14 1191

count10_mean_ordered <- count10[order(count10$mean_sd, decreasing=T),]
count10_mean_ordered_sd0.15 <- count10_mean_ordered %>% filter(mean_sd>0.15) #1,027

godmc_candidate <- count10_mean_ordered_sd0.15 %>% filter(cpg %in% count10_weighted_ordered_sd0.15$cpg) #844

seale <- read_delim("../../03_FeatureCheck/data/VMP/Seale_Meta-analysis_results_VMPs_blood.txt") %>% filter(`VMP.class` != "non-VMP") %>% mutate(study="Seale") %>% dplyr::select(CpG, study)
slieker <- read_delim("../../03_FeatureCheck/data/VMP/Slieker_VMP.csv") %>% mutate(study="Slieker") %>% dplyr::select(CpG, study)
grant_all <- read_delim("../../03_FeatureCheck/data/VMP/Grant_allVMP.csv")
grant_lowICC <- read_delim("../../03_FeatureCheck/data/VMP/Grant_VMP_with_low_ICC.csv")
grant <- grant_all %>% filter((Name %in% grant_lowICC$Name)==F) %>% mutate(CpG=Name, study="Grant") %>% dplyr::select(CpG, study)

publication_candidate <- seale %>% filter(CpG %in% slieker$CpG) %>% filter(CpG %in% grant$CpG) %>% mutate(cpg=CpG)
vmeQTL_list5 <- rbind(publication_candidate %>% select(cpg) %>% mutate(dataset="publications"), 
                      godmc_candidate %>% select(cpg) %>% mutate(dataset="godmc"))
colnames(vmeQTL_list5) <- c("ID","dataset")
write.table(vmeQTL_list5, file="../data/list5_candidateCpGs.txt", col=T, row=F, sep="\t", quote=F)

godmc <- count10_weighted_ordered_sd0.15 %>% mutate(CpG=cpg, study="GoDMC") %>% dplyr::select(CpG, study)

cpginput <- list(Seale=seale$CpG, Slieker=slieker$CpG, Grant=grant$CpG, GoDMC=godmc$CpG)

library(VennDiagram)
library(grid)
venn_obj <- venn.diagram(
  x = cpginput,
  category.names = c("Seale", "Slieker", "Grant", "GoDMC"),
  filename = NULL,                      # NULL prints to display instead of writing a PNG
  fill = c("#264653", "#2A9D8F", "#E9C46A", "#E76F51"),       # Fill colors
  alpha = 0.5,                          # Transparency
  col = "white",                        # Border line color
  lwd = 2,                              # Border line width
  cat.cex = 1.1,                        # Label font size
  cat.fontface = "bold"
)

pdf("/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/03_FeatureCheck/data/VMP/GoDMCsd_VMP_venn_plot.pdf", width = 8, height = 8)
grid.draw(venn_obj)
dev.off()

count10_weighted_ordered_sd0.2 %>% filter(cpg %in% grant$Name) %>% nrow() # 175
count10_weighted_ordered_sd0.2 %>% filter(cpg %in% seale$CpG) %>% nrow() # 13
count10_weighted_ordered_sd0.2 %>% filter(cpg %in% slieker$CpG) %>% nrow() # 0
count10_weighted_ordered_sd0.2 %>% filter(cpg %in% slieker$CpG) %>% filter(cpg %in% seale$CpG) %>% nrow() # 0

count10_mean_ordered_sd0.2 %>% filter(cpg %in% grant$Name) %>% nrow() # 204
count10_mean_ordered_sd0.2 %>% filter(cpg %in% seale$CpG) %>% nrow() # 17
count10_mean_ordered_sd0.2 %>% filter(cpg %in% slieker$CpG) %>% nrow() # 0
count10_mean_ordered_sd0.2 %>% filter(cpg %in% slieker$CpG) %>% filter(cpg %in% seale$CpG) %>% nrow() # 0

seale %>% filter(CpG %in% slieker$CpG) %>% nrow() #5,668
ageVMP_overlap <- seale %>% filter(CpG %in% slieker$CpG)
ageVMP_overlap_ordered <- ageVMP_overlap[order(ageVMP_overlap$FDR),]
ageVMP_overlap %>% filter(CpG %in% grant$Name) #1433
