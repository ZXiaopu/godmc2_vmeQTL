library(data.table)

args <- commandArgs(T)
chr <- args[1]
index <- args[2]

setwd(paste0("../data/Meta_Results/chr",chr))

drm_markers  <- fread(paste0("DRM_random_meta_results_cpg_index", index, "1.txt"), select = "MarkerName")[[1]]
svlm_markers <- fread(paste0("SVLM_random_meta_results_cpg_index", index, "1.txt"), select = "MarkerName")[[1]]

exclude_dt <- data.table(MarkerName = unique(c(drm_markers, svlm_markers)))
rm(drm_markers, svlm_markers); gc() # Free memory immediately

bf_dt <- fread(paste0("BF_random_meta_results_cpg_index", index, "1.txt"), drop = "valid_indices")

BF_only <- bf_dt[!exclude_dt, on = "MarkerName"]
rm(bf_dt, exclude_dt); gc()

BF_only[, SNP_Probe := MarkerName]
BF_only[, c("SNP", "A1", "A2", "CpG") := tstrsplit(MarkerName, "_", fixed = TRUE)]
n_asso <- uniqueN(BF_only$SNP_Probe)
n_cpg  <- uniqueN(BF_only$CpG)
message(sprintf("Results detected by BF but not DRM/SVLM are: n(asso) = %d n(cpg) = %d", n_asso, n_cpg))
fwrite(BF_only, file = paste0("Missing_asso_run_by_BF_only_index", index, "_260904.txt"), sep = "\t", quote = FALSE)
