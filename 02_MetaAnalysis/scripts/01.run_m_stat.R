library(getmstatistic)
library(gridExtra)       
library(ggplot2)
library(tidyr)
library(dplyr)
library(readr)
library(data.table)

input_path <- "../data"

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

    # Count once
    counts_all <- o[, .N, by = .(SNP, Probe, SNP_Probe)]
    counts_ref <- ref[, .N, by = .(SNP, Probe, SNP_Probe)]

    keep <- unique(c(counts_all[N >= 10, SNP_Probe],counts_ref[N >= 10, SNP_Probe]))
    o3 <- o[J(keep), on = "SNP_Probe"]

    setDT(o3)
    asso <- o3[
                o3[
                    order(Probe, file != "bib_eur_mother"),
                    .SD[1],
                    by = Probe
                    ],
            on = "SNP_Probe"
            ]

    removeSNP <- asso %>% group_by(SNP) %>% tally() %>% filter(n>=30)
    asso <- asso %>% filter((SNP %in% removeSNP$SNP)==F)

 #   if(include_two_cohorts == T){
 #       o %>% filter(file %in% c("bib_eur_mother")) %>% group_by(Probe) %>% tally() %>% nrow()
 #       ref <- o %>% filter(file %in% c("bib_eur_mother"))
 #       o1 <- o %>% filter(SNP_Probe %in% ref$SNP_Probe) %>% group_by(SNP, Probe, SNP_Probe) %>% tally() %>% filter(n>=10)
 #       o2 <- o %>% group_by(SNP, Probe, SNP_Probe) %>% tally() %>% filter(n>=10)
 #       o3 <- o %>% filter(SNP_Probe %in% c(o1$SNP_Probe, o2$SNP_Probe))
 #   } else{
 #       o <- o %>% filter((file %in% c("bib_eur_mother","GLAKU"))==F)
 #       o3 <- o %>% group_by(SNP, Probe, SNP_Probe) %>% tally() %>% filter(n>=10)
 #   }

 #   asso <- data.frame()
 #   for (c in unique(o3$Probe)){
 #       tmp <- o3 %>% filter(Probe == c)
 #       
 #       if ("bib_eur_mother" %in% tmp$file){
 #           tmp <- tmp %>% filter(file == "bib_eur_mother")
 #       }
 #
 #        asso0 <- o3 %>% filter(SNP_Probe %in% tmp$SNP_Probe[1])
 #        asso <- rbind(asso, asso0)
 #   }

    ## bib_eur_mother cohort does not have as many association available as others
    ## run m-statistics for two datasets
    ## ds1: all 21 cohorts
    ## ds2: excluding bib_eur_mother
    
    ds1 <- asso
    ds1m <- getmstatistic(ds1$b, ds1$SE, ds1$SNP_Probe, ds1$file, save_dir=paste0(input_path,"/Mstat_P",pval, "/",method))
#    ds2m <- getmstatistic(ds2$b, ds2$SE, ds2$SNP_Probe, ds2$file, save_dir=paste0(input_path,"/Mstat_P5.8e-14/", method, "_chr",chr,"_20cohorts"))
    save(ds1m, ds1, file = paste0(input_path, "/Mstat_P", pval, "/", method, "_allchr_10qtls_Mtat.RData"))
}

run_m(BFfiles, "BF", 5.8e-14, TRUE)
run_m(DRMfiles, "DRM", 5.8e-14, TRUE)
run_m(SVLMfiles, "SVLM", 5.8e-14, TRUE)
