library(tidyr)
library(readr)
library(dplyr)

arguments <- commandArgs(T)
chr <- arguments[1]
vQTL_path <- arguments[2]

info <- read_delim(paste0(vQTL_path, "/chr", chr, "_info"),delim="/",col_names=c("cohort","chr","output"))

info1 <- info %>% separate(output, sep="_genetic", into=c("remove1","output1")) %>%
        separate(output1, sep="_cpg", into=c("genetic_chunk","output2")) %>%
        select(cohort, chr, genetic_chunk)

cohort <- read_delim(paste0(vQTL_path, "/Cohort_Information_short.csv")) %>% mutate(cohort = Study) %>%  select(cohort, Nsamples_3a)

df <- merge(info1, cohort, by.x="cohort")

write.table(df, paste0(vQTL_path, "/chr",chr,"_info_full"),col=F,row=F, sep=",", quote=F)
