library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)

p <- args[1]
print(p)
jma_files <- list.files(path = p, pattern = "jma.cojo$")

d <- lapply(paste0(p,"/", jma_files), function(x) {
                        message("reading ", x) 
                        read_delim(x, col_types=cols(refA=col_character())) %>% mutate(file=gsub(paste0(p,"/"),"",x))
                        }) %>% bind_rows()
d1 <- d %>%
        mutate(file1 = gsub(".cojo.jma.cojo","",file)) %>% 
        separate(file1, sep="_", into=c("cpg","method")) %>%
        select(SNP, cpg, method, b, se, p)

snplist <- d1$SNP %>% unique() %>% data.frame()

write.table(d1, file = paste0(p,"/vmeQTL_vCpG_pair_after_LD_COJO.txt"), col=T, row=F, sep="\t", quote=F)
write.table(snplist, file = paste0(p, "/vmeQTL_vCpG_pair_after_LD_COJO.snplist"), col=F, row=F, sep="\t", quote=F)
