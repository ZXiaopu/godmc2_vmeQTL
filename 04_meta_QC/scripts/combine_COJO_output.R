library(readr)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)

p <- args[1]
setwd(p)

## COJO have been conducted on these CpGs
jma_files <- list.files(path = ".", pattern = "jma.cojo$")

## all CpGs in this chr
cpglist_files <- list.files(path = ".", pattern = "cpglist")

c <- lapply(cpglist_files, function(x){
                        read_delim(x, delim=" ",col_names=F) %>% 
                            mutate(CpG = X1, file=x) %>% 
                            separate(file, sep="_", into=c("cpglist","method")) %>%
                            select(CpG, method)}) %>% bind_rows()

c <- c %>% mutate(cpg_method = paste(CpG, method, sep="_"))
c1 <- c %>% filter(cpg_method %in% gsub(".cojo.jma.cojo","",jma_files)==F)

tmp <- lapply(c1$cpg_method, function(x){read_delim(paste0(x,".snplist"), delim="\t") %>% 
                                            mutate(file=x,file1=x)}) %>% 
                                            bind_rows() %>% separate(file1, sep="_", into=c("cpg","method"))

tmp1 <- lapply(c(1:nrow(tmp)), function(x){read_delim(paste0(tmp$file[x], ".cojo.ma")) %>% 
                                            filter(SNP == tmp$SNP[x]) %>%
                                            mutate(SNP = tmp$SNP[x], cpg=tmp$cpg[x], method = tmp$method[x]) %>%
                                            select(SNP, cpg, method, b, se, p)
                                            }) %>% bind_rows()

d <- lapply(jma_files, function(x) {
                        message("reading ", x) 
                        read_delim(x, col_types=cols(refA=col_character())) %>% mutate(file=x)
                        }) %>% bind_rows()
d1 <- d %>% 
        mutate(file1 = gsub(".cojo.jma.cojo","",file)) %>% 
        separate(file1, sep="_", into=c("cpg","method")) %>%
        select(SNP, cpg, method, b, se, p)

out <- rbind(tmp1, d1)
snplist <- out$SNP %>% unique() %>% data.frame()

write.table(out, file = paste0(p,"/vmeQTL_vCpG_pair_after_LD_COJO.txt"), col=T, row=F, sep="\t", quote=F)
write.table(snplist, file = paste0(p, "/vmeQTL_vCpG_pair_after_LD_COJO.snplist"), col=F, row=F, sep="\t", quote=F)
