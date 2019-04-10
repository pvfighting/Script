# significant variant-gene pairs for meta-eQTL analysis result
# April. 09, 2019
library(tidyverse)

meta_eqtl_file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.rep.MetaSoft.output.tsv"
data = read_tsv(meta_eqtl_file, col_names=F,  skip =1)
names(data)[1:16] = c("RSID","STUDY","PVALUE_FE","BETA_FE","STD_FE","PVALUE_RE","BETA_RE","STD_RE","PVALUE_RE2","STAT1_RE2","STAT2_RE2","PVALUE_BE","I_SQUARE","Q","PVALUE_Q","TAU_SQUARE")
data = data %>% select(1:16)
data = subset(data, PVALUE_RE2 != "NA")
data$FDR_FE = p.adjust(data$PVALUE_FE)
data$FDR_RE = p.adjust(data$PVALUE_RE)
data$FDR_RE2 = p.adjust(data$PVALUE_RE2)
RE2_FDR = data %>% filter(PVALUE_RE2 < 0.05)

write.table(RE2_FDR, file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.rep.MetaSoft.output.P05.txt", row.names=F, col.names=T, sep = "\t", quote=F )
