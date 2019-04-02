library(tidyverse)
ADGWAS_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADLoci.genes.rmPos.txt"
r2_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.r2.rmPos.txt"
r2 = read.table(r2_file,  header=T, sep = "\t", stringsAsFactors=F)
ADGWAS = read.delim(ADGWAS_file, header=T, sep = "\t", stringsAsFactors=F) 
ADGWAS = ADGWAS %>% filter(USE=="YES") # 2339 AD GWAS SNPs out of 2339
r2=r2%>%filter(CHR_A!="CHR_A")
#
r2$region = ADGWAS$GENE[match(r2$SNP_A, ADGWAS$SNP)]  #2024 snps
#r22=r2%>%filter(is.na(region)) #remove the SNP not in ADLoci.genes.txt
#r2=r2%>%filter(!is.na(region)) #remove the SNP not in ADLoci.genes.txt
SNP_region = r2 %>% group_by(SNP_B) %>% summarise(region = MaxTable(region)) %>% as.data.frame
df = data.frame(SNP = unique(c(ADGWAS$SNP, SNP_region$SNP_B)), GWAS_SNP = "NO", region = "", stringsAsFactors=F)
rownames(df) = df$SNP
df[ADGWAS$SNP, "GWAS_SNP"] = "YES"
df[SNP_region$SNP_B, "region"] = SNP_region$region
df[ADGWAS$SNP, "region"] = ADGWAS$GENE
df = df %>% filter(SNP != "SNP_B") # 10142 variants #need to check with Tao
write.table(df, file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.txt", row.names=F, col.names=T, sep = "\t", quote=F )

MaxTable <- function(InVec, mult = FALSE) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  if (isTRUE(mult)) {
    levels(InVec)[A == max(A)]
  } 
  else levels(InVec)[which.max(A)]
}


####################################################
#generate the file ADGWAS.SNPs.Proxy.region.SNPsPos.txt
####################################################

ADGWAS_Proxy_region_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.txt"
EUR_rim_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/1000G/Chr1_22.EUR.rd.bim"
Proxy_region = read.table(ADGWAS_Proxy_region_file,header=T,sep="\t",stringsAsFactors = F) #10142 SNPs
EUR_rim = read.table(EUR_rim_file,header = F,sep = "\t",stringsAsFactors = F)
names(EUR_rim) = c("CHR", "SNP", "GeneticDist", "Pos", "A1", "A2")
Proxy_region_pos = merge(Proxy_region,EUR_rim,by = "SNP",all.x=FALSE) #10061 some rs ID cannot be found in EUR 1000G
Proxy_region_pos = Proxy_region_pos%>%filter(!grepl("INS",A1)) 
Proxy_region_pos = Proxy_region_pos%>%filter(!grepl("CN",A2)) #10051

Proxy_region_pos = Proxy_region_pos[,c(1,2,3,4,6,7,8)] #10051
write.table(Proxy_region_pos, file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt", row.names=F, col.names=T, sep = "\t", quote=F )

fortabix = Proxy_region_pos[,c(4,5)]
write.table(fortabix, file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt.fortabix", row.names=F, col.names=F, sep = "\t", quote=F )

#then go to python.py  dbSNP.aano.py

###############################################
#rm redundancy
############################
# R code
library(tidyverse)
gwas_proxy_variants = read.table("/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final", header=T, sep = "\t", stringsAsFactors=F)
# read in 10051 variants
gwas_proxy_variants = gwas_proxy_variants %>% arrange(Chr, Pos, Ref, Alt, desc(GWAS_SNP)) %>% distinct(Chr, Pos, Ref, Alt, .keep_all=T)
# after remove one of the duplicates, 10051 variants left
write.table(gwas_proxy_variants, file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.rd", col.names=T, row.names=F, sep = "\t", quote=F)
# 3537 GWAS SNPs , gwas_proxy_variants %>% filter(GWAS_SNP=="YES") %>% dim
# 18330 proxy SNPs, gwas_proxy_variants %>% filter(GWAS_SNP=="NO") %>% filter(nchar(Ref)==1 & nchar(Alt)==1) %>%dim
# 2036 proxy indels, gwas_proxy_variants %>% filter(GWAS_SNP=="NO") %>% filter(nchar(Ref)!=1 | nchar(Alt)!=1) %>%dim
