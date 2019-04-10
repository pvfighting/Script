# overlap variant-gene pairs for meta-eQTL analysis result
# April. 09, 2019
library(tidyverse)

dis_eqtl_file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.dis.MetaSoft.output.FDR05.txt"
rep_eqtl_file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.rep.MetaSoft.output.P05.txt"

data_dis = read.table(dis_eqtl_file, header=T, sep = "\t", stringsAsFactors=F)
data_rep = read.table(rep_eqtl_file, header=T, sep = "\t", stringsAsFactors=F)
#names(data)[1:16] = c("RSID","STUDY","PVALUE_FE","BETA_FE","STD_FE","PVALUE_RE","BETA_RE","STD_RE","PVALUE_RE2","STAT1_RE2","STAT2_RE2","PVALUE_BE","I_SQUARE","Q","PVALUE_Q","TAU_SQUARE")
data_dis = data_dis %>% select(1:16)
data_dis = subset(data_dis, PVALUE_RE2 != "NA")
data_rep = data_rep %>% select(1:16)
data_rep = subset(data_rep, PVALUE_RE2 != "NA")

overlap = intersect(data_dis$RSID,data_rep$RSID)
data_overlap = data_dis[which(data_dis$RSID %in% overlap),]
#data_overlap is the data in the discovery now

data_overlap = data_overlap %>% separate(RSID, into = c("GENE", "SNP"), sep = "__", remove=F)
SNP_anno_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.rd"
SNP_anno = read_tsv(SNP_anno_file)
SNP_anno$SNP = paste(SNP_anno$SNP, SNP_anno$Ref, SNP_anno$Alt, sep = ":")
SNP_anno$SNP_chrpos = paste(SNP_anno$Chr, SNP_anno$Pos,SNP_anno$Ref, SNP_anno$Alt, sep = ":")
data_overlap$SNP = do.call(rbind, strsplit(data_overlap$SNP, '_', fixed=T))[,1]
data_overlap$SNP_rs = SNP_anno$SNP[match(data_overlap$SNP, SNP_anno$SNP_chrpos)]
data_overlap$SNPPos = SNP_anno$Pos[match(data_overlap$SNP_rs, SNP_anno$SNP)]
data_overlap$SNPChr = SNP_anno$Chr[match(data_overlap$SNP_rs, SNP_anno$SNP)]
data_overlap$Region = SNP_anno$Region[match(data_overlap$SNP_rs, SNP_anno$SNP)]

gtf_file = "/Users/jiajiepeng/workspace/p1/reference/hg19/gencode.v24lift37.annotation.gtf.gene.txt"
gtf = read.table(gtf_file, header=T, stringsAsFactors=F, sep = "\t")
gtf$geneID = do.call(rbind, strsplit(gtf$gene_id, '.', fixed=T))[,1]

data_overlap$Gene_name = gtf$gene_name[match(data_overlap$GENE, gtf$geneID)]
data_overlap$PVALUE_RE2_rep = data_rep$PVALUE_RE2[match(data_overlap$RSID,data_rep$RSID)]

selectOutput = data_overlap%>%select(RSID,SNP_rs,Region,Gene_name,PVALUE_RE2,PVALUE_RE2_rep)%>% rename(PVALUE_RE2_dis = PVALUE_RE2)
write.table(selectOutput, file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.dis_rep.output.txt", row.names=F, col.names=T, sep = "\t", quote=F )

