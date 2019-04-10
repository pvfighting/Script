# Manhattan plot for meta-eQTL analysis of 8 studies
# Jan. 12, 2018
library(tidyverse)

#meta_eqtl_file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.rep.MetaSoft.output.tsv"
meta_eqtl_file = "/Users/jiajiepeng/workspace/p1/MetaAna/AD.meta_eqtl.cis.7.MetaSoft.output.tsv"
data = read_tsv(meta_eqtl_file, col_names=F,  skip =1)
names(data)[1:16] = c("RSID","STUDY","PVALUE_FE","BETA_FE","STD_FE","PVALUE_RE","BETA_RE","STD_RE","PVALUE_RE2","STAT1_RE2","STAT2_RE2","PVALUE_BE","I_SQUARE","Q","PVALUE_Q","TAU_SQUARE")
data = data %>% select(1:16)
data = subset(data, PVALUE_FE != "NA")
data = data %>% separate(RSID, into = c("GENE", "SNP"), sep = "__", remove=F)
data$FDR_FE = p.adjust(data$PVALUE_FE)
data$FDR_RE = p.adjust(data$PVALUE_RE)
data$FDR_RE2 = p.adjust(data$PVALUE_RE2)
SNP_anno_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.rd"
SNP_anno = read_tsv(SNP_anno_file)
SNP_anno$SNP = paste(SNP_anno$SNP, SNP_anno$Ref, SNP_anno$Alt, sep = ":")
SNP_anno$SNP_chrpos = paste(SNP_anno$Chr, SNP_anno$Pos,SNP_anno$Ref, SNP_anno$Alt, sep = ":")
data$SNP = do.call(rbind, strsplit(data$SNP, '_', fixed=T))[,1]
data$SNP_rs = SNP_anno$SNP[match(data$SNP, SNP_anno$SNP_chrpos)]
data$SNPPos = SNP_anno$Pos[match(data$SNP_rs, SNP_anno$SNP)]
data$SNPChr = SNP_anno$Chr[match(data$SNP_rs, SNP_anno$SNP)]
data$Region = SNP_anno$Region[match(data$SNP_rs, SNP_anno$SNP)]

# annotate FE_FDR, RE_FDR, RE2_FDR
gtf_file = "/Users/jiajiepeng/workspace/p1/reference/hg19/gencode.v24lift37.annotation.gtf.gene.txt"
gtf = read.table(gtf_file, header=T, stringsAsFactors=F, sep = "\t")
gtf$geneID = do.call(rbind, strsplit(gtf$gene_id, '.', fixed=T))[,1]
#FE_FDR = data %>% filter(FDR_FE < 0.05) %>% filter(PVALUE_FE !=0) %>% select(GENE, SNP_rs, SNPPos, SNPChr, Region, PVALUE_FE, FDR_FE, BETA_FE, STD_RE) %>% rename(SNP = SNP_rs, FDR=FDR_FE)# 156 unique genes
#RE_FDR = data %>% filter(FDR_RE < 0.05) %>% filter(PVALUE_RE !=0) %>% select(GENE, SNP_rs, SNPPos, SNPChr, Region, PVALUE_RE, FDR_RE, BETA_RE, STD_RE) %>% rename(SNP = SNP_rs, FDR=FDR_RE)# 90 unique genes
RE2_FDR = data %>% filter(FDR_RE2 < 0.05) %>% filter(PVALUE_RE2 !=0) %>% select(GENE, SNP_rs, SNPPos, SNPChr, Region, PVALUE_RE2, FDR_RE2, STAT1_RE2, STAT2_RE2) %>% rename(SNP = SNP_rs, FDR=FDR_RE2) # 167 unique genes, Mar.6 2018
Annotate <- function(df){
  eqtl = df
  # select the most significant SNP for each gene
  eqtl = eqtl %>% rename(gene = GENE)
  by_gene = eqtl %>% arrange(FDR) %>% group_by(gene) %>% slice(1:1) %>% as.data.frame %>% arrange(FDR)
  by_gene_count = eqtl %>% arrange(FDR) %>% group_by(gene) %>% summarise(count = length(gene)) %>% as.data.frame 
  by_gene$Number_of_variants = by_gene_count$count[match(by_gene$gene, by_gene_count$gene)]
  by_gene$Gene_name = gtf$gene_name[match(by_gene$gene, gtf$geneID)]
  by_gene$CHR = gtf$CHR[match(by_gene$gene, gtf$geneID)]
  by_gene$CHR = sub('chr','',by_gene$CHR)
  by_gene$Start = gtf$Start[match(by_gene$gene, gtf$geneID)]
  by_gene$End = gtf$End[match(by_gene$gene, gtf$geneID)]
  by_gene$SNP_is_GWAS = SNP_anno$GWAS_SNP[match(by_gene$SNP, SNP_anno$SNP)]
  by_gene$SNP_in_Region = SNP_anno$Region[match(by_gene$SNP, SNP_anno$SNP)] 
  
  by_gene[is.na(by_gene)] = ""
  return(by_gene)
}
meta_egenes = Annotate(RE2_FDR)


RE2_all = data %>% filter(PVALUE_RE2 !=0) %>% select(GENE, SNP_rs, SNPPos, SNPChr, Region, PVALUE_RE2, FDR_RE2, STAT1_RE2, STAT2_RE2) %>% rename(SNP = SNP_rs, FDR=FDR_RE2) %>%filter(!((SNPChr==17)&(SNPPos >=56384022) ))
 #remove some variants which may not corrcet labeled. related with ABI3 loci in Chr 17


data_uniq = RE2_all # show snp-gene pairs
region_freq = as.data.frame(table(data_uniq$Region))
region_to_remove = as.character(subset(region_freq, Freq < 50)$Var1)
data_uniq = subset(data_uniq, !(Region %in% region_to_remove))
data_uniq = data_uniq[order(data_uniq$SNPChr, data_uniq$SNPPos),]
min_region_pos = aggregate(SNPPos ~ Region, data_uniq, min) #left most position for SNPs in each region
index = match(data_uniq$Region, min_region_pos$Region)
data_uniq$newPos = data_uniq$SNPPos - min_region_pos$SNPPos[index] # for each region align to [0-X]
the_lag = aggregate(newPos ~ Region, data_uniq, max)
the_lag = the_lag[match(unique(data_uniq$Region), the_lag$Region),] # sort the_lag by region name
the_lag$addValue = Reduce(`+`, the_lag$newPos, accumulate=T) # set coordinates shift
the_lag$addValue = c(0, the_lag$addValue[1:nrow(the_lag)-1] + 150000*(1:(nrow(the_lag)-1)) )
data_uniq$newPos = data_uniq$newPos + the_lag$addValue[match(data_uniq$Region, the_lag$Region)]

if(nrow(the_lag)%%2 ==1){
  the_lag$class = c(rep(1:2,nrow(the_lag)/2),1)
}else{
  the_lag$class = rep(1:2,nrow(the_lag)/2)
}

data_uniq$class = as.factor(the_lag$class[match(data_uniq$Region, the_lag$Region)])
labelPos = aggregate(newPos ~ Region, data_uniq, range) # position for labeling text
labelPos$newPos = apply(labelPos$newPos, 1, mean)
the_lag$labelPos = labelPos$newPos[match(the_lag$Region, labelPos$Region)]
the_lag$Chr = SNP_anno$Chr[match(the_lag$Region, SNP_anno$Region)]
#Add vertical lines
tmp = aggregate(newPos~Region,data_uniq, range) # pos range in new coordinating system
vline_offset = 150000/2
the_lag$region_vline_start =  tmp$newPos[match(the_lag$Region, tmp$Region),1] - vline_offset
the_lag$region_vline_end   =  tmp$newPos[match(the_lag$Region, tmp$Region),2] + vline_offset

# real/true median positon (in Mbp) for each region
tmp = aggregate(SNPPos~Region, data_uniq, range)
the_lag$realPos_start = tmp$SNPPos[match(the_lag$Region, tmp$Region),1]
the_lag$realPos_start = paste(round(the_lag$realPos_start/1000000,1),'',sep ="")
the_lag$realPos_end = tmp$SNPPos[match(the_lag$Region, tmp$Region),2]
the_lag$realPos_end = paste(round(the_lag$realPos_end/1000000,1),'',sep ="")

# mark egenes
meta_egenes$SNPnewPos = data_uniq$newPos[match(meta_egenes$SNP, data_uniq$SNP)]
meta_egenes = meta_egenes %>% rename(Pvalue = PVALUE_RE2)
meta_egenes_top = meta_egenes %>% arrange(Pvalue) %>% group_by(SNP_in_Region) %>% slice(1:1)
egene_count_by_region = meta_egenes %>% group_by(SNP_in_Region) %>% summarise(egene_count = length(gene)) # show how many egenes per locus
meta_egenes_top = meta_egenes_top %>% mutate(egene_count = egene_count_by_region$egene_count[match(SNP_in_Region, egene_count_by_region$SNP_in_Region)])
# addon_genes = c("GBAP1", "CCNT2", "MCCC1", "SNCA-AS1", "HLA-DQB1", "HLA-DQB1-AS1", "KANSL1", "LRRC37A", "LRRC37A2", "LRRC37A4P", "MAPT", "MAPT-AS1", "HSD17B1")
# meta_egenes_top_addon = meta_egenes %>% filter(Gene_name %in% addon_genes)
# meta_egenes_top = rbind(meta_egenes_top%>%ungroup, meta_egenes_top_addon)

# threshold line: max p value with FDR<0.05
P_threshold = data_uniq %>% filter(FDR<0.05) %>% .$PVALUE_RE2 %>% max
P_threshold_log = -log10(P_threshold)
the_lag$SRegion = sub('/.*$', '', the_lag$Region)  #simplified region names, using one gene


######### show -log10(Pvlaues) > FDR 0.05 in log2 scale #########
thin_data_uniq = data_uniq %>% mutate( logP = -log10(PVALUE_RE2) )
thin_data_uniq_up = thin_data_uniq %>% filter( logP >= P_threshold_log)
thin_data_uniq_down = thin_data_uniq %>% filter( logP < P_threshold_log)
offset = 0.5
thin_data_uniq_down = thin_data_uniq_down %>% mutate(logP_scaled = log2(P_threshold_log)-offset + (logP - min(.$logP))/(max(.$logP)-min(.$logP)) * (offset) ) # normalize the -log10(pvalues) to (log2(P_threshold_log), log2(P_threshold_log)-offset)
thin_data_uniq_up = thin_data_uniq_up %>% mutate(logP_scaled = log2(logP))
thin_data_uniq = rbind(thin_data_uniq_up, thin_data_uniq_down)

# thin the graph
DIST = 5000000
thin_data_uniq$pos_thin = round(thin_data_uniq$newPos/DIST, 2)
thin_data_uniq$logp_thin = round(thin_data_uniq$logP_scaled/3, 2) 
thin_data_uniq = thin_data_uniq %>% arrange(pos_thin, logp_thin) %>% group_by(pos_thin, logp_thin) %>% sample_n(1)

# manhattan plot for CZI grant
pdf("/Users/jiajiepeng/workspace/p1/plot/Thin_manhattan.AD.meta.cis7.pdf", width=5, height=2)
g = ggplot(thin_data_uniq, aes(x=newPos, logP_scaled))+ geom_point(aes(color=class), size=0.8, stroke=0) 
g = g + geom_hline(yintercept = log2(-log10(P_threshold)), color="red")
g = g + scale_x_continuous(breaks = the_lag$labelPos, labels = the_lag$SRegion, expand = c(0.025,0))
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=6), axis.text.y = element_text(hjust = 1, size=6))
g = g + theme(panel.grid.major = element_blank(), panel.grid.minor.x= element_blank(), axis.ticks.length =  unit(2, "pt") ) 
g = g + geom_vline(xintercept= unique(c(the_lag$region_vline_start[2:nrow(the_lag)-1], the_lag$region_vline_end[2:nrow(the_lag)-1])), color="white", alpha = 1, size=0.3, linetype =6)
g = g + geom_text(data = the_lag, aes(x = labelPos, y = log2(P_threshold_log)-offset-0.3, label= Chr), size = 2, check_overlap = F, angle=0)
g = g + geom_point(data = meta_egenes_top, aes(x = SNPnewPos, y = log2(-log10(Pvalue))), size=1.2, color = "goldenrod1",shape=18)
g = g + geom_text(data = meta_egenes_top, aes(x = SNPnewPos, y = log2(-log10(Pvalue)), label = Gene_name), nudge_y=0.2, size = 1.5, check_overlap = F, angle= 45, hjust =0)
g = g + scale_colour_discrete(guide=F)
# g = g + xlab("Regions")
# g = g + ylab("-log10(Pvalue)")
g = g + scale_y_continuous(breaks = log2( c(10, 20, 40, 80, 150, 250) ), labels =  c(10, 20, 40, 80, 150, 250), minor_breaks = log2(seq(10,250,10)), limits = c(log2(P_threshold_log)-offset-0.35,log2(400)) )
#g = g + theme(axis.title = element_text(size = 8))
g = g + theme(axis.title = element_blank())
g 
dev.off()