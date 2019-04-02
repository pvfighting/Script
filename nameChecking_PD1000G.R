library(tidyverse)
ADGWAS_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADLoci.genes.txt"
oneKG_bim_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/Chr1_22.EUR.rd.bim"
ADGWAS = read.delim(ADGWAS_file, header=T, sep = "\t",stringsAsFactors=F) 
oneKG = read.table(oneKG_bim_file, header=F, sep = "\t", stringsAsFactors=F)
names(oneKG) = c("CHR", "Marker", "GeneticDist", "Pos", "A1", "A2")
#ADGWAS = ADGWAS %>% filter(Use=="YES") # 895 AD GWAS SNPs out of 895
#ADGWAS$CHR_BP = paste(ADGWAS$CHR_ID, ADGWAS$CHR_POS, sep=":")
oneKG$CHR_Pos = paste(paste(paste(oneKG$CHR, oneKG$Pos,sep=":"),oneKG$A1,sep="_"),oneKG$A2,sep = "_")
comm_pos = intersect(ADGWAS$uniqID.a1a2, oneKG$CHR_Pos)
ADGWAS_sub = ADGWAS %>% filter(uniqID.a1a2 %in% comm_pos) # 1997 /2357
oneKG_sub = oneKG %>% filter(CHR_Pos %in% comm_pos) # 1997
oneKG_sub$ADGWAS_marker = ADGWAS_sub$uniqID.a1a2[match(oneKG_sub$CHR_Pos, ADGWAS_sub$uniqID.a1a2)]
diff_marker = oneKG_sub[oneKG_sub$Marker!=oneKG_sub$ADGWAS_marker,] %>% select(Marker, A1, A2, CHR_Pos, ADGWAS_marker)
comm_markers = intersect(ADGWAS$SNP, oneKG$Marker) # 508 common markers