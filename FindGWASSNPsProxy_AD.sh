cd /PHShome/jq34/MEGA/hg19/1000G

# Extract EUR population
awk 'BEGIN{OFS="\t"} {if($3=="EUR"){print $1,$1}}' /PHShome/jq34/neurogen/Tao/hg19/1000G/integrated_call_samples_v3.20130502.ALL.panel > EUR/EUR.list.txt
cd ./EUR

# Extract EUR population and filter variants with MAF < 0.05。# In Tao's ms, this threshold is 0.01.  We set it as 0.05
for i in $(seq 1 22);do
plink --bfile \
/PHShome/jq34/neurogen/referenceGenome/Homo_sapiens/UCSC/hg19/Annotation/Variation/1000G/Chr${i} \
--keep EUR.list.txt --maf 0.05 --make-bed --out byChr/Chr${i}.EUR
done

cd ./byChr
#remove duplicated markers (non-biallelic markers)

for i in $(seq 1 22);do
cut -f 2 Chr${i}.EUR.bim | sort | uniq -d > Chr${i}.EUR.bim.dup
plink --bfile Chr${i}.EUR --exclude Chr${i}.EUR.bim.dup --make-bed --out Chr${i}.EUR.rd
done

cat Chr*.EUR.rd.bim > Chr1_22.EUR.rd.bim
# Copy EUR folder to Orchestra cluster
# Working on Orchestra cluster
cd /home/jp360/MEGA/1000G/EUR
# R code to check the marker name consistance between AD GWAS and 1000G
# For double-checking
#R code, on the laptop nameChecking_PD1000G.R

#Calculate r2
cd /home/jp360/MEGA/ADGWAS
awk 'BEGIN{FS="\t"}$16 == "YES"' ADLoci.genes.txt | cut -f1 > ADLoci.genes.txt.col1 # 

cd /home/jp360/MEGA/1000G/EUR/byChr_rd
for i in $(seq 1 22);do
plink --bfile Chr${i}.EUR.rd --r2 dprime --ld-snp-list /home/jp360/MEGA/ADGWAS/ADLoci.genes.txt.col1.rsID.rmPos --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0.4 --out Chr${i}.EUR.rd
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' Chr${i}.EUR.rd.ld > Chr${i}.EUR.rd.ld.txt
done
cat Chr*.EUR.rd.ld.txt > Chr1_22.rd.ld.txt

#plink --bfile Chr21.EUR.rd --r2 dprime --ld-snp-list /home/jp360/MEGA/ADGWAS/ADLoci.genes.txt.col1.rsID.rmPos --ld-window-kb 250 --ld-window 99999 --ld-window-r2 0.4 --out Chr21.EUR.rd
#awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8}' Chr21.EUR.rd.ld > Chr21.EUR.rd.ld.txt


cd /home/jp360/MEGA/ADGWAS

####ADGWAS_proxy_region.R#####
library(tidyverse)
PDGWAS_file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS27old+35novelLoci.txt"
r2_file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.r2.txt"
r2 = read.table(r2_file,  header=T, sep = "\t", stringsAsFactors=F)
PDGWAS = read.table(PDGWAS_file, header=T, sep = "\t", stringsAsFactors=F) 
PDGWAS = PDGWAS %>% filter(Use=="YES") # 3537 PD GWAS SNPs out of 3723
r2$region = PDGWAS$Candidate.Gene[match(r2$SNP_A, PDGWAS$q)]
SNP_region = r2 %>% group_by(SNP_B) %>% summarise(region = MaxTable(region)) %>% as.data.frame
df = data.frame(SNP = unique(c(PDGWAS$q, SNP_region$SNP_B)), GWAS_SNP = "NO", region = "", stringsAsFactors=F)
rownames(df) = df$SNP
df[PDGWAS$q, "GWAS_SNP"] = "YES"
df[SNP_region$SNP_B, "region"] = SNP_region$region
df[PDGWAS$q, "region"] = PDGWAS$Candidate.Gene
df = df %>% filter(SNP != "SNP_B") # 23949 variants
write.table(df, file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.rmPos.txt", row.names=F, col.names=T, sep = "\t", quote=F )

MaxTable <- function(InVec, mult = FALSE) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  if (isTRUE(mult)) {
    levels(InVec)[A == max(A)]
  } 
  else levels(InVec)[which.max(A)]
}


#######generate ADGWAS.SNPs.Proxy.region.SNPsPos.txt#############
###ADGWAS_proxy_region.R####

ADGWAS_Proxy_region_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.txt"
EUR_rim_file = "/Users/jiajiepeng/workspace/p1/ADGWAS/1000G/Chr1_22.EUR.rd.bim"
Proxy_region = read.table(ADGWAS_Proxy_region_file,header=T,sep="\t",stringsAsFactors = F) #10142 SNPs
EUR_rim = read.table(EUR_rim_file,header = F,sep = "\t",stringsAsFactors = F)
names(EUR_rim) = c("CHR", "SNP", "GeneticDist", "Pos", "A1", "A2")
Proxy_region_pos = merge(Proxy_region,EUR_rim,by = "SNP",all.x=FALSE) #10061
Proxy_region_pos = Proxy_region_pos%>%filter(!grepl("INS",A1)) 
Proxy_region_pos = Proxy_region_pos%>%filter(!grepl("CN",A2)) #10051

Proxy_region_pos = Proxy_region_pos[,c(1,2,3,4,6,7,8)] #10051
write.table(Proxy_region_pos, file = "/Users/jiajiepeng/workspace/p1/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt", row.names=F, col.names=T, sep = "\t", quote=F )

##################
tabix -R ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt.fortabix  /home/tw83/twang/1000G/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz > ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt.1000G

## oneKG.anno.py
## oneKG.anno.py identify the ref allele and alt allele
class BIM:
	def __init__(self, CHR, pos, marker, A1, A2):
		self.CHR = CHR
		self.pos = pos
		self.marker = marker
		self.A1 = A1
		self.A2 = A2
oneKG_anno_file = "/home/jp360/MEGA/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt.1000G"
PD_risk_SNP_file = "/home/jp360/MEGA/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.txt"
#PD_risk_SNP_file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.txt"
output_file = "/home/jp360/MEGA/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.txt"
#output_file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.anno.txt"
out = open(output_file, 'w')
oneKG_rsID_dic = {}
oneKG_chrpos_dic = {}
for line in open(oneKG_anno_file):
	linesplit = line.strip().split("\t")
	CHR = linesplit[0]
	pos = linesplit[1]
	rsID = linesplit[2]
	refAllele = linesplit[3]
	altAllele = linesplit[4]
	for altAllele in altAllele.split(','):
		value = BIM(CHR, pos, rsID, refAllele, altAllele)
		if rsID in oneKG_rsID_dic:
			oneKG_rsID_dic[rsID].append(value)
		else:
			oneKG_rsID_dic[rsID] = [value]

		if CHR+":"+pos in oneKG_chrpos_dic:
			oneKG_chrpos_dic[CHR+":"+pos].append(value)
		else:
			oneKG_chrpos_dic[CHR+":"+pos] = [value]

for line in open(PD_risk_SNP_file):
	linesplit = line.strip().split("\t")
	if not line.startswith('rs'):
		continue
	CHR = linesplit[3]
	pos = linesplit[4]
	rsID = linesplit[0]
	A1 = linesplit[5]
	A2 = linesplit[6]

	if rsID in oneKG_rsID_dic:
		tmp = []
		marker_in_kg = rsID
		for bim in oneKG_rsID_dic[rsID]:
			tmp.append((bim.A1, bim.A2))
		if (A1, A2) in tmp:
			refAllele = A1
			altAllele = A2
		elif (A2, A1) in tmp:
			refAllele = A2
			altAllele = A1
		else:
			refAllele = "mismatch"
			altAllele = "mismatch"
	elif CHR+":"+pos in oneKG_chrpos_dic:
		tmp = {}
		for bim in oneKG_chrpos_dic[CHR+":"+pos]:
			tmp[bim.A1 + ':' + bim.A2] = bim.marker
		if A1 + ':' + A2 in tmp:
			refAllele = A1
			altAllele = A2
			marker_in_kg = tmp[A1 + ':' + A2]
		elif A2 + ':' + A1 in tmp:
			refAllele = A2
			altAllele = A1
			marker_in_kg = tmp[A2 + ':' + A1]
		else:
			refAllele = "mismatch2"
			altAllele = "mismatch2"
			marker_in_kg = "NA"
	else:
		refAllele = "mismatch3"
		altAllele = "mismatch3"
		marker_in_kg = "NA"

	outline = line.strip() + "\t" + "\t".join([refAllele, altAllele, marker_in_kg]) + '\n'
	out.write(outline)
out.close()
###################################
# search dbSNP and add ref, alt alleles: PDGWAS.SNPs.Proxy.region.SNPsPos.anno.txt.manuallycorrect
grep mismatch /home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.anno.txt
#didnot find mismatch

#awk 'BEGIN{OFS="\t"; print "SNP\tGWAS_SNP\tRegion\tChr\tPos\tRef\tAlt"}{if(NR == FNR){array[$1]=$11"\t"$12} if(NR > FNR) {if(array[$1]!=""){print $1,$2,$3,$4,$5,array[$1]}else{print $1,$2,$3,$4,$5,$8,$9} }}' \
#PDGWAS.SNPs.Proxy.region.SNPsPos.anno.txt.manuallycorrect PDGWAS.SNPs.Proxy.region.SNPsPos.anno.txt > PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final

awk 'BEGIN{OFS="\t"; print "SNP\tGWAS_SNP\tRegion\tChr\tPos\tRef\tAlt"}{print $1,$2,$3,$4,$5,$8,$9}' \
ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.txt > PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final

awk 'BEGIN{OFS=":"} {if(NR !=1) {print $4,$5,$6,$7}}' /home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final > PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final.PlinkExtraction

######################
# Feb. 2, 2018 fixed the duplicate variants
# Due to marker name disconcordance between GWAS SNPs and 1000G variants, there was ~20 variants are duplicate in 
# PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final, we need to remove them
# R code
library(tidyverse)
gwas_proxy_variants = read.table("/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final", header=T, sep = "\t", stringsAsFactors=F)
# read in 23931 variants
gwas_proxy_variants = gwas_proxy_variants %>% arrange(Chr, Pos, Ref, Alt, desc(GWAS_SNP)) %>% distinct(Chr, Pos, Ref, Alt, .keep_all=T)
# after remove one of the duplicates, 23,903 variants left
write.table(gwas_proxy_variants, file = "/home/tw83/twang/AMP/Rerun/PDGWAS/PDGWAS.SNPs.Proxy.region.SNPsPos.anno.filtered.final.rd", col.names=T, row.names=F, sep = "\t", quote=F)
##############

awk 'BEGIN{OFS=":"} {if(NR !=1) {print $4,$5,$6,$7}}' /home/jp360/MEGA/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.rd > ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.PlinkExtraction.rd








