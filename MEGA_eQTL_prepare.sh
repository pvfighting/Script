cd /home/jp360/MEGA/P60_ADHC/eQTL_prepare
ln -s /home/jp360/MEGA/P60_ADHC/Imputation_1000G/plinkout/imputepostqc.bed  merged.bed
ln -s /home/jp360/MEGA/P60_ADHC/Imputation_1000G/plinkout/imputepostqc.bim  merged.bim
ln -s /home/jp360/MEGA/P60_ADHC/Imputation_1000G/plinkout/imputepostqc.fam  merged.fam

plink --bfile merged --maf 0.05 --make-bed --out DATA
#6709258 variants and 76 people pass filters and QC.

bsub -q priority -W 120:00 -M 90000 python _ReNameSNPs.1000G.py   DATA.bim   DATA.renamed.bim
ln -s /home/jp360/MEGA/P60_ADHC/eQTL_prepare/DATA.bed  DATA.renamed.bed

###########generate DATA.renamed.fam#########

#################################################
## extract PD risk SNPs for eqtl

bsub -q priority -W 120:00 -M 90000 plink --bfile DATA.renamed --extract /home/jp360/MEGA/ADGWAS/ADGWAS.SNPs.Proxy.region.rmPos.SNPsPos.anno.filtered.final.PlinkExtraction.rd \
--make-bed --out DATA.eQTL
# 76 people (30 males, 46 females) loaded from .fam.--extract: 8602 variants remaining.
cut -f 2 DATA.eQTL.bim | awk 'BEGIN{FS=":"}{print $1":"$2":"$3":"$4,$4}'  > DATA.eQTL.bim.recode_allele
plink --bfile DATA.eQTL --recode A --recode-allele DATA.eQTL.bim.recode_allele  --out DATA.eQTL.recode
rowsToCols DATA.eQTL.recode.raw  DATA.eQTL.recode.raw.T

sed -i -e '1d' -e '3,6d' DATA.eQTL.recode.raw.T
# Done genotype for eQTL
# SNP loci for eQTL

cut -f 2 DATA.eQTL.bim | awk 'BEGIN{FS=":"}{print $1":"$2":"$3":"$4"_"$4, "chr"$1, $2}' > SNPs.loci.txt
sed -i '1isnp\tchr\tpos' SNPs.loci.txt