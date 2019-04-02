#!/usr/local/bin/perl  -w
# If the dataset contains heterozygous hapliod genotypes problem, so first clean data using --set-hh-missing
#/home/jp360/MEGA/P60_ADHC
plink --bfile MEGAP60.GC025.ADHC.zero --set-hh-missing --make-bed --out DATA.CleanSNPs

#Total genotyping rate is 0.995205.
#1755868 variants and 81 people pass filters and QC.
####### Step1 remove subjects with call rates < 95%
plink --bfile DATA.CleanSNPs --mind 0.05 --out DATA.filterSubject --make-bed
#1755868 variants and 78 people pass filters and QC.

####### Step2 remove Sex unmatched ## 
plink --bfile DATA.filterSubject --check-sex --out DATA.filterSubject.Sex 
#1765741 variants and 78 people pass filters and QC.
#1755868 variants and 78 people pass filters and QC.
####### Step3 remove SNP calling rate less than 95%
plink --bfile DATA.filterSubject --geno 0.05 --make-bed --out DATA.fileterSNPs 
#1734976 variants and 78 people pass filters and QC.

####### Step4 Remove Hardy-Weinberg equilibrium SNPs
plink --bfile DATA.fileterSNPs --hwe 1e-6 --out DATA.Hardy --make-bed 
#--hwe: 10 variants removed due to Hardy-Weinberg exact test.
#1734967 variants and 78 people pass filters and QC.

####### Step5 Remove Test-mishap SNPs
#Total genotyping rate is 0.998947.
#1744658 variants and 78 people pass filters and QC.
#ote: No phenotypes present.
#--test-mishap: 0 loci checked (1744658 skipped).
plink --bfile DATA.Hardy --test-mishap --out DATA.Mis 
awk 'BEGIN{OFS="\t"}{if($8 < 1e-9){print $1,$2,$3,$4,$5,$6,$7,$8,$9 }}' DATA.Mis.missing.hap |  grep 'HETERO' | cut -f1 | sort | uniq > DATA.mishap.txt
plink --bfile DATA.Hardy --exclude DATA.mishap.txt --out DATA.Mis --make-bed 
#1734967 variants and 78 people pass filters and QC.
####### Step6 Remove SNPs with MAF < 0.01
plink --bfile DATA.Mis --maf 0.05 --out DATA.MAF --make-bed
#933568 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#801399 variants and 78 people pass filters and QC.

###### Step7, Step8
plink --bfile DATA.MAF --indep 50 5 2 --out DATA.MAF  
#Pruning complete.  588915 of 801399 variants removed.
plink --bfile DATA.MAF --extract DATA.MAF.prune.in --make-bed --out DATA.MAF.subset
#remove linkage disequilibrium
#212484 variants and 78 people pass filters and QC

####### Step7 Remove # PLINK Heterozygosity subjects
##### remove F outlier
plink --bfile DATA.MAF.subset --het --out DATA.Het
##analysis F column of DATA.Het.het using R, generated subjects file call DATA.QC.outlier.txt with +/- 4sd of F value;
# 2 individuals were removed

awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' DATA.Het.het > DATA.Het.het1 ## reformat
##R code
data = read.table('DATA.Het.het1', header = T, sep = "\t", stringsAsFactors = F)
sd = sd(data$F)
mean = mean(data$F)
data_normal = subset(data, F >= mean-4*sd); data_normal = subset(data_normal, F <= mean+4*sd)
outlier = setdiff(data$IID, data_normal$IID)
data_outlier = subset(data, IID %in% outlier)
write.table(data_outlier, "DATA.QC.outlier.txt", quote = F, sep = "\t", row.names = F)

awk 'NR !=1 {print $1, $2}' DATA.QC.outlier.txt > DATA.QC.outlier.txt1
mv DATA.QC.outlier.txt1 DATA.QC.outlier.txt

######## Step8 IBS/IBD Filtering subjects

plink --bfile DATA.MAF.subset --genome --out DATA.MAF.subset.relatedness
#212484 variants and 78 people pass filters and QC.
##may have problem in this step
awk '{if($10>0.1875){print $1,$2}}' DATA.MAF.subset.relatedness.genome | sed 1d >> DATA.QC.outlier.txt

awk '{if($10>0.98){print $3,$4}}' DATA.MAF.subset.relatedness.genome | sed 1d >> DATA.QC.outlier.txt

####### removed subjects in Step8 and Step9 together
plink --bfile DATA.MAF --remove DATA.QC.outlier.txt --make-bed --out DATA.HC.QCed
#801399 variants and 67 people pass filters and QC.

# DATA.HC.QCed.bed/bim/fam # 335 people, have genotype
plink -bfile DATA.HC.QCed --indep 50 5 2 --out DATA.HC.QCed.prune
plink -bfile DATA.HC.QCed --extract DATA.HC.QCed.prune.prune.in --make-bed --out DATA.HC.QCed.pruned
#204796 variants and 67 people pass filters and QC.
cp DATA.HC.QCed.pruned.bim DATA.HC.QCed.pruned.pedsnp
cp DATA.HC.QCed.pruned.fam DATA.HC.QCed.pruned.pedind
###add a step to remove the snp with longname
perl /home/jp360/EIG-master/bin/smartpca.perl -i DATA.HC.QCed.pruned.bed -a DATA.HC.QCed.pruned.pedsnp -b DATA.HC.QCed.pruned.pedind -o DATA.HC.QCed.pruned.pca -p DATA.HC.QCed.pruned.plot -e DATA.HC.QCed.pruned.eval -l DATA.HC.QCed.pruned.smartpca.log -m 5 -k 3 -t 2 
#perl /home/jp360/EIG-master/bin/smartpca.perl -i DATA.HC.QCed.pruned.bed -a DATA.HC.QCed.pruned.pedsnp -b DATA.HC.QCed.pruned.pedind -o DATA.HC.QCed.pruned.pca -p DATA.HC.QCed.pruned.plot -e DATA.HC.QCed.pruned.eval -l DATA.HC.QCed.pruned.smartpca.log -m 5 -k 3 -t 2 
 # Top 3 pca is in DATA.HC.QCed.pruned.pca.evec
awk '/REMOVED/ {print $3}' DATA.HC.QCed.pruned.smartpca.log | awk 'BEGIN{FS=":";OFS="\t"}{print $1,$2}' > DATA.HC.QCed.pruned.smartpca.outlier.txt
plink --bfile DATA.HC.QCed --remove DATA.HC.QCed.pruned.smartpca.outlier.txt --make-bed --out DATA.HC.QCed.Pass
#801399 variants and 67 people pass filters and QC.
 ## R code: format top 3 genotyping principle components
setwd("/home/jp360/MEGAtest")
data = read.table("DATA.HC.QCed.longname.pruned.pca.evec", header = F, stringsAsFactors = F)
colnames(data) = c("ID", "PCA_1", "PCA_2", "PCA_3","V5")
data = subset(data, select=c("ID", "PCA_1", "PCA_2", "PCA_3"))
data$ID = do.call(rbind, strsplit(data$ID, ":", fixed=T))[,1]
#WGS_key = read.csv("/home/tw83/twang/AMP/ROSMAP/WGS/VCF/AMP-AD_rosmap_WGS_id_key.csv", header=T, stringsAsFactors=F, check.names=F)
#data$ID = WGS_key$projid[match(data$ID, WGS_key$WGS_id)]
write.table(t(data), file = "DATA.HC.QCed.longname.pruned.pca.evec.cov", quote=F, sep = "\t", row.names=T, col.names=F)
## End R code

#here here here...
############ Run 1000G checking tool################
plink --bfile DATA.HC.QCed.Pass --freq --out DATA.HC.QCed.Pass

bsub -q priority -W 120:00 -M 90000 perl /home/tw83/twang/AMP/tools/HRC-1000G-check-bim.pl \
 -b DATA.HC.QCed.Pass.bim \
 -f DATA.HC.QCed.Pass.frq \
 -r /home/tw83/twang/AMP/tools/1000GP_Phase3_combined.legend \
 -g -p 'EUR'

sh Run-plink.sh #generate by the last step
mkdir byChr
mv DATA.HC.QCed.Pass-updated-* byChr/
cd byChr
## Recode Plink to vcf
#Plink2VCF.sh
#!/bin/bash
#4.2 Recode Plink to vcf
for (( i = 1; i <= 22; i++ )); do
	plink --bfile DATA.HC.QCed.Pass-updated-chr$i --recode vcf --out DATA.HC.QCed.Pass-updated-chr$i
done
	#4.3 Sort and Compress vcf file
[[ ! -d ./vcfgz ]] && mkdir vcfgz
for (( i = 1; i <= 22; i++ )); do
	vcf-sort DATA.HC.QCed.Pass-updated-chr${i}.vcf | bgzip -c > ./vcfgz/chr${i}.vcf.gz
done

# Job submission to Michigan imputation server
##########################
# unzip.sh
#! /usr/bin/bash
password=WAa3teHMFx/3pD
for i in $(seq 1 22);do
	7za e chr_${i}.zip -p$password
done

cd /home/jp360/MEGA/P60_ADHC/Imputation_1000G
mkdir plinkout
##########################
# VCFQC.sh
#!/bin/bash
chnum=$1
plink --vcf chr${chnum}.dose.vcf.gz --make-bed --out s1_${chnum} --const-fid
plink --bfile s1_${chnum} --bmerge s1_${chnum} --merge-mode 6 --out missnp_${chnum}
plink --bfile s1_${chnum} --exclude missnp_${chnum}.missnp --make-bed --out s2_${chnum}
plink --bfile s2_${chnum} --list-duplicate-vars --out dupvar_${chnum}
plink --bfile s2_${chnum} --exclude dupvar_${chnum}.dupvar --make-bed --out s3_${chnum}
gunzip -d chr${chnum}.info.gz
plink --bfile s3_${chnum} --qual-scores chr${chnum}.info 7 1 1 --qual-threshold 0.3 --make-bed --out ./plinkout/QCed_chr$chnum
rm missnp_${chnum}.*
rm dupvar_${chnum}.*
rm s1_${chnum}.*
rm s2_${chnum}.*
rm s3_${chnum}.*
#############
bsub -q short -W 1:00 sh VCFQC.sh 1
bsub -q short -W 1:00 sh VCFQC.sh 2
for chnum in $(seq 2 22);
  do
  bsub -q short -W 1:00 sh VCFQC.sh $chnum
done
cd ./plinkout
echo QCed_chr2.bed QCed_chr2.bim QCed_chr2.fam > merge.list
for chnum in $(seq 3 22);do
	echo QCed_chr${chnum}.bed QCed_chr${chnum}.bim QCed_chr${chnum}.fam >> merge.list
done
plink --bfile QCed_chr1 --merge-list merge.list --make-bed --out imputepostqc
###########################



######### Genotyping PCA figure with hapmap ######## 
ln -s DATA.HC.QCed.longname.pruned.bim DATA.HC.QCed.longname.pruned.pca.bim
ln -s DATA.HC.QCed.longname.pruned.bed DATA.HC.QCed.longname.pruned.pca.bed
ln -s DATA.HC.QCed.longname.pruned.fam DATA.HC.QCed.longname.pruned.pca.fam
ln -s /home/tw83/twang/AMP/hapmap3/hg19/Hapmap3.QC.bim  Hapmap3.QC.bim
ln -s /home/tw83/twang/AMP/hapmap3/hg19/Hapmap3.QC.bed  Hapmap3.QC.bed
ln -s /home/tw83/twang/AMP/hapmap3/hg19/Hapmap3.QC.fam  Hapmap3.QC.fam
awk 'BEGIN{OFS="\t"} {print $1, $1":"$4":"$5":"$6, $3,$4,$5,$6}' Hapmap3.QC.bim > Hapmap3.QC.bim.newname
mv Hapmap3.QC.bim Hapmap3.QC.bim.source
mv Hapmap3.QC.bim.newname Hapmap3.QC.bim
awk 'BEGIN{OFS="\t"} {print $1, $1":"$4":"$5":"$6, $3,$4,$5,$6}' DATA.HC.QCed.longname.pruned.pca.bim > DATA.HC.QCed.longname.pruned.pca.bim.newname
mv DATA.HC.QCed.longname.pruned.pca.bim DATA.HC.QCed.longname.pruned.pca.bim.source
mv DATA.HC.QCed.longname.pruned.pca.bim.newname DATA.HC.QCed.longname.pruned.pca.bim
awk 'BEGIN{OFS="\t"} NR == FNR {array[$2"\t"$5"\t"$6]=1} NR > FNR {if(array[$2"\t"$5"\t"$6]==1){print $2}}' Hapmap3.QC.bim DATA.HC.QCed.longname.pruned.pca.bim > ShareSNP.txt
plink -bfile DATA.HC.QCed.longname.pruned.pca --extract ShareSNP.txt --make-bed --out DATA.HC.QCed.longname.pruned.pca.Share
####### merge DATA.QC and HapMAP3
plink --bfile Hapmap3.QC --extract ShareSNP.txt --make-bed --out Hapmap3.QC.Share 
plink --bfile DATA.HC.QCed.longname.pruned.pca.Share --bmerge Hapmap3.QC.Share.bed Hapmap3.QC.Share.bim Hapmap3.QC.Share.fam --make-bed --out DATA.Hapmap3
plink --bfile DATA.Hapmap3 --indep 50 5 2 --out DATA.Hapmap3.prune 
plink -bfile DATA.Hapmap3 --extract DATA.Hapmap3.prune.prune.in --make-bed --out DATA.Hapmap3.pruned 
cp DATA.Hapmap3.pruned.bim DATA.Hapmap3.pruned.pedsnp
cp DATA.Hapmap3.pruned.fam DATA.Hapmap3.pruned.pedind
awk '{print $1,$2,$3,$4,$5,-9}' DATA.Hapmap3.pruned.pedind > DATA.Hapmap3.pruned.pedind9
perl /home/tw83/bin/EIG6.1.1/bin/smartpca.pl -i DATA.Hapmap3.pruned.bed -a DATA.Hapmap3.pruned.pedsnp -b DATA.Hapmap3.pruned.pedind9 -o DATA.Hapmap3.pruned.pca -p DATA.Hapmap3.pruned.plot -e DATA.Hapmap3.pruned.eval -l DATA.Hapmap3.smartpca.log -m 5 -k 2 -t 2

