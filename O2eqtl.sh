cd /home/jp360/ADEQTL

ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.fam PDMAP.ALL.QC.Done.fam
ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.bim PDMAP.ALL.QC.Done.bim
ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.bed PDMAP.ALL.QC.Done.bed

cut -d ' ' -f 2 PDMAP.ALL.QC.Done.fam > PDMAP.ALL.QC.Done_sujectID.txt
python ./sepFirstline.py -i ./genes.fpkm.cufflinks.TCPY.uniq.mergeSameSample.txt -o genes.fpkm.cufflinks.TCPY.uniq.names.txt
cut -d '_' -f 2 genes.fpkm.cufflinks.TCPY.uniq.names.txt > genes.fpkm.cufflinks.TCPY.uniq.subjectID.txt
sort PDMAP.ALL.QC.Done_subjectID.txt genes.fpkm.cufflinks.TCPY.uniq.subjectID.txt | uniq -d > OverlapSubjectID_PDMAP-TCPY.txt
cat genes.fpkm.cufflinks.TCPY.uniq.names.txt | grep -Ff OverlapSubjectID_PDMAP-TCPY.txt | cut -d '_' --output-delimiter=' '  -f 1,2 > OverlapSubjectID_PDMAP-TCPY.keep.txt

#OverlapSubjectID_PDMAP-TCPY.keep.txt extract the overlap subjects in PDMAP dataset and TCPY dataset
plink --bfile PDMAP.All.QC.Done --keep OverlapSubjectID_PDMAP-TCPY.keep.txt --make-bed --out PDMAP.TCPY 
## R code to rename fam file
fam = read.table("DATA.fam", header=F, sep = " ", stringsAsFactors=F)
names_map = read.table("name.change.tab", header= F, sep = "\t", stringsAsFactors=F)
names_map$V3 = do.call(rbind, strsplit(names_map$V1, '_', fixed=T))[,2]
fam$V2 = names_map$V2[match(fam$V2, names_map$V3)]
write.table(fam, file = "DATA.fam.new", col.names=F, row.names=F, sep = " ", quote=F)
## End R
#
plink --bfile PDMAP.TCPY --freq --out PDMAP.TCPY
#By itself, --freq writes a minor allele frequency report to plink.frq. 

######### Run 1000G checking tool ####  Imputation preparation and checking #this command could generate the Run-plink.sh file for 
bsub -q priority -W 120:00 -M 90000 perl /home/tw83/twang/AMP/tools/HRC-1000G-check-bim.pl \
 -b PDMAP.TCPY.bim \
 -f PDMAP.TCPY.frq \
 -r /home/tw83/twang/AMP/tools/1000GP_Phase3_combined.legend \
 -g -p 'EUR'

sh Run-plink.sh 
#this separate the SNP by chromosomes~  Force a particular reference (A1) allele
#where are the files like Exclude-PDMAP.All.QC.Done-HRC.txt Chromosome-PDMAP.All.QC.Done-HRC.txt 
#Position-PDMAP.All.QC.Done-HRC.txt Strand-Flip-PDMAP.All.QC.Done-HRC.txt Force-Allele1-PDMAP.All.QC.Done-HRC.txt from?
# These files are the output of last command. What is the difference between -HRC and -1000G ? how to set it ?  Ask tao about this. 

mkdir byChr
mv  mv DATA-updated-* byChr/
cd byChr
## Recode Plink to vcf
#Plink2VCF.sh
#!/bin/bash
#4.2 Recode Plink to vcf
for (( i = 1; i <= 23; i++ )); do
	plink --bfile DATA-updated-chr$i --recode vcf --out DATA-updated-chr$i
done
	#4.3 Sort and Compress vcf file
[[ ! -d ./vcfgz ]] && mkdir vcfgz
for (( i = 1; i <= 23; i++ )); do
	vcf-sort DATA-updated-chr${i}.vcf | bgzip -c > ./vcfgz/chr${i}.vcf.gz
done

# Job submission
sftp://orchestra.med.harvard.edu/home/tw83/twang/AMP/BrainCode/1000G/byChr/vcfgz

########################## prepare the vcf files
# unzip.sh
#! /usr/bin/bash
password=yQ7wZ/8SswWyIE
for i in $(seq 1 22);do
	7za e chr_${i}.zip -p$password
done

cd /home/tw83/twang/AMP/BrainCode/1000G/Imputation_1000G/
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

