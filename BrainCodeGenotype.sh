cd /home/jp360/BrainCodeEQTL

ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.fam PDMAP.ALL.QC.Done.fam
ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.bim PDMAP.ALL.QC.Done.bim
ln -s /home/tw83/twang/AMP/BrainCode/PDMAP.All.QC.Done.bed PDMAP.ALL.QC.Done.bed
#PDMAP.ALL.QC.Done is the brainCode QCed data. In the TCPY project, we should got this data after QC the Mega chip TCPY data.
mkdir 1000G
cd 1000G
ln -s /home/tw83/twang/AMP/BrainCode/Imputation/QC/84samples/name.change.tab name.change.tab
cut -f 1 name.change.tab | cut -d '_' --output-delimiter=' '  -f 1,2 > 84samples.keep.txt
# 84 samples.keep.txt. "HC BN10-39" This file is generated for selecting data from PDMAP.ALL.QC.Done (plink format)
plink --bfile ../PDMAP.ALL.QC.Done --keep 84samples.keep.txt --make-bed --out DATA
## R code to rename fam file
fam = read.table("DATA.fam", header=F, sep = " ", stringsAsFactors=F)
names_map = read.table("name.change.tab", header= F, sep = "\t", stringsAsFactors=F)
names_map$V3 = do.call(rbind, strsplit(names_map$V1, '_', fixed=T))[,2]
fam$V2 = names_map$V2[match(fam$V2, names_map$V3)]
write.table(fam, file = "DATA.fam.new", col.names=F, row.names=F, sep = " ", quote=F)
## End R. What is this step for ?
mv DATA.fam DATA.fam.source
mv DATA.fam.new DATA.fam

#--freq writes a minor allele frequency report to plink.frq.
######### Run 1000G checking tool ####
bsub -q priority -W 2:00 -M 90000 perl /home/tw83/twang/AMP/tools/HRC-1000G-check-bim.pl \
 -b ./DATA.bim \
 -f ./DATA.frq \
 -r /home/tw83/twang/AMP/tools/1000GP_Phase3_combined.legend \
 -g -p 'EUR'

