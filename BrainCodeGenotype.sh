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