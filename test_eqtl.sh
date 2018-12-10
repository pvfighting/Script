cd /PHShome/jq34/code

ln -s /data/neurogen/MEGA/Plate44-46-May26-2017/MEGAP44-46.GC025.fam /PHShome/jq34/MEGAP44-46.GC025.fam
cut -d ' ' -f 2 /PHShome/jq34/MEGAP44-46.GC025.fam > /PHShome/jq34/MEGAP44-46_sampleID.GC025.txt
python sepFirstline.py -i ../genes.fpkm.cufflinks.TCPY.uniq.xls -o ../genes.fpkm.cufflinks.TCPY.uniq.names.txt

cut -d '_' -f 2 /PHShome/jq34/genes.fpkm.cufflinks.TCPY.uniq.names.txt > /PHShome/jq34/temp_sampleID1.txt
cut -d 'P' -f1 /PHShome/jq34/MEGAP44-46_sampleID.GC025.txt | sed 's/.$//' > /PHShome/jq34/temp_sampleID2.txt


