cd /PHShome/jq34/

ln -s /data/neurogen/MEGA/Plate44-46-May26-2017/MEGAP44-46.GC015.fam /PHShome/jq34/MEGAP44-46.GC015.fam
ln -s /data/neurogen/MEGA/Plate44-46-May26-2017/MEGAP44-46.GC025.fam /PHShome/jq34/MEGAP44-46.GC025.fam
cut -d ' ' -f 2 /PHShome/jq34/MEGAP44-46.GC025.fam > /PHShome/jq34/MEGAP44-46_sampleID.GC025.txt
python ./code/sepFirstline.py -i ./genes.fpkm.cufflinks.TCPY.uniq.mergeSameSample.txt -o genes.fpkm.cufflinks.TCPY.uniq.names.txt

cut -d '_' -f 2 /PHShome/jq34/genes.fpkm.cufflinks.TCPY.uniq.names.txt > /PHShome/jq34/temp_sampleID1.txt
cut -d 'P' -f1 /PHShome/jq34/MEGAP44-46_sampleID.GC025.txt | sed 's/.$//' > /PHShome/jq34/temp_sampleID2.txt
cut -d 'P' -f2 /PHShome/jq34/MEGAP44-46_sampleID.GC025.txt | cut -d '_' -f1 > /PHShome/jq34/temp_sampleID2_plate.txt


sort /PHShome/jq34/temp_sampleID1.txt /PHShome/jq34/temp_sampleID2.txt | uniq -d > /PHShome/jq34/SamplesGenoExpre.txt
cat /PHShome/jq34/genes.fpkm.cufflinks.TCPY.uniq.names.txt | grep -Ff /PHShome/jq34/SamplesGenoExpre.txt | cut -d '_' --output-delimiter=' '  -f 1,2 > samples.keep.txt