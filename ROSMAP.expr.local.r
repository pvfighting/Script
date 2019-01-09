library(tidyverse)
counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

setwd("/Users/jiajiepeng/Documents/GitHub/Script")

############### Transform from counts to TPMs ##############
gene_length_file = "data/reference/hg38/gencode.v24.annotation.gtf.gene_length_unionexion.txt" # union exon length
gene_length =read.table(gene_length_file, header=F, sep = "\t", stringsAsFactors=F)
names(gene_length) = c("ENSGID", "length")
# remove gene name version
gene_length$ENSGID= do.call(rbind, strsplit(gene_length$ENSGID, '.', fixed=T))[,1]

ROSMAP_count_file = "data/ROSMAP/ROSMAP_all_counts_matrix.txt"
ROSMAP_count = read.table(ROSMAP_count_file, header=T, sep = "\t", stringsAsFactors=F, check.names=F)
ROSMAP_count = ROSMAP_count[-c(1:4),]
ROSMAP_count$feature = do.call(rbind, strsplit(ROSMAP_count$feature, '.', fixed=T))[,1]
comm = intersect(gene_length$ENSGID, ROSMAP_count$feature)

rownames(ROSMAP_count) = ROSMAP_count$feature; ROSMAP_count = ROSMAP_count[,-1]
rownames(gene_length) = gene_length$ENSGID; 
ROSMAP_count = ROSMAP_count[comm,]
gene_length = gene_length[comm,]
ROSMAP_tpm = counts_to_tpm(ROSMAP_count, gene_length$length)



############### Sample QC based on alignment summary metrics ###############
# 1) Keep samples with > 10 million mapped reads and 70% mappability with mapping quality of Q20 or higher signifying that the aligner estimates a 1/100 (or smaller) chance that the alignment is wrong. 
# remove 5 sample out of 639. 634 samples left 
map_summary = read.table("data/ROSMAP/ROSMAP_all_metrics_matrix.txt", header=T, sep = "\t", stringsAsFactors=F)
map_summary$HQ_mappability = map_summary$AlignmentSummaryMetrics__PF_HQ_ALIGNED_READS / map_summary$AlignmentSummaryMetrics__PF_READS
map_summary_keep = map_summary %>% filter(HQ_mappability>0.7 & AlignmentSummaryMetrics__PF_HQ_ALIGNED_READS > 10000000)

############### Select HC individuals #################
# 367 people left
ROSMAP_clinical = read.table("data/ROSMAP/ROSMAP_clinical_rename.txt", header=T, stringsAsFactors=F)
ROSMAP_clinical_HC = ROSMAP_clinical %>% filter(cogdx <=3)
ROSMAP_clinical_HC = ROSMAP_clinical_HC %>%rename(mrna_id = rename)%>% mutate(mrna_id = sub('X','', mrna_id))
tmp = do.call(rbind, strsplit(ROSMAP_clinical_HC$mrna_id, '_', fixed=T))
ROSMAP_clinical_HC$mrna_id_new = paste(tmp[,1], tmp[,2], sep = "_")
ROSMAP_clinical_HC$mrna_id_batch = tmp[,3]
comm_mrna_id = intersect(colnames(ROSMAP_tpm), ROSMAP_clinical_HC$mrna_id_new)# 370 people
comm_sample = intersect(map_summary_keep$sample, comm_mrna_id) # 367 people
ROSMAP_tpm = ROSMAP_tpm[, sort(comm_sample)] # Matrix, 60554*367, colnames:"03_120405"  "04_120405"


################ Fiter low-expressed genes ###############
# 1) keep genes with > 0.1 TPM in at least 20% of samples and >= 6 reads in at least 20% samples
ROSMAP_tpm = as.data.frame(ROSMAP_tpm)
keep_genes_idx <- (rowMeans(ROSMAP_tpm>0.1)>0.2) 
ROSMAP_tpm = ROSMAP_tpm[keep_genes_idx,] #38471 * 367
ROSMAP_count = ROSMAP_count[rownames(ROSMAP_tpm), colnames(ROSMAP_tpm)]
keep_genes_idx <- (rowMeans(ROSMAP_count>=6)>0.2)
ROSMAP_tpm = ROSMAP_tpm[keep_genes_idx,] 
############### Sample QC based on TPM table (before QC) ################
# RLE plot
library(ape)
library(reshape2)
message("generating RLE plot...")
logtpm = log10(ROSMAP_tpm + 1e-4)# so row value of 0 will be -4 in the transformed value
rle=logtpm-apply(logtpm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="TPM", id="ID")
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
outputfile="data/ROSMAP/QCFigures/samples.QC.plot.RLE.pdf"
pdf(outputfile, width=8, height=4)
par(mar=c(3,3,3,3))
boxplot(TPM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.3, main="Relative Log Expression", xlab="", ylab="RLE", frame=F)
abline(h=0, col='red',lwd=1)
dev.off()

# clustering
sampleDists = 1 - cor(ROSMAP_tpm, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
hcphy = as.phylo(hc)
co = ROSMAP_clinical_HC$mrna_id_batch[match(hcphy$tip.label, ROSMAP_clinical_HC$mrna_id_new)]
co = as.factor(co)
levels(co) = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
co = as.character(co)

pdf("data/ROSMAP/QCFigures/samples.QC.plot.cluster-hcp.pdf",width = 5, height = 5 )
par(mar=c(3,3,3,3))
plot(hcphy, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.)")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
dev.off()
# D statistics
pdf("data/ROSMAP/QCFigures/samples.QC.plot.Dstat.pdf")
D = apply(1-sampleDists, 1, median)
hist(D, breaks=100, ylab="Number of samples", xlab="D-statistic", main="Histogram of D-statistic")
#legend("topleft", paste(names(sort(D[which(D<0.9)])), round(sort(D[which(D<0.9)]),2)), bty='n')
dev.off()

# Gender mismatch
message("generating gender-match plot...")
pdf("data/ROSMAP/QCFigures/samples.QC.plot.Gender.pdf", width = 5, height = 5)
chrX='ENSG00000229807'  # XIST
chrY='ENSG00000129824'  # RPS4Y1
d=as.data.frame(t(logtpm[c(chrX,chrY),])); colnames(d)=c("chrX","chrY")
d$Sex = ROSMAP_clinical_HC$msex[match(rownames(d), ROSMAP_clinical_HC$mrna_id_new)]
plot(d$chrX, d$chrY, xlab="Expression of XIST", ylab="Expression of RPS4Y1", col= 'white',bg=ifelse(d$Sex==0,'red','blue'), pch=21, bty="n", main="Gender-specific expression")
#text(subset(d, chrX>0 & chrX<1.5 & chrY<0.2), rownames(subset(d, chrX>0 & chrX<1.5 & chrY<0.2)),pos=2, cex=0.5)
legend('bottomleft',pch=21,c("Female","Male"), col='white',pt.bg=c("red","blue"), bty='n', cex=.5)
dev.off()

################################### Adjust covariates###########
RIN_path = "data/ROSMAP/AMP_AD_ROSMAP_Broad-Rush_RNA_Seq_RIN.txt"
RIN = read.table(RIN_path, header = T, sep = "\t")
keep_genes_idx <- (rowSums(ROSMAP_tpm>0.1)>=10) 
expr= ROSMAP_tpm[keep_genes_idx,]
# logorithm
expr=log10(expr+1e-4)  # so row value of 0 will be -4 in the transformed value
# outlier correction: quantile normalization with order preserved. Now TPM is changed to rank normalized gene expression.
m=apply(expr, 1, mean); sd=apply(expr, 1, sd)
expr = t(apply(expr, 1, rank, ties.method = "average"));
#expr = qnorm(expr / (ncol(expr)+1));  # to standard normalization
expr = qnorm(expr / (ncol(expr)+1), mean=m, sd=sd)  # or, to preserve the mean and sd of each gene
rm(m,sd)

expr = as.data.frame(expr)

covs = ROSMAP_clinical_HC
rownames(covs) = covs$mrna_id_new
covs = covs[names(expr),]

covs$RIN = RIN$RINcontinuous[match(covs$mrna_id_new, RIN$Sampleid)]
covs$batch = as.factor(covs$mrna_id_batch)
covs$msex = as.factor(covs$msex)
covs$age_death = as.numeric(gsub('\\+', "", covs$age_death))
covs$RIN = as.numeric(covs$RIN)
covs$pmi = as.numeric(covs$pmi)
library(sva)

message("# adjusting expression with covariates...")
######################
#Combat remove batch effects

modcombat = model.matrix(~1, data = covs)
#01.09.2019 expr --- data.matrix(expr)   revised by Jiajie
combat_expr = ComBat(dat=data.matrix(expr), batch=covs$mrna_id_batch, mod = modcombat)
sampleDists = 1 - cor(expr, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
hcphy = as.phylo(hc)
co = ROSMAP_clinical_HC$mrna_id_batch[match(hcphy$tip.label, ROSMAP_clinical_HC$mrna_id_new)]
co = as.factor(co)
levels(co) = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
co = as.character(co)
pdf("samples.QC.plot.cluster-hcp.combat.pdf",width = 5, height = 5 )
par(mar=c(3,3,3,3))
# before combat 
plot(hcphy, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.) before combat")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
# after combat
sampleDists_combat = 1 - cor(combat_expr, method='spearman')
hc_combat = hclust(as.dist(sampleDists_combat),method = "complete")
hcphy_combat = as.phylo(hc_combat)
plot(hcphy_combat, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.) post combat")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
dev.off()

######################
# sva adjusted
Mod = model.matrix(~msex+RIN+age_death+pmi, data=covs) # full model(adjustment variables + variables of interest)
#Mod = model.matrix(~batch, data=covs)
Mod0 = model.matrix(~1,data=covs) # Null model (vairables of interest)

svaobj = sva(as.matrix(combat_expr),Mod, Mod0)
fsvaobj = fsva(dbdat=as.matrix(combat_expr),mod=Mod,sv=svaobj, newdat=as.matrix(combat_expr))
residuals = fsvaobj$db

message("# run RLE on SVA normalized quantification data ...")
######################

## RLE before and after SVA
#setwd("/home/tw83/twang/AMP/Rerun/RNAseq_Reprocessesing/ROSMAP/QCFigures/")
pdf("RLE.plot.sva.pdf", width=8, height=4)
res=data.frame(expr, check.names = F)
rle1=res-apply(res, 1, median) # before SVA
res=data.frame(residuals, check.names = F)
rle2=res-apply(res, 1, median) # after SVA
rle=melt(cbind(ID=rownames(rle1), rle1), variable.name = "Sample",value.name ="TPM", id="ID")
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(TPM ~ bymedian, data=rle, ylim=c(-4,4), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot before adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)
rle=melt(cbind(ID=rownames(rle2), rle2), variable.name = "Sample",value.name ="TPM", id="ID")
bymedian <- with(rle, reorder(Sample, TPM, IQR))  # sort by IQR
par(mar=c(7,3,3,1))
boxplot(TPM ~ bymedian, data=rle, ylim=c(-4,4), outline=F, las=2, boxwex=1, col='gray', cex.axis=0.5, main="RLE plot after adjustment", xlab="", ylab="Relative log expression (RLE)")
abline(h=0, col='red',lwd=1)
dev.off()

## mean gene expression distribution before and after SVA
pdf("expression.hist.sva.plot.pdf", width=8, height=8)
par(mfrow=c(2,1))
hist(apply(expr,1,mean), breaks=100, xlab="Rank normalized expression log10(TPM)", main="Expression distribution before adjustment")
hist(apply(residuals,1,mean), breaks=100, xlab="Rank normalized expression log10(TPM)", main="Expression distribution after adjustment")
dev.off()

## clustering before and after sva
sampleDists = 1 - cor(expr, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
hcphy = as.phylo(hc)
co = ROSMAP_clinical_HC$mrna_id_batch[match(hcphy$tip.label, ROSMAP_clinical_HC$mrna_id_new)]
co = as.factor(co)
levels(co) = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
co = as.character(co)
pdf("cluster-hcp.sva.plot.pdf",width = 8, height = 8 )
par(mar=c(3,3,3,3))
# before sva 
plot(hcphy, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.) before SVA")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
# after sva
sampleDists_sva = 1 - cor(residuals, method='spearman')
hc_sva = hclust(as.dist(sampleDists_sva),method = "complete")
hcphy_sva = as.phylo(hc_sva)
plot(hcphy_sva, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.) post SVA")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
dev.off()


message("# save final quantification data into file")
######################
write.table(format(residuals, digits=4,nsmall=4), file = "data/ROSMAP/expression.postSVA.xls", sep="\t", col.names = NA, quote=F,row.names = TRUE)
rownames(ROSMAP_clinical_HC) = ROSMAP_clinical_HC$mrna_id_new
ROSMAP_clinical_RNAseqSamples = ROSMAP_clinical_HC[colnames(residuals),]
write.table(ROSMAP_clinical_RNAseqSamples, file = "data/ROSMAP/ROSMAP_clinical_RNAseqQCedSamples.txt", sep="\t", col.names=TRUE, quote=F, row.names=FALSE)
