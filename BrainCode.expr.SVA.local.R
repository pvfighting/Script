setwd("/Users/jiajiepeng/Documents/GitHub/Script")

library(reshape2)
library(sva)
library(ape)
SAMPLE_GROUP="HCILBPD_SNDA"
samplelist=read.table(paste0("data/samplelist.",SAMPLE_GROUP), header = F, stringsAsFactors = FALSE)[,1]
#list the samples we used 
expr_file = "data/genes.fpkm.cuffnorm.allSamples.uniq.xls"
covariance_file = "data/BRAINCODE_Sequencing_Log_covariates_table.tsv"
expr = read.table(expr_file, header=T, check.names=F)
rownames(expr) = expr[,1]; expr = expr[,-1];
#read the names of transcripts
#colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr)) 
covs = read.delim(covariance_file, stringsAsFactors =F)

covs=covs[match(samplelist, covs$sampleName), ]
#select the covs of samles we used in the following analysis
expr=subset(expr, select = samplelist)
 # change the sample ID to subject ID
colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr)) #this step may update based on different name format
  # convert to factor
covs$batch = as.factor(covs$batch)
covs$sex = as.factor(covs$sex)
covs$readsLength = as.factor(covs$readsLength)
message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))
# filtering expression data... keep the genes with expression values larger than 0.1 in at least 10 smaples.
expr=expr[rowSums(expr>0.1)>=10,]
message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))

############### Sample QC based on FPKM table (before QC) ################
# RLE plot
library(ape)
library(reshape2)
message("generating RLE plot...")
logfpkm = log10(expr + 1e-4)# so row value of 0 will be -4 in the transformed value
rle=logfpkm-apply(logfpkm, 1, median) # change "/" to "-" so that we got log(fold-change) which centered on 0 on the RLE plot.
rle=melt(cbind(ID=rownames(rle), rle), variable.name = "Sample",value.name ="FPKM", id="ID")
bymedian <- with(rle, reorder(Sample, FPKM, IQR))  # sort by IQR
outputfile="QCFigures/BrainCode/samples.QC.plot.RLE.pdf"
pdf(outputfile, width=8, height=4)
par(mar=c(3,3,3,3))
boxplot(FPKM ~ bymedian, data=rle, outline=F, las=2, boxwex=1, col='gray', cex.axis=0.3, main="Relative Log Expression", xlab="", ylab="RLE", frame=F)
abline(h=0, col='red',lwd=1)
dev.off()

## clustering
sampleDists = 1 - cor(expr, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
hcphy = as.phylo(hc)
co = covs$batch[match(hcphy$tip.label, covs$subjectID)]
co = as.factor(co)
levels(co) = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
co = as.character(co)

pdf("QCFigures/BrainCode/samples.QC.plot.cluster-hcp.pdf",width = 5, height = 5 )
par(mar=c(3,3,3,3))
plot(hcphy, tip.col = co, type = "unrooted", cex=.2, lab4ut='axial',underscore = T, main="Clustering of samples (Spearman - Cor.)")
Xcol = c("gray0", "red", "orange", "yellow", "green", "blue", "purple", "wheat", "violet")
Xtext = c("batch 0", "batch 1","batch 2","batch 3","batch 4","batch 5","batch 6","batch 7","batch 8")
legend('bottomleft',pch=21,Xtext, col='white',pt.bg=Xcol, cex=.5)
dev.off()

message(" # transforming RPKM to rank normalized gene expression ...")
  ######################
  # logorithm
  expr=log10(expr+1e-4)  # so row value of 0 will be -4 in the transformed value
  # outlier correction: quantile normalization with order preserved. Now RPKM is changed to rank normalized gene expression.
  m=apply(expr, 1, mean); sd=apply(expr, 1, sd)
  expr = t(apply(expr, 1, rank, ties.method = "average"));
  #expr = qnorm(expr / (ncol(expr)+1));  # to standard normalization
  expr = qnorm(expr / (ncol(expr)+1), mean=m, sd=sd)  # or, to preserve the mean and sd of each gene
  
  rm(m,sd)

# sva adjusted
Mod = model.matrix(~batch+sex+RIN+age+PMI, data=covs)
#Mod0 = model.matrix(~1,data=covs)

## estimate the latent factors
#n.sv = num.sv(expr, Mod, method="leek")
#svaobj = sva(as.matrix(expr),Mod, Mod0, n.sv=n.sv)
svaobj = sva(as.matrix(expr),Mod)
fsvaobj = fsva(dbdat=as.matrix(expr),mod=Mod,sv=svaobj, newdat=as.matrix(expr))

residuals = fsvaobj$db
message("# save final quantification data into file")
######################
write.table(format(residuals, digits=4,nsmall=4), file = "data/HCILBPD_SNDA.expression.postSVA.xls", sep="\t", col.names = NA, quote=F,row.names = TRUE)

## mean gene expression distribution before and after SVA
pdf("expression.hist.sva.plot.pdf", width=8, height=8)
par(mfrow=c(2,1))
hist(apply(expr,1,mean), breaks=100, xlab="Rank normalized expression FPKM", main="Expression distribution before adjustment")
hist(apply(residuals,1,mean), breaks=100, xlab="Rank normalized expression FPKM", main="Expression distribution after adjustment")
dev.off()

## RLE before and after SVA
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

## clustering before and after sva
sampleDists = 1 - cor(expr, method='spearman')
hc=hclust(as.dist(sampleDists),method = "complete")
hcphy = as.phylo(hc)
co = covs$batch[match(hcphy$tip.label, covs$subjectID)]
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
