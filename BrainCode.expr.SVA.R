library(reshape2)
library(sva)

SAMPLE_GROUP="HCILBPD_SNDA"
samplelist=read.table(paste0("/home/tw83/twang/AMP/BrainCode/geneExpr/samplelist.",SAMPLE_GROUP), header = F, stringsAsFactors = FALSE)[,1]
#list the samples we used 
expr_file = "/home/tw83/twang/AMP/BrainCode/geneExpr/genes.fpkm.cuffnorm.allSamples.uniq.xls"
covariance_file = "/home/tw83/twang/AMP/BrainCode/geneExpr/BRAINCODE_Sequencing_Log_covariates_table.tsv"
expr = read.table(expr_file, header=T, check.names=F)
rownames(expr) = expr[,1]; expr = expr[,-1];
#read the names of transcripts

covs = read.delim(covariance_file, stringsAsFactors =F)

covs=covs[match(samplelist, covs$sampleName), ]
#select the covs of samles we used in the following analysis
expr=subset(expr, select = samplelist)
 # change the sample ID to subject ID
colnames(expr)=gsub(".*_(.*)_.*_.*_rep.*", "\\1", colnames(expr))
  # convert to factor
covs$batch = as.factor(covs$batch)
covs$sex = as.factor(covs$sex)
covs$readsLength = as.factor(covs$readsLength)
message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))
# filtering expression data...
expr=expr[rowSums(expr>0.1)>=10,]
message(paste(" -- now expression matrix has",nrow(expr),"rows and",ncol(expr),"columns"))

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
write.table(format(residuals, digits=4,nsmall=4), file = "/home/tw83/twang/AMP/BrainCode/geneExpr/HCILBPD_SNDA.expression.postSVA.xls", sep="\t", col.names = NA, quote=F,row.names = TRUE)

