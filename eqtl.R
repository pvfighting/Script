library(MatrixEQTL)
wd = "/Users/jiajiepeng/workspace/p1/eqtl_data"
setwd(wd)
gene_loci = "/Users/jiajiepeng/workspace/p1/eqtl_data/gencode.v24lift37.annotation.gtf.gene.loci"
gene_expr = "/Users/jiajiepeng/workspace/p1/eqtl_data/HCAD_TCPY.expression.postSVA.eqtl.xls"
snps_loci = "/Users/jiajiepeng/workspace/p1/eqtl_data/SNPs.loci.txt"
snps_file = "/Users/jiajiepeng/workspace/p1/eqtl_data/DATA.eQTL.recode.raw.T"  
pca_cov_file = "/Users/jiajiepeng/workspace/p1/eqtl_data/DATA.HC.QCed.pruned.pca.evec.diag.cov" 


message("## loading SNP data...")
######################
snps = SlicedData$new();  # # GxN matrix
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;      # read file in slices of 5,000 rows
snps$LoadFile(snps_file);

## Load gene expression data
message("## loading expression data...")
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(gene_expr);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";
cvrt$fileOmitCharacters = "NA";
cvrt$fileSkipRows = 1;
cvrt$fileSkipColumns = 1;
if(length(pca_cov_file > 0))
{
  cvrt$LoadFile(pca_cov_file)
}

common = intersect(colnames(gene), colnames(snps))

snps$ColumnSubsample(match(common, colnames(snps)))
gene$ColumnSubsample(match(common, colnames(gene)))
cvrt$ColumnSubsample(match(common, colnames(cvrt)))

snpspos = read.table(snps_loci, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_loci, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = 'Final.trans.eQTL.P1e4.test.xls' , 
  pvOutputThreshold = 1e-4, #1e-8,
  useModel = modelLINEAR, 
  #useModel = modelLINEAR_CROSS, 
  errorCovariance = numeric(),
  verbose = FALSE,
  output_file_name.cis = "Final.cis.eQTL.P1.test.xls",
  pvOutputThreshold.cis = 1,  # no effect when min.pv.by.genesnp = TRUE
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = 1e6, #1e6
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

pdf(file="diagnostics_matrixeqlt.all.nodia.pdf", paper="usr")
plot(me, pch = 16, cex = 0.7);
dev.off()


message(paste("eQTL analysis is done and found", me$trans$neqtls,"trans SNP-gene pairs with significance p < 1e-4."))

message(paste("eQTL analysis is done and found", me$cis$neqtls,"cis SNP-gene pairs with significance p < 1."))