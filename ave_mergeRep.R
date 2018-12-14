#check and merge replicates of the same sample
expr_file = "/Users/jiajiepeng/workspace/p1/data/genes.fpkm.cufflinks.TCPY.uniq.xls"
expr = read.table(expr_file, header=T, check.names=F)
rownames(expr) = expr[,1]
expr = expr[,-1]
tempCol = colnames(expr)
tempCol = sort(tempCol)

subjectIDCol=gsub(".*_(.*)_.*_.*_rep.*", "\\1", tempCol)
tempSubjectID = ""
tempSampleRepSet = c()
exprNEW = expr
for (v in 1:length(tempCol)){
  if (v==1){
    tempSubjectID = subjectIDCol[v]
    tempSampleRepSet =c(tempCol[v])
   # print(tempSampleRepSet)
  }
  else{
    if(tempSubjectID == subjectIDCol[v]){
      tempSampleRepSet = append(tempSampleRepSet,tempCol[v])
      #print(tempCol[v])
      #print(tempSampleRepSet)
    }
    else{
      if(length(tempSampleRepSet)>1){
        combineSampleID = combineSampleIDs(tempSampleRepSet)
       # print(combineSampleID)
        tempMeans =  rowMeans(expr[,tempSampleRepSet])
        print(combineSampleID)
        exprNEW = cbind(exprNEW,tempMeans)
        exprNEW = exprNEW[,!names(exprNEW) %in% tempSampleRepSet]
        colnames(exprNEW)[length(exprNEW)]<-combineSampleID
        tempSubjectID = subjectIDCol[v]
        tempSampleRepSet=c(tempCol[v])
        
      }
      else{
        tempSubjectID = subjectIDCol[v]
        tempSampleRepSet =c(tempCol[v])
    #    print(tempSampleRepSet)
      }
    }
  }
  #print(tempCol[v])
}
#write.csv(exprNEW,file="/Users/jiajiepeng/workspace/p1/data/genes.fpkm.cufflinks.TCPY.uniq.mergeSameSample.csv",sep="\t")
write.table (data.frame("tracking_id"=rownames(exprNEW),exprNEW,check.names=FALSE), file ="/Users/jiajiepeng/workspace/p1/data/genes.fpkm.cufflinks.TCPY.uniq.mergeSameSample.txt", sep ="\t", row.names =FALSE, col.names =TRUE, quote =FALSE)
combineSampleIDs <- function(tempSampleRepSet){
  temp1 = gsub("(.*_.*_.*_.*)_rep.*", "\\1", tempSampleRepSet[1])
  temp2 = gsub(".*_.*_.*_.*(_rep.*)", "\\1", tempSampleRepSet[1])
  for(i in 2:length(tempSampleRepSet)){
    temp3 = gsub(".*_.*_.*_(.*)_rep.*", "\\1", tempSampleRepSet[i])
    temp1 = paste(temp1,temp3,sep="-")
  }
  temp1 = paste(temp1,temp2,sep="")
  #print(temp1)
  return(temp1)
}

