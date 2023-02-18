
#' Title Filter out gene symbols that don't exist in the DEG (weight) gene symbols
#'
#' @param x pre-compiled molecular concepts. Each element in the molecular concepot list has a couple of dozens genes.
#' @param WeightGeneSymbol Gene symbols of weight genes
#' @param inlist Default is TRUE.
#'
#' @return Moleculare concept lists that contains gene symbols common with weight genes
#' @export
#'
#' @examples filterGeneSymb(x=feature.list, WeightGeneSymbol=names(weight),inlist=TRUE)
filterGeneSymb <- function (x,WeightGeneSymbol,inlist) {
  if (inlist){
    x <- x[x %in% WeightGeneSymbol]
  }
  else{
    x <- x[!x %in% WeightGeneSymbol]
  }
  return (x)
}


#' Title normalize the weight of genes
#'
#' @param x A vector of gene weight
#'
#' @return normalized weight
#' @export
#'
#' @examples normalize(weight)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


#' Title calculate ES based on gene weight (sorted in decreasing order)
#'
#' @param weight A vector of gene weight
#' @param posList index of genes in the concept list that are common with weighted genes.
#'
#' @return Calculated enrichment score
#' @export
#'
#' @examples cal_ES(weight, posList=VectorGeneSymbol)
cal_ES<-function(weight,posList){
  if (sum(weight<0)>0) {
    stop (paste("negative values found in weight, aborting calculation of ES",weight,sep=" "))#there cannot be nagative values in weight
  }
  weight<-sort(weight,decreasing=TRUE)
  onIndice<-which(names(weight) %in% posList)
  stepDown<-1/(length(weight)-length(onIndice))
  norWeight<-rep(-stepDown,length(weight))
  names(norWeight)<-names(weight)
  norWeight[onIndice]<-weight[onIndice]/sum(weight[onIndice]) #corrected: should not be abs(weight[onIndice]/sum(weight[onIndice])) and there cannot be negative values
  cumNorWeight<-cumsum(norWeight)
  minES<-min(cumNorWeight)
  maxES<-max(cumNorWeight)
  if (abs(minES)>=abs(maxES)){
    return(0)
  }else{
    return(maxES)
  }
}


#' Title calculate permutated ES
#'
#' @param weight A vector of gene weight
#' @param numOnList index of genes in the concept list that are common with weighted genes.
#' @param nPermu  Default is 1000
#'
#' @return Permuted enrichment score of weight genes
#' @export
#'
#' @examples perm_weightedKS(weight, numOnList=10, nPermu=1000)
perm_weightedKS<-function(weight,numOnList,nPermu=1000){
  permu=sapply(1:nPermu,function(x) cal_ES(weight,sample(names(weight),numOnList)))
  permu<-as.numeric(permu)
  return(permu)
}


#' Title transform the negative tail of the weight into small positive values
#'
#' @param weight signed q-values calculated from limmaDEG or scdeDEG
#'
#' @return positive or negative weight.
#' @export
#'
#' @examples weight.transform(weight) or weight.transform(-weight)
weight.transform<-function(weight){
  weight=ifelse(weight<=0,min(weight[weight>0])*((weight-min(weight[weight<=0]))/(max(weight[weight<=0])-min(weight[weight<=0]))),weight)
  return(weight)
}


#' Title Run K-S test for mcolecular concepts
#'
#' @param weight signed q-values calculated from limmaDEG or scdeDEG
#' @param compare.list pre-compiled molecular concepts
#' @param myPermu Permuted calculation of enrichment scores
#' @param p.cut To cut off the pathways by p-value threshold
#' @param min.numOnList The minimum numbeer of common genes between gene signatures and uniConSig score genes
#' @param cal.qValue From pathway enrichment
#'
#' @return Weighted K-S test output for molecular concepts.
#' @export
#'
#' @examples weightedKSV2(weight,compare.list,myPermu,p.cut=0.05,min.numOnList=5,cal.qValue=T)
weightedKSV2<-function(weight,compare.list,myPermu,p.cut=0.05,min.numOnList=5,cal.qValue=T){
  mytest = lapply(compare.list, function (x) {
    ES = cal_ES(weight, x)
    onIndice = which(names(weight) %in% x)
    numOnList = length(onIndice)
    if(numOnList < min.numOnList|ES==0) {
      return(c(0, 1, ""))
    } else {
      permu = myPermu[[as.character(numOnList)]]
      NES = ES/mean(permu)
      pValue = length(permu[which(permu >= ES)])/length(permu)
      return(c(NES,pValue,"NA"))
    }
  })
  mytest.out =t(bind_rows(mytest))
  colnames(mytest.out) = c("NES","pValue","qValue")
  mytest.out = data.frame(Compare.List = rownames(mytest.out), mytest.out,stringsAsFactors = F)
  rownames(mytest.out) = NULL
  mytest.out = mytest.out[mytest.out[, 'qValue'] == "NA", ]
  if (cal.qValue==T){
    mytest.out[, 'qValue']=qvalue(as.numeric(mytest.out[, 'pValue']), pi0 = 1)$qvalues
  }
  mytest.out=mytest.out[as.numeric(mytest.out[, 'pValue'])<p.cut,]
  mytest.out[2:4] <- mutate_all(mytest.out[2:4], function(x) as.numeric(as.character(x)))
  return(mytest.out)
}


#' Title calculated weighted KS using weightedKSV2 based on a gene list with t-values as weight (such as the gene list sorted by t-values)
#'
#' @param weight signed q-values of DEG genes calculated from limmaDEG or scdeDEG
#' @param feature.list pre-compiled molecular concepts
#' @param high.enrich default is TRUE. Then, it uses 'weight' to rank subjects
#' @param p.cut Default 0.5
#' @param correct.outlier Default is FALSE. User can choose to correct outliers
#' @param minsize Filter out concepts that have genes less than minsize
#'
#' @return Data.frame of K-S test results. This K-S test calculate the enrichment of molecular concepts in feature.list
#' @export
#'
#' @examples weightedKS.batch(weight, feature.list=feature.list,high.enrich=TRUE,p.cut=0.05,correct.outlier=FALSE,minsize=5)
weightedKS.batch<-function(weight,feature.list, high.enrich=TRUE,p.cut=0.05,correct.outlier=FALSE,minsize=5){
  feature.list<-lapply(feature.list, filterGeneSymb, WeightGeneSymbol=names(weight),inlist=TRUE)
  feature.list<-feature.list[which(sapply(feature.list,function(x) length(x)>=minsize))]
  if(high.enrich){
    print("Using Weight to rank subjects")
    tmp.weight<-sort(weight) ##use 1-AUC to favor the sensitive cell lines
  }else{
    print("Using 1-weight to rank subjects")
    tmp.weight=sort(max(weight)-weight)
  }
  if (correct.outlier){
    tmp.weight=tmp.weight-min(tmp.weight)
    mySlope=max(tmp.weight)/which.max(tmp.weight)
    max.index=which.max((mySlope*1:length(tmp.weight)-tmp.weight)/sqrt(1+mySlope^2))
    cut=tmp.weight[max.index]
    outlier.slope=((max(tmp.weight)-cut)/max(tmp.weight))/((length(tmp.weight)-max.index)/length(tmp.weight))
    if (outlier.slope>1){
      print (paste("outlier cut is ",outlier.slope,": >1, correcting outliers of weight",paste=""))
      tmp.weight=ifelse(tmp.weight>cut,nthroot(tmp.weight,2), tmp.weight) #weight transform nthroot = root of 1/n, here we used square root
    }
  }
  tmp.weight=sort(normalize(tmp.weight),decreasing=TRUE)
  matchn=sapply(feature.list,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(matchn))
  print(paste("Pre-calculating 1000 ES permutations",sep=" ")) #corrected: precalculate all possible positive numbers less than the max number
  myPermu=lapply(matchn,function(j) perm_weightedKS(tmp.weight,j));names(myPermu)=as.numeric(matchn)
  print("Pre-calculation finished")
  ks.result<-weightedKSV2(weight=tmp.weight,compare.list=feature.list,myPermu=myPermu,p.cut=p.cut,min.numOnList=minsize)
  return(ks.result)
}


#' Title run weighted KS for a weighted subject list using weightedKSV2 module
#'
#' @param weight signed q-values calculated from limmaDEG or scdeDEG
#' @param signed TRUE is default. This means that you will use signed q-value to run weighted k-s test.
#' @param feature.list pre-compiled molecular concepts
#' @param minsize Filter out concepts that have genes less than minsize
#' @param correct.outlier Default is FALSE. User can choose to correct outliers
#' @param transformNegWeight Default is FALSE. It transforms the negative tail of the weight into small positive values
#'
#' @return K-S test result.
#' @export
#'
#' @examples run.weightedKS(weight,signed=TRUE,feature.list=feature.list,minsize=10,correct.outlier=FALSE,transformNegWeight=FALSE )
run.weightedKS<-function(weight,signed=TRUE,feature.list,minsize=10,correct.outlier=FALSE,transformNegWeight=FALSE){
  feature.select<-lapply(feature.list,filterGeneSymb ,WeightGeneSymbol=names(weight),inlist=TRUE)
  feature.select<-feature.select[which(sapply(feature.select,function(x) length(x)>=minsize))]
  if (signed==TRUE){
    if (transformNegWeight==TRUE){
      posWeight=weight.transform(weight)
      negWeight=weight.transform(-weight)
    }else{
      posWeight=weight
      negWeight=-weight
    }
    ks.descending<-weightedKS.batch(weight=posWeight,feature.select, high.enrich=TRUE,correct.outlier=correct.outlier,minsize=minsize)
    ks.ascending<-weightedKS.batch(weight=negWeight,feature.select, high.enrich=TRUE,correct.outlier=correct.outlier,minsize=minsize)
  }else{
    ks.descending<-weightedKS.batch(weight,feature.select,high.enrich=TRUE,correct.outlier=correct.outlier,minsize=minsize)
    ks.ascending<-weightedKS.batch(weight,feature.select,high.enrich=FALSE,correct.outlier=correct.outlier,minsize=minsize)
  }
  ks.all=list(ks.descending=ks.descending,ks.ascending=ks.ascending)
  return(ks.all)
}









