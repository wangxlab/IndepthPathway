#' Install and load necessary packages for IndepthPathway
#'
#' @return install and load packages for IndepthPathway
#' @export
#'
#' @examples LoadPackage()
LoadPackage <- function() {
  packages=c("BiocManager", "stringr","qvalue","limma","dplyr","pheatmap","igraph","otuSummary", "data.table",  "influential")## "gplots", moments, "vegan", "tidyr", "pROC",  "corrplot", "pacman","matrixStats", "msigdbr", "RColorBrewer","fgsea",
  #if (!require("BiocManager")) install.packages("BiocManager")
  #Checking if the package belongs to CRAN or to Bioconductor and installing them accordingly.
  for(lib in packages){
    if(!lib %in% installed.packages()){
      if(lib %in% available.packages()[,1]){
        install.packages(lib,dependencies=TRUE)
      }else {(BiocManager::install(lib))
      }}
  }
  #Loading the libraries
  sapply(packages,require,character=TRUE)
}


#' Read molecular concepts or gene signature from gmt files
#'
#' @param GeneSigFile molecular concept file in gmt format
#' @param min to filter out concepts that have genes less than 'min' value
#'
#' @return molecular concepts in list
#' @export
#'
#' @examples read_concepts(GeneSigFile)
read_concepts <- function(GeneSigFile, min = 10) { # min: set the cutoff for minimal number of genes in the concepts or pathways
  Lines <- data.table::fread(GeneSigFile, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  read.list <- lapply (Lines, function(x) {
    genotype=strsplit(x, "\t")[[1]]
    gname=genotype[1]
    return(list(gname,genotype[-c(1:2)]))
  })
  genotype.list=lapply(read.list, `[[`, 2)
  names(genotype.list)= lapply(read.list, `[[`, 1)
  genotype.list=genotype.list[lengths(genotype.list)>=min]
  return(genotype.list)
}


#' Title Read gct and cls files
#'
#' @param gctFile my .gct file name
#' @param clsFile my .cls file name
#'
#' @return List of expression data and class data
#' @export
#'
#' @examples read.gct(gctFile=gctFile,clsFile=clsFile)
read.gct<-function(gctFile,clsFile){
  expData=read.table(gctFile, stringsAsFactors=FALSE, skip=2, header=TRUE, row.names=NULL, check.names=FALSE, fill=TRUE,sep="\t")
  expUnique<-expData[which(!duplicated(expData[,1])),]
  rownames(expUnique)<-expUnique[,1]
  expUnique=expUnique[,-c(1,2)]
  expUnique=data.frame(lapply(expUnique, function(x) as.numeric(as.character(x))),
                       check.names=FALSE, row.names = rownames(expUnique))
  expUnique<-expUnique[(which(!apply(expUnique,1,var)==0)),]  ## Remove genes (rows) that have zero variance
  myClass=read.table(clsFile, sep="\t", stringsAsFactors=FALSE, header=FALSE, skip=1, row.names=NULL, check.names=FALSE, fill=TRUE)
  colnames(myClass)=colnames(expUnique)
  myClass=t(myClass)
  colnames(myClass)=c("Class")
  myClass=as.data.frame(myClass,stringsAsFactors = FALSE,check.names=FALSE)
  return(list(expData=expUnique,class=myClass))
}


#' Limma to find differentially expressed genes
#'
#' @param gctFile Gene expression data in GCT file format
#' @param clsFile Class data in cls format.
#' @param groups.order two classes
#'
#' @return Limma input data in list
#' @export
#'
#' @examples limmaDEG(gctFile,clsFile,groups.order=c("ClassA","ClassB"))
limmaDEG<-function(gctFile,clsFile,groups.order=c("ClassA","ClassB")){ # the comparison is based on increasingly sorted orders of group labels with positive t value showing upregulated in class 1
  expData=read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
  expUnique<-expData[which(!duplicated(expData[,1])),]
  rownames(expUnique)<-expUnique[,1]
  expUnique=expUnique[,-c(1,2)]
  expUnique=data.frame(lapply(expUnique, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(expUnique))
  expUnique<-expUnique[(which(!apply(expUnique,1,var)==0)),]  ## Remove genes (rows) that have zero variance
  myClass=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, comment.char="",row.names=NULL, check.names=F, fill=TRUE)
  if (typeof(groups.order)=="NULL"){
    myClass<-factor(as.vector(myClass[2,]),levels=myClass[1,2:3])
  }else{
    myClass<-factor(as.vector(myClass[2,]),levels=groups.order)
  }
  label=sort(unique(as.character(myClass)))

  LimmaInput <- list(ExpData=expUnique,myClass=myClass)
  return(LimmaInput)
}


#' Limma to find differentially expressed genes
#'
#' @param LimmaFit Gene expression data in gct format. But first two rows were removed.
#'
#' @return data.frame of Limma output
#' @export
#'
#' @examples LimmaWeight(LimmaFit)
LimmaWeight <- function(LimmaFit=LimmaFit) {
      fit.eBayes <- limma::eBayes(LimmaFit)
      LimmaOut<-limma::topTable(fit.eBayes,number=100000000, adjust.method="BH")  # BH, Benjamini and Hochberg FDR adjusted p-value.
      colnames(LimmaOut)<-c("Log2FC","AveExpr", "T.Value","P.Value","Q.Value","DGE.odds")
      LimmaOut$Signed.Q.Value=sign(LimmaOut$T.Value)*(-log10(LimmaOut$Q.Value))
      LimmaOut<-LimmaOut[order(LimmaOut$Signed.Q.Value,decreasing = TRUE),]
      return(LimmaOut)
}


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
#' @examples normalizeWeight(weight)
normalizeWeight <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}


#' Title calculate ES based on gene weight (sorted in decreasing order)
#'
#' @param tmp.weight A vector of gene weight
#' @param posList index of genes in the concept list that are common with weighted genes.
#'
#' @return Calculated enrichment score
#' @export
#'
#' @examples cal_ES(tmp.weight, posList=compare.list[[1]])
cal_ES<-function(tmp.weight,posList=compare.list[[1]]){
  if (sum(tmp.weight<0)>0) {
    stop (paste("negative values found in weight, aborting calculation of ES",tmp.weight,sep=" "))#there cannot be negative values in weight
  }
  weight<-sort(tmp.weight,decreasing=TRUE)
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
#' @param tmp.weight A vector of gene weight
#' @param numOnList index of genes in the concept list that are common with weighted genes.
#' @param nPermu  Default is 1000
#'
#' @return Permuted enrichment score of weight genes
#' @export
#'
#' @examples perm_weightedKS(tmp.weight, numOnList=10, nPermu=1000)
perm_weightedKS<-function(tmp.weight,numOnList=10,nPermu=1000){
  permu=sapply(1:nPermu,function(x) cal_ES(tmp.weight,sample(names(tmp.weight),numOnList)))
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
#' @examples weight.transform(weight)
#' @examples weight.transform(-weight)
weight.transform<-function(weight){
  weight=ifelse(weight<=0,min(weight[weight>0])*((weight-min(weight[weight<=0]))/(max(weight[weight<=0])-min(weight[weight<=0]))),weight)
  return(weight)
}


#' Title Run K-S test for mcolecular concepts
#'
#' @param tmp.weight signed q-values calculated from limmaDEG or scdeDEG
#' @param compare.list pre-compiled molecular concepts
#' @param myPermu Permuted calculation of enrichment scores
#' @param p.cut To cut off the pathways by p-value threshold
#' @param min.numOnList The minimum numbeer of common genes between gene signatures and uniConSig score genes
#' @param cal.qValue From pathway enrichment
#'
#' @return Weighted K-S test output for molecular concepts.
#' @export
#'
#' @examples weightedKSV2(tmp.weight,compare.list,myPermu=myPermu,p.cut=0.05,min.numOnList=5,cal.qValue=TRUE)
weightedKSV2<-function(tmp.weight,compare.list,myPermu=myPermu,p.cut=0.05,min.numOnList=5,cal.qValue=TRUE){
  mytest = lapply(compare.list, function (x) {
    ES = cal_ES(tmp.weight, x)
    onIndice = which(names(tmp.weight) %in% x)
    numOnList = length(onIndice)
    if(numOnList < min.numOnList | ES==0) {
      return(c(0, 1, ""))
    } else {
      permu = myPermu[[as.character(numOnList)]]
      NES = ES/mean(permu)
      pValue = length(permu[which(permu >= ES)])/length(permu)
      return(c(NES,pValue,"NA"))
    }
  })
  mytest.out =t(dplyr::bind_rows(mytest))
  colnames(mytest.out) = c("NES","pValue","qValue")
  mytest.out = data.frame(Compare.List = rownames(mytest.out), mytest.out,stringsAsFactors = FALSE)
  rownames(mytest.out) = NULL
  mytest.out = mytest.out[mytest.out[, 'qValue'] == "NA", ]
  if (cal.qValue==TRUE){
    mytest.out[, 'qValue']=qvalue::qvalue(as.numeric(mytest.out[, 'pValue']), pi0 = 1)$qvalues
    mytest.out=mytest.out[as.numeric(mytest.out[, 'pValue'])<p.cut,]
    mytest.out[2:4] <- dplyr::mutate_all(mytest.out[2:4], function(x) as.numeric(as.character(x)))
  } else{
    mytest.out=mytest.out[as.numeric(mytest.out[, 'pValue'])<p.cut,]
    mytest.out[2:3] <- dplyr::mutate_all(mytest.out[2:3], function(x) as.numeric(as.character(x)))
  }

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
  tmp.weight=sort(normalizeWeight(tmp.weight),decreasing=TRUE)
  matchn=sapply(feature.list,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(matchn))
  print(paste("Pre-calculating 1000 ES permutations",sep=" ")) #corrected: precalculate all possible positive numbers less than the max number
  myPermu=lapply(matchn,function(j) perm_weightedKS(tmp.weight,j));names(myPermu)=as.numeric(matchn)
  print("Pre-calculation finished")
  ks.result<-weightedKSV2(tmp.weight,compare.list=feature.list,myPermu=myPermu,p.cut=p.cut,min.numOnList=minsize)
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


#' Title Calculate weight based on Ochiai Index
#'
#' @param targetList Target gene
#' @param feature.list Genes symbols in Gene Signatures (MSigDB)
#' @param method "Ochiai" or "Jaccard"
#'
#' @return Calculated weight of target genes
#' @export
#'
#' @examples CalWeight(targetList,feature.list,method="Ochiai")
CalWeight<-function(targetList,feature.list,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  tmp.intersect<-intersect(targetList,feature.list)
  tmp.union<-base::union(targetList,feature.list)
  if (method=="Ochiai"){
    tmp.weight=length(tmp.intersect)/sqrt(length(targetList)*length(feature.list))
    tmp.weight.1=(length(tmp.intersect)-1)/sqrt(length(targetList)*length(feature.list))
  }else if (method=="Jaccard"){
    tmp.weight=length(tmp.intersect)/length(tmp.union)
    tmp.weight.1=(length(tmp.intersect)-1)/length(tmp.union)
  }
  return(c(tmp.weight,tmp.weight.1,length(tmp.intersect),length(targetList),length(feature.list)))
}


#' Title batch calculate weight for a target list and a compendia of comparing lists
#'
#' @param targetList A vector of gene symbols
#' @param feature.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param method "Ochiai" or "Jaccard"
#'
#' @return Weight of target genes
#' @export
#'
#' @examples batch_CalWeight(targetList=targetList,feature.list=feature.list,method="Ochiai")
batch_CalWeight<-function(targetList,feature.list,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  targetList=as.character(targetList)
  target.result<-data.frame(matrix(ncol = 6, nrow = 0))
  for (i in 1:length(compare.list)){
    target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],CalWeight(targetList,feature.list[[i]],method=method))
  }
  colnames(target.result)<-c("Feature.List","Weight","Weight.1","Intersect","Target.Size","Compare.Size")
  return(target.result)
}


#' Title Calculate uniConSig scores based on ks test results
#'
#' @param up.ks KS test results by weight genes
#' @param down.ks KS test reults by -weight genes
#' @param preCalmatrix Precalculated matrix to penalize the redundancies in molecular concepts
#' @param feature.list pre-compiled molecular concepts
#' @param p.cut Default is 0.01
#' @param q.cut Default is 0.25
#' @param NES.cut Default is 0
#' @param power Default is 1
#' @param root Default is 1
#' @param ECNpenalty Default is 0.5
#'
#' @return uniConSig output. Data.frame of up.uniConSig and down.uniConSig scores for whole genes
#' @export
#'
#' @examples cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list)
cal.uniConSig.ks<-function(up.ks,down.ks,preCalmatrix,feature.list,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5){
  up.ks[is.na(up.ks)] <- 0
  down.ks[is.na(down.ks)] <- 0
  up.ks.sub<-up.ks[up.ks$NES>NES.cut & up.ks$pValue<p.cut & up.ks$qValue<q.cut,-c(3,4)]
  up.ks.mat <- as.matrix(up.ks.sub[-1])
  row.names(up.ks.mat) <- up.ks.sub$Compare.List
  down.ks.sub<-down.ks[down.ks$NES>NES.cut & down.ks$pValue<p.cut & down.ks$qValue<q.cut,-c(3,4)]
  down.ks.mat <- as.matrix(down.ks.sub[-1])
  row.names(down.ks.mat) <- down.ks.sub$Compare.List
  up.list=list()
  down.list=list()
  result<-data.frame(matrix(ncol = 3, nrow = 0))
  for (i in 1:nrow(preCalmatrix)){
    tmp.line<-unlist(strsplit(preCalmatrix[i,],"\t"))
    subjectID=as.character(tmp.line[1])
    if(length(tmp.line)==2){
      next
    }
    tmp.epsilon<-as.data.frame(do.call(rbind, strsplit(tmp.line[3:length(tmp.line)],"@")))
    colnames(tmp.epsilon)<-c("Compare.List","Epsilon")
    ECN=sum(1/(as.numeric(tmp.epsilon[,"Epsilon"])^root))
    up.epsilon<-tmp.epsilon[tmp.epsilon$Compare.List %in% row.names(up.ks.mat),]
    down.epsilon<-tmp.epsilon[tmp.epsilon$Compare.List %in% row.names(down.ks.mat),]
    up.weight<- as.data.frame(up.ks.mat[which(row.names(up.ks.mat) %in% up.epsilon$Compare.List),"NES"])
    up.weight$Compare.List <- as.character(row.names(up.weight))
    colnames(up.weight)=c("NES","Compare.List")
    up.data<-as.matrix(merge(up.epsilon,up.weight,by.x="Compare.List",by.y="Compare.List"))
    up.data=cbind(up.data,normWeight=(as.numeric(up.data[,"NES"])^power)/(as.numeric(up.data[,"Epsilon"])^root))
    score.sensitive<-sum(as.numeric(up.data[,"normWeight"]))/(ECN^ECNpenalty)
    up.data=as.data.frame(up.data)
    colnames(up.data)[-1]=c(paste(subjectID,colnames(up.data[-1]),sep=":"))
    up.list[[paste(subjectID,sep="")]]<-up.data
    down.weight<- as.data.frame(down.ks.mat[which(row.names(down.ks.mat) %in% down.epsilon$Compare.List),"NES"])
    down.weight$Compare.List <- as.character(row.names(down.weight))
    colnames(down.weight)=c("NES","Compare.List")
    res.data<-as.matrix(merge(down.epsilon,down.weight,by.x="Compare.List",by.y="Compare.List"))
    res.data=cbind(res.data,normWeight=(as.numeric(res.data[,"NES"])^power)/(as.numeric(res.data[,"Epsilon"])^root))
    score.resistant<-sum(as.numeric(res.data[,"normWeight"]))/(ECN^ECNpenalty)
    res.data=as.data.frame(res.data)
    colnames(res.data)[-1]=c(paste(subjectID,colnames(res.data[-1]),sep=":"))
    down.list[[paste(subjectID,sep="")]]<-res.data
    result[nrow(result) + 1,]<-c(tmp.line[1],score.sensitive,score.resistant)
    if(i %% 2000==0){
      print(paste("Processed ",i," subjects",sep=""))
    }
  }
  colnames(result)<-c("subjectID","up.uniConSig","down.uniConSig")
  return(result)
}


#' Title Calculate uniConSig scores for genes.
#'
#' @param targetList A vector of gene symbols
#' @param feature.list pre-compiled molecular concepts
#' @param preCalmatrix Precalculated matrix to penalize the redundancies in molecular concepts
#' @param minsize Default 10
#' @param weight.cut Default 0.05
#' @param power Default 1
#' @param root Default 1
#' @param ECNpenalty Default 0.5
#' @param method Default 'Ochiai
#' @param rm.overfit Remove overfitted weight genes
#'
#' @return uniConSig output by non-weighted method
#' @export
#'
#' @examples cal.uniConSig(targetList=targetList_Up,feature.list=feature.list,preCalmatrix,rm.overfit=FALSE)
cal.uniConSig<-function(targetList,feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai",rm.overfit=FALSE){ #method=="Ochiai" or "Jaccard"
  if (!exists("feature.list")){
    stop("please provide the list of concepts to variable: feature.list")
  }else if (!exists("preCalmatrix")){
    stop("please provide the preCalmatrix for concepts to variable: preCalmatrix")
  }
  target.weight=batch_CalWeight(targetList=targetList,feature.list=feature.list,method=method)
  target.weight=target.weight[as.numeric(target.weight$Compare.Size)>minsize&as.numeric(target.weight$Weight)>weight.cut,]
  if (nrow(target.weight)<10){
    print(paste("There are only",nrow(target.weight),"signature features with weights, this means that the target list is functionally heterogeneous, please lower the cutoff for weight.cut"))
    return(NULL)
  }else{
    result<-data.frame(matrix(ncol = 3, nrow = 0))
    for (i in 1:nrow(preCalmatrix)){
      tmp.line<-unlist(strsplit(preCalmatrix[i,],"\t"))
      subjectID=as.character(tmp.line[1])
      if(length(tmp.line)==2){
        result[nrow(result) + 1,]<-c(tmp.line[1],0)
      }
      tmp.strsplit=strsplit(tmp.line[3:length(tmp.line)],"@")
      tmp.strsplit=tmp.strsplit[which(lengths(tmp.strsplit)==2)]
      tmp.epsilon<-as.data.frame(do.call(rbind, tmp.strsplit),stringsAsFactors=FALSE)
      colnames(tmp.epsilon)<-c("Feature.List","Epsilon")
      tmp.epsilon[,"Epsilon"]=(as.numeric(tmp.epsilon[,"Epsilon"]))^root

      #number of concepts/ mean of epsilon
      ECN=nrow(tmp.epsilon)/mean(as.numeric(tmp.epsilon[,"Epsilon"]),trim=0.3)
      tmp.data<-merge(tmp.epsilon,target.weight,by.x="Feature.List",by.y="Feature.List")
      tmp.data$Epsilon=as.character(tmp.data$Epsilon)
      tmp.data$Feature.List=as.character(tmp.data$Feature.List)
      if (nrow(tmp.data)==0){
        next
      }
      if (rm.overfit==TRUE){
        if (subjectID %in% targetList){
          tmp.data$Weight4Cal=tmp.data$Weight.1
        }else{
          tmp.data$Weight4Cal=tmp.data$Weight
        }
      }else{
        tmp.data$Weight4Cal=tmp.data$Weight
      }
      tmp.data$Epsilon=as.numeric(as.character(tmp.data$Epsilon))
      tmp.data$Weight4Cal=as.numeric(as.character(tmp.data$Weight4Cal))
      #ECN=sum(1/(tmp.data$Epsilon^root))
      uniConSig<-sum((tmp.data$Weight4Cal^power)/(tmp.data$Epsilon))/(ECN^ECNpenalty)
      result[nrow(result) + 1,]<-c(tmp.line[1],uniConSig,ifelse(subjectID %in% targetList,1,0))
      if(i %% 5000==0){
        print(paste("Processed ",i," subject IDs",sep=""))
      }
    }
    colnames(result)<-c("subjectID","uniConSig","Target.List")
    result$uniConSig=normalizeWeight(as.numeric(as.character(result$uniConSig)))
    result=result[order(result$uniConSig, decreasing = TRUE),]
    rownames(result)=1:nrow(result)
    return(result)
  }
}

#' Title perform CSEA2 for pathway enrichment analysis based on a dichotomous target gene list
#'
#' @param targetScore Up or Down.uniConSig scores of whole genes.
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param p.cut To cut off the pathways by p-value threshold
#' @param minsize To cut off gene signatures (MSigDB) by the minimun number of genes
#' @param min.numOnList The minimum numbeer of common genes between gene signatures and uniConSig score genes
#' @param transformNegWeight To transform uniConSig score or normalize it.
#' @param cal.qValue From pathway enrichment
#'
#' @return CSEA or WCSEA output. Data.frame of pathway enrichment.
#' @export
#'
#' @examples CSEA2(targetScore=setNames(as.numeric(uniConSigResult$up.uniConSig), uniConSigResult$subjectID),compare.list,p.cut=1)
CSEA2<-function(targetScore=setNames(as.numeric(uniConSigResult$up.uniConSig), uniConSigResult$subjectID), compare.list,p.cut=0.05,minsize=5,min.numOnList=5,transformNegWeight=FALSE,cal.qValue=TRUE){
  if (!exists("compare.list")){
    stop("please provide the list of concepts to variable: feature.list")
  }
  if (any (duplicated(names(targetScore))) ){
    stop("please make sure that there is no duplicated gene names in your weights")
  }
  #tmp.weight=targetScore$uniConSig
  #names(tmp.weight)=targetScore$subjectID
  if (transformNegWeight==TRUE){
    targetScore=weight.transform(targetScore)
  }else{
    targetScore=normalizeWeight(targetScore)
  }
  tmp.weight=targetScore
  tmp.weight=tmp.weight[which(tmp.weight!=0)]
  if (length(tmp.weight)<50){
    stop(paste("There are only",length(tmp.weight),"genes have uniConSig scores, cannot perform CSEA"))
  }else if(length(tmp.weight)<1000){
    print (paste("There are only",length(tmp.weight),"genes have uniConSig scores, CSEA results may not be reliable"))
  }else {
    print (paste("There are",length(tmp.weight),"genes have uniConSig scores"))
  }
  compare.list=lapply(compare.list,filterGeneSymb, WeightGeneSymbol=names(tmp.weight),inlist=TRUE)
  compare.list<-compare.list[which(sapply(compare.list,function(x) length(x)>=minsize))]
  matchn=sapply(compare.list,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(matchn))
  matchn=matchn[matchn>=min.numOnList]
  print(paste("Pre-calculating 1000 ES permutations",sep=" ")) #corrected: precalculate all possible positive numbers less than the max number
  myPermu=lapply(matchn,function(j) perm_weightedKS(tmp.weight,j));names(myPermu)=as.numeric(matchn)
  print("Pre-calculation finished")
  ks.result<-weightedKSV2(tmp.weight,compare.list,myPermu=myPermu,p.cut=p.cut,min.numOnList=min.numOnList,cal.qValue=cal.qValue)
  ks.result=ks.result[order(ks.result$NES,decreasing = TRUE),]
  if (nrow(ks.result)==0){
    print("no pathway is found to be significant by KS test")
    return(NULL)
  }else{
    row.names(ks.result)=1:nrow(ks.result)
    return(ks.result)
  }
}


#' Title disambiguate pathway result by CSEA or WCSEA
#'
#' @param GEA.result Pathway enrichment results
#' @param uniConSigResult uniConSig result
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param upPathways TRUE or FALSE. TRUE will use up.uniConSig scores
#' @param topn The number of pathways to disambiguate.
#' @param p.cut Default is 0.05
#' @param p.adjust.method A method to adjust p-values; "BH","bonferroni" or NULL
#'
#' @return disambiguated pathway enrichment result
#' @export
#'
#' @examples disambiguationCSEA(GEA.result=down.CSEA.result,uniConSigResult=uniConSigResult,compare.list=compare.list,upPathways=FALSE,topn=30,p.cut=0.01,p.adjust.method="bonferroni")
disambiguationCSEA<-function(GEA.result,uniConSigResult,compare.list,upPathways=NULL,topn=30,p.cut=0.05,p.adjust.method="BH"){ #upPathways TRUE or FALSE (for W-CSEA) or NULL (NULL is for D-CSEA). p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
  topPathway=GEA.result$Compare.List[1:topn]
  topPathway.list=compare.list[topPathway]
  topPathway.del=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=topPathway.list[[i]]
    for (j in 1:length(topPathway.list)){
      pathwayy=topPathway.list[[j]]
      if (i==j){
        next
      }
      list=list(a=pathwayx[which(!pathwayx %in% pathwayy)])
      names(list)=c(paste(names(topPathway.list)[i],":",names(topPathway.list)[j],sep=""))
      topPathway.del=c(topPathway.del,list)
    }
  }
  if (is.null(upPathways)){
    CSEA.del<-CSEA2(targetScore=setNames(as.numeric(uniConSigResult$uniConSig), uniConSigResult$subjectID),compare.list=topPathway.del,p.cut=2)
  }else if (upPathways==TRUE){
    CSEA.del<-CSEA2(targetScore=setNames(as.numeric(uniConSigResult$up.uniConSig), uniConSigResult$subjectID),compare.list=topPathway.del,p.cut=2)
  }else{
    CSEA.del<-CSEA2(targetScore=setNames(as.numeric(uniConSigResult$down.uniConSig), uniConSigResult$subjectID),compare.list=topPathway.del,p.cut=2)
  }
  if(typeof(p.adjust.method)!="NULL"){#http://rcompanion.org/rcompanion/f_01.html
    CSEA.del$p.adjust=p.adjust(CSEA.del$pValue,method = p.adjust.method)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = TRUE)
  for (k in 1:nrow(ij)){
    delname=paste(rownames(tmp.bin)[ij[k,"row"]],colnames(tmp.bin)[ij[k,"col"]],sep=":")
    if (delname %in% CSEA.del$Compare.List){
      if(typeof(p.adjust.method)!="NULL"){
        tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(CSEA.del$p.adjust[CSEA.del$Compare.List==delname])
      }else{
        tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(CSEA.del$pValue[CSEA.del$Compare.List==delname])
      }
    }else{
      tmp.bin[ij[k,"row"],ij[k,"col"]]=1
    }
  }
  index=which(tmp.bin>p.cut,arr.ind = TRUE)
  remove=c()
  for (x in 1:nrow(index)){
    rowname=rownames(tmp.bin)[index[x,"row"]]
    colname=colnames(tmp.bin)[index[x,"col"]]
    if (rowname==colname | rowname %in%remove |colname %in% remove ){
      next
    }
    remove=c(remove,rowname)
  }
  filter.result=GEA.result[1:topn,]
  filter.result=filter.result[which(!filter.result$Compare.List %in% remove),]
  return(list(pathway.filtered=filter.result,pathway.todel=remove,disambiguate.matrix=tmp.bin))
}



#' Title disambiguate pathway result by GSEA
#'
#' @param GEA.result Pathway enrichment results
#' @param weight target gene weight
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param topn The number of pathways to disambiguate.
#' @param p.cut Default is 0.05
#' @param transformNegWeight To transform uniConSig score or normalize it.
#' @param p.adjust.method A method to adjust p-values; "BH","bonferroni" or NULL
#'
#' @return disambiguated pathway enrichment result by GSEA
#' @export
#'
#' @examples disambiguationGSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=30,p.cut=0.8,p.adjust.method="bonferroni")
disambiguationGSEA<-function(GEA.result,weight,compare.list,topn=30,p.cut=0.05,transformNegWeight=FALSE,p.adjust.method="BH"){ #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
  topPathway=GEA.result$Compare.List[1:topn]
  topPathway.list=compare.list[topPathway]
  topPathway.del=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=topPathway.list[[i]]
    for (j in 1:length(topPathway.list)){
      pathwayy=topPathway.list[[j]]
      if (i==j){
        next
      }
      list=list(a=pathwayx[which(!pathwayx %in% pathwayy)])
      names(list)=c(paste(names(topPathway.list)[i],":",names(topPathway.list)[j],sep=""))
      topPathway.del=c(topPathway.del,list)
    }
  }
  GSEA.del<-CSEA2(targetScore=weight,compare.list=topPathway.del,p.cut=2,minsize=5,transformNegWeight=transformNegWeight)
  if(typeof(p.adjust.method)!="NULL"){#http://rcompanion.org/rcompanion/f_01.html
    GSEA.del$p.adjust=p.adjust(GSEA.del$pValue,method = p.adjust.method)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = TRUE)
  for (k in 1:nrow(ij)){
    delname=paste(rownames(tmp.bin)[ij[k,"row"]],colnames(tmp.bin)[ij[k,"col"]],sep=":")
    if (delname %in% GSEA.del$Compare.List){
      if(typeof(p.adjust.method)!="NULL"){
        tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(GSEA.del$p.adjust[GSEA.del$Compare.List==delname])
      }else{
        tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(GSEA.del$pValue[GSEA.del$Compare.List==delname])
      }
    }else{
      tmp.bin[ij[k,"row"],ij[k,"col"]]=1
    }
  }
  index=which(tmp.bin>p.cut,arr.ind = TRUE)
  remove=c()
  for (x in 1:nrow(index)){
    rowname=rownames(tmp.bin)[index[x,"row"]]
    colname=colnames(tmp.bin)[index[x,"col"]]
    if (rowname==colname | rowname %in%remove |colname %in% remove ){
      next
    }
    remove=c(remove,rowname)
  }
  filter.result=GEA.result[1:topn,]
  filter.result=filter.result[which(!filter.result$Compare.List %in% remove),]
  #rownames(filter.result)=1:nrow(filter.result)
  return(list(pathway.filtered=filter.result,pathway.todel=remove,disambiguate.matrix=tmp.bin))
}



#' Title Calculation the functional associations of disambiguated pathway enrichment result
#'
#' @param topPathway up or down.disambiguation data
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param feature.list pre-compiled molecular concepts
#' @param preCalmatrix Precalculated matrix to penalize the redundancies in molecular concepts
#' @param minsize Default 10. Filter out concepts that have genes less than minsize
#' @param rm.overfit Default FALSE.
#'
#' @return A similarity matrix between top N pathways
#' @export
#'
#' @examples pathwayAssociation(topPathway=topPathwayUp, compare.list, feature.list, preCalmatrix,  minsize=10, rm.overfit=FALSE)
pathwayAssociation<-function(topPathway=topPathwayUp,compare.list,feature.list,preCalmatrix,minsize=10,rm.overfit=FALSE){ #upPathways TRUE or FALSE
  #topPathwayLarge = topPathway[which(sapply(topPathway,function(x) length(unlist(compare.list[x]) )>=minsize))]
  topPathway.list=compare.list[topPathway]
  uniConSig.list=lapply(topPathway.list,function(targetList){
        cal.uniConSig(targetList,feature.list=feature.list,preCalmatrix=preCalmatrix,minsize=minsize,rm.overfit=rm.overfit)
  })
  topPathway.CSEA=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=names(topPathway.list)[i]
    tmp.uniConSig=uniConSig.list[[pathwayx]]
    tmp.weight=setNames(as.numeric(tmp.uniConSig$uniConSig),tmp.uniConSig$subjectID)
    if (length(tmp.weight[tmp.weight!=0])<50){
      print(paste(pathwayx,": There are only",length(tmp.weight),"genes have uniConSig scores, cannot perform CSEA. Assuming there is no functional relations with other pathways."))
      tmp.CSEA=data.frame(Compare.List=names(topPathway.list),NES=ifelse(names(topPathway.list)==pathwayx,5,0),pValue=ifelse(names(topPathway.list)==pathwayx,0,1),qValue=ifelse(names(topPathway.list)==pathwayx,0,1))
    }else{
      tmp.CSEA=CSEA2(targetScore=tmp.weight, compare.list=topPathway.list,p.cut=2,minsize=0, cal.qValue=FALSE) #make sure to use p.cut=2 and minsize=0 to obtain CSEA results for all pathways.
      topPathway.missed=names(topPathway.list)[!names(topPathway.list) %in% tmp.CSEA$Compare.List]
      if (length(topPathway.missed)>0){
        for (pathway in topPathway.missed){
          tmp.CSEA=rbind(tmp.CSEA,c(pathway,0,1,1))
        }
      }
    }
    list=list(a=tmp.CSEA)
    names(list)=c(paste(names(topPathway.list)[i],sep=""))
    topPathway.CSEA=c(topPathway.CSEA,list)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = TRUE)
  for (k in 1:nrow(ij)){
    pathwayi=rownames(tmp.bin)[ij[k,"row"]]
    pathwayj=colnames(tmp.bin)[ij[k,"col"]]
    tmp.CSEA=topPathway.CSEA[[pathwayi]]
    if (pathwayi != pathwayj){
      tmp.bin[ij[k,"row"],ij[k,"col"]]=tmp.CSEA$NES[tmp.CSEA$Compare.List==pathwayj]
    }
    #print (k)
  }
  mode(tmp.bin)<-"numeric"
  diag(tmp.bin) <- max(tmp.bin)
  return(tmp.bin)
}


#' Title Option1.Draw heatmaps showing the functional associations between selected top pathways
#'
#' @param matrixData PathwayAssociation data
#' @param clustering Default is TRUE
#' @param distanceMatric Default is 'euclidean'
#' @param fontSize Default is 0.7
#' @param sepWidth Default is 0.1
#'
#' @return heatmaps showing the functional associations between selected top pathways
#' @export
#'
#' @examples pathway.heatmap(matrixData=up.assoc,clustering = TRUE,fontSize=5)
pathway.heatmap<-function(matrixData,clustering=TRUE,distanceMatric="euclidean",fontSize=0.7,sepWidth=0.1){
  rownames(matrixData)=gsub("_"," ",rownames(matrixData))
  colnames(matrixData)=gsub("_"," ",colnames(matrixData))
  my_palette <- colorRampPalette(c("dodger blue", "white", "red"))(n = 1000)
  pheatmap::pheatmap(matrixData, cluster_rows=clustering, cluster_cols=clustering,color=my_palette, show_rownames = TRUE, show_colnames = TRUE,fontsize=fontSize,legend=TRUE)
  #return(pathway.heatmap)
}


#' Title Merge upregulated and downregulated pathways. For pathways that are both up or down regulated, the most significant rank will be retained
#'
#' @param uppathway filtered up_disambiguate
#' @param downpathway filtered down_disambiguate
#'
#' @return Merged pathways. The pathways of the most significant rank will be retained
#' @export
#'
#' @examples mergePathway(uppathway,downpathway)
mergePathway<-function(uppathway,downpathway){
  pathway.ambiguous=intersect(uppathway$Compare.List,downpathway$Compare.List)
  uppathway$Type="UP"
  downpathway$Type="DOWN"
  if (length(pathway.ambiguous)==0){
    pathway.combine=rbind(uppathway,downpathway)
  }else{
    ambiguous.uprank=setNames(rownames(uppathway)[match(pathway.ambiguous,uppathway$Compare.List)],pathway.ambiguous)
    ambiguous.downrank=setNames(rownames(downpathway)[match(pathway.ambiguous,downpathway$Compare.List)],pathway.ambiguous)
    pathway.ambiguous.type=setNames(ifelse(as.numeric(ambiguous.uprank[pathway.ambiguous])<=as.numeric(ambiguous.downrank[pathway.ambiguous]),"UP","DOWN"),pathway.ambiguous)
    uppathway=uppathway[!uppathway$Compare.List %in% names(pathway.ambiguous.type)[pathway.ambiguous.type=="DOWN"],]
    downpathway=downpathway[!downpathway$Compare.List %in% names(pathway.ambiguous.type)[pathway.ambiguous.type=="UP"],]
    pathway.combine=rbind(uppathway,downpathway)
  }
  #print(pathway.combine[1:3,1:3])
  return(pathway.combine)
}


#' Title Option2.Draw network figure for significant pathways based on association file
#'
#' @param pathway.out Merged pathway data
#' @param assoc pathwayMergeAssoc
#' @param NES.cut Default is 2
#' @param node.size Default is 2
#' @param line.thickness Default is 1
#' @param font.size Default is 0.3
#' @param wrap.string.width Default is 10
#' @param FirstCharacterUpperCase Default is FALSE
#'
#' @return network of the functional associations between selected top pathways
#' @export
#'
#' @examples draw.network(pathway.out=pathway.merge,assoc=pathwayMergeAssoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
draw.network<-function(pathway.out,assoc,NES.cut=2,node.size=2,line.thickness=1,font.size=0.3,wrap.string.width=10,FirstCharacterUpperCase=FALSE){
  links=otuSummary::matrixConvert(assoc, colname = c("from", "to", "similarity"))
  links=links[links$similarity>NES.cut,]
  if ("Type" %in% colnames(pathway.out)){
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type=pathway.out$Type[match(row.names(assoc),pathway.out$Compare.List)],stringsAsFactors = F)
  }else{
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type="NA",stringsAsFactors = F)
  }
  #https://kateto.net/netscix2016.html
  net <- influential::graph_from_data_frame(d=links, vertices=nodes, directed=F)
  colrs <- c("deepskyblue","tomato1")
  igraph::V(net)$color <- colrs[as.factor(igraph::V(net)$type)]#set node color
  igraph::V(net)$size <- igraph::V(net)$size #set node size
  igraph::V(net)$label <- gsub("HALLMARK_|REACTOME_|KEGG_|WP_|BIOCARTA_|PID_","",names(igraph::V(net)))
  igraph::V(net)$label <- gsub("_"," ",igraph::V(net)$label)
  igraph::V(net)$label <- stringr::str_wrap(igraph::V(net)$label,wrap.string.width)
  if(FirstCharacterUpperCase==T) igraph::V(net)$label <- stringr::str_to_title(igraph::V(net)$label)
  igraph::E(net)$width <- igraph::E(net)$similarity*line.thickness
  igraph::E(net)$arrow.size <-0
  igraph::E(net)$color <- 'cornsilk2'
  #dev.new(width=12, height=5, unit="in")
  #windowsFonts("Arial" = windowsFont("Arial"))
  network=plot(net,vertex.label.color="black",vertex.label.cex = font.size,vertex.label.family="Times", vertex.label.font=2,
               vertex.label.degree = 20,vertex.frame.color="white",asp=0.5,layout=igraph::layout.fruchterman.reingold(net, niter=10000))
  return(network)
}




