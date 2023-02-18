#' Limma to find differentially expressed genes
#'
#' @param ExpDataGCT Gene expression data in gct format. But first two rows were removed.
#' @param ClassData Class data in cls format. But first two rows were removed.
#' @param groups.order two classes
#'
#' @return data.frame of Limma output
#' @export
#'
#' @examples limmaDEG(ExpDataGCT,ClassData,groups.order=c("ClassA","ClassB"))
limmaDEG <- function(ExpDataGCT,ClassData,groups.order=NULL) { # the comparison is based on increasingly sorted orders of group labels with positive t value showing upregulated in class 1
  #expData=read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
  expUnique<-ExpDataGCT[which(!duplicated(ExpDataGCT[,1])),]
  rownames(expUnique)<-expUnique[,1]
  expUnique=expUnique[,-c(1,2)]
  expUnique=data.frame(lapply(expUnique, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(expUnique))
  expUnique<-expUnique[(which(!apply(expUnique,1,var)==0)),]  ## Remove genes (rows) that have zero variance
  #myClass=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, comment.char="",row.names=NULL, check.names=F, fill=TRUE)
  if (typeof(groups.order)=="NULL"){
    ClassData<-factor(ClassData[2,],levels=ClassData[1,2:3])
  }else{
    ClassData<-factor(myClass[2,],levels=groups.order)
  }
  label=sort(unique(as.character(ClassData)))
  fit<-lmFit(expUnique,design=model.matrix(~ClassData))
  fit.eBayes <- eBayes(fit)
  limmaOut<-topTable(fit.eBayes,number=100000000, adjust.method="BH")  # BH, Benjamini and Hochberg FDR adjusted p-value.
  colnames(limmaOut)<-c("Log2FC","AveExpr", "T.Value","P.Value","Q.Value","DGE.odds")
  limmaOut$Signed.Q.Value=sign(limmaOut$T.Value)*(-log10(limmaOut$Q.Value))
  limmaOut<-limmaOut[order(limmaOut$Signed.Q.Value,decreasing = TRUE),]
  return(limmaOut)
}


#' SCDE to find differentially expressed genes
#'
#' @param ExpDataGCT Gene expression data in gct format. But first two rows were removed.
#' @param ClassData Class data in cls format. But first two rows were removed.
#' @param groups.order two classes
#' @param min.lib.size default 100
#' @param min.reads default 1
#' @param min.detected default 1
#' @param SignPzscore default 0.75
#'
#' @return data.frame of SCDE output
#' @export
#'
#' @examples scdeDEG(ExpDataGCT, ClassData, groups.order=c("ClassA","ClassB"))
scdeDEG<-function(ExpDataGCT,ClassData, groups.order=NULL, min.lib.size=100, min.reads=1,min.detected=1,SignPzscore=0.75)  {
  #countData=read.table(readCountGCT, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
  countUnique<-ExpDataGCT[which(!duplicated(ExpDataGCT[,1])),]
  rownames(countUnique)<-countUnique[,1]
  countUnique=countUnique[,-c(1,2)]
  countUnique.df=data.frame(lapply(countUnique, function(x) as.numeric(as.character(x))),
                            check.names=F, row.names = rownames(countUnique))    # 7744 177
  cleanData=countUnique.df
  # set the sample group (sg)
  #myClass=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, comment.char="",row.names=NULL, check.names=F, fill=TRUE)
  if (typeof(groups.order)=="NULL"){
    ClassData<-factor(ClassData[2,],levels=ClassData[1,2:3])
  }else{
    ClassData<-factor(ClassData[2,],levels=groups.order)
  }
  ClassData_TP <- t(ClassData)
  myClass<-as.factor(ClassData_TP[1,])
  names(myClass)<-colnames(countUnique.df)
  myClass.selc<-myClass[names(myClass) %in% colnames(cleanData)]   # 160

  ## Fitting error models
  cleanData<-apply(cleanData,2,function(x) {storage.mode(x) <- 'integer'; x})

  o.ifm <- scde.error.models(counts=cleanData, groups=myClass.selc, n.cores=1, threshold.segmentation=TRUE,
                             save.crossfit.plots=FALSE, save.model.plots=FALSE, verbose=1)    # 160 6
  # the expected expression magnitudes (very poor fits)
  valid.cells <- o.ifm$corr.a > 0
  #table(valid.cells)
  o.ifm <- o.ifm[valid.cells, ]

  # estimate gene expression prior
  o.prior <- scde.expression.prior(models=o.ifm, counts=cleanData, length.out=400, show.plot=FALSE)  # 401 4

  # define two groups of cells: because one sample was removed from "o.ifm"
  groups <-myClass.selc[names(myClass.selc) %in% rownames(o.ifm)]
  names(groups) <- row.names(o.ifm)

  expDiff <- scde.expression.difference(o.ifm, cleanData, o.prior, groups=groups, n.randomizations=150, n.cores=10, verbose=0)  #3700 6
  colnames(expDiff)<-c("conf.int.low","foldChange.MLE","conf.int.up","conservative.foldChange","Signed.P.Zscore", "Signed.Q.Value")
  expDiff$Signed.Q.Value <- (expDiff$Signed.Q.Value)*(-1)
  expDiff.order<-expDiff[order(expDiff$Signed.P.Zscore,decreasing=T), ]

  return(expDiff.order)
}

