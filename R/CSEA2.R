#' Title perform CSEA2 for pathway enrichment analysis based on a dichotomous target gene list
#'
#' @param target.score Up or Down.uniConSig scores of whole genes.
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
#' @examples CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=1)
CSEA2<-function(target.score,compare.list,p.cut=0.05,minsize=5,min.numOnList=5,transformNegWeight=FALSE,cal.qValue=T){
  if (!exists("compare.list")){
    stop("please provide the list of concepts to variable: feature.list")
  }
  if (length(Filter(any,duplicated(names(target.score)))  ) ){
    stop("please make sure that there is no duplicated gene names in your weights")
  }
  #tmp.weight=target.score$uniConSig
  #names(tmp.weight)=target.score$subjectID
  if (transformNegWeight==TRUE){
    target.score=weight.transform(target.score)
  }else{
    target.score=normalize(target.score)
  }
  tmp.weight=target.score
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

