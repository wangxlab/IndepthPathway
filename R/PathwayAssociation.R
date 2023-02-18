#' Title Calculation the functional associations of deambiguated pathway enrichment result
#'
#' @param DeambiguatePathway up or down.deambiguation data
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param feature.list pre-compiled molecular concepts
#' @param preCalmatrix Precalculated matrix to penalize the redundancies in molecular concepts
#' @param selectn The maximum number of pathwasy for similarity matrix.
#' @param minsize Default 10. Filter out concepts that have genes less than minsize
#' @param rm.overfit Default FALSE.
#'
#' @return A similarity matrix between top N pathways
#' @export
#'
#' @examples pathwayAssociation(DeambiguatePathway=up.deambiguate, compare.list, feature.list, preCalmatrix, selectn=30, minsize=10)
pathwayAssociation<-function(DeambiguatePathway, compare.list,feature.list,preCalmatrix, selectn=30, minsize=10,rm.overfit=FALSE) {
  topPathway=DeambiguatePathway[[1]]$Compare.List[1:min(c(selectn,nrow(DeambiguatePathway[[1]])))]
  topPathway=topPathway[which(sapply(topPathway,function(x) length(compare.list[[x]])>=minsize))]
  topPathway.list=compare.list[topPathway]
  uniConSig.list=lapply(topPathway.list,function(pathwayx){
    cal.uniConSig(target.list=pathwayx,feature.list=feature.list,preCalmatrix=preCalmatrix,minsize=minsize,rm.overfit=rm.overfit)
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
      tmp.CSEA=CSEA2(target.score=tmp.weight, compare.list=topPathway.list,p.cut=2,minsize=0,cal.qValue=F) #make sure to use p.cut=2 and minsize=0 to obtain CSEA results for all pathways.
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
  ij=which(tmp.bin==0, arr.ind = T)
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


