#' Title Deambiguate pathway result by CSEA or WCSEA
#'
#' @param GEA.result Pathway enrichment results
#' @param uniConSig.result uniConSig result
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param upPathways TRUE or FALSE. TRUE will use up.uniConSig scores
#' @param topn The number of pathways to deambiguate.
#' @param p.cut Default is 0.05
#' @param p.adjust.method A method to adjust p-values; "BH","bonferroni" or NULL
#'
#' @return Deambiguated pathway enrichment result
#' @export
#'
#' @examples disambiguation.CSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=30,p.cut=0.01,p.adjust.method="bonferroni")
deambiguation.CSEA<-function(GEA.result,uniConSig.result,compare.list,upPathways=NULL,topn=30,p.cut=0.05,p.adjust.method="BH"){ #upPathways TRUE or FALSE (for W-CSEA) or NULL (NULL is for D-CSEA). p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
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
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=2)
  }else if (upPathways==TRUE){
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=2)
  }else{
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=2)
  }
  if(typeof(p.adjust.method)!="NULL"){#http://rcompanion.org/rcompanion/f_01.html
    CSEA.del$p.adjust=p.adjust(CSEA.del$pValue,method = p.adjust.method)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = T)
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



#' Title Deambiguate pathway result by GSEA
#'
#' @param GEA.result Pathway enrichment results
#' @param weight target gene weight
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param topn The number of pathways to deambiguate.
#' @param p.cut Default is 0.05
#' @param transformNegWeight To transform uniConSig score or normalize it.
#' @param p.adjust.method A method to adjust p-values; "BH","bonferroni" or NULL
#'
#' @return Deambiguated pathway enrichment result by GSEA
#' @export
#'
#' @examples deambiguation.GSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=min(c(topn,nrow(down.GSEA.result))),p.cut=0.8,p.adjust.method="bonferroni")
deambiguation.GSEA<-function(GEA.result,weight,compare.list,topn=30,p.cut=0.05,transformNegWeight=FALSE,p.adjust.method="BH"){ #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
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
  GSEA.del<-CSEA2(target.score=weight,compare.list=topPathway.del,p.cut=2,minsize=5,transformNegWeight=transformNegWeight)
  if(typeof(p.adjust.method)!="NULL"){#http://rcompanion.org/rcompanion/f_01.html
    GSEA.del$p.adjust=p.adjust(GSEA.del$pValue,method = p.adjust.method)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = T)
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



