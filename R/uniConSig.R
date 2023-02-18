#' Title Calculate weight based on Ochiai Index
#'
#' @param list1 Target gene
#' @param list2 Genes symbols in Gene Signatures (MSigDB)
#' @param method "Ochiai" or "Jaccard"
#'
#' @return Calculated weight of target genes
#' @export
#'
#' @examples
CalWeight<-function(list1,list2,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  tmp.intersect<-intersect(list1,list2)
  tmp.union<-union(list1,list2)
  if (method=="Ochiai"){
    tmp.weight=length(tmp.intersect)/sqrt(length(list1)*length(list2))
    tmp.weight.1=(length(tmp.intersect)-1)/sqrt(length(list1)*length(list2))
  }else if (method=="Jaccard"){
    tmp.weight=length(tmp.intersect)/length(tmp.union)
    tmp.weight.1=(length(tmp.intersect)-1)/length(tmp.union)
  }
  return(c(tmp.weight,tmp.weight.1,length(tmp.intersect),length(list1),length(list2)))
}





#' Title batch calculate weight for a target list and a compendia of comparing lists
#'
#' @param target.list A vector of gene symbols
#' @param compare.list Gene Signatures provided by Molecular Signature Database (MSigDB)
#' @param method "Ochiai" or "Jaccard"
#'
#' @return Weight of target genes
#' @export
#'
#' @examples batch_CalWeight(target.list=target.list,compare.list=feature.list,method="Ochiai")
batch_CalWeight<-function(target.list,compare.list,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  target.list=as.character(target.list)
  target.result<-data.frame(matrix(ncol = 6, nrow = 0))
  for (i in 1:length(compare.list)){
    target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],CalWeight(target.list,compare.list[[i]],method=method))
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
#' @param target.list A vector of gene symbols
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
#' @examples cal.uniConSig(target.list=target.list_Up,feature.list=feature.list,preCalmatrix,rm.overfit=F)
cal.uniConSig<-function(target.list,feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai",rm.overfit){ #method=="Ochiai" or "Jaccard"
  if (!exists("feature.list")){
    stop("please provide the list of concepts to variable: feature.list")
  }else if (!exists("preCalmatrix")){
    stop("please provide the preCalmatrix for concepts to variable: preCalmatrix")
  }
  target.weight=batch_CalWeight(target.list=target.list,compare.list=feature.list,method=method)
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
        if (subjectID %in% target.list){
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
      result[nrow(result) + 1,]<-c(tmp.line[1],uniConSig,ifelse(subjectID %in% target.list,1,0))
      if(i %% 5000==0){
        print(paste("Processed ",i," subject IDs",sep=""))
      }
    }
    colnames(result)<-c("subjectID","uniConSig","Target.List")
    result$uniConSig=normalize(as.numeric(as.character(result$uniConSig)))
    result=result[order(result$uniConSig, decreasing = TRUE),]
    rownames(result)=1:nrow(result)
    return(result)
  }
}





