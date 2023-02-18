#' Read molecular concepts or gene signature from gmt files
#'
#' @param fileName molecular concept file in gmt format
#' @param min to filter out concepts that have genes less than 'min' value
#'
#' @return molecular concepts in list
#' @export
#'
#' @examples read_concepts("h.all.v7.5.1.symbols.gmt")
read_concepts <- function(fileName, min = 10) { # min: set the cutoff for minimal number of genes in the concepts or pathways
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
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
#' @examples read.gct(gctFile="scData_quiescentVsActive.gct",clsFile="scData_quiescentVsActive.cls")
read.gct<-function(gctFile,clsFile){
  expData=read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
  expUnique<-expData[which(!duplicated(expData[,1])),]
  rownames(expUnique)<-expUnique[,1]
  expUnique=expUnique[,-c(1,2)]
  expUnique=data.frame(lapply(expUnique, function(x) as.numeric(as.character(x))),
                       check.names=F, row.names = rownames(expUnique))
  expUnique<-expUnique[(which(!apply(expUnique,1,var)==0)),]  ## Remove genes (rows) that have zero variance
  myClass=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, row.names=NULL, check.names=F, fill=TRUE)
  colnames(myClass)=colnames(expUnique)
  myClass=t(myClass)
  colnames(myClass)=c("Class")
  myClass=as.data.frame(myClass,stringsAsFactors = F,check.names=F)
  return(list(expData=expUnique,class=myClass))
}
