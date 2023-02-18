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
#' @examples pathway.heatmap(matrixData=up.assoc,clustering = T,fontSize=5)
pathway.heatmap<-function(matrixData,clustering=TRUE,distanceMatric="euclidean",fontSize=0.7,sepWidth=0.1){
  rownames(matrixData)=gsub("_"," ",rownames(matrixData))
  colnames(matrixData)=gsub("_"," ",colnames(matrixData))
  my_palette <- colorRampPalette(c("dodger blue", "white", "red"))(n = 1000)
  pheatmap(matrixData, cluster_rows=clustering, cluster_cols=clustering,color=my_palette, show_rownames = T, show_colnames = T,fontsize=fontSize,legend=T)
  #return(pathway.heatmap)
}


#' Title Option2.Draw network figure for significant pathways based on association file
#'
#' @param pathway.out Merged pathway data
#' @param assoc pathway.merge.assoc
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
#' @examples draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
draw.network<-function(pathway.out,assoc,NES.cut=2,node.size=2,line.thickness=1,font.size=0.3,wrap.string.width=10,FirstCharacterUpperCase=F){
  links=matrixConvert(assoc, colname = c("from", "to", "similarity"))
  links=links[links$similarity>NES.cut,]
  if ("Type" %in% colnames(pathway.out)){
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type=pathway.out$Type[match(row.names(assoc),pathway.out$Compare.List)],stringsAsFactors = F)
  }else{
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type="NA",stringsAsFactors = F)
  }
  #https://kateto.net/netscix2016.html
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)
  colrs <- c("deepskyblue","tomato1")
  V(net)$color <- colrs[as.factor(V(net)$type)]#set node color
  V(net)$size <- V(net)$size #set node size
  V(net)$label <- gsub("HALLMARK_|REACTOME_|KEGG_|WP_|BIOCARTA_|PID_","",names(V(net)))
  V(net)$label <- gsub("_"," ",V(net)$label)
  V(net)$label <- stringr::str_wrap(V(net)$label,wrap.string.width)
  if(FirstCharacterUpperCase==T) V(net)$label <- stringr::str_to_title(V(net)$label)
  E(net)$width <- E(net)$similarity*line.thickness
  E(net)$arrow.size <-0
  E(net)$color <- 'cornsilk2'
    #dev.new(width=12, height=5, unit="in")
  #windowsFonts("Arial" = windowsFont("Arial"))
  network=plot(net,vertex.label.color="black",vertex.label.cex = font.size,vertex.label.family="Times", vertex.label.font=2,vertex.label.degree = 20,vertex.frame.color="white",asp=0.5,layout=layout.fruchterman.reingold(net, niter=10000))
  return(network)
}


#' Title Option4.Draw heatmap for selected pathway
#'
#' @param genomatrix Gene expression data
#' @param geneset A specific gene signature
#' @param weights weight gene
#' @param fontsize_row Default is 6
#' @param fontsize_col Default is 6
#' @param cluster_rows Default is FALSE
#' @param annCol Default is NULL
#'
#' @return  heatmap of gnes for selected gene signature
#' @export
#'
#' @examples plot.geneset.matrix(genomatrix=expData,geneset=compare.list[["HALLMARK_E2F_TARGETS"]],weights=weight,annCol = sampleClass)
plot.geneset.matrix<-function(genomatrix,geneset,weights,fontsize_row=6,fontsize_col=6,cluster_rows=F,annCol=NULL){
  genomatrix=genomatrix[rownames(genomatrix) %in% geneset,]
  genomatrix=genomatrix[rowSums(genomatrix)!=0,]
  genomatrix=t(scale(t(genomatrix),center=TRUE,scale=TRUE))
  quantile.range <- quantile(genomatrix, probs = seq(0, 1, 0.01),na.rm = T)
  myBreaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  myColor  <- colorRampPalette(c("skyblue", "white", "red"))(length(myBreaks) - 1)
  weights.geneset=weights[names(weights) %in% rownames(genomatrix)]
  weights.geneset=weights.geneset[names(weights.geneset) %in% row.names(genomatrix)]
  genomatrix=genomatrix[names(weights.geneset)[order(weights.geneset,decreasing=T)],]
  if (!is.vector(weights, mode = "numeric")){
    stop("weights must be named numeric")
  }
  annRow=data.frame(row.names = names(weights),statistical.value=weights,stringsAsFactors = F,check.names = F)
  ann_colors = list(
    statistical.value=colorRampPalette(c("deepskyblue","white","red"))(n=500)
  )
  if (typeof(annCol)=="NULL"){
    pheatmap(genomatrix, cluster_rows=cluster_rows, cluster_cols=F,color=myColor, breaks=myBreaks,show_rownames = T, show_colnames = T,annotation_row=annRow,annotation_colors = ann_colors,fontsize_row=fontsize_row,fontsize_col=fontsize_col,legend=T)
  }else{
    pheatmap(genomatrix, cluster_rows=cluster_rows, cluster_cols=F,color=myColor, breaks=myBreaks,show_rownames = T, show_colnames = T,annotation_row=annRow,annotation_col=annCol,annotation_colors = ann_colors,fontsize_row=fontsize_row,fontsize_col=fontsize_col,legend=T)
  }
}


#' Title Draw interactome network for a selected pathway based on DEG weight
#'
#' @param links Interactome file, for example "interactome_all.hs.sym"
#' @param genelist Gene symbols in a specific gene signature
#' @param geneweight Weight gene
#' @param target Default is NULL
#' @param pdf.file Output file name.
#' @param self.interaction Default is FALSE
#' @param filter.coreEnrich Default is TRUE
#'
#' @return Interactome network for a selected pathway based on DEG weight
#' @export
#'
#' @examples draw.interactome.network(links=HuRI.interactions,genelist=genelist,geneweight=geneweight,target=NULL,pdf.file="My.HuRI.network.pdf",self.interaction=F,filter.coreEnrich=T)
draw.interactome.network<-function(links,genelist,geneweight,target=NULL,pdf.file,self.interaction=F,filter.coreEnrich=T){
  library(igraph);library(RColorBrewer);library(fgsea)
  geneweight=geneweight[!is.na(geneweight)]
  genelist=genelist[genelist %in% names(geneweight)]
  if (filter.coreEnrich==T){
    enrich.plot=plotEnrichment(genelist, normalize(geneweight))
    weight.order=setNames(c("zero",geneweight),order(c(0,geneweight),decreasing=T))
    if(abs(max(enrich.plot$data$y))>abs(min(enrich.plot$data$y))){
      cut.off.max=as.numeric(weight.order[enrich.plot$data$x[which.max(enrich.plot$data$y)]])
      genelist=genelist[genelist %in% names(geneweight)[geneweight>cut.off.max]]
    }else{
      cut.off.min=as.numeric(weight.order[enrich.plot$data$x[which.min(enrich.plot$data$y)]])
      genelist=genelist[genelist %in% names(geneweight)[geneweight<cut.off.min]]
    }
  }
  if(typeof(target)!="NULL") genelist=c(genelist,target)
  if(ncol(links)!=4) stop("links must have four columns: from, to,PMID,Source")
  colnames(links)=c("from", "to","PMID","Source")
  links=links[links$from %in% genelist & links$to %in% genelist,]
  if (nrow(links)<5){
    stop("there are too little genes that have links, please increase the gene numbers in the geneset")
  }
  links=links[!duplicated(apply(links[,1:2],1,function(x) paste(sort(x),collapse=''))),]
  if (self.interaction==F){
    links=links[links$from!=links$to,]
  }
  if(typeof(target)!="NULL"){
    links$target.links=ifelse(links$from %in% target | links$to %in% target,1,0)
    links=links[order(links$target.links),]
  }
  genes=unique(c(links$from,links$to))
  genes=genes[genes %in% names(geneweight)]
  if (typeof(names(genelist))!="NULL"){
    if (typeof(target)!="NULL"){
      nodes=data.frame(id=genes,size=1,type=ifelse(genes %in% target,"Target",names(genelist)[match(genes,genelist)]),stringsAsFactors = F)
    }else{
      nodes=data.frame(id=genes,size=1,type=names(genelist)[match(genes,genelist)],stringsAsFactors = F)
    }
  }else{
    if (typeof(target)!="NULL"){
      nodes=data.frame(id=genes,size=1,type=ifelse(genes %in% target,"Target","Pathway"),stringsAsFactors = F)
    }else{
      nodes=data.frame(id=genes,size=1,type="Pathway",stringsAsFactors = F)
    }
  }
  #https://kateto.net/netscix2016.html
  net <- graph_from_data_frame(d=links[,1:2], vertices=nodes, directed=F)
  V(net)$size <- V(net)$size*3 #set node size
  V(net)$label <- names(V(net))
  pal <- brewer.pal(length(unique(nodes$type))+2, "Dark2")
  V(net)$color <- as.factor(nodes$type[match(names(V(net)), nodes$id)])
  if (typeof(target)!="NULL"){
    E(net)$color <-ifelse(links$from %in% target | links$to %in% target,'red',ifelse(grepl("Hi-union|Lit-BM",links$Source),'grey60','cornsilk2'))
  }else{
    E(net)$color <- ifelse(grepl("Hi-union|Lit-BM",links$Source),'grey60','cornsilk2')
  }
  my_palette1 <- colorRampPalette(c("dodgerblue3", "white"))(n = 10000)
  my_palette2 <- colorRampPalette(c("white", "red"))(n = 10000)
  geneweight=geneweight[names(geneweight) %in% V(net)$label]
  label.colors<- c(setNames(my_palette1[cut(c(0,geneweight[geneweight<0]),10000)][-1],names(geneweight)[geneweight<0]),
                   setNames(my_palette2[cut(c(0,geneweight[geneweight>0]),10000)][-1],names(geneweight)[geneweight>0]))
  node.colors=label.colors[match(V(net)$label,names(label.colors))]
  E(net)$width <-  9
  E(net)$arrow.size <-0
  par(mar=c(1,1,1,1))
  mycircle <- function(coords, v=NULL, params) {
    vertex.color <- params("vertex", "color")
    if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
    }
    vertex.size  <- 1/200 * params("vertex", "size")
    if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
    }
    vertex.frame.color <- params("vertex", "frame.color")
    if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
    }
    vertex.frame.width <- params("vertex", "frame.width")
    if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
    }
    mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
           vertex.size, vertex.frame.width,
           FUN=function(x, y, bg, fg, size, lwd) {
             symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                     circles=size, add=TRUE, inches=FALSE)
           })
  }
  add.vertex.shape("fcircle", clip=igraph.shape.noclip,
                   plot=mycircle, parameters=list(vertex.frame.color=1,
                                                  vertex.frame.width=1))
  pdf(file=pdf.file,width=100, height=50)
  plot(net,vertex.color=node.colors,vertex.shape="fcircle", vertex.frame.color=V(net)$color,
       vertex.frame.width=15,vertex.label.color=V(net)$color,vertex.size=5,vertex.label.cex=5,edge.width=9,vertex.label.color="black",vertex.label.degree = 6,vertex.label.dist=1,asp=1,layout=layout_nicely)
  graphics.off()
}



