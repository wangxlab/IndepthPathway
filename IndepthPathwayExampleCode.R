
###########################################################################################################
##IMPORTANT: this code block is required for all calculations including both uniConSig and CSEA.
setwd("YourWorkingDirectory")

library(IndepthPathway)
LoadPackage() #Must run: this step automatically install and load all packages required for DGE and IndepthPathway analysis.
# data(feature.list)
# data(preCalmatrix)
# data(compare.list)
# data(ExpData) # Example gct data
# data(ClassData) # Example cls data
####################################################################################################################
##Option I. Perform WCSEA for pathway analysis based on a weighted gene list (for example, the statistical values from DGE analysis).
##This method is more applicable to scRNAseq.
####################################################################################################################

##### ==== If you use your own .gct and .cls files ##### ====
#provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
#please refer to: https://www.genepattern.org/file-formats-guide

# gctFile="/YourDirectory/scData_quiescentVsActive.gct"
# ExpData <- read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")

# clsFile="/YourDirectory/scData_quiescentVsActive.cls"
# ClassData <- read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, comment.char="",row.names=NULL, check.names=F, fill=TRUE)


##### ==== Users can choose SCDE or Limma to calculate signed q-value. ##### ====
## DGE Option 1 The below code perform DGE analysis based on Limma for bulk and single cell gene expression data. This is the profered method.
## provide .gct file that contains bulk or single cell gene expression data and .cls file that defines the cell groups, as in the below example.
## for the format of .gct and .cls files, please refer to: https://www.genepattern.org/file-formats-guide
LimmaOut=limmaDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active"))
weight=setNames(LimmaOut$Signed.Q.Value,row.names(LimmaOut))#signed q values will be used as weights for WCSEA analysis

## DGE Option 2 The below code perform DGE analysis based on SCDE for single cell gene expression data.
##This step takes about 10~30 min according to the gct size.
# SCDEOut <- scdeDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active"))
# weight=setNames(SCDEOut$Signed.Q.Value,row.names(SCDEOut))#signed q values will be used as weights for WCSEA analysis

## To run SCDE for Differential Expression Gene (DEG) calculation, please install 'flexmix ver.2.3-13" first and install 'scde ver.1.99.1' or 'ver.1.99.2'
## 1) Install 'flexmix ver.2.3-13"
# require(devtools)
# install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")
## 2) Install 'scde v1.99.1'
## a. Go to https://hms-dbmi.github.io/scde/package.html and download 'scde-1.99.1.tar'
## b. Go to terminal and run the line below
## $ R CMD INSTALL /YourOwnDirectlry/scde-1.99.1.tar
## R CMD INSTALL /Users/lees18/Library/CloudStorage/OneDrive-UPMC/CCBR-XWANGLAB10/19_SymSim_wCSEA2/10_0_indepthPathway_Code/scde-1.99.1.tar
## 3) After installing flexmix and scde, restart R. After restarting, Check flexmix and scde versions.
#library(flexmix)
#library(scde)
# sessionInfo()

#calculate uniConSig cores that compute functional relations of human genes underlying the highly weighted genes.
#This step took 20 mins when testing, it is likely to take a long time in large dataset
ks.result=run.weightedKS(weight,signed=T,feature.list)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list)

# perform pathway enrichment analysis
# The pathway enrighment analysis took 10 mins each when testing
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)

#disambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to disambiguate
up.deambiguate<-deambiguation.CSEA(GEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,
                                     topn=min(c(topn,nrow(up.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni")  #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
down.deambiguate<-deambiguation.CSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,
                                       topn=min(c(topn,nrow(down.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni") #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
#View filtered top pathways
View(up.deambiguate$pathway.filtered)
View(down.deambiguate$pathway.filtered)

####################################################################################################################
#VISULIZATION OF PATHWAYS
#option 1-2 can be used for WCSEA and GSEA. option 3-5 can be used for CSEA, WCSEA, and GSEA
####################################################################################################################
####Option 1: draw heatmaps showing the functional associations between selected top pathways___
## Calculation the functional associations of the selected top pathways revealed by WCSEA analysis and return a similarity matrix between top pathways
## The PathwayAssociation took 30~45 mins according to your computer specs.
## Minsize specifies the cutoff of min number of genes required the provided pathways. To compute pathway association it is recommended that the selected pathways have at least 10 genes.
selectn=30 #specify the number of top pathways to compute functional associations.
up.assoc <- pathwayAssociation(topPathway=up.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(up.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(down.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
pathway.heatmap(matrixData=up.assoc,clustering=T,fontSize=5)
pathway.heatmap(matrixData=down.assoc,clustering=T,fontSize=5)

####Option 2: draw network showing the functional associations between selected top pathways___
## It merges up-regulated and down-regulated pathways. In the rare case that a pathway is ranked to the top in both up and downregulated pathway lists, the most significant rank will be retained.
## From the merged pathways, pathwayAssociation() will calculate the functional associations
## In draw.network(), the parameter, NES.cut specifies the cut off to filter the edges in they pathway association network. The higher the NES.cut, the less edges will be shown. It is recommended to set NES.cut between 1.5 and 2.
#node.size specifies the size of the nodes, line.thickness specifies the thickness of the lines,font.size specifies the font sizes of the pathway names,wrap.string.width specifies the width to wrap the pathway names
selectn=30
pathway.merge=merge.pathway(up.pathway=up.deambiguate[[1]][1:min(c(selectn,nrow(up.deambiguate[[1]]))),],down.pathway=down.deambiguate[[1]][1:min(c(selectn,nrow(down.deambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="WCSEA_PathwayAssocHeatmapNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
graphics.off()

####Option 3: Draw enrichment plot for a selected pathway___
## This generates an enrichment plot of a specific pathway as in GSEA.
Pathway="HALLMARK_E2F_TARGETS"
plotEnrichment(compare.list[[Pathway]], normalize(weight))
weight.order=setNames(weight,order(weight,decreasing=T))

####Option 4: Draw heat map for a selected pathway. annCol is optional.___
## The parameter, filter.coreEnrich=TRUE, will filter core enriched genes in the pathway geneset for heatmap
Pathway="HALLMARK_E2F_TARGETS"
exp.list=read.gct(gctFile="./scDataset/scData_quiescentVsActive.gct",clsFile="./scDataset/scData_quiescentVsActive.cls")
expData=exp.list$expData
sampleClass=exp.list$class
expData=expData[,rownames(sampleClass)[order(sampleClass$Class)]]
plot.geneset.matrix(genomatrix=expData,geneset=compare.list[[Pathway]],weights=weight,annCol = sampleClass,filter.coreEnrich=T)

####Option 5: Draw interactome network for a selected pathway and your target gene of interests___
## This is a powerful tool for exploring novel mechanisms for your target gene.
## It allows you to specify the target gene of interest, so you can analyze the interactions of your gene of interests with the genes in the selected enriched pathways
## Parameters: links, input interaction pairs. genelist, input genes to be ploted. geneweight, input the differential expression ratios or
## statistics, target, input the target gene of interest (i.e. target="MYC"). self.interaction, set to false to filter homodimerizations,
## filter.coreEnrich, only plot core enriched genes in the selected pathways. Set to True to reduce the complexity of the network.
interactions=read.delim("./PathwayDb/interactome_all.hs.sym",stringsAsFactors = F,check.names = F,header = T,sep="\t") # interactome_all.hs.sym contains interactome data from HuRI and Entrez Gene
geneset.select=compare.list[c("REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX")] # provide the standard names of the selected pathways to plot. The names need to be exactly the same as in the pathway enrichment result table.
pdf.file="./Results/interactome.network.pdf"
genelist = setNames(unlist(geneset.select),rep(names(geneset.select), times = sapply(geneset.select, length))) #generate the gene list based on the selected pathways to plot
genelist=genelist[!duplicated(genelist)] #remove duplicated genes
geneweight=setNames(limma$T.Value,row.names(limma)) #set the weight of the genes based on signed statistical value from DGE analysis
draw.interactome.network(links=interactions,genelist=genelist,geneweight=geneweight,target="MYC",pdf.file=pdf.file,self.interaction=F,filter.coreEnrich=T)

####################################################################################################################
##Option II. Perform CSEA for Pathway Enrichment Analysis based on an experimentally defined gene set (i.e., upregulated or downregulated genes)
####################################################################################################################

#Perform CSEA for dichotomous experimental gene list from scDataset
#compare.list=c(read_concepts("PathwayDb/h.all.v7.5.1.symbols.gmt"),read_concepts("PathwayDb/c2.cp.v7.5.1.symbols.gmt"),read_concepts("PathwayDb/c2.cgp.v2022.1.Hs.symbols.gmt"))
#target.list<-read.table ("scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1

#here we use FDR q value cut off of 0.1 to select up or down-regulated genes
target.list_Up <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value > log10(0.1)]; length(target.list_Up)
target.list_Dw <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value < -log10(0.1)]; length(target.list_Dw)

#perform deep functional interpretation of the target gene list and calculate uniConSig scores.The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig_Up=cal.uniConSig(target.list=target.list_Up,feature.list=feature.list,preCalmatrix,rm.overfit=F)
#The CSEA2 function took 10 mins when testing
CSEA.result_Up<-CSEA2(setNames(as.numeric(uniConSig_Up$uniConSig), uniConSig_Up$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways

uniConSig_Dw=cal.uniConSig(target.list=target.list_Dw,feature.list=feature.list,preCalmatrix,rm.overfit=F)
#The CSEA2 function took 10 mins when testing
CSEA.result_Dw<-CSEA2(setNames(as.numeric(uniConSig_Dw$uniConSig), uniConSig_Dw$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways

#disambiguate top enriched pathways.
#This step took 20 mins when testing
topn=100 #specify the number of top pathways to disambiguate
deambiguate_CSEA_Up<-deambiguation.CSEA(GEA.result=CSEA.result_Up,uniConSig.result=uniConSig_Up,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result_Up))),
                                  p.cut=0.01,p.adjust.method="bonferroni")#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
topn=100 #specify the number of top pathways to disambiguate
deambiguate_CSEA_Dw<-deambiguation.CSEA(GEA.result=CSEA.result_Dw,uniConSig.result=uniConSig_Dw,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result_Dw))),
                                        p.cut=0.01,p.adjust.method="bonferroni")#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
#compute functional associations between selected top pathways
#assoc <- pathwayAssociation(topPathway=deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(disambiguate[[1]])))],compare.list,feature.list,preCalmatrix)
CSEA_up.assoc <- pathwayAssociation(DeambiguatePathway=deambiguate_CSEA_Up,  compare.list, feature.list, preCalmatrix,selectn=30, minsize=10)
CSEA_down.assoc <- pathwayAssociation(DeambiguatePathway=deambiguate_CSEA_Dw, compare.list,feature.list,preCalmatrix,selectn=30, minsize=10)

#draw heatmaps showing the functional associations between selected top pathways
pdf(file="CSEA_Up_PathwayAssocHeatmapNetwork_20230225.pdf",width=10, height=10)
pathway.heatmap(matrixData=CSEA_up.assoc,clustering = TRUE,fontSize=8)

#draw network figure for top pathways. NES.cut determines the levels of significance for the pathway similarity to be shown as edges. The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=deambiguate_CSEA_Up[[1]],assoc=CSEA_up.assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
graphics.off()

####################################################################################################################
##Option III. Perform GSEA for pathway analysis
##This method is more applicable to pathway analysis of bulk gene expression data.
####################################################################################################################
# perform pathway enrichment analysis (!Please make sure to input gene names with standard HGNC symbols)
up.GSEA.result<-CSEA2(target.score=weight,compare.list,p.cut=0.05)

down.GSEA.result<-CSEA2(target.score=-weight,compare.list,p.cut=0.05)

#disambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to disambiguate
up.deambiguate<-deambiguation.GSEA(GEA.result=up.GSEA.result,weight=weight,compare.list=compare.list,topn=min(c(topn,nrow(up.GSEA.result))),
                                     p.cut=0.05,p.adjust.method="bonferroni")#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
down.deambiguate<-deambiguation.GSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=min(c(topn,nrow(down.GSEA.result))),
                                       p.cut=0.05,p.adjust.method="bonferroni")
