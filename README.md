# IndepthPathway
 
# How to use IndepthPathway to perform WCSEA, CSEA, and GSEA pathway analysis

__IndepthPathway: an integrated tool for in-depth pathway enrichment analysis based on bulk and single cell sequencing data__

The current version of IndepthPathway is 2.8.6. This was tested on R v4.2.1

## A. Introduction

Single-cell sequencing (SCS) enables exploring the pathways and processes of cells and cell populations, however, there is a paucity of pathway enrichment methods designed to tolerate the high noise and low gene coverage of SCS technology. When gene expression data are noisy and signals are sparse, testing pathway enrichment based on the genes measured may not yield statistically significant results which is particularly problematic when detecting the pathways enriched in less abundant cells that are vulnerable to disturbances.

Here we developed a Weighted Concept Signature Enrichment Analysis (WCSEA) algorithm specialized for pathway enrichment analysis from single cell transcriptomics (scRNA-seq), taking account of the levels of differential expressions to detect different magnitudes of pathway alterations, and substantially improve its noise resistance. WCSEA took a broader approach for assessing the functional relations of pathway gene sets to a target gene list, and leverage the universal concept signature of the target gene list (the cumulative signature of molecular concepts characteristic of the target gene list), to tolerate the high noise and low coverage of this technology.


Leveraging its unique strength, IndepthPathway will promote the application of bulk and single cell sequencing technologies to explore the cellular pathway mechanisms at more precise resolution.

## B.	Basic requirements

The uniConSig and CSEA modules are compiled in an R package "IndepthPathway" held at https://github.com/wangxlab/IndepthPathway 

To install the package, first install R from CRAN: https://cran.r-project.org/

For more user-friendly interface, R-Studio can be installed from here: https://www.rstudio.com/products/rstudio/download/


## C. Installing R Packages

```
install.packages("devtools")
library(devtools)

devtools::install_github("wangxlab/IndepthPathway")
library(IndepthPathway)

IndepthPathway::LoadPackage() ##MUST RUN: this step automatically install and load all packages required for DGE and IndepthPathway analysis.
```

## D.	How to run IndepthPathway for pathway enrichment analysis

### Example codes for IndepthPathway are available at IndepthPathwayExampleCode.R

Check the molecule concept dataset, preCalmatrix, and compare.list we compiled.
Check the example gct and cls data

```
setwd("your working directory")
data(feature.list)
data(preCalmatrix)
data(compare.list)
data(ExpData) # Example gct data
data(ClassData) # Example cls data
```

### Option I. Perform WCSEA for pathway analysis based on a weighted gene list (for example, the signed statistical values from DGE analysis). 
##This method generate more reproducible and indepth enrichment results for functional pathways based on scRNAseq data than GSEA. 
___The W-CSEA analysis took 50 mins in total when testing, it might vary depends on the dataset___
```

##### Users can choose SCDE or Limma to calculate signed q-value. 
## DGE Option 1 The below code perform DGE analysis based on Limma for bulk and single cell gene expression data. This is the profered method.
## provide .gct file that contains bulk or single cell gene expression data and .cls file that defines the cell groups, as in the below example.
## for the format of .gct and .cls files, please refer to: https://www.genepattern.org/file-formats-guide
limmaOut=limmaDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active"))
weight=setNames(limmaOut$Signed.Q.Value,row.names(limmaOut))#signed q values will be used as weights for WCSEA analysis

## DGE Option 2 The below code perform DGE analysis based on SCDE for single cell gene expression data. 
## Please install 'flexmix v2.3-13' and 'scde v1.99.1' first, which is explained in the beginning of "IndepthPathway_ExampleCode_v2.8.5.R" code
##This step takes about 10~30 min according to the gct size. 
# SCDEOut <- scdeDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active")) 
# weight=setNames(SCDEOut$Signed.Q.Value,row.names(SCDEOut))#signed q values will be used as weights for WCSEA analysis

## calculate uniConSig cores that compute functional relations of human genes underlying the highly weighted genes.
## This step took 20 mins when testing, it is likely to take a long time in large dataset
ks.result=run.weightedKS(weight,signed=T,feature.list)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list)

## perform pathway enrichment analysis
## The pathway enrighment analysis took 10 mins each when testing
## Parameter, p.cut, is to get significant pathways with the p.cut threshold. To get the unfiltered results for all pathways, set p.cut=2.
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)

## Disambiguate up and down-regulated pathways
## Parameter, p.cut, is to get significant pathways with the p.cut threshold. 
## upPathways set to TRUE or FALSE (for W-CSEA).
## Users can choose p.adjust.method: NULL, "bonferroni" or "BH"
topn=30 #specify the number of top pathways to disambiguate
up.disambiguate<-disambiguationCSEA(GEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,
                               topn=min(c(topn,nrow(up.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni")  #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
down.disambiguate<-disambiguationCSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,
                               topn=min(c(topn,nrow(down.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni") #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"

```
### VISULIZATION OF PATHWAY ENRICHMENT ANALYSIS RESULTS
___Option 1: draw heatmaps showing the functional associations between selected top pathways___
```
## Calculation the functional associations of the selected top pathways revealed by WCSEA analysis and return a similarity matrix between top pathways
## The PathwayAssociation took 30~45 mins according to your computer specs.
## Minsize specifies the cutoff of min number of genes required the provided pathways. To compute pathway association it is recommended that the selected pathways have at least 10 genes.
selectn=30 #specify the number of top pathways to compute functional associations. 
up.assoc <- pathwayAssociation(topPathway=up.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(up.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(down.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
pathway.heatmap(matrixData=up.assoc,clustering=T,fontSize=5)
pathway.heatmap(matrixData=down.assoc,clustering=T,fontSize=5)
```

___Option 2: draw network showing the functional associations between selected top pathways___
```
## It merges up-regulated and down-regulated pathways. In the rare case that a pathway is ranked to the top in both up and downregulated pathway lists, the most significant rank will be retained.
## From the merged pathways, pathwayAssociation() will calculate the functional associations
## In draw.network(), the parameter, NES.cut specifies the cut off to filter the edges in they pathway association network. The higher the NES.cut, the less edges will be shown. It is recommended to set NES.cut between 1.5 and 2.
#node.size specifies the size of the nodes, line.thickness specifies the thickness of the lines,font.size specifies the font sizes of the pathway names,wrap.string.width specifies the width to wrap the pathway names
selectn=30
pathway.merge=mergePathway(up.pathway=up.deambiguate[[1]][1:min(c(selectn,nrow(up.deambiguate[[1]]))),],down.pathway=down.deambiguate[[1]][1:min(c(selectn,nrow(down.deambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="WCSEA_PathwayAssocHeatmapNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
graphics.off()
```

___Option 3: Draw enrichment plot for a selected pathway___
```
## This generates an enrichment plot of a specific pathway as in GSEA.
Pathway="HALLMARK_E2F_TARGETS"
plotEnrichment(compare.list[[Pathway]], normalize(weight))
weight.order=setNames(weight,order(weight,decreasing=T))
```

___Option 4: Draw heat map for a selected pathway. annCol is optional.___
```
## The parameter, filter.coreEnrich=TRUE, will filter core enriched genes in the pathway geneset for heatmap
Pathway="HALLMARK_E2F_TARGETS"
exp.list=read.gct(gctFile="./scDataset/scData_quiescentVsActive.gct",clsFile="./scDataset/scData_quiescentVsActive.cls")
expData=exp.list$expData
sampleClass=exp.list$class
expData=expData[,rownames(sampleClass)[order(sampleClass$Class)]]
plot.geneset.matrix(genomatrix=expData,geneset=compare.list[[Pathway]],weights=weight,annCol = sampleClass,filter.coreEnrich=T)
```

___Option 5: Draw interactome network for a selected pathway and your target gene of interests___
```
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
```

## Option II. Perform CSEA for Pathway Enrichment Analysis based on an experimentally defined gene set (i.e., upregulated or downregulated genes)
## Perform CSEA for dichotomous experimental gene list (i.e., upregulated or downregulated genes)
___The CSEA analysis took 50 mins in total when testing, it might vary depends on the dataset___

```
#here we use FDR q value cut off of 0.1 to select up or down-regulated genes
target.list <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value > log10(0.1)]; length(target.list_Up) # select upregulated genes
#target.list <- rownames(LimmaOut)[LimmaOut$Signed.Q.Value < -log10(0.1)]; length(target.list_Dw) # select downregulated genes

## perform deep functional interpretation of the target gene list and calculate uniConSig scores.The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,rm.overfit=F)
## The CSEA2 function took 10 mins when testing
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways

## Disambiguate top enriched pathways.
## This step took 20 mins when testing
topn=100 #specify the number of top pathways to deambiguate
disambiguate<-disambiguationCSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result))),p.cut=0.01)

## Compute functional associations between selected top pathways
selectn=30 #specify the number of top pathways to compute associations
assoc <- pathwayAssociation(topPathway=deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(deambiguate[[1]])))],compare.list,feature.list,preCalmatrix)

## Draw heatmaps showing the functional associations between selected top pathways 
pdf(file="CSEA_PathwayAssocHeatmapNetwork.pdf",width=10, height=10)
pathway.heatmap(matrixData=assoc,clustering = TRUE,fontSize=8)

## Draw network figure. NES.cut determines the levels of significance for the pathway similarity to be shown as edges. The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=deambiguate[[1]],assoc=assoc,NES.cut=2)
graphics.off()
```

## Option III. Perform GSEA for pathway analysis
```
## This method is more appliable to pathway analysis of bulk gene expression data.
## Perform limma DGE analysis for single cell gene expression data. 
## Provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
## Please refer to: https://www.genepattern.org/file-formats-guide
weight=setNames(limma$Signed.Q.Value,row.names(limma))#signed q values will be used as weights for GSEA analysis

## Perform pathway enrichment analysis
up.GSEA.result<-CSEA2(target.score=weight,compare.list,p.cut=0.05)
down.GSEA.result<-CSEA2(target.score=-weight,compare.list,p.cut=0.05)

## Disambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to deambiguate
up.disambiguate<-disambiguationGSEA(GEA.result=up.GSEA.result,weight=weight,compare.list=compare.list,topn=min(c(topn,nrow(up.GSEA.result))),p.cut=0.01)
down.disambiguate<-disambiguationGSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=min(c(topn,nrow(down.GSEA.result))),p.cut=0.01)

###VISULIZATION OF PATHWAY ASSOCIATION NETWORK (Option 2).
## Draw network showing the functional associations between selected top pathways 
selectn=30
pathway.merge=mergePathway(up.pathway=up.deambiguate[[1]][1:min(c(selectn,nrow(up.deambiguate[[1]]))),],down.pathway=down.deambiguate[[1]][1:min(c(selectn,nrow(down.deambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="GSEA_PathwayAssocHeatmapNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2)
graphics.off()
```


```
