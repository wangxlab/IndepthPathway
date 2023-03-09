# IndepthPathway
 
# How to use IndepthPathway to perform WCSEA, CSEA, and GSEA pathway analysis

__IndepthPathway: an integrated tool for in-depth pathway enrichment analysis based on bulk and single cell sequencing data__

The current version of IndepthPathway is 2.8.5. This was tested on R v4.2.1

## A. Introduction

â€¢	Single-cell sequencing (SCS) enables exploring the pathways and processes of cells and cell populations, however, there is a paucity of pathway enrichment methods designed to tolerate the high noise and low gene coverage of SCS technology. When gene expression data are noisy and signals are sparse, testing pathway enrichment based on the genes measured may not yield statistically significant results which is particularly problematic when detecting the pathways enriched in less abundant cells that are vulnerable to disturbances.

â€¢	Here we developed a Weighted Concept Signature Enrichment Analysis (WCSEA) algorithm specialized for pathway enrichment analysis from single cell transcriptomics (scRNA-seq), taking account of the levels of differential expressions to detect different magnitudes of pathway alterations, and substantially improve its noise resistance. WCSEA took a broader approach for assessing the functional relations of pathway gene sets to a target gene list, and leverage the universal concept signature of the target gene list (the cumulative signature of molecular concepts characteristic of the target gene list), to tolerate the high noise and low coverage of this technology.

â€¢	We then incorporated WCSEA into a R package called â€œIndepthPathwayâ€ for biologists to broadly leverage this method for pathway analysis based on bulk and single cell sequencing data. Through simulating the technical variability and dropouts in gene expression characteristic of scRNA-seq, WCSEA yielded overall low deviations in pathway enrichment results. This could be attributed to the computation of the universal concept signature prior to pathway enrichment analysis, which make WCSEA more resistant to noise and missing values of individual gene expressions.

â€¢	Leveraging its unique strength, IndepthPathway will promote the application of bulk and single cell sequencing technologies to explore the cellular pathway mechanisms at more precise resolution.

## B.	Basic requirements

â€¢	The uniConSig and CSEA modules are compiled in an R package â€œIndepthPathwayâ€ held at https://github.com/wangxlab/IndepthPathway 

â€¢	To install the package, first install R from CRAN: https://cran.r-project.org/

â€¢	For more user-friendly interface, R-Studio can be installed from here: https://www.rstudio.com/products/rstudio/download/


## C. Installing R Packages

```
install.packages("devtools")
library(devtools)

devtools::install_github("wangxlab/IndepthPathway")
library(IndepthPathway)

LoadPackage() ## Load necessary packages to run IndepthPathway
```

## D.	How to run IndepthPathway for pathway enrichment analysis
# Example codes for IndepthPathway are available at IndepthPathwayExampleCode.R

â€¢	Check the molecule concept dataset, preCalmatrix, and compare.list we compiled.
â€¢ Check the example gct and cls data
```
setwd("your working directory")
data(feature.list)
data(preCalmatrix)
data(compare.list)
data(ExpData) # Example gct data
data(ClassData) # Example cls data
```

### Option I. Perform WCSEA for pathway analysis based on a weighted gene list (for example, the statistical values from DGE analysis). This method is more appliable to scRNAseq. 
___The W-CSEA analysis took 50 mins in total when testing, it might vary depends on the dataset___
```

##### Users can choose SCDE or Limma to calculate signed q-value. 
#You need to install 'flexmix v2.3-13' and 'scde v1.99.1' first. How to install thme is explained in the beginning of "IndepthPathway_ExampleCode_v2.8.5.R" code
#This is by SCDE.  This takes about 10~30 min according to the gct size. 
SCDEOut <- scdeDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active")) 
weight=setNames(SCDEOut$Signed.Q.Value,row.names(SCDEOut))#signed q values will be used as weights for WCSEA analysis

## This is by Limma. 
#perform limma DGE analysis for single cell gene expression data. 
#provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
#please refer to: https://www.genepattern.org/file-formats-guide
# limmaOut=limmaDEG(ExpDataGCT=ExpData, ClassData=ClassData, groups.order=c("Quiescent","Active"))
# weight=setNames(limmaOut$Signed.Q.Value,row.names(limmaOut))#signed q values will be used as weights for WCSEA analysis

#calculate uniConSig cores that compute functional relations of human genes underlying the highly weighted genes. correct.overfit should be set to False for WCSEA analysis
#This step took 20 mins when testing, it is likely to take a long time in large dataset
ks.result=run.weightedKS(weight,signed=T,feature.list,correct.overfit=F)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,correct.overfit=FALSE)

# perform pathway enrichment analysis
# The pathway enrighment analysis took 10 mins each when testing
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)

#deambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to deambiguate
up.deambiguate<-deambiguation.CSEA(GEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,
                               topn=min(c(topn,nrow(up.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni")  #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
down.deambiguate<-deambiguation.CSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,
                               topn=min(c(topn,nrow(down.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni") #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
									 

___upPathways use TRUE/FALSE for W-CSEA or use NULL for CSEA___

```
### VISULIZATION OF PATHWAYS
___Option 1: draw heatmaps showing the functional associations between selected top pathways___
```
selectn=30 #specify the number of top pathways to compute associations. Minsize specifies the cutoff of min number of genes required the provided pathway
#The PathwayAssociation took 30 mins when testing
up.assoc <- pathwayAssociation(topPathway=up.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(up.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(down.deambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
pathway.heatmap(matrixData=up.assoc,clustering = T,fontSize=5)
pathway.heatmap(matrixData=down.assoc,clustering = T,fontSize=5)
```

___Option 2: draw network showing the functional associations between selected top pathways___
```
selectn=30
pathway.merge=merge.pathway(up.pathway=up.deambiguate[[1]][1:min(c(selectn,nrow(up.deambiguate[[1]]))),],down.pathway=down.deambiguate[[1]][1:min(c(selectn,nrow(down.deambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="WCSEA_PathwayAssocHeatmapNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2)
graphics.off()
```


## Option II. Perform CSEA for Pathway Enrichment Analysis based on an experimentally defined gene set (i.e., upregulated or downregulated genes)
### Perform CSEA for dichotomous experimental gene list from scDataset
___The CSEA analysis took 50 mins in total when testing, it might vary depends on the dataset___

```
target.list <- rownames(SCDEOut)[SCDEOut$Signed.Q.Value > 0]; length(target.list_Up)

#perform deep functional interpretation of the target gene list and calculate uniConSig scores.The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,rm.overfit=F)
#The CSEA2 function took 10 mins when testing
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways

#deambiguate top enriched pathways.
#This step took 20 mins when testing
topn=100 #specify the number of top pathways to deambiguate
deambiguate<-deambiguation.CSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result))),p.cut=0.01)

#compute functional associations between selected top pathways
selectn=30 #specify the number of top pathways to compute associations
assoc <- pathwayAssociation(topPathway=deambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(deambiguate[[1]])))],compare.list,feature.list,preCalmatrix)

#draw heatmaps showing the functional associations between selected top pathways 
pdf(file="CSEA_PathwayAssocHeatmapNetwork.pdf",width=10, height=10)
pathway.heatmap(matrixData=assoc,clustering = TRUE,fontSize=8)

#draw network figure. NES.cut determines the levels of significance for the pathway similarity to be shown as edges. The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=deambiguate[[1]],assoc=assoc,NES.cut=2)
graphics.off()
```

## Option III. Perform GSEA for pathway analysis
```
#This method is more appliable to pathway analysis of bulk gene expression data.
#perform limma DGE analysis for single cell gene expression data. 
#provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
#please refer to: https://www.genepattern.org/file-formats-guide
weight=setNames(limma$Signed.Q.Value,row.names(limma))#signed q values will be used as weights for GSEA analysis

# perform pathway enrichment analysis
up.GSEA.result<-CSEA2(target.score=weight,compare.list,p.cut=0.05)
down.GSEA.result<-CSEA2(target.score=-weight,compare.list,p.cut=0.05)

#deambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to deambiguate
up.deambiguate<-deambiguation.GSEA(GEA.result=up.GSEA.result,weight=weight,compare.list=compare.list,topn=min(c(topn,nrow(up.GSEA.result))),p.cut=0.01)
down.deambiguate<-deambiguation.GSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=min(c(topn,nrow(down.GSEA.result))),p.cut=0.01)

###VISULIZATION OF PATHWAYS
#draw network showing the functional associations between selected top pathways 
selectn=30
pathway.merge=merge.pathway(up.pathway=up.deambiguate[[1]][1:min(c(selectn,nrow(up.deambiguate[[1]]))),],down.pathway=down.deambiguate[[1]][1:min(c(selectn,nrow(down.deambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="GSEA_PathwayAssocHeatmapNetwork.pdf",width=20, height=20)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2)
graphics.off()
```
