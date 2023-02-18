
## To run SCDE for Differential Expression Gene (DEG) calculation, you need to to install 'flexmix ver.2.3-13" first and install 'scde ver.1.99.1' or 'ver.1.99.2'

## 1) Install 'flexmix ver.2.3-13"
require(devtools)
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")

## 2) Install 'scde v1.99.1'
## a. Go to https://hms-dbmi.github.io/scde/package.html and download 'scde-1.99.1.tar'
## b. Go to terminal and run the line below
## $ R CMD INSTALL /YourOwnDirectlry/scde-1.99.1.tar
## R CMD INSTALL /Users/lees18/Library/CloudStorage/OneDrive-UPMC/CCBR-XWANGLAB10/19_SymSim_wCSEA2/10_0_indepthPathway_Code/scde-1.99.1.tar

## 3) After installing flexmix and scde, restart R. After restarting, Check flexmix and scde versions.
library(flexmix)
library(scde)
sessionInfo()

###########################################################################################################
##IMPORTANT: this code block is required for all calculations including both uniConSig and CSEA.
setwd("/Users/lees18/Library/CloudStorage/OneDrive-UPMC/CCBR-XWANGLAB10/19_SymSim_wCSEA2/52_1_GithubRcode_IndepthPathway_v2.8.5")

data(feature.list)
data(preCalmatrix)
data(compare.list)
data(ExpData)
data(ClassData)
####################################################################################################################
##Option I. Perform W-CSEA for pathway analysis based on a weighted gene list (for example, the statistical values from DGE analysis).
##This method is more appliable to scRNAseq.
####################################################################################################################
gctFile="/Users/lees18/Library/CloudStorage/OneDrive-UPMC/CCBR-XWANGLAB10/19_SymSim_wCSEA2/10_0_indepthPathway_Code_BySCDE_HSC/indepthPathway_Original/scDataset_GSE68981/scData_quiescentVsActive.gct"
clsFile="/Users/lees18/Library/CloudStorage/OneDrive-UPMC/CCBR-XWANGLAB10/19_SymSim_wCSEA2/10_0_indepthPathway_Code_BySCDE_HSC/indepthPathway_Original/scDataset_GSE68981/scData_quiescentVsActive.cls"

ExpData <- read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
ClassData <- read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, comment.char="",row.names=NULL, check.names=F, fill=TRUE)

use_data_raw()
usethis::use_data(ExpData, compress="xz")
usethis::use_data(ClassData, compress="xz")


##### Users can choose SCDE or Limma to calculate signed q-value.
## This is by SCDE.  This takes about 10~30 min according to the gct size.
SCDEOut <- scdeDEG(readCountGCT=gctFile, clsFile=clsFile, groups.order=c("Quiescent","Active"))
saveRDS(SCDEOut, "SCDEOut_quiescentVsActive.rds")
weight=setNames(SCDEOut$Signed.Q.Value,row.names(SCDEOut))#signed q values will be used as weights for WCSEA analysis

## This is by Limma.
#perform limma DGE analysis for single cell gene expression data.
#provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
#please refer to: https://www.genepattern.org/file-formats-guide
# limma=limmaDGE(gctFile="./indepthPathway_Original/scDataset_GSE68981/scData_quiescentVsActive.gct", clsFile="./indepthPathway_Original/scDataset_GSE68981/scData_quiescentVsActive.cls",groups.order=c("Quiescent","Active"))
# saveRDS(limma, file="LimmaOut_SignedQVal_quiescentVsActive.rds")
# weight=setNames(limma$Signed.Q.Value,row.names(limma))#signed q values will be used as weights for WCSEA analysis

#calculate uniConSig cores that compute functional relations of human genes underlying the highly weighted genes.
#This step took 20 mins when testing, it is likely to take a long time in large dataset
ks.result=run.weightedKS(weight,signed=T,feature.list)
saveRDS(ks.result, file="WeightedKSResult_BySCDE_quiescentVsActive.rds")
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile)
saveRDS(uniConSig.result, file="uniConSig.result_BySCDE_quiescentVsActive.rds")

# perform pathway enrichment analysis
# The pathway enrighment analysis took 10 mins each when testing
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
saveRDS(up.CSEA.result, file="up.WCSEA.result_BySCDE_quiescentVsActive.rds")
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
saveRDS(down.CSEA.result, file="down.WCSEA.result_BySCDE_quiescentVsActive.rds")

#disambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to disambiguate
up.deambiguate<-disambiguation.CSEA(GEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,
                                     topn=min(c(topn,nrow(up.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni")  #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(up.disambiguate, file="up.disambiguate_BySCDE_quiescentVsActive.rds")

down.deambiguate<-disambiguation.CSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,
                                       topn=min(c(topn,nrow(down.CSEA.result))),p.cut=0.01,p.adjust.method="bonferroni") #p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
saveRDS(down.disambiguate, file="down.disambiguate_BySCDE_quiescentVsActive.rds")

#View filtered top pathways
View(up.disambiguate$pathway.filtered)
View(down.disambiguate$pathway.filtered)

####################################################################################################################
#VISULIZATION OF PATHWAYS
#option 1-2 can be used for WCSEA and GSEA. option 3-5 can be used for CSEA, WCSEA, and GSEA
####################################################################################################################
#Option 1: draw heatmaps showing the functional associations between selected top pathways
#selectn=30 #specify the number of top pathways to compute associations. Minsize specifies the cutoff of min number of genes required the provided pathway
#The PathwayAssociation took 30 mins when testing

#up.assoc <- pathwayAssociation(topPathway=up.disambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(up.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)

up.assoc <- pathwayAssociation(DeambiguatePathway=up.deambiguate,  compare.list, feature.list, preCalmatrix,selectn=30, minsize=10)
saveRDS(up.assoc, file="up.assoc_BySCDE_Pair1A_EC_H1.rds")
down.assoc <- pathwayAssociation(DeambiguatePathway=down.deambiguate, compare.list,feature.list,preCalmatrix,minsize=10)
saveRDS(down.assoc, file="down.assoc_BySCDE_Pair1A_EC_H1.rds")

pdf(file="WCSEA_PathwayUpAssocHeatmap.pdf",width=10, height=10)
pathway.heatmap(matrixData=up.assoc,clustering = T,fontSize=5)
graphics.off()

pdf(file="WCSEA_PathwayDownAssocHeatmap.pdf",width=10, height=10)
pathway.heatmap(matrixData=down.assoc,clustering = T,fontSize=5)
graphics.off()

#Option 2: draw network showing the functional associations between selected top pathways
#when making pathway association networks, it requires CSEA to compute the functional associations between pathways in the step that calls pathwayAssociation.
## If you encounter the warning message "there are 0 genes have uniconsig scores" when calculating functional associations for a specific "pathwayx",
## it means that the program cannot find sufficient functional associations within the "pathwayx" , so that it cannot compute the associations of that pathway with other pathways.
selectn=30
pathway.merge=merge.pathway(up.pathway=up.disambiguate[[1]][1:min(c(selectn,nrow(up.disambiguate[[1]]))),],down.pathway=down.disambiguate[[1]][1:min(c(selectn,nrow(down.disambiguate[[1]]))),])
pathway.merge.assoc <- pathwayAssociation(topPathway=pathway.merge$Compare.List,compare.list,feature.list,preCalmatrix,minsize=10)
pdf(file="WCSEA_PathwayAssocHeatmapNetworkPosSignQ.pdf",width=10, height=10)
draw.network(pathway.out=pathway.merge,assoc=pathway.merge.assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
graphics.off()

#Option 3: draw enrichment plot for selected pathway
Pathway="HALLMARK_E2F_TARGETS"
pdf(file="WCSEA_Enrichmentplot.pdf",width=10, height=10)
plotEnrichment(compare.list[[Pathway]], normalize(weight))
weight.order=setNames(weight,order(weight,decreasing=T))
graphics.off()

#Option 4: draw heatmap for selected pathway. annCol is optional
Pathway="HALLMARK_E2F_TARGETS"
exp.list=read.gct(gctFile=gctFile,clsFile=clsFile)
expData=exp.list$expData
sampleClass=exp.list$class
expData=expData[,rownames(sampleClass)[order(sampleClass$Class)]]
pdf(file="WCSEA_PathwayHeatmap.pdf",width=10, height=10)
plot.geneset.matrix(genomatrix=expData,geneset=compare.list[[Pathway]],weights=weight,annCol = sampleClass)
graphics.off()

#Option 5: draw interactome network for a selected pathway and your target gene of interests. This is a powerful tool for exploring novel mechanisms for your target gene.
## It allows you to specify the target gene of interests, so you can analyze the interactions of your gene of interests with the pathways genes.
#Parameters: links, input interaction paris. genelist, input genes to be ploted. geneweight, input the differential expression ratios or statistics, target,
## input the target gene of interest (i.e. target="MYC"). self.interaction, set to false to filter homodimerizations, filter.coreEnrich, only plot core enriched genes
HuRI.interactions=read.delim("./PathwayDb/interactome_all.hs.sym",stringsAsFactors = F,check.names = F,header = T,sep="\t")
geneset.select=compare.list[c("REACTOME_RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY")]
Pathway.name="REACTOME_RESOLUTION_OF_AP_SITES_VIA_THE_MULTIPLE_NUCLEOTIDE_PATCH_REPLACEMENT_PATHWAY"
genelist = setNames(unlist(geneset.select),rep(names(geneset.select), times = sapply(geneset.select, length)))
genelist=genelist[!duplicated(genelist)]
#geneweight=setNames(limma$T.Value,row.names(limma))
geneweight=setNames(SCDEOut$Signed.P.Zscore,row.names(SCDEOut))
#draw.interactome.network(links=HuRI.interactions,genelist=genelist,geneweight=geneweight,target=NULL,pdf.file=paste0("./Results/",Pathway.name,".HuRI.network.pdf"),self.interaction=F,filter.coreEnrich=T)
InteractouOutFile <- paste0(Pathway.name,".HuRI.network.pdf")
draw.interactome.network(links=HuRI.interactions,genelist=genelist,geneweight=geneweight,target=NULL,pdf.file=InteractouOutFile,self.interaction=F,filter.coreEnrich=T)

####################################################################################################################
##Option II. Perform CSEA for Pathway Enrichment Analysis based on an
##experimentally defined gene set (i.e., upregulated or downregulated genes)
####################################################################################################################

#Perform CSEA for dichotomous experimental gene list from scDataset
#compare.list=c(read_concepts("PathwayDb/h.all.v7.5.1.symbols.gmt"),read_concepts("PathwayDb/c2.cp.v7.5.1.symbols.gmt"),read_concepts("PathwayDb/c2.cgp.v2022.1.Hs.symbols.gmt"))
target.list<-read.table ("scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1

#perform deep functional interpretation of the target gene list and calculate uniConSig scores.The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,rm.overfit=F)
#The CSEA2 function took 10 mins when testing
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways

#disambiguate top enriched pathways.
#This step took 20 mins when testing
topn=100 #specify the number of top pathways to disambiguate
disambiguate<-disambiguation.CSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result))),
                                  p.cut=0.01,p.adjust.method="bonferroni")#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"

#compute functional associations between selected top pathways
selectn=30 #specify the number of top pathways to compute associations
assoc <- pathwayAssociation(topPathway=disambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(disambiguate[[1]])))],compare.list,feature.list,preCalmatrix)

#draw heatmaps showing the functional associations between selected top pathways
pdf(file="CSEA-PathwayAssocHeatmapNetwork_20220627.pdf",width=10, height=10)
pathway.heatmap(matrixData=assoc,clustering = TRUE,fontSize=8)

#draw network figure for top pathways. NES.cut determines the levels of significance for the pathway similarity to be shown as edges. The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=disambiguate[[1]],assoc=assoc,NES.cut=2,node.size=1,line.thickness=1,font.size=0.4,wrap.string.width=15)
graphics.off()

####################################################################################################################
##Option III. Perform GSEA for pathway analysis
####################################################################################################################
#This method is more appliable to pathway analysis of bulk gene expression data.
compare.list=c(read_concepts("./PathwayDb/h.all.v7.5.1.symbols.gmt"),read_concepts("./PathwayDb/c2.cp.v7.5.1.symbols.gmt"),read_concepts("./PathwayDb/c2.cgp.v2022.1.Hs.symbols.gmt"))

#perform limma DGE analysis for single cell gene expression data.
#provide .gct file that contains single cell gene expression data and .cls file that defines the cell groups, as in the example.
#please refer to: https://www.genepattern.org/file-formats-guide
limma=limmaDGE(gctFile=gctFile,clsFile=clsFile)
weight=setNames(limma$Signed.Q.Value,row.names(limma))#signed q values will be used as weights for GSEA analysis

# perform pathway enrichment analysis (!Please make sure to input gene names with standard HGNC symbols)
up.GSEA.result<-CSEA2(target.score=weight,compare.list,p.cut=0.05)
down.GSEA.result<-CSEA2(target.score=-weight,compare.list,p.cut=0.05)

#disambiguate up and downregulated pathways
topn=30 #specify the number of top pathways to disambiguate
up.disambiguate<-disambiguation.GSEA(GEA.result=up.GSEA.result,weight=weight,compare.list=compare.list,topn=min(c(topn,nrow(up.GSEA.result))),
                                     p.cut=0.05,p.adjust.method="bonferroni")#p.adjust.method-- chose one: NULL or "bonferroni" or "BH"
down.disambiguate<-disambiguation.GSEA(GEA.result=down.GSEA.result,weight=-weight,compare.list=compare.list,topn=min(c(topn,nrow(down.GSEA.result))),
                                       p.cut=0.05,p.adjust.method="bonferroni")

