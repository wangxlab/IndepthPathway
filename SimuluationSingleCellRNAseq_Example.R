############# =============== ############# =============== ############# =============== ############# =============== ############# ===============
### Install required R packages for SymSim
############# =============== ############# =============== ############# =============== ############# =============== ############# ===============
# The required Bioconductor packages can be installed as follows in R:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IRanges")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("S4Vectors")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")

## To install SymSim
library("devtools")
devtools::install_github("YosefLab/SymSim")

############# =============== SymSim (Synthetic model of multiple variability factors for Simulation) example code ============== #############
############# =============== ############# =============== ############# =============== ############# ===============
library(IndepthPathway)
library(dplyr)
library(data.table)

setwd("YourWorkingDirectory")

# data(ExpData) # gct data
# data(GenelengthData)
# data(len2nfrag)
############# =============== ############# =============== ############# =============== ############# ===============

## Check ExpData
ExpData[1:3,1:5]
# HGNC.symbol    MGI.symbol D11_HSC B11_HSC C11_HSC
# 1     TRAPPC2 0610009B22Rik       0       0       0
# 2    KIAA1841 0610010F05Rik       0       0       0
# 3    C17orf49 0610010K14Rik       0       0       0

## Delete the "MGI.symbol" column and change 'HGNC.symbol' column to rownames.  Duplicated rownames are not allowed, so I removed duplicated rows in advance.
ExpDataProc <- ExpData %>% dplyr::select(-MGI.symbol) %>% dplyr::filter(!duplicated(HGNC.symbol)) %>% tibble::column_to_rownames("HGNC.symbol"); ExpDataProc[1:3,1:3]
#           D11_HSC B11_HSC C11_HSC
# TRAPPC2        0       0       0
# KIAA1841       0       0       0
# C17orf49       0       0       0
ExpDataProc[] <- lapply(ExpDataProc, as.integer)

GenelengthDataMatch <- GenelengthData$V2[match(rownames(ExpDataProc), GenelengthData$V1)]; length(GenelengthDataMatch)  # 11506
GenelengthDataMatchProc <-as.numeric(GenelengthDataMatch[!is.na(GenelengthDataMatch)]); length(GenelengthDataMatchProc)  # 11427
ExpDataProc_RmvNA <-  ExpDataProc[!is.na(GenelengthDataMatch), ]; dim(ExpDataProc_RmvNA)  # 11427 177


protocol<-"nonUMI"; MyAlpha <- 0.04;  MyP <- 18; permutationNumb <- 2
for (permu in seq(1:permutationNumb) ) {
  Observed_counts <- True2ObservedCounts(true_counts=ExpDataProc_RmvNA, protocol=protocol, alpha_mean=MyAlpha,
                                         alpha_sd=0.02, gene_len=GenelengthDataMatchProc, nPCR1= MyP, depth_mean=5e4, depth_sd=3e3)
  Observed_count.matrix<-Observed_counts[[1]]
  rownames(Observed_count.matrix)<-rownames(ExpDataProc_RmvNA)
  colnames(Observed_count.matrix)<-colnames(ExpDataProc_RmvNA)

  Observed.FeatureCount.gct <- MakeGCT(Observed_count.matrix)

  GctOutFile<-paste0("Simulated_",protocol,"_a",MyAlpha,"p", MyP,"_permu",permu,".gct")
  fwrite(Observed.FeatureCount.gct, file=GctOutFile,col.names=F,row.names=F,sep="\t",quote=F)
}



