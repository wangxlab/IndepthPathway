#' Install and load necessary packages
#'
#' @return install and load packages
#' @export
#'
#' @examples LoadPackage()
LoadPackage <- function() {
        packages=c("stringr","qvalue","vegan","tidyr","moments","limma","dplyr","gplots",
                   "RColorBrewer","corrplot","pheatmap","igraph","otuSummary","pROC",
                   "matrixStats","pacman","fgsea","RColorBrewer","devtools","data.table","scde")
        if (!require("BiocManager")) install.packages("BiocManager")
        #Checking if the package belongs to CRAN or to Bioconductor and installing them accordingly.
        for(lib in packages){
          if(!lib %in% installed.packages()){
            if(lib %in% available.packages()[,1]){
              install.packages(lib,dependencies=TRUE)
            }else {(BiocManager::install(lib))
            }}
        }
        #Loading the libraries
        sapply(packages,require,character=TRUE)
}


