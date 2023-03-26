#' @title Pre-compiled molecular concepts
#'
#' @description A data set of pre-compiled molecular concepts
#'
#' @format A list format. List names are concepts and elements are gene symbols
#' \describe{
#'   \item{KEGG_CITRATE_CYCLE_TCA_CYCLE}{Gene symbols in the concept of KEGG_CITRATE_CYCLE_TCA_CYCLE}
#'   \item{KEGG_PENTOSE_PHOSPHATE_PATHWAY}{Gene symbols in the concept of KEGG_PENTOSE_PHOSPHATE_PATHWAY}
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"feature.list"


#' @title Gene signature data obtained from MSigDB
#'
#' @description Gene signature data provided by Molecular Signature Database (MSigDB)
#'
#' @format A list format. List names are gene signatures and elements are gene symbols
#' \describe{
#'   \item{KEGG_CITRATE_CYCLE_TCA_CYCLE}{Gene symbols in the concept of KEGG_CITRATE_CYCLE_TCA_CYCLE}
#'   \item{KEGG_PENTOSE_PHOSPHATE_PATHWAY}{Gene symbols in the concept of KEGG_PENTOSE_PHOSPHATE_PATHWAY}
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"compare.list"


#' @title Precalculated matrix to penalize the redundancies in molecular concepts
#'
#' @description Precalculated matrix to penalize the redundancies in molecular concepts
#'
#' @format A gmt file format.
#' \describe{
#'   \item{ACSS2}{KEGG_GLYCOLYSIS_GLUCONEOGENESIS@6.40917382433143\tKEGG_PYRUVATE_METABOLISM@6.01715179470204}
#'   \item{GCK}{KEGG_GLYCOLYSIS_GLUCONEOGENESIS@16.9306840571511\tKEGG_GALACTOSE_METABOLISM@13.6040224314813}
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"preCalmatrix"


#' @title Gene expression data in gct file format
#'
#' @description Gene expression data in gct file format
#'
#' @format Data.frame of gene expression data. The first two rows were removed.
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"ExpData"


#' @title Class data in cls file format
#'
#' @description Class data in cls file format
#'
#' @format Data.frame of class data. The first two rows were removed.
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"ClassData"

#' @title Gene symbol and gene length data
#'
#' @description Gene length data for 11,506 genes
#'
#' @format Data.frame of gene length data
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"GenelengthData"


#' @title Transcript length to number of fragments (for the nonUMI protocol)
#'
#' @description From transcript length to number of fragments (for the nonUMI protocol)
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"len2nfrag"


