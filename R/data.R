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
#'   \item{ACSS2}{KEGG_GLYCOLYSIS_GLUCONEOGENESIS}{KEGG_PYRUVATE_METABOLISM}
#'   \item{GCK}{KEGG_GLYCOLYSIS_GLUCONEOGENESIS}{KEGG_GALACTOSE_METABOLISM}
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"preCalmatrix"


#' @title gct file
#'
#' @description gct file
#'
#' @format txt file
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"gctFile"


#' @title cls file
#'
#' @description cls file
#'
#' @format txt file
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"clsFile"


#' @title Molecular pathway or gene signature file
#'
#' @description Molecular pathway or gene signature file to make compare.list
#'
#' @format file name
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"GeneSigFile"



#' @title Weight of genes: Signed.Q.Value
#'
#' @description Limma or SCDE can calculate the weight of genes.
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"weight"


#' @title LimmaFit data
#'
#' @description LimmaFit data
#'
#' @format limma
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"LimmaFit"


#' @title LimmaOut
#'
#' @description Limma Output
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"LimmaOut"


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


#' @title uniConSigResult example
#'
#' @description Calculate uniConSig scores for all genes
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"uniConSigResult"


#' @title targetScore example
#'
#' @description targetScore obtained by uniConSig scores for all genes
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"targetScore"


#' @title tmp.weight example
#'
#' @description tmp.weight is the same as targetScore which is obtained by uniConSig scores for all genes
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"tmp.weight"


#' @title targetList example
#'
#' @description targetList obtained by Limma or SCDE.
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"targetList"


#' @title up.CSEA.result example
#'
#' @description up.CSEA.result is obtained by CSEA2() when you input up.uniConSig scores
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"up.CSEA.result"


#' @title down.CSEA.result example
#'
#' @description down.CSEA.result is obtained by CSEA2() when you input down.uniConSig scores
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"down.CSEA.result"


#' @title up.GSEA.result example
#'
#' @description up.GSEA.result is obtained by CSEA2() when you input gene weight positive
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"up.GSEA.result"


#' @title down.CSEA.result example
#'
#' @description down.GSEA.result is obtained by CSEA2() when input gene weight negative
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"down.GSEA.result"


#' @title genelist example to draw interactome network for a selected pathway and your target gene of interests.
#'
#' @description You can analyze the interactions of your gene of interests with the pathways genes.
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"genelist"


#' @title Disambiguation of CSEA2 result
#'
#' @description Disambiguate upregulated pathways
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"up.disambiguate"


#' @title Disambiguation of CSEA2 result
#'
#' @description Disambiguate downregulated pathways
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"down.disambiguate"


#' @title Calculate associations between top pathways.
#'
#' @description Calculate associations between selected top up-regulated pathways.
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"up.assoc"


#' @title Calculate associations between top pathways.
#'
#' @description Calculate associations between selected top down-regulated pathways.
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"down.assoc"


#' @title Calculate associations between merged pathways of up- and down-
#'
#' @description Calculate associations between merged up- and down-regulated pathways.
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"pathwayMergeAssoc"


#' @title Merged pathways of top up- and down-regulated
#'
#' @description Merged pathways of top up- and down-regulated
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"pathway.merge"


#' @title Gene interaction pair data.
#'
#' @description This gene data provide information of gene interactions to make a network.
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"HuRI.interactions"


#' @title Expression count matrix example. You can make GCT from this matrix.
#'
#' @description Expression count matrix  obtained from True2ObservedCounts.
#'
#' @format matrix
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"ExpMatrix"


#' @title KS test to calculate enrichment of feature list
#'
#' @description The weight of genes are used to compute functional relations of feature list
#'
#' @format list
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"ks.result"


#' @title Genes of Signed.Q.Value > 0 from Limma or SCDE out
#'
#' @description Limma or SCDE can calculate the weight of genes. The genes of positive weight will be used to calculate uniConSig score
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"targetList_Up"


#' @title Genes of Signed.Q.Value < 0 from Limma or SCDE out
#'
#' @description Limma or SCDE can calculate the weight of genes. The genes of negative weight will be used to calculate uniConSig score
#'
#' @format vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"targetList_Dw"



#' @title Normalized ExpressionData of scRNA-seq.
#'
#' @description ExpressionData can be obtained from gct file.
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"ExpressionData"


#' @title Processed true_counts ExpData of scRNA-seq.
#'
#' @description ExpressionData can be obtained from gct file. This is before simulation.
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"true_counts"


#' @title The length data for the genes in the ExpressionData
#'
#' @description This gene length data is from SymSim
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"gene_len"


#' @title Top pathway names from disambiguated CSEA2 out
#'
#' @description Top pathway names of up-regulated.
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"topPathwayUp"


#' @title Top pathway names from disambiguated CSEA2 out
#'
#' @description Top pathway names of down-regulated.
#'
#' @format a vector
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"topPathwayDw"


#' @title Permuted calculation of enrichment scores
#'
#' @description Permuted calculation of enrichment scores
#'
#' @format List
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"myPermu"


#' @title Top 30 pathways from up.disambituated pathways.
#'
#' @description Top 30 up.pathways.
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"uppathway"


#' @title Top 30 pathways from down.disambituated pathways.
#'
#' @description Top 30 down.pathways.
#'
#' @format data.frame
#' \describe{
#' }
#' @source <https://github.com/wangxlab/IndepthPathway>
"downpathway"
