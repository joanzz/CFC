#' Dataset: mRNA expression
#'
#'A glioblastoma (GBM) protein coding gene expression dataset downloaded from TCGA.
#'
#'
#'\itemize{
#'  \item Rows are genes
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name mRNAexp
#'@examples
#' data(mRNAexp)
NULL


#' Dataset: miRNA expression
#'
#'A glioblastoma (GBM) miRNA expression dataset downloaded from TCGA.
#'
#'\itemize{
#'  \item Rows are miRNAs
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name miRNAexp
#'@examples
#' data(miRNAexp)
NULL

#' Dataset: lncRNA expression
#'
#'A glioblastoma (GBM) lncRNA expression dataset downloaded from TCGA.
#'
#'\itemize{
#'  \item Rows are lncRNAs
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name lncRNAexp
#'@examples
#' data(lncRNAexp)
NULL

#' Dataset: methylation
#'
#'A methylation dataset downloaded from TCGA.
#'
#'\itemize{
#'  \item Rows are features
#'  \item Columns are cancer samples
#'}
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name methylation
#'@examples
#' data(methylation)
NULL

#' Dataset: All metabolism geneset
#'
#'A metabolism dataset downloaded from the literature.
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name pathway_all
#'@examples
#' data(pathway_all)
#' @references Peng X, Chen Z, Farshidfar F, Xu X, Lorenzi PL, Wang Y, Cheng F, Tan L, Mojumdar K, Du D, Ge Z, Li J, Thomas GV, Birsoy K, Liu L, Zhang H, Zhao Z, Marchand C, Weinstein JN; Cancer Genome Atlas Research Network, Bathe OF, Liang H. Molecular Characterization and Clinical Relevance of Metabolic Expression Subtypes in Human Cancers. Cell Rep. 2018 Apr 3;23(1):255-269.e4. doi: 10.1016/j.celrep.2018.03.077. PMID: 29617665; PMCID: PMC5916795.
NULL

#' Dataset: The results of cluster
#'
#'A dataset from Clusterofcluster
#'
#'@docType data
#'@keywords datasets
#'@format A list
#'@name results
#'@examples
#' data(results)
NULL

#' Dataset: clinical data
#'
#'The clinical data downloaded from TCGA
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name meta
#'@examples
#' data(meta)
NULL

#' Dataset: Data PAC
#'
#'The data for plot PAC
#'
#'@docType data
#'@keywords datasets
#'@format A data matrix
#'@name dot_df
#'@examples
#' data(dot_df)
NULL

#' Dataset: exprData
#'
#'The data for drug sensitivity analysis
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name exprData
#'@examples
#' data(exprData)
NULL

#' Dataset: studyResponse
#'
#'The data for drug sensitivity analysis
#'
#'@docType data
#'@keywords datasets
#'@format A list
#'@name studyResponse
#'@examples
#' data(studyResponse)
NULL

#' Dataset: bortIndex
#'
#'The data for drug sensitivity analysis
#'
#'@docType data
#'@keywords datasets
#'@format A list
#'@name bortIndex
#'@examples
#' data(bortIndex)
NULL

#' Dataset: CNV
#'
#'Raw CNV data for bladder cancer downloaded in XENA
#'
#'@docType data
#'@keywords datasets
#'@format A dataset
#'@name CNV
#'@examples
#' data(CNV)
NULL

#' Dataset: CNVmatrix
#'
#'A matrix transformed by CNV data
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name CNVmatrix
#'@examples
#' data(CNVmatrix)
NULL

#' Dataset: entr_expr
#'
#'a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name entr_expr
#'@examples
#' data(entr_expr)
NULL

#' Dataset: entr_adja
#'
#'an adjacency matrix representing a connected PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name entr_adja
#'@examples
#' data(entr_adja)
NULL

#' Dataset: entr_labels
#'
#'Classification (label) of the sample
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name entr_labels
#'@examples
#' data(entr_labels)
NULL

#' Dataset: drugToCellLineDataCgp
#'
#'To map cell lines to .CEL file names
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name drugToCellLineDataCgp
#'@examples
#' data(drugToCellLineDataCgp)
NULL

#' Dataset: gdsc_brainarray_syms
#'
#'The gene expression data
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name gdsc_brainarray_syms
#'@examples
#' data(gdsc_brainarray_syms)
NULL

#' Dataset: drugSensitivityDataCgp
#'
#'The drug sensitivity data
#'
#'@docType data
#'@keywords datasets
#'@format A matrix
#'@name drugSensitivityDataCgp
#'@examples
#' data(drugSensitivityDataCgp)
NULL
