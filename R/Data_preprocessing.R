#' Data imputation
#'
#' Data imputation for features with missing values
#'@importFrom impute impute.knn
#'@param Data A matrix representing the genomic data such as gene expression data, miRNA expression data. \cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#'@param fun A character value representing the imputation type. The optional values are shown below:
#'\itemize{
#' \item "median". The NAs will be replaced by the median of the existing values of this feature in all samples.
#' \item "mean". The NAs will be replaced by the mean of the existing values of this feature in all samples.
#' \item "knn". It will apply the "impute" package to impute the missing values. This is a common way to process
#' the missing observation for MicroArray dataset(Default value).
#'
#'}
#'@return
#' The data matrix after imputation (without NAs).
#'@examples
#' Data=matrix(runif(1000),nrow = 50,ncol = 20)
#' geneName=paste("Gene", 1:50, sep = " ")
#' sampleName=paste("Sample", 1:20, sep = " ")
#' rownames(Data)=geneName
#' colnames(Data)=sampleName
#' index=sample(c(1:1000),60)
#' Data[index]=NA
#' result=data.imputation(Data,fun="knn")
#'
#' @references
#' Xu T, Le TD, Liu L, Su N, Wang R, Sun B, Colaprico A, Bontempi G, Li J. CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation and visualization. Bioinformatics. 2017 Oct 1;33(19):3131-3133. doi: 10.1093/bioinformatics/btx378. PMID: 28605519.
#'@export
#'
data.imputation<-function(Data,fun="knn")
{
  if(fun=="median")
  {
    result=apply(Data,1,function(x){
      x<-as.numeric(x)
      x[is.na(x)] =median(x, na.rm=TRUE)
      x
    })
    result=t(result)
  } else if(fun=="mean")
  {
    result=apply(Data,1,function(x){
      x<-as.numeric(x)
      x[is.na(x)] =mean(x, na.rm=TRUE)
      x
    })
    result=t(result)
  } else if(fun=="knn")
  {
    result=impute.knn(Data)$data
  }
  colnames(result)=colnames(Data)
  result
}


#' Data filter
#'
#' Data filter for features with missing values
#' @param Data A matrix representing the genomic data such as gene expression data, miRNA expression data. \cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#' @param percentage The maximum proportion in the sample with missing values
#'
#' @return
#' The data matrix after filter.
#' @export
#'
#' @examples
#' Data=matrix(runif(1000),nrow = 50,ncol = 20)
#' geneName=paste("Gene", 1:50, sep = " ")
#' sampleName=paste("Sample", 1:20, sep = " ")
#' rownames(Data)=geneName
#' colnames(Data)=sampleName
#' index=sample(c(1:1000),300)
#' Data[index]=NA
#' result=data.filter(Data,percentage=0.6)
data.filter=function(Data,percentage=0.8){
  Data[is.na(Data)]=0
  sum=as.integer(length(colnames(Data))*percentage)
  flag <- apply(Data, 1, function(x) sum(x==0) < sum)
  Data <- Data[which(flag),]
  Data
}

#'Data normalization
#'
#'Conduct normalization for dataset.
#'
#'@param Data A matrix representing the genomic data such as gene expression data, miRNA expression data. \cr
#' For the matrix, the rows represent the genomic features, and the columns represent the samples.
#'@param type A character value representing the normalization type. The optional values are shown below:
#'\itemize{
#' \item "feature_Median". The default value. Normalize dataset by sweeping the median values of each feature.
#' \item "feature_Mean". Normalize dataset by sweeping the mean values of each feature.
#' \item "feature_zscore". Conduct z_score normalization for each feature.
#' \item "sample_zscore".  Conduct z_score normalization for each samples.
#'
#'}
#'@param log2 A logical value. If TRUE, the data is transform as log2(x+1). This is commonly used for RNAseq data.
#'@return
#' The normalized data matrix.
#'@examples
#' data(mRNAexp)
#' result=data.normalization(mRNAexp,type="feature_Median",log2=FALSE)
#'
#'@references
#'Xu T, Le TD, Liu L, Su N, Wang R, Sun B, Colaprico A, Bontempi G, Li J. CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation and visualization. Bioinformatics. 2017 Oct 1;33(19):3131-3133. doi: 10.1093/bioinformatics/btx378. PMID: 28605519.
#'@export
#'
data.normalization<-function(Data,type="feature_Median",log2=FALSE)
{
  if(log2)
  {
    data=log2(Data+1)
  } else{
    data=Data
  }

  if(type=="feature_Median")
  {
    result=sweep(data,1,apply(data,1,function(x) median(x, na.rm = TRUE)))
  } else if(type=="feature_Mean")
  {
    result=sweep(data,1,apply(data,1,function(x) mean(x, na.rm = TRUE)))
  } else if(type=="feature_zscore")
  {
    var_row=apply(data,1,var)
    index=which(var_row<1e-10)
    if(length(index)>0)
    {
      data=data[-index,]
      cat("The features with the zero variance have been removed.")
    }
    result=t(scale(t(data)))
  } else if(type=="sample_zscore")
  {
    var_col=apply(data,2,var)
    index=which(var_col<1e-10)
    if(length(index)>0)
    {
      data=data[-index,]
      cat("The samples with the zero variance have been removed.")
    }
    result=scale(data)
  }
  result
}

#' Quantile Normalize By Feature
#'
#' A cross-platform normalization method using FSQN to improve comparability of DNA microarray and RNA-seq datasets
#' @importFrom preprocessCore normalize.quantiles.use.target
#' @param matrix_to_normalize A matrix (m x n) with m samples as rows, and n features as columns.
#' @param target_distribution_matrix matrix (m2 x n) with m2 samples as rows, and n features as columns to use as the target distribution.
#'
#' @return A normalized matrix
#' @export
#'
#'@references Franks JM, Cai G, Whitfield ML. Feature specific quantile normalization enables cross-platform classification of molecular subtypes using gene expression data. Bioinformatics. 2018 Jun 1;34(11):1868-1874. doi: 10.1093/bioinformatics/bty026. PMID: 29360996; PMCID: PMC5972664.
#'
#' @examples
#' set.seed(7)
#' target <- matrix(rnorm(100*150, mean = 1, sd = 1), nrow = 100, ncol = 150)
#' test <- matrix(rnorm(30*150, mean = 2, sd = 2), nrow = 30, ncol = 150)
#' normalized_test <- quantileNormalizeByFeature(test, target)
#'
quantileNormalizeByFeature <- function(matrix_to_normalize,
                                       target_distribution_matrix){

  if (ncol(matrix_to_normalize) != ncol(target_distribution_matrix)){
    cat("ERROR: Data matrices are not compatible - column lengths differ!")
  }
  else{

    data.qn <- matrix(0, nrow = nrow(matrix_to_normalize),
                      ncol = ncol(matrix_to_normalize))

    for (i in 1:ncol(matrix_to_normalize)){
      feature.to.normalize <- matrix_to_normalize[,i]
      target.feature.dist <- target_distribution_matrix[,i]
      result <- normalize.quantiles.use.target(
        x = as.matrix(feature.to.normalize),
        target = target.feature.dist,
        copy = TRUE)
      data.qn[,i] <- result
    }
    rownames(data.qn) = rownames(matrix_to_normalize)
    colnames(data.qn) = colnames(matrix_to_normalize)
    return(data.qn)
  }
}


