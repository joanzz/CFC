#' Select Metabolisms
#'
#' Select the type of metabolism you are interested in for downstream data analysis.
#'
#' @param Data A matrix representing the genomic data such as gene expression data.
#' @param type Type of metabolism of interest(Such as Lipid Metabolism,Carbohydrate Metabolism,Amino_acid Metabolism)
#'
#' @return A set of genes containing specific metabolic pathways
#' @export
#'
#' @examples
#' data(mRNAexp)
#' data_lipid=metabolism(mRNAexp,type="Lipid")
#'

metabolism <- function(Data,type="Lipid"){
  if(type=="Lipid"){
    Data=Data[row.names(Data) %in% pathway_all[,"Lipid"],]
  }else if(type=="Amino_acid"){
    Data=Data[row.names(Data) %in% pathway_all[,"Amino_acid"],]
  }else if(type=="Carbohydrate"){
    Data=Data[row.names(Data) %in% pathway_all[,"Carbohydrate"],]
  }else if(type=="Energy"){
    Data=Data[row.names(Data) %in% pathway_all[,"Energy"],]
  }else if(type=="Nucleotide"){
    Data=Data[row.names(Data) %in% pathway_all[,"Nucleotide"],]
  }else if(type=="TCA_cycle"){
    Data=Data[row.names(Data) %in% pathway_all[,"TCA_cycle"],]
  }else if(type=="Vitamin_cofactor"){
    Data=Data[row.names(Data) %in% pathway_all[,"Vitamin_cofactor"],]
  }
  Data
}

#' SbyPCA
#'
#' This function is based on the prcomp(), we write a shell for it and make it easy to use on genomic data.
#'
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param PC_percent A numeric values in [0,1]  representing the ratio  of principal component is seclected.
#' @param scale A bool variable, If true, the Data is normalized before PCA.
#'
#' @return
#' A new matrix with full or part Principal Component in new projection space.
#' @references
#'  Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'
#'
#' @examples
#' data(mRNAexp)
#' data1=SbyPCA(mRNAexp, PC_percent=0.9,scale = TRUE)
#'
#' @export
#'
SbyPCA<-function(Data,PC_percent=1,scale = TRUE)
{
  newData=t(Data)
  aa=prcomp(newData,scale = scale)
  vars <- apply(aa$x, 2, var)
  props <- vars / sum(vars)
  xx=as.vector(cumsum(props))
  num=which(xx>PC_percent)[1]
  coeff=aa$ro
  score=(newData%*%coeff)[,1:num]
  result=t(score)
  result
}

#' SbyVar
#'
#' Biological feature (such as gene) selection based on the most variance.
#'
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param cut.type A character value representing the selection type. The optional values are shown below:
#' \itemize{
#' \item "topk"
#' \item "cutoff"
#' }
#' @param value A numeric value.\cr
#' If the cut.type="topk", the top number of value features are selected.\cr
#' If the cut.type="cutoff", the features with (var>value) are selected.
#' @return An extracted subset data matrix with most variance features from the input data matrix.
#' @references
#'  Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'
#'
#' @examples
#' data(mRNAexp)
#' data1=SbyVar(mRNAexp, cut.type="topk",value=1000)
#'
#' @export
SbyVar<-function(Data, cut.type="topk", value)
{
  vars=apply(Data,1,var)
  feature_num=length(vars)
  hist(vars, breaks=feature_num*0.1, col="red",
       main="Expression (Variance) distribution",
       xlab="The Variance of feature")
  if(cut.type=="topk")
  {
    index= sort(vars,decreasing = TRUE,index.return=TRUE)
    if(value>nrow(Data))
    {
      value=nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff=index$x[value]
    abline(v=cutoff,col = "blue",lty = 5,lwd=1.5)
    index=index$ix[1:value]
    selectData=Data[index,]
  }
  if(cut.type=="cutoff")
  {
    abline(v=value,col = "blue",lty = 5,lwd=1.5)
    index=which(vars>value)
    selectData=Data[index,]
  }
  selectData
}

#' SbyMAD
#'
#'Biological feature (such as gene) selection based on the MAD.
#'
#' @param Data A data matrix representing the genomic data measured in a set of samples.
#' For the matrix, the rows represent the genomic features, and the columns represents the samples.
#' @param cut.type A character value representing the selection type. The optional values are shown below:
#' \itemize{
#' \item "topk"
#' \item "cutoff"
#' }
#' @param value A numeric value.\cr
#' If the cut.type="topk", the top number of value features are selected.\cr
#' If the cut.type="cutoff", the features with (MAD>value) are selected.
#' @return An extracted subset data matrix with the most variant MAD features from the input data matrix.
#' @references
#'  Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'
#'
#' @examples
#' data(mRNAexp)
#' data1=SbyMAD(mRNAexp, cut.type="topk",value=1000)
#'
#' @export
SbyMAD<-function(Data, cut.type="topk", value)
{
  mads=apply(Data,1,mad)
  feature_num=length(mads)
  hist(mads, breaks=feature_num*0.1, col="red",
       main="Expression (MAD) distribution",
       xlab="The MAD of feature")
  if(cut.type=="topk")
  {
    index= sort(mads,decreasing = TRUE,index.return=TRUE)
    if(value>nrow(Data))
    {
      value=nrow(Data)
      cat("Warning: the feature selection number is beyond the original feature numnber")
    }
    cutoff=index$x[value]
    abline(v=cutoff,col = "blue",lty = 5,lwd=1.5)
    index=index$ix[1:value]
    selectData=Data[index,]
  }
  if(cut.type=="cutoff")
  {
    abline(v=value,col = "blue",lty = 5,lwd=1.5)
    index=which(mads>value)
    selectData=Data[index,]
  }
  selectData


}







