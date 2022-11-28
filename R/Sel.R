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
