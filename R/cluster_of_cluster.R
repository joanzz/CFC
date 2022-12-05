#' Cluster of cluster
#'
#' ClusterofCluster is a multi-method based on the R package "ConsensusClusterPlus".
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @param mRNAexp A matrix of mRNA expression
#' @param miRNAexp A matrix of miRNA expression
#' @param lncRNAexp A matrix of lncRNA expression
#' @param methylation A matrix of methylation
#' @param CNVmatrix A matrix of Copy number variation
#' @param maxK integer value. maximum cluster number  for Consensus Clustering Algorithm to evaluate.
#' @param reps  integer value. number of subsamples(in other words, The iteration number of each cluster number)
#' @param clusterAlg character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
#' @param distance character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param pItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param pFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param plot Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param innerLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param seed optional numerical value. sets random seed for reproducible results.
#'
#' @return A list of length maxK. Each element is a list containing consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class asssignments).
#'
#' @seealso \code{ConsensusClusterPlus}
#'
#' @export
#'
#' @examples
#'data(mRNAexp)
#'data(miRNAexp)
#'data(lncRNAexp)
#'data(methylation)
#'data(CNVmatrix)
#'mRNAexp=data.filter(mRNAexp,percentage=0.6)
#'miRNAexp=data.filter(miRNAexp,percentage=0.6)
#'lncRNAexp=data.filter(lncRNAexp,percentage=0.6)
#'mRNAexp[mRNAexp==0]<-NA
#'miRNAexp[miRNAexp==0]<-NA
#'lncRNAexp[lncRNAexp==0]<-NA
#'mRNAexp<-as.matrix(mRNAexp)
#'miRNAexp<-as.matrix(miRNAexp)
#'lncRNAexp<-as.matrix(lncRNAexp)
#'mRNAexp=data.imputation(mRNAexp,fun="knn")
#'miRNAexp=data.imputation(miRNAexp,fun="knn")
#'lncRNAexp=data.imputation(lncRNAexp,fun="knn")
#'results = Cluster_of_cluster(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,maxK=6,reps=50,pItem=0.8,pFeature=1,clusterAlg="hc",distance="pearson",innerLinkage="complete",seed=1262118388.71279,plot="pdf")
#'
Cluster_of_cluster <- function(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,
                               maxK=6,reps=50,pItem=0.8,pFeature=1,# notice which maxK you set
                               clusterAlg="hc", ### "pam","hc","km"
                               distance="pearson",
                               innerLinkage="complete",#implement consensus clustering with innerLinkage="complete".
                               seed=1262118388.71279,
                               plot="pdf"){
  #mRNA
  rt_mRNA <- mRNAexp
  title_mRNA=tempdir()
  title_mRNA="mRNA_cluster"
  rt_nor_mRNA <- as.matrix(rt_mRNA)
  results_mRNA = ConsensusClusterPlus(rt_nor_mRNA,
                                      maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                      clusterAlg=clusterAlg, ### "pam","hc","km"
                                      distance=distance,
                                      innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                      seed=seed,
                                      plot=plot)   # plot="png" will keep figure separately
  maxK_mRNA = maxK
  Kvec_mRNA = 2:maxK_mRNA
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC_mRNA = rep(NA,length(Kvec_mRNA))
  names(PAC_mRNA) = paste("K=",Kvec_mRNA,sep="") # from 2 to maxK
  for(i in Kvec_mRNA){
    M = results_mRNA[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC_mRNA[i-1] = Fn(x2) - Fn(x1)
  }
  optK_mRNA = Kvec_mRNA[which.min(PAC_mRNA)]
  ####
  ####
  #miRNA
  rt_miRNA <- miRNAexp
  title_miRNA=tempdir()
  title_miRNA="miRNA"
  rt_nor_miRNA <- as.matrix(rt_miRNA)
  results_miRNA = ConsensusClusterPlus(rt_nor_miRNA,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                       clusterAlg=clusterAlg, ### "pam","hc","km"
                                       distance=distance,
                                       innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                       seed=seed,
                                       plot=plot)   # plot="png" will keep figure separately
  maxK_miRNA = maxK  #
  Kvec_miRNA = 2:maxK_miRNA
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC_miRNA = rep(NA,length(Kvec_miRNA))
  names(PAC_miRNA) = paste("K=",Kvec_miRNA,sep="") # from 2 to maxK
  for(i in Kvec_miRNA){
    M = results_miRNA[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])# 计算累积概率密度， 概率密度累积之和等于1
    PAC_miRNA[i-1] = Fn(x2) - Fn(x1)#
  }
  optK_miRNA = Kvec_miRNA[which.min(PAC_miRNA)]
  ##lncRNA
  rt_lncRNA <- lncRNAexp
  title_lncRNA=tempdir()
  title_lncRNA="lncRNA"
  rt_nor_lncRNA <- as.matrix(rt_lncRNA)
  results_lncRNA = ConsensusClusterPlus(rt_nor_lncRNA,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                        clusterAlg=clusterAlg, ### "pam","hc","km"
                                        distance=distance,
                                        innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                        seed=seed,
                                        plot=plot)   # plot="png" will keep figure separately
  ### chunk 3   PAC implementation
  maxK_lncRNA = maxK  #
  Kvec_lncRNA = 2:maxK_lncRNA
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC_lncRNA = rep(NA,length(Kvec_lncRNA))
  names(PAC_lncRNA) = paste("K=",Kvec_lncRNA,sep="") # from 2 to maxK
  for(i in Kvec_lncRNA){
    M = results_lncRNA[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC_lncRNA[i-1] = Fn(x2) - Fn(x1)
  }
  #The optimal K
  optK_lncRNA = Kvec_lncRNA[which.min(PAC_lncRNA)]
  ###
  ###methylation
  rt_met <- methylation
  title_met=tempdir()
  title_met="methylation"
  rt_nor_met <- as.matrix(rt_met)
  results_met = ConsensusClusterPlus(rt_nor_met,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                     clusterAlg=clusterAlg, ### "pam","hc","km"
                                     distance=distance,
                                     innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                     seed=seed,
                                     plot=plot)   # plot="png" will keep figure separately

  ### chunk 3   PAC implementation
  maxK_met = maxK  #
  Kvec_met = 2:maxK_met
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC_met = rep(NA,length(Kvec_met))
  names(PAC_met) = paste("K=",Kvec_met,sep="") # from 2 to maxK
  for(i in Kvec_met){
    M = results_met[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC_met[i-1] = Fn(x2) - Fn(x1)
  }
  #The optimal K
  optK_met = Kvec_met[which.min(PAC_met)]
  ###cnv
  rt_cnv <- CNVmatrix
  title_cnv=tempdir()
  title_cnv="CNV"
  rt_nor_cnv <- as.matrix(rt_cnv)
  results_cnv = ConsensusClusterPlus(rt_nor_cnv,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                     clusterAlg=clusterAlg, ### "pam","hc","km"
                                     distance=distance,
                                     innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                     seed=seed,
                                     plot=plot)   # plot="png" will keep figure separately

  ### chunk 3   PAC implementation
  maxK_cnv = maxK  #
  Kvec_cnv = 2:maxK_cnv
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC_cnv = rep(NA,length(Kvec_cnv))
  names(PAC_cnv) = paste("K=",Kvec_cnv,sep="") # from 2 to maxK
  for(i in Kvec_cnv){
    M = results_cnv[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC_cnv[i-1] = Fn(x2) - Fn(x1)
  }
  #The optimal K
  optK_cnv = Kvec_cnv[which.min(PAC_cnv)]
  ####
  opt_mRNA=results_mRNA[[optK_mRNA]]["consensusClass"]
  opt_miRNA=results_miRNA[[optK_miRNA]]["consensusClass"]
  opt_lncRNA=results_lncRNA[[optK_lncRNA]]["consensusClass"]
  opt_met=results_met[[optK_met]]["consensusClass"]
  opt_cnv=results_cnv[[optK_cnv]]["consensusClass"]

  #mRNA
  matrix_mRNA=as.data.frame(opt_mRNA)
  matrix_mRNA$rowname=substr(rownames(matrix_mRNA),1,15)
  matrix_mRNA$rowname=gsub("[.]","-",matrix_mRNA$rowname)
  matrix_mRNA <- matrix_mRNA[!duplicated(matrix_mRNA$rowname),]
  rownames(matrix_mRNA)=matrix_mRNA$rowname
  cluster_mRNA <- as.data.frame(matrix(0,length(row.names(matrix_mRNA)),optK_mRNA))
  rownames(cluster_mRNA)<-rownames(matrix_mRNA)
  for (i in 1:length(row.names(cluster_mRNA))) {
    cluster_mRNA[i,matrix_mRNA[i,1]]=1
  }
  c_mRNA = 1:optK_mRNA
  colnames(cluster_mRNA) = paste("mClass",c_mRNA,sep="")

  #miRNA
  matrix_miRNA=as.data.frame(opt_miRNA)
  matrix_miRNA$rowname=substr(rownames(matrix_miRNA),1,15)
  matrix_miRNA$rowname=gsub("[.]","-",matrix_miRNA$rowname)
  matrix_miRNA <- matrix_miRNA[!duplicated(matrix_miRNA$rowname),]
  rownames(matrix_miRNA)=matrix_miRNA$rowname
  cluster_miRNA <- as.data.frame(matrix(0,length(row.names(matrix_miRNA)),optK_miRNA))
  rownames(cluster_miRNA)<-rownames(matrix_miRNA)
  for (i in 1:length(row.names(cluster_miRNA))) {
    cluster_miRNA[i,matrix_miRNA[i,1]]=1
  }
  c_miRNA = 1:optK_miRNA
  colnames(cluster_miRNA) = paste("miClass",c_miRNA,sep="")

  #lncRNA
  matrix_lncRNA=as.data.frame(opt_lncRNA)
  matrix_lncRNA$rowname=substr(rownames(matrix_lncRNA),1,15)
  matrix_lncRNA$rowname=gsub("[.]","-",matrix_lncRNA$rowname)
  matrix_lncRNA <- matrix_lncRNA[!duplicated(matrix_lncRNA$rowname),]
  rownames(matrix_lncRNA)=matrix_lncRNA$rowname
  cluster_lncRNA <- as.data.frame(matrix(0,length(row.names(matrix_lncRNA)),optK_lncRNA))
  rownames(cluster_lncRNA)<-rownames(matrix_lncRNA)
  for (i in 1:length(row.names(cluster_lncRNA))) {
    cluster_lncRNA[i,matrix_lncRNA[i,1]]=1
  }
  c_lncRNA = 1:optK_lncRNA
  colnames(cluster_lncRNA) = paste("lncClass",c_lncRNA,sep="")
  #methylation
  matrix_met=as.data.frame(opt_met)
  matrix_met$rowname=substr(rownames(matrix_met),1,15)
  matrix_met$rowname=gsub("[.]","-",matrix_met$rowname)
  matrix_met <- matrix_met[!duplicated(matrix_met$rowname),]
  rownames(matrix_met)=matrix_met$rowname
  cluster_met <- as.data.frame(matrix(0,length(row.names(matrix_met)),optK_met))
  rownames(cluster_met)<-rownames(matrix_met)
  for (i in 1:length(row.names(cluster_met))) {
    cluster_met[i,matrix_met[i,1]]=1
  }
  c_met = 1:optK_met
  colnames(cluster_met) = paste("methylation",c_met,sep="")
  #cnv
  matrix_cnv=as.data.frame(opt_cnv)
  matrix_cnv$rowname=substr(rownames(matrix_cnv),1,15)
  matrix_cnv$rowname=gsub("[.]","-",matrix_cnv$rowname)
  matrix_cnv <- matrix_cnv[!duplicated(matrix_cnv$rowname),]
  rownames(matrix_cnv)=matrix_cnv$rowname
  cluster_cnv <- as.data.frame(matrix(0,length(row.names(matrix_cnv)),optK_cnv))
  rownames(cluster_cnv)<-rownames(matrix_cnv)
  for (i in 1:length(row.names(cluster_cnv))) {
    cluster_cnv[i,matrix_met[i,1]]=1
  }
  c_cnv = 1:optK_cnv
  colnames(cluster_cnv) = paste("methylation",c_cnv,sep="")

  ###
  cluster_mRNA$rownames=row.names(cluster_mRNA)
  cluster_miRNA$rownames=row.names(cluster_miRNA)
  cluster_lncRNA$rownames=row.names(cluster_lncRNA)
  cluster_met$rownames=row.names(cluster_met)
  cluster_cnv$rownames=row.names(cluster_cnv)

  m1 <- merge(cluster_mRNA, cluster_miRNA, by.x = "rownames", by.y = "rownames")
  m2<- merge(cluster_lncRNA, m1, by.x = "rownames", by.y = "rownames")
  m3<-merge(cluster_met, m2, by.x = "rownames", by.y = "rownames")
  m4<-merge(cluster_cnv, m3, by.x = "rownames", by.y = "rownames")
  row.names(m4)=m4$rownames
  m4=m4[,-1]
  ####
  rt_nor <- t(m4)
  title=tempdir()
  title="Cluster_of_Cluster"
  rt_nor <- as.matrix(rt_nor)
  results = ConsensusClusterPlus(rt_nor,maxK=maxK,reps=reps,pItem=pItem,pFeature=pFeature,# notice which maxK you set
                                 clusterAlg=clusterAlg, ### "pam","hc","km"
                                 distance=distance,
                                 innerLinkage=innerLinkage,#implement consensus clustering with innerLinkage="complete".
                                 seed=seed,
                                 plot=plot)   # plot="png" will keep figure separately
  return(results)

}
