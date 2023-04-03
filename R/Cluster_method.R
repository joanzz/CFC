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
#'methylation=data.filter(methylation,percentage=0.6)
#'CNVmatrix=data.filter(CNVmatrix,percentage=0.6)
#'mRNAexp[mRNAexp==0]<-NA
#'miRNAexp[miRNAexp==0]<-NA
#'lncRNAexp[lncRNAexp==0]<-NA
#'methylation[methylation==0]<-NA
#'CNVmatrix[CNVmatrix==0]<-NA
#'mRNAexp<-as.matrix(mRNAexp)
#'miRNAexp<-as.matrix(miRNAexp)
#'lncRNAexp<-as.matrix(lncRNAexp)
#'methylation<-as.matrix(methylation)
#'CNVmatrix<-as.matrix(CNVmatrix)
#'mRNAexp=data.imputation(mRNAexp,fun="mean")
#'miRNAexp=data.imputation(miRNAexp,fun="mean")
#'lncRNAexp=data.imputation(lncRNAexp,fun="mean")
#'methylation=data.imputation(methylation,fun="mean")
#'CNVmatrix=data.imputation(CNVmatrix,fun="mean")
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

#' Execute SNF(Similarity Network Fusion ) based on five datasets.
#'
#  SNF is a multi-omics data processing method that constructs a fusion patient similarity network.
#'
#' Based on the existing SNF algorithm, we improve it, increase the dataset type, and add the sample preprocessing function, because different sources of datasets, the sample format is not the same, researchers usually spend a lot of time in the unified sample format part, using our improved clustering algorithm,
#' users do not need to unify the sample format, directly input into the clustering algorithm, greatly saving the user's time.
#'
#' @param mRNAexp A matrix of mRNA expression
#' @param miRNAexp A matrix of miRNA expression
#' @param lncRNAexp A matrix of lncRNA expression
#' @param methylation A matrix of methylation
#' @param CNVmatrix A matrix of Copy number variation
#' @param clusterNum A integer representing the return cluster number
#' @param K Number of nearest neighbors
#' @param alpha Variance for local model
#' @param t Number of iterations for the diffusion process
#' @param plot Logical value. If true, draw the heatmap for the distance matrix with samples ordered to form clusters.
#'
#' @return A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'
#'  \item \strong{originalResult} : The clustering result of the original SNF algorithm"
#'
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  }
#' @references
#' B Wang, A Mezlini, F Demir, M Fiume, T Zu, M Brudno, B Haibe-Kains, A Goldenberg (2014) Similarity Network Fusion: a fast and effective method to aggregate multiple data types on a genome wide scale. Nature Methods. Online. Jan 26, 2014
#' Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'
#' @examples
#' data(mRNAexp)
#' data(miRNAexp)
#' data(lncRNAexp)
#' data(methylation)
#' data(CNVmatrix)
#' result=ExecuteSNF(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix, clusterNum=3, K=20, alpha=0.5, t=20)
#' result$group
#'
ExecuteSNF<-function(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,clusterNum, K=20, alpha=0.5, t=20,plot=TRUE)
{
  colnames(mRNAexp)=substr(colnames(mRNAexp),1,15)
  colnames(mRNAexp)=gsub("[.]","-",colnames(mRNAexp))
  colnames(miRNAexp)=substr(colnames(miRNAexp),1,15)
  colnames(miRNAexp)=gsub("[.]","-",colnames(miRNAexp))
  colnames(lncRNAexp)=substr(colnames(lncRNAexp),1,15)
  colnames(lncRNAexp)=gsub("[.]","-",colnames(lncRNAexp))
  colnames(methylation)=substr(colnames(methylation),1,15)
  colnames(methylation)=gsub("[.]","-",colnames(methylation))
  colnames(CNVmatrix)=substr(colnames(CNVmatrix),1,15)
  colnames(CNVmatrix)=gsub("[.]","-",colnames(CNVmatrix))
  a=intersect(x=colnames(miRNAexp),y=colnames(mRNAexp))
  b=intersect(x=colnames(lncRNAexp),y=a)
  c=intersect(x=colnames(methylation),y=b)
  d=intersect(x=colnames(CNVmatrix),y=c)
  mRNAexp1=mRNAexp[,match(d,colnames(mRNAexp))]
  miRNAexp1=miRNAexp[,match(d,colnames(miRNAexp))]
  lncRNAexp1=lncRNAexp[,match(d,colnames(lncRNAexp))]
  methylation1=methylation[,match(d,colnames(methylation))]
  CNVmatrix1=CNVmatrix[,match(d,colnames(CNVmatrix))]
  datasets=list(mRNAexp=mRNAexp1,miRNAexp=miRNAexp1,lncRNAexp=lncRNAexp1,methylation=methylation1,CNVmatrix=CNVmatrix1)
  W_temp=list()
  for(i in 1:length(datasets))
  {
    distance=dist2(as.matrix(t(datasets[[i]])), as.matrix(t(datasets[[i]])))
    W_temp[[i]] = affinityMatrix(distance, K, alpha)
  }
  W = SNF(W_temp, K=K, t=t)
  group =spectralClustering(W,clusterNum)

  diag(W)=0
  diag(W)=max(W)
  distanceMatrix=W
  attr(distanceMatrix,'class')="Similarity"

  if(plot)
    displayClusters(W, group)
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=group)
  result
}

#' Execute iCluster based on five datasets.
#'
#' Shen (2009) proposed a latent variable regression with a lasso constraint for joint modeling of multiple omics
#' data types to identify common latent variables that can be used to cluster patient samples into biologically and clinically relevant disease subtypes.\cr
#' This function is based on the R package "iCluster".\cr
#' Based on the existing SNF algorithm, we improve it, increase the dataset type, and add the sample preprocessing function, because different sources of datasets, the sample format is not the same, researchers usually spend a lot of time in the unified sample format part, using our improved clustering algorithm,
#' users do not need to unify the sample format, directly input into the clustering algorithm, greatly saving the user's time.
#'
#' @importFrom iCluster iCluster2 plotiCluster
#' @param mRNAexp A matrix of mRNA expression
#' @param miRNAexp A matrix of miRNA expression
#' @param lncRNAexp A matrix of lncRNA expression
#' @param methylation A matrix of methylation
#' @param CNVmatrix A matrix of Copy number variation
#' @param k Number of subtypes for the samples
#' @param lambda Penalty term for the coefficient matrix of the iCluster model
#' @param scalar Logical value. If true, a degenerate version assuming scalar covariance matrix is used.
#' @param max.iter  maximum iteration for the EM algorithm
#' @param scale Logical value. If true, the genomic features in the matrix is centered.
#'
#' @return A list with the following elements.
#'\itemize{
#'   \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'
#'  \item \strong{originalResult} : The clustering result of the original function "iCluster2()".
#'
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  }
#'
#'
#' @details
#'  For iCluster algorithm, it cannot process high-dimensional data, otherwise it is very very time-consuming or reports a mistake.
#'  Based on test, it could smoothly run for the matrix with around 1500 features. Normally it need feature selection step first to reduce feature number.
#' @references
#' Ronglai Shen, Adam Olshen, Marc Ladanyi. (2009). Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. Bioinformatics 25, 2906-2912.\cr
#' Ronglai Shen, Qianxing Mo, Nikolaus Schultz, Venkatraman E. Seshan, Adam B. Olshen, Jason Huse, Marc Ladanyi, Chris Sander. (2012). Integrative Subtype Discovery in Glioblastoma Using iCluster. PLoS ONE 7, e35236
#' Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#' @seealso \code{\link{iCluster2}}
#' @examples
#' data(mRNAexp)
#' data(miRNAexp)
#' data(lncRNAexp)
#' data(methylation)
#' data(CNVmatrix)
#' result=ExecuteiCluster(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix, k=3, lambda=list(0.44,0.33,0.28))
#' result$group
#'
ExecuteiCluster<-function(mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix, k, lambda=NULL, scale=TRUE, scalar=FALSE, max.iter=10)
{
  colnames(mRNAexp)=substr(colnames(mRNAexp),1,15)
  colnames(mRNAexp)=gsub("[.]","-",colnames(mRNAexp))
  colnames(miRNAexp)=substr(colnames(miRNAexp),1,15)
  colnames(miRNAexp)=gsub("[.]","-",colnames(miRNAexp))
  colnames(lncRNAexp)=substr(colnames(lncRNAexp),1,15)
  colnames(lncRNAexp)=gsub("[.]","-",colnames(lncRNAexp))
  colnames(methylation)=substr(colnames(methylation),1,15)
  colnames(methylation)=gsub("[.]","-",colnames(methylation))
  colnames(CNVmatrix)=substr(colnames(CNVmatrix),1,15)
  colnames(CNVmatrix)=gsub("[.]","-",colnames(CNVmatrix))
  a=intersect(x=colnames(miRNAexp),y=colnames(mRNAexp))
  b=intersect(x=colnames(lncRNAexp),y=a)
  c=intersect(x=colnames(methylation),y=b)
  d=intersect(x=colnames(CNVmatrix),y=c)
  mRNAexp1=mRNAexp[,match(d,colnames(mRNAexp))]
  miRNAexp1=miRNAexp[,match(d,colnames(miRNAexp))]
  lncRNAexp1=lncRNAexp[,match(d,colnames(lncRNAexp))]
  methylation1=methylation[,match(d,colnames(methylation))]
  CNVmatrix1=CNVmatrix[,match(d,colnames(CNVmatrix))]
  datasets=list(mRNAexp=mRNAexp1,miRNAexp=miRNAexp1,lncRNAexp=lncRNAexp1,methylation=methylation1,CNVmatrix=CNVmatrix1)
  data1=list()
  for(i in 1:length(datasets))
  {
    data1[[i]]=t(datasets[[i]])
  }

  fit=iCluster2(datasets=data1, k=k, lambda=lambda, scale=scale, scalar=scalar, max.iter=10)

  plotiCluster(fit=fit, label=rownames(data1[[1]]))
  group=fit$clusters
  result=list(group=group,originalResult=fit)
  result
}

#' Execute Consensus Clustering based on five datasets.
#'
#' This function is based on the R package "ConsensusClusterPlus".\cr
#' Based on the existing SNF algorithm, we improve it, increase the dataset type, and add the sample preprocessing function, because different sources of datasets, the sample format is not the same, researchers usually spend a lot of time in the unified sample format part, using our improved clustering algorithm,
#' users do not need to unify the sample format, directly input into the clustering algorithm, greatly saving the user's time.
#'
#' @importFrom ConsensusClusterPlus ConsensusClusterPlus
#' @param clusterNum A integer representing the return cluster number, this value should be less
#' than maxClusterNum(maxK). This is the only additional parameter in our function compared to the original
#' R package "ConsensusClusterPlus". All the other parameters are compatible to the function "ConsensusClusterPlus().
#' @param mRNAexp A matrix of mRNA expression
#' @param miRNAexp A matrix of miRNA expression
#' @param lncRNAexp A matrix of lncRNA expression
#' @param methylation A matrix of methylation
#' @param CNVmatrix A matrix of Copy number variation
#' @param maxK integer value. maximum cluster number  for Consensus Clustering Algorithm to evaluate.
#' @param reps  integer value. number of subsamples(in other words, The iteration number of each cluster number)
#' @param clusterAlg character value. cluster algorithm. 'hc' heirarchical (hclust), 'pam' for paritioning around medoids, 'km' for k-means upon data matrix, 'kmdist' for k-means upon distance matrices (former km option), or a function that returns a clustering.
#' @param distance character value. 'pearson': (1 - Pearson correlation), 'spearman' (1 - Spearman correlation), 'euclidean', 'binary', 'maximum', 'canberra', 'minkowski" or custom distance function.
#' @param  title character value for output directory. This title can be an absolute or relative path
#' @param pItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param pFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param plot Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param innerLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param finalLinkage Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param writeTable Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsItem Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param weightsFeature Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param verbose Please refer to the "ConsensusClusterPlus" package for detailed information.
#' @param corUse Please refer to the "ConsensusClusterPlus" package for detailed information.

#' @return A list with the following elements.
#'\itemize{
#'  \item \strong{group} : A vector represent the group of cancer subtypes. The order is corresponding to the the samples in the data matrix.\cr
#'
#'   This is the most important result for all clustering methods, so we place it as the first component. The format of group
#'   is consistent across different algorithms and therefore makes it convenient for downstream analyses. Moreover, the format
#'   of group is also compatible with the K-means result and the hclust (after using the cutree() function).
#'
#'  \item \strong{distanceMatrix} : It is a sample similarity matrix. The more large value between samples in the matrix, the more similarity the samples are.
#'
#'   We extracted this matrix from the algorithmic procedure because it is useful for similarity analysis among the samples based on the clustering results.
#'
#'  \item \strong{originalResult} : The clustering result of the original function "ConsensusClusterPlus()"
#'
#'  Different clustering algorithms have different output formats. Although we have the group component which has consistent format for all of the algorithms (making it easy for downstream analyses), we still keep the output from the original algorithms.
#'  }
#'
#' @details
#'  If the data is a list containing the matched mutli-genomics  data matrices like the input as "ExecuteiCluster()" and "ExecuteSNF()",
#'   we use "z-score" to normalize features for each data matrix first. Then all the normalized data matrices from the data list are concatenated
#'   according to samples. The concatenated data matrix is the samples with a long features (all features in the data list).
#'   Our purpose is to make convenient comparing the different method with same dataset format. See examples.
#'
#' @seealso \code{ConsensusClusterPlus}
#'
#' @references
#' Monti, S., Tamayo, P., Mesirov, J., Golub, T. (2003) Consensus Clustering: A Resampling-Based Method for Class Discovery and Visualization of Gene Expression Microarray Data. Machine Learning, 52, 91-118.
#' Xu,Taosheng \email{taosheng.x@@gmail.com}, Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'
#' @examples
#' data(mRNAexp)
#' data(miRNAexp)
#' data(lncRNAexp)
#' data(methylation)
#' data(CNVmatrix)
#' result=ExecuteCC(clusterNum=3,mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,maxK=5,clusterAlg="hc",distance="pearson")
#' result$group
#'
ExecuteCC<-function(clusterNum,
                    mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,maxK=10,clusterAlg="hc",
                    distance="pearson",title="ConsensusClusterResult",
                    reps=500, pItem=0.8, pFeature=1,plot="png",
                    innerLinkage="average", finalLinkage="average",
                    writeTable=FALSE,weightsItem=NULL,weightsFeature=NULL,
                    verbose=FALSE,corUse="everything")
{
  colnames(mRNAexp)=substr(colnames(mRNAexp),1,15)
  colnames(mRNAexp)=gsub("[.]","-",colnames(mRNAexp))
  colnames(miRNAexp)=substr(colnames(miRNAexp),1,15)
  colnames(miRNAexp)=gsub("[.]","-",colnames(miRNAexp))
  colnames(lncRNAexp)=substr(colnames(lncRNAexp),1,15)
  colnames(lncRNAexp)=gsub("[.]","-",colnames(lncRNAexp))
  colnames(methylation)=substr(colnames(methylation),1,15)
  colnames(methylation)=gsub("[.]","-",colnames(methylation))
  colnames(CNVmatrix)=substr(colnames(CNVmatrix),1,15)
  colnames(CNVmatrix)=gsub("[.]","-",colnames(CNVmatrix))
  a=intersect(x=colnames(miRNAexp),y=colnames(mRNAexp))
  b=intersect(x=colnames(lncRNAexp),y=a)
  c=intersect(x=colnames(methylation),y=b)
  e=intersect(x=colnames(CNVmatrix),y=c)
  mRNAexp1=mRNAexp[,match(e,colnames(mRNAexp))]
  miRNAexp1=miRNAexp[,match(e,colnames(miRNAexp))]
  lncRNAexp1=lncRNAexp[,match(e,colnames(lncRNAexp))]
  methylation1=methylation[,match(e,colnames(methylation))]
  CNVmatrix1=CNVmatrix[,match(e,colnames(CNVmatrix))]
  CNVmatrix1=CNVmatrix1[which(rowSums(CNVmatrix1) > 0),]
  d=list(mRNAexp=mRNAexp1,miRNAexp=miRNAexp1,lncRNAexp=lncRNAexp1,methylation=methylation1,CNVmatrix=CNVmatrix1)
  if(is.list(d))
  {
    temp=NULL
    for(i in 1: length(d))
    {
      temp=rbind(temp,d[[i]])
    }
    temp=t(scale(t(temp)))
  }
  else
    temp=d
  originalResult=ConsensusClusterPlus(
    temp, maxK=maxK,clusterAlg=clusterAlg,
    distance=distance,title=title,
    reps=reps, pItem=pItem, pFeature=pFeature,plot=plot,
    innerLinkage=innerLinkage, finalLinkage=finalLinkage,
    writeTable=writeTable,weightsItem=weightsItem,weightsFeature=weightsFeature,
    verbose=verbose,corUse=corUse)
  group=originalResult[[clusterNum]][["consensusClass"]]
  distanceMatrix=originalResult[[clusterNum]][["consensusMatrix"]]
  attr(distanceMatrix,'class')="Similarity"
  #icl=calcICL(result,title =fileName,plot="png" )
  result=list(group=group,distanceMatrix=distanceMatrix,originalResult=originalResult)
  result
}


