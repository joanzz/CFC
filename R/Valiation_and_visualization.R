#' Surplot
#'
#' Draw a KM survival map
#'
#' @importFrom survival survfit Surv
#' @importFrom survminer ggsurvplot
#' @param results The results after clustering
#' @param meta The clinical data downloaded from TCGA
#' @param optK The optimal K
#'
#' @return A figure
#' @export
#'
#' @examples
#'data(results)
#'data(meta)
#'maxK = 6
#'Kvec = 2:maxK
#'x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
#'PAC = rep(NA,length(Kvec))
#'names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
#'for(i in Kvec){
#'  M = results[[i]]$consensusMatrix
#'  Fn = ecdf(M[lower.tri(M)])
#'  PAC[i-1] = Fn(x2) - Fn(x1)
#'}
#'optK = Kvec[which.min(PAC)]
#'Surplot(results,meta,optK=3)
#'
Surplot <-function(results,meta,optK=2){
  Cluster <- results[[optK]]$consensusClass
  names(Cluster)=gsub('[.]', '-', names(Cluster))
  names(Cluster)=substr(names(Cluster),1,15)
  meta=meta[match(names(Cluster),rownames(meta)),]
  meta$Cluster = Cluster
  sfit <- survfit(Surv(OS.time, OS) ~ Cluster,data = meta)
  ggsurvplot(sfit,data=meta,pval = T,palette = "jco",conf.int = TRUE)
}


#' Dif.limma
#'
#' Differently Expression Analysis for genomic data. We apply limma package to conduct the analysis.
#' @importFrom limma lmFit eBayes makeContrasts contrasts.fit voom topTable
#' @param Tumor_Data A matrix representing the genomic data of cancer samples such as gene expression data, miRNA expression data.\cr
#' For the matrix, the rows represent the genomic features, and the columns represent the cancer samples.
#' @param Normal_Data A matrix representing the genomic data of Normal samples.\cr
#' For the matrix, the rows represent the genomic features corresponding to the Tumor_Data, and the columns represent the normal samples.
#' @param group A vector representing the subtype of each tumor sample in the Tumor_Data. The length of group is equal to the column number of Tumor_Data.
#' @param topk The top number of different expression features that we want to extract in the return result.
#' @param sort.by This is a parmeter of "topTable() in limma pacakge". "Character string specifying statistic to rank genes by. Possible values for topTable and toptable are "logFC", "AveExpr", "t", "P", "p", "B" or "none". (Permitted synonyms are "M" for "logFC", "A" or "Amean" for "AveExpr", "T" for "t" and "p" for "P".) Possibilities for topTableF are "F" or "none". Possibilities for topTreat are as for topTable except for "B"."
#' @param adjust.method This is a parmeter of "topTable() in limma pacakge".  Refer to the "method used to adjust the p-values for multiple testing. Options, in increasing conservatism, include "none", "BH", "BY" and "holm". See p.adjust for the complete list of options. A NULL value will result in the default adjustment method, which is "BH"."
#' @param RNAseq A bool type representing the input datatype is a RNASeq or not. Default is FALSE for microarray data.
#' @return
#' A list representing the differently expression for each subtype comparing to the Normal group.
#' @examples
#' data("mRNAexp")
#' mRNAexp=as.matrix(mRNAexp)
#' Normal_Data=as.matrix(mRNAexp[,sample(1:100,20)])
#' result=Dif.limma(Tumor_Data=mRNAexp,Normal_Data=Normal_Data,group=NULL,topk=NULL,RNAseq=FALSE)
#'
#' @references
#' Smyth, Gordon K. "Limma: linear models for microarray data." Bioinformatics and computational biology solutions using R and Bioconductor.
#' Springer New York, 2005. 397-420.
#' Xu T, Le TD, Liu L, Su N, Wang R, Sun B, Colaprico A, Bontempi G, Li J. CancerSubtypes: an R/Bioconductor package for molecular cancer subtype identification, validation and visualization. Bioinformatics. 2017 Oct 1;33(19):3131-3133. doi: 10.1093/bioinformatics/btx378. PMID: 28605519.
#' @export

Dif.limma<-function(Tumor_Data,Normal_Data,group=NULL,topk=NULL,sort.by="p", adjust.method="BH",RNAseq=FALSE)
{
  if(is.null(group))
    group=rep(1,ncol(Tumor_Data))
  groupN=sort(unique(group))
  len=length(groupN)

  mylist.names <- paste("Subtype",order(groupN))
  result <- vector("list", length(mylist.names))
  names(result) <- mylist.names

  num=ncol(Normal_Data)
  for(i in 1:len)
  {
    index=which(group==groupN[i])
    #label=c(rep(0,num),rep(1,length(index)))
    #label=as.factor(label)
    #design <- model.matrix(~-1+label)
    #colnames(design)=c("Normal","Cancer")

    Normal=NULL
    Cancer=NULL
    design=cbind(Normal=c(rep(1,num), rep(0,length(index))), Cancer=c(rep(0,num), rep(1,length(index))))

    Data=cbind(Normal_Data,Tumor_Data[,index])
    if(RNAseq){
      mR <- voom(Data, design, plot=TRUE)
      mR <-as.matrix(mR)}
    else
      mR=Data

    ###In order to return the index of features, set the same name for two feautre
    ### then restore the name
    name1=rownames(mR)[1]
    name2=rownames(mR)[2]
    rownames(mR)[1]="repeat"
    rownames(mR)[2]="repeat"
    ########## mR ############
    mRfit=lmFit(mR, design)
    mRfit=eBayes(mRfit)
    contrast.matrix=makeContrasts(CancervNormal=Cancer - Normal, levels=design)
    mRfit2=contrasts.fit(mRfit, contrast.matrix)
    mRfit2=eBayes(mRfit2)

    if(is.null(topk))
    {
      topk=nrow(mR)
    }
    mRresults=topTable(mRfit2, number= topk, sort.by=sort.by, adjust.method=adjust.method)
    ######restore feature name
    mRresults[which(rownames(mRresults)=="1"),1]=name1
    mRresults[which(rownames(mRresults)=="2"),1]=name2

    result[[i]]=mRresults
  }
  result
}

#' Plot optimal K
#'
#' Plot a figure for optimal K.
#'
#'@importFrom ggplot2 theme_classic geom_line scale_y_continuous theme xlab aes
#'@importFrom RColorBrewer brewer.pal
#' @param dot_df A dataset containing K,DeltaArea,PAC columns.
#'
#' @return A figure.
#' @export
#'
#' @examples
#' data(dot_df)
#' PACplot(dot_df)
#'
PACplot <- function(dot_df){
  set1 <- c(brewer.pal(9,"Set1"))
  p1 <- ggplot() + theme_classic()+
    geom_line(mapping = aes(x = dot_df$k,
                            y = dot_df$DeltaArea),
              group=1 , size=2 ,
              color=set1[2])  +

    geom_line(mapping = aes(x = dot_df$k,
                            y = dot_df$PAC*5 ),
              group = 1, size = 2, color = set1[1]) +


    scale_y_continuous(breaks = seq(0,1,0.2), "Delta Area",
                       sec.axis = sec_axis(~./5,
                                           name = "PAC"))

  p1 <- p1 + theme(axis.text.x = element_text(size = 16,
                                              vjust = 0.5,
                                              hjust = 0.5
  ))+
    theme(axis.text.y = element_text(size = 16,
                                     vjust = 0.5,
                                     hjust = 0.5))+
    theme(axis.text.y.right = element_text(size = 16,
                                           vjust = 0.5,
                                           hjust = 0.5))



  p1 <- p1 + xlab("K Number") +
    theme(axis.title.x = element_text(size = 16,

                                      face = "bold",
                                      vjust = 0.5,
                                      hjust = 0.5))+


    theme(axis.title.y = element_text(size = 16,

                                      face = "bold",
                                      vjust = 0.5,
                                      hjust = 0.5))

  p1


}

#' Drug senstivity
#'
#' Predict drug senstivity for a drug in CGP,and visualize with ggplot.
#'
#' @importFrom ggplot2 ggplot geom_boxplot theme_bw ylab
#' @import pRRophetic
#' @param testMatrix a gene expression matrix with gene names as row ids and sample names as column ids
#' @param studyResponse Drug handling information
#' @param bortIndex Drug handling information
#' @param drug the name of the drug for which you would like to predict sensitivity, one of A.443654, A.770041, ABT.263, ABT.888, AG.014699, AICAR, AKT.inhibitor.VIII, AMG.706, AP.24534, AS601245, ATRA, AUY922, Axitinib, AZ628, AZD.0530, AZD.2281, AZD6244, AZD6482, AZD7762, AZD8055, BAY.61.3606, Bexarotene, BI.2536, BIBW2992, Bicalutamide, BI.D1870, BIRB.0796, Bleomycin, BMS.509744, BMS.536924, BMS.708163, BMS.754807, Bortezomib, Bosutinib, Bryostatin.1, BX.795, Camptothecin, CCT007093, CCT018159, CEP.701, CGP.082996, CGP.60474, CHIR.99021, CI.1040, Cisplatin, CMK, Cyclopamine, Cytarabine, Dasatinib, DMOG, Docetaxel, Doxorubicin, EHT.1864, Elesclomol, Embelin, Epothilone.B, Erlotinib, Etoposide, FH535, FTI.277, GDC.0449, GDC0941, Gefitinib, Gemcitabine, GNF.2, GSK269962A, GSK.650394, GW.441756, GW843682X, Imatinib, IPA.3, JNJ.26854165, JNK.9L, JNK.Inhibitor.VIII, JW.7.52.1, KIN001.135, KU.55933, Lapatinib, Lenalidomide, LFM.A13, Metformin, Methotrexate, MG.132, Midostaurin, Mitomycin.C, MK.2206, MS.275, Nilotinib, NSC.87877, NU.7441, Nutlin.3a, NVP.BEZ235, NVP.TAE684, Obatoclax.Mesylate, OSI.906, PAC.1, Paclitaxel, Parthenolide, Pazopanib, PD.0325901, PD.0332991, PD.173074, PF.02341066, PF.4708671, PF.562271, PHA.665752, PLX4720, Pyrimethamine, QS11, Rapamycin, RDEA119, RO.3306, Roscovitine, Salubrinal, SB.216763, SB590885, Shikonin, SL.0101.1, Sorafenib, S.Trityl.L.cysteine, Sunitinib, Temsirolimus, Thapsigargin, Tipifarnib, TW.37, Vinblastine, Vinorelbine, Vorinostat, VX.680, VX.702, WH.4.023, WO2009093972, WZ.1.84, X17.AAG, X681640, XMD8.85, Z.LLNle.CHO, ZM.447439.
#' @param tissueType specify if you would like to traing the models on only a subset of the CGP cell lines (based on the tissue type from which the cell lines originated). This be one any of "all" (for everything, default option), "allSolidTumors" (everything except for blood), "blood", "breast", "CNS", "GI tract" ,"lung", "skin", "upper aerodigestive"
#' @param batchCorrect How should training and test data matrices be homogenized. Choices are "eb" (default) for ComBat, "qn" for quantiles normalization or "none" for no homogenization.
#' @param selection How should duplicate gene ids be handled. Default is -1 which asks the user. 1 to summarize by their or 2 to disguard all duplicates.
#' @param dataset The datasets you want to choose, including 2014 and 2016
#' @param colours The color you are interested
#'
#' @return A figure.
#' @export
#'
#' @references Geeleher P, Cox N, Huang RS. pRRophetic: an R package for prediction of clinical chemotherapeutic response from tumor gene expression levels. PLoS One. 2014 Sep 17;9(9):e107468. doi: 10.1371/journal.pone.0107468. PMID: 25229481; PMCID: PMC4167990.
#' @examples
#' data(exprData)
#' data(studyResponse)
#' drugsensitivity(testMatrix=exprData,
#' studyResponse=studyResponse,
#' bortIndex=bortIndex,
#' drug="Bosutinib",
#' tissueType = "all",
#' batchCorrect = "eb",
#' selection=1,
#' dataset = "cgp2014",colours=c('#e94753','#47a4e9'))
#'
drugsensitivity <-function(testMatrix=exprData,
                           studyResponse=studyResponse,
                           bortIndex=bortIndex,
                           drug="Bosutinib",
                           tissueType = "all",
                           batchCorrect = "eb",
                           selection=1,
                           dataset = "cgp2014",colours=c('#e94753','#47a4e9')){
  predictedPtype <- suppressWarnings(pRRopheticPredict(testMatrix=testMatrix,
                                                       drug=drug,
                                                       tissueType = tissueType,
                                                       batchCorrect = batchCorrect,
                                                       selection=selection,
                                                       dataset = dataset))
  df <- stack(list(NR = predictedPtype[((studyResponse == 'PGx_Responder = NR') & bortIndex)],
                   R = predictedPtype[((studyResponse == 'PGx_Responder = R') & bortIndex)]))
  ggplot(data = df,
         aes(y = values,
             x = ind))+
    geom_boxplot(alpha = 0.3,
                 fill = colours)+
    theme_bw()+
    ylab(paste("Predicted", drug,"Sensitivity")) +
    xlab('Clinical Response')
}

#' Entropy rate
#'
#' This function computes the entropy rate (SR) and visualization.
#'
#' @importFrom ggpubr ggboxplot stat_compare_means
#' @importFrom igraph arpack
#' @param entr_expr a genome-wide expression vector of a sample with names labeling the genes (also annotated to entrez gene IDs) of same length as rows of adj.m.
#' @param entr_adja an adjacency matrix representing a connected PPI network, with rownames/colnames of the matrix annotated to Entrez gene IDs.
#' @param entr_labels Classification (label) of the sample
#' @param local logical, if TRUE also compute the local Shannon (normalised) entropies.
#' @param method scalar specifying method to compute stochastic matrix: 1=stochastic vector of a node is independent of the node's gene expression value, 2=stochastic vector of a node is dependent on the node's gene expression value through use of the betaFn. Default is 1.
#' @param quants.v optionally, a vector of length 2 specifying the quantiles for defining low and high expression (only needed if method=2 is used).
#'
#' @return A boxplot about the entropy.
#' @export
#'
#' @references CompSR.R DoIntegPIN.R  Author: Andrew Teschendorff Date: 31 Dec 2013.
#' @examples
#' data("entr_adja")
#' data("entr_expr")
#' data("entr_labels")
#' ComENTR(entr_expr,entr_adja,entr_labels,local=TRUE,method=1)
#'
ComENTR <-function(entr_expr,entr_adja,entr_labels,local=TRUE,method=1,quants.v=c(0.1,0.9)){
  pin.o <- DoIntegPIN(entr_adja,entr_expr)
  #require(igraph)
  fa <- function(x,extra=NULL) {
    as.vector(pin.o$a %*% x)
  }
  ap.o <- arpack(fa, options=list(n=nrow(pin.o$a),nev=1,which="LM",maxiter=30120),sym=TRUE);
  v <- ap.o$vectors
  lambda <- ap.o$values
  maxSR <- log(lambda)

  ####
  selS.idx = c(1:ncol(pin.o$e))
  tmpE.m <- pin.o$e[,selS.idx]
  sr.lo <- list()
  k.v <- rowSums(pin.o$a)
  for(s in 1:ncol(tmpE.m)){
    p.m <- matrix(0,nrow=length(tmpE.m[,s]),ncol=length(tmpE.m[,s]))
    rownames(p.m) <- rownames(pin.o$a)
    colnames(p.m) <- rownames(pin.o$a)

    if(method==1){
      for(g in 1:nrow(pin.o$a)){
        nn.idx <- which(pin.o$a[g,]==1)
        term2.v <- tmpE.m[,s][nn.idx]/sum(tmpE.m[,s][nn.idx])
        p.m[g,nn.idx] <- term2.v
      }
    }else if(method==2){
      alpha.v <- quantile(as.vector(tmpE.m[,s]),probs=quants.v)
      beta.v <- betaFn(tmpE.m[,s],alpha.v)
      for(g in 1:nrow(pin.o$a)){
        nn.idx <- which(pin.o$a[g,]==1)
        term2.v <- tmpE.m[,s][nn.idx]/sum(tmpE.m[,s][nn.idx])
        p.m[g,nn.idx] <- (1-beta.v[g])/k.v[g] + beta.v[g]*term2.v
      }
    }


    fp <- function(x,extra=NULL) {
      as.vector(p.m %*% x)
    }
    fpt <- function(x,extra=NULL) {
      as.vector(t(p.m) %*% x)
    }
    ap.o <- arpack(fpt, options=list(n=nrow(p.m),nev=1,which="LM"),sym=FALSE);
    invP.v <- abs(as.numeric(ap.o$vectors));
    invP.v <- invP.v/sum(invP.v);
    S.v <- apply(p.m,1,CompS);
    SR <- sum(invP.v*S.v);
    if(is.null(maxSR)==FALSE){## if provided then normalise relative to maxSR
      SR <- SR/maxSR;
    }
    if(local){
      CompNS <- function(p.v){
        tmp.idx <- which(p.v>0);
        if(length(tmp.idx)>1){
          S <-  - sum( p.v[tmp.idx]*log(p.v[tmp.idx]) )/log(length(tmp.idx));
        }
        else { ### one degree nodes have zero entropy, avoid singularity.
          S <- 0;
        }
        return(S);
      }
      NS.v <- apply(p.m,1,CompNS);
    }else{
      NS.v <- NULL;
    }
    sr.lo[[s]] <- list(sr=SR,inv=invP.v,s=S.v,ns=NS.v)

    print(s)}
  sr.v <- vector();
  for(s in 1:length(sr.lo)){
    sr.v[s] <- sr.lo[[s]]$sr;
  }
  sr.v <- as.data.frame(sr.v)
  sr.v$sample <- colnames(tmpE.m)
  names(sr.v)[1] <- 'entr'
  entr=merge(sr.v, entr_labels, by = "sample")
  t.test(entr~REC,data = entr)
  p <- ggboxplot(entr, x="REC", y="entr", color = "REC")
  p <- p+stat_compare_means(method = "t.test")
  p
}

#'silhouette_SimilarityMatrix
#'
#'Silhouette refers to a method of interpretation and validation of consistency within clusters of data.
#'The technique provides a succinct graphical representation of how well each object lies within its cluster (From Wiki).\cr
#'Note that: This function is a rewriting version of the function "silhouette()" in R package cluster.
#'   The original function "silhouette()" is to compute the silhouette information based on a dissimilarity matrix.
#'   Here the silhouette_SimilarityMatrix() is to solve the computation based on the similarity matrix.
#'   The result of the silhouette_SimilarityMatrix() is compatible to the function "Silhouette()".
#'
#'@param group A vector represent the cluster label for a set of samples.
#'@param similarity_matrix A similarity matrix between samples
#'@return
#' An object, sil, of class silhouette which is an [n x 3] matrix with
#' attributes. The colnames correspondingly are c("cluster", "neighbor", "sil_width").
#' @details
#' For each observation i, the return sil[i,] contains the cluster to which i belongs as well as the neighbor
#' cluster of i (the cluster, not containing i, for which the average
#' dissimilarity between its observations and i is minimal),
#' and the silhouette width s(i) of the observation.
#' @examples
#' data(mRNAexp)
#' data(miRNAexp)
#' data(lncRNAexp)
#' data(methylation)
#' data(CNVmatrix)
#' result=ExecuteCC(clusterNum=3,mRNAexp,miRNAexp,lncRNAexp,methylation,CNVmatrix,maxK=5,clusterAlg="hc",distance="pearson")
#' sil=silhouette_SimilarityMatrix(result$group, result$distanceMatrix)
#' plot(sil)
#'
#' @seealso \code{\link{silhouette}}
#'
#'
#'@references
#' Rousseeuw, P.J. (1987) Silhouettes: A graphical aid to the interpretation and validation of cluster analysis. J. Comput. Appl. Math., 20, 53-65.
#' Xu,Taosheng \email{taosheng.x@@gmail.com},Thuc Le \email{Thuc.Le@@unisa.edu.au}
#'@export
#'
silhouette_SimilarityMatrix<-function(group, similarity_matrix)
{
  similarity_matrix=as.matrix(similarity_matrix)
  similarity_matrix<-(similarity_matrix+t(similarity_matrix))/2
  diag(similarity_matrix)=0
  normalize <- function(X) X / rowSums(X)
  similarity_matrix<-normalize(similarity_matrix)

  n <- length(group)
  if(!all(group == round(group))) stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if(k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if(doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))

  wds <- matrix(NA, n,3, dimnames =list(names(group), c("cluster","neighbor","sil_width")))
  for(j in 1:k)
  {
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index, drop = FALSE], 2,
                          function(r) tapply(r, group[!index], mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index,"neighbor"] <- cluster_id[-j][maxC]
    s.i <- if(Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj - 1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i) / pmax(b.i, a.i), 0)
    } else 0
    wds[index,"sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}

