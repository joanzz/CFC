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
#' Plot a figure for optimal K
#'@importFrom ggplot2 theme_classic geom_line scale_y_continuous theme xlab
#'@importFrom RColorBrewer brewer.pal
#' @param dot_df A dataset containing K,DeltaArea,PAC columns
#'
#' @return A figure
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

#' Predict drug senstivity for a drug in CGP,and visualize with ggplot.
#'
#' @importFrom ggplot2 ggplot geom_boxplot theme_bw
#' @importFrom pRRophetic pRRopheticPredict
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
#' @return A figure
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





