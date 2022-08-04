# Author: Rafael Bargiela, PhD. Bangor University (UK). 2021

# Full analysis of PSMs abundance based on NSAF (Normalized Spectra Abundance factor) and lnNSAF
################################################################################################
  # M: matrix with PSMs for A and B samples or datasets
  # length.vector: Vector with the same length as nrow of M with the length of the proteins to be analyse
  # several.samples: Logical, indicating if A and B are datasets with more than one sample (TRUE), or if M is just a two columns matrix (FALSE)
  # g1 and g2: When several.samples is TRUE, g1 and g2 indicates which columns belong to A (g1) and B (g2)
  # group.analysis: method for group analysis. If 'mean', average among groups is calculated, if 'sum' analysis is based on summation of PSMs of each group.
  # min.count.sample: Minimum number of PSMs on any of the samples for each of the proteins to analysed
  # min.total.count: Minimum total number of PSMs for each protein to be analysed
  
PSManalysis<-function(M,length.vector,several.samples=FALSE,g1=NA,g2=NA,group.analysis="sum",min.count.sample=2,min.total.count=5){

  if(nrow(M)!=length(length.vector)){stop("length.vector is of different size than number of rows in M")}
  if(several.samples==TRUE){
    if(is.numeric(g1)==FALSE){stop("g1 can't be NA when several.samples is TRUE")}
    if(is.numeric(g2)==FALSE){stop("g2 can't be NA when several.samples is TRUE")}
    if(grepl("^mean|sum$",group.analysis,perl=TRUE)==FALSE){stop("group.analysis method must be 'mean' or 'sum' if several.samples are set")}    
  }else{
    if(ncol(M)>2){warning("Input matrix have more than 2 columns, only first and second will be used for analysis")}
  }  
  
  v<-vector("numeric")
  for(r in 1:nrow(M)){
    if(several.samples==TRUE){
      if(max(M[r,])>=min.count.sample & sum(M[r,])>=min.total.count){
        v<-append(v,r,length(v))        
      }else{
      }
    }else{
      if(max(M[r,])>=min.count.sample){         # Getting rid of low abundance proteins
        v<-append(v,r,length(v))                # In case A and B are result of merging several dataset samples
      }else{
      }
    }
  }
  
  FM<-{}
  if(several.samples==TRUE){
    if(group.analysis=="sum"){
      A<-rowSums(as.matrix(M[v,g1]))
      B<-rowSums(as.matrix(M[v,g2]))
    }
    if(group.analysis=="mean"){
      A<-rowMeans(as.matrix(M[v,g1]))
      B<-rowMeans(as.matrix(M[v,g2]))
    }
    filtered.proteins<-cbind(A,B) # low abundance proteins discarded    
  }else{
    filtered.proteins<-as.matrix(M[v,]) # low abundance proteins discarded 
  }
  L<-length.vector[v]      # AA length of filtered.proteins
  PSMadj<-filtered.proteins+1 # PSM+1 matrix
  # Corrected relative abundance (PSM+1 normalized)
  # PSMs are normalized over the samples with less total number of PSM+1
  PSMsum<-colSums(PSMadj)
  if(PSMsum[1]<PSMsum[2]){
    Ca<-PSMadj[,1]
    Cb<-PSMadj[,2]*(PSMsum[1]/PSMsum[2])
  }else{
    Cb<-PSMadj[,2]
    Ca<-PSMadj[,1]*(PSMsum[2]/PSMsum[1])
  }
  PSM<-matrix(c(filtered.proteins,PSMadj,Ca,Cb),nc=6,dimnames = list(rownames(filtered.proteins),c("PSM of filtered proteins A","PSM of filtered proteins B","PSM+1 A","PSM+1 B","norm PSM+1 A","norm PSM+1 B")))
  
  # Statistics
  Rab<-Ca/Cb # Abundance ratio Sample A:B
  Rba<-Cb/Ca # Abundance ratio Sample B:A
  Gt<-as.numeric(2*(Ca*log(Ca/((Ca+Cb)/2))+Cb*log(Cb/((Ca+Cb)/2))) ) # G-test
  Gp<-as.numeric(pchisq(Gt,df=1,lower.tail = FALSE))                 # G-tet p-value by chi-squared
  GpFDR<-as.numeric(p.adjust(Gp,method="BH"))                        # False Discovery Rate (FDR) of the G-test p-value by Benjamin-Hochberg method
  # Complementary error function over adjusted normal distribution for log2 protein abundance ratio (R)
  lRab<-log(Rab,base=2) # Getting corrected Abundance ratio A:B over normal distribution mean
  fitAB<-fitdistr(lRab,densfun="normal")
  mAB<-fitAB$estimate[1]
  sdAB<-fitAB$estimate[2]
  Rab.corrected<-(lRab-mAB)/(sdAB/sqrt(2))
  lRba<-log(Rba,base=2) # Getting corrected Abundance Ratio B:A over normal distribution mean
  fitBA<-fitdistr(lRba,densfun="normal")  
  mBA<-fitBA$estimate[1]
  sdBA<-fitBA$estimate[2]
  Rba.corrected<-(lRba-mBA)/(sdBA/sqrt(2))
  
  # erfc<-function(x) 2*pnorm(-sqrt(2)*x,lower.tail=FALSE) # Complementary error function (over A:B values)
  # RabERFC<-erfc(abs(Rab.corrected))
  erfc<-function(x,mean,sd) 2*pnorm(-sqrt(2)*x,mean=mean,sd=sd,lower.tail=TRUE) # Complementary error function (over A:B values)
  RabERFC<-erfc(abs(Rab.corrected),mAB,sdAB)
  erfcFDR<-as.numeric(p.adjust(RabERFC,method="BH")) 
  
  G<-matrix(c(Rab,Rba,Rab.corrected,Rba.corrected,Gt,Gp,GpFDR,RabERFC,erfcFDR),nc=9,dimnames = list(rownames(filtered.proteins),c("Abundance ratio (R A:B)","Abundance ratio (R B:A)","R A:B corrected","R B:A corrected","G-test","G-test p-value","adjusted G-test (FDR)","ERFC","adjusted ERFC (FDR)")))
  
  # NSAF values
  SAF1<-as.numeric(Ca/L)
  SAF2<-as.numeric(Cb/L)
  NSAF1<-as.numeric(SAF1/sum(SAF1))
  NSAF2<-as.numeric(SAF2/sum(SAF2))
  lnNSAF1<-as.numeric(log(NSAF1))
  lnNSAF2<-as.numeric(log(NSAF2))
  lnFC<-lnNSAF1-lnNSAF2 # fold change (on sustraction format because the use of logarithms)
  fitlnNSAF1<-fitdistr(lnNSAF1,densfun = "normal")
  fitlnNSAF2<-fitdistr(lnNSAF2,densfun="normal")
  fitlnFC<-fitdistr(lnFC,densfun="normal")
  AF<-matrix(c(SAF1,SAF2,NSAF1,NSAF2,lnNSAF1,lnNSAF2,lnFC),nc=7,dimnames = list(rownames(filtered.proteins),c("SAF A","SAF B","NSAF A","NSAF B","ln(NSAF) A","ln(NSAF) B","Fold change (ln(NSAFa)-ln(NSAFb))")))
  
  #
  Gauss.Rab<-c(mAB,sdAB)
  Gauss.Rba<-c(mBA,sdBA)
  Gauss.lnNSAF1<-c(fitlnNSAF1$estimate[1],fitlnNSAF1$estimate[2])
  Gauss.lnNSAF2<-c(fitlnNSAF2$estimate[1],fitlnNSAF2$estimate[2])  
  Gauss.lnFC.NSAF<-c(fitlnFC$estimate[1],fitlnFC$estimate[2])
  NDE<-matrix(c(Gauss.Rab,Gauss.Rba,Gauss.lnNSAF1,Gauss.lnNSAF2,Gauss.lnFC.NSAF),nc=2,byrow=TRUE,dimnames=list(c("log2 Ratio A:B","log2 Ratio B:A","lnNSAF A","lnNSAF B","Fold change ln(NSAF)"),c("Mean","SD")))
  # Assembling final object
  FM$PSM.adjustment<-PSM
  FM$Statistics<-G
  FM$Abundance.factor<-AF
  FM$Normal.Dist.Estimates<-NDE
  return(FM)
}
