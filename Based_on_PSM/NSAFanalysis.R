# Author: Rafael Bargiela, PhD. Bangor University (UK). 2021

# Full analysis of PSMs abundance based on NSAF (Normalized Spectra Abundance factor) and lnNSAF
################################################################################################
  # M: Matrix with PSMs for each sample or datasets
  # Metadata: Matrix containing samples metadata per row, with same rows than columns on M. It must have a column named "Group" or defined by Metadata.group.col
  # Metadata.group.col: Character if length 1, definining Metadata column were groups are defined, in case column Group is not included in the Metadata matrix. If Group column exsists, Metadata.group.col is ignored.
  # length.vector: Vector with the same length as nrow of M with the length of the proteins to be analyzed. It could be a column inside M, in that case name it as 'length.vector.column'.
  # groups.vector : Numerical or Character vector of length 2, defining the groups to be compare if group.comparison TRUE.
  # filtering: Filter proteins by min.count.per and min.total.count. Set FALSE in case using a set of previously filtered proteins.
  # group.comparison: Comparing all groups against selected control group
  # control: Position of the group considered as Control or initial to compare the rest of groups, if group.comparison is TRUE
  # min.count.per: Minimum percentage of samples within a group where a protein must be present. If a protein is present in the minimum percentage selected in any of the groups, it is selected, otherwise is discarded for further analysis. (default 60%, set 0 to include all proteins)
  # include.small.groups: logical. If TRUE, groups with n samples < 3 are included in the filtering, so when the expression of a protein meets the criteria within this groups, they are included in the final table, as it would happens with any other group. Default FALSE.
  # LongAnalysis : Logical, if TRUE, methods to adapt the data for Longitudinal Studies are applied to imputed data previously to NSAF and lnNSAF calculation. Method is default set to "limma", unique method currently, to applied Mixed Models to fit the data and remove random effect due to individuals variable.
  # method : Method to applied for Longitudinal analysis. Default, Mixed Models from limma package.

PSManalysis<-function(M,length.vector=NULL,Metadata,Metadata.group.col=NA,filtering=TRUE,min.count.per=60,min.total.count=2,include.small.groups=FALSE,group.comparison=FALSE,groups.vector=NULL,control=1,group.analysis="sum",LongAnalysis=FALSE,method="limma"){
  require(MASS)
  require(limma)
  # require(insight)
  # CHECKINGS -----------------------------------------------------------------------
  if(is.null(length.vector)==TRUE){
    if(sum(grepl("length.vector.column",colnames(M)))<1){
      stop("ERROR: Protein length needs to be provided:
          length.vector is NA and length.vector.column has not been detected as any column name at M")
    }else{
      lv<-grep("length.vector.column",colnames(M),value=FALSE,perl=TRUE)
      message(paste("Column",lv,"will be used as protein length vector"))
      length.vector<-M[,lv]
      names(length.vector)<-rownames(M)
      M<-M[,-lv]
    }
  }else{
    if(nrow(M)!=length(length.vector)){stop("length.vector is of different size than number of rows in M")}else{
        names(length.vector)<-rownames(M)
    }
  }
  if(nrow(Metadata)!=ncol(M)){
    stop("Number of rows on Metadata must be equal to number of columns on M")
  }else{
    if(sum(grepl("^Group$",colnames(Metadata),perl=TRUE))>0){
      message(c("Following groups were found:",paste(unique(Metadata$Group),collapse=" ")))
    }else{
        if(sum(grepl(paste("^",Metadata.group.col,"$",sep=""),colnames(Metadata),perl=TRUE))>0){
            message(c("Following groups were found:",paste(unique(Metadata[,Metadata.group.col]),collapse=" ")))
              Metadata$Group<-Metadata[,Metadata.group.col]
        }else{
          stop(paste("ERROR: Couldn't find Metatadata.group.col on Metadata:",Metadata.group.col,sep=" "))
        }
    }
  }

  if(group.comparison==TRUE & length(groups.vector)!=2){
    stop("ERROR: groups.vector has the wrong length (it must be two)")
  }else{
    if(sum(grepl(paste("^",paste(groups.vector,collapse="|"),"$",sep=""),Metadata$Group,perl=TRUE))==2){
        message(c("Groups to analyze:",paste(unique(groups.vector),collapse=" ")))         
    }
  }

  if(grepl("^mean|sum$",group.analysis,perl=TRUE)==FALSE){stop("group.analysis method must be 'mean' or 'sum' if several.samples are set")}
  # END CHECKINGS --------------------------------------------------------------------
  FM<-{} # Final Matrix
  G<-unique(Metadata$Group) #  # Getting the different groups
  print(G)
  if(filtering==TRUE){
    # Analysing the abundance of each protein on each group
    message("- Filtering proteins")
    Gprots<-vector("character")
    for(g in 1:length(G)){
      group.samples<-grep(G[g],Metadata$Group,value=FALSE,perl=TRUE) # Getting samples belonging to each group
      if(include.small.groups==FALSE & length(group.samples)<3){
          message(paste("Group",G[g],"won't be considered for filtering due to n<3",sep=" "))
      }else{
        gM<-as.matrix(M[,group.samples]) # Specific Matrix for each group
        if(min.count.per==0){
          min.count<-1 
        }else{
          min.count<-ceiling(ncol(gM)*min.count.per/100) 
        }
        for(r in 1:nrow(gM)){
          if(sum(gM[r,]!=0,na.rm=TRUE)>=min.count & sum(gM[r,])>=min.total.count){
            Gprots<-append(Gprots,rownames(gM)[r],length(Gprots))
          }
        }
      }
    }
    filtered.prots<-unique(Gprots)
    filt.length<-length.vector[filtered.prots]
    message(paste("-- ",length(filtered.prots)," Proteins filtered"))
    FM$PSMs.filtered<-M[filtered.prots,]
  }else{
    message("- Skipping filtration of proteins")
    message("-- WARNING: Ensure your proteins were previously filtered")
    FM$PSMs.filtered<-M
    filt.length<-length.vector
  }
  # NORMALIZATIONS of PSMs values ----------------------------------------------
  message("- Normalization of filtered proteins")
    ## Getting Adjusted PSM (PSM+1)-------------------------------------------
    PSMadj<-FM$PSMs.filtered+1 #  Abundance Adjusted Matrix. Corrected relative abundance (PSM+1)
    PSMadjTot<-colSums(PSMadj)
    smallest<-which.min(PSMadjTot) # Sample with smallest number of PSM
    PSMnorm<-sapply(1:ncol(PSMadj),function(x){
      return(PSMadj[,x]*(PSMadjTot[smallest]/PSMadjTot[x]))
    },simplify=TRUE)
    message(paste("-- Normalization done multiplying by the ratio between each sample and the sample with smallest total number of PSMs:",colnames(PSMadj)[smallest]))
    colnames(PSMnorm)<-paste(colnames(PSMadj),"Norm PSM+1")
    rownames(PSMnorm)<-rownames(FM$PSMs.filtered)
    ## In case of Longitudinal (time) analysis
    if(LongAnalysis==TRUE){
      message("NOTE: For Longitudinal analysis, metadata must have 'Time', 'Individuals' and 'Group' columns for each sample")
        if(method=="limma"){
          ### Limma method with Humanzee package
          require(limma)
          require(Humanzee)
          message("Using Mixed Models (MM), with limma and Humanzee package")
          message("Current Linear Mixed model is: Intensity~1+Individuals; with Group as block")   
          design<-model.matrix(~1+Individuals,data=Metadata)
          block<-Metadata$Group
          dup_corrs <- duplicateCorrelation(PSMnorm,design = design, block = block)
          mdl.fit<-Humanzee::ruv_mixed_model(PSMnorm,ndups = 1,design = design, block = block,correlation = dup_corrs$cons)
          M.fit <- t( design %*% t(mdl.fit$coef) ) + mdl.fit$resid
          FM$PSMs.Norm.Mixed.Model.Fit.limma<-M.fit
          message("Using Imputed data after random effect removal")
          FM$PSMs.Norm.NoMixedModel<-FM$PSMs.Norm 
          FM$PSMs.Norm<-FM$PSMs.Norm.Mixed.Model.Fit.limma
        }else{
          stop("Wrong method selection for Longitudinal Analysis")
        }
    }else{
      message("Using Imputed data without previous random effect removal")
          FM$PSMs.Norm<-PSMnorm
    }
    ## Getting means of the log2 Norm Adjusted PSM for each group ------------
      log2.PSMs.Norm.means<-sapply(1:length(G),function(x){
        g.samples<-grep(G[x],Metadata$Group,value=FALSE,perl=TRUE)
        g.log2.PSMs.Norm.mean<-rowMeans(as.matrix(log(FM$PSMs.Norm[,g.samples],base=2)))
        return(g.log2.PSMs.Norm.mean)
      },simplify=TRUE)
      colnames(log2.PSMs.Norm.means)<-paste(G,"mean[log2(Norm PSM+1)]")
      rownames(log2.PSMs.Norm.means)<-rownames(FM$PSMs.filtered)
      FM$log2.PSMs.Norm.means<-log2.PSMs.Norm.means

  # NORMALIZED SPECTRA ABUNDANCE FACTOR (NSAF) Analysis ----------------------- 
  message("- Calculating Normalized Spectral Abundance Factor (NSAF)")
    SAF<-FM$PSMs.Norm/filt.length
    NSAF<-sapply(1:ncol(SAF),function(x){
        return(SAF[,x]/sum(SAF[,x]))
      },simplify=TRUE)
    colnames(NSAF)<-paste(colnames(SAF),"NSAF")
    lnNSAF<-log(NSAF)
    colnames(lnNSAF)<-paste(colnames(SAF),"lnNSAF")
    NSAFfull<-cbind(NSAF,lnNSAF)
    rownames(NSAFfull)<-rownames(FM$PSMs.filtered)
    FM$NSAF<-NSAFfull
    ## Fitting to Normal distribution ----------------------- 
    message("- Fitting NSAF values to Normal distribution")
      NDE<-sapply(1:length(G),function(x){
          g.samples<-grep(G[x],Metadata$Group,value=FALSE,perl=TRUE)
          gLnNSAF<-lnNSAF[,g.samples]
          fitlnNSAFx<-fitdistr(gLnNSAF,densfun = "normal")
          return(c(fitlnNSAFx$estimate[1],fitlnNSAFx$estimate[2]))
        },simplify=TRUE)
      NDE<-t(NDE)
      rownames(NDE)<-G
      colnames(NDE)<-c("Mean NSAF fitted to Norm. Dist.","SD NSAF fitted to Norm. Dist.")
      FM$NDE<-NDE
    ## Getting means of the lnNSAF for each group ------------
      lnNSAF.means<-sapply(1:length(G),function(x){
        g.samples<-grep(G[x],Metadata$Group,value=FALSE,perl=TRUE)
        g.lnNSAF.mean<-rowMeans(as.matrix(lnNSAF[,g.samples]))
        return(g.lnNSAF.mean)
      },simplify=TRUE)
      colnames(lnNSAF.means)<-paste(G,"mean(lnNSAF)")
      rownames(lnNSAF.means)<-rownames(FM$PSMs.filtered)
      FM$lnNSAF.means<-lnNSAF.means
    ## Final NSAF table --------------------------------------


  # COMPARATIVE ANALYSIS AMONG GROUPS ------------------------------------------------
  # Groups will be compared against selected control/reference group (first group by default)
  ## Three different Abundance analysis to perform: 
    #### based on log2(Norm adjusted PSM)
    #### based on Abundance Ratio (R)
    #### based on lnNSAF
  # ----------------------------------------------------------------------------------   
  if(group.comparison=="TRUE"){
    message("- Comparison Analysis Between Groups")
    ## Getting samples and tables for control group
    control.name<-groups.vector[control]
    c.samples<-grep(control.name,Metadata$Group,value=FALSE,perl=TRUE)
    ## NOTE: Below use of as.matrix is to avoid error when ncol = 1
    c.NSAFtable<-as.matrix(NSAF[,c.samples])
    c.lnNSAFtable<-as.matrix(lnNSAF[,c.samples])
    c.log2.norm.PSM<-as.matrix(log(FM$PSMs.Norm[,c.samples],base=2))
    c.log2.norm.PSM.mean<-rowMeans(c.log2.norm.PSM)

    aG<-groups.vector[-control] # Groups to analyze against control
    message(paste("-- Control group: ",groups.vector[control]))
    message(c("-- Comparing groups: ",paste(aG,collpase=" "))) 
    FCEstimates<-data.frame(matrix(ncol=2,nrow=length(aG)))
    FCEstimates[control,1]<-NA
    FCEstimates[control,2]<-NA

    for(g in 1:length(aG)){
      message(paste("-- Analyzing group: ",aG[g],"against ",control.name))
      ## Getting samples and tables for each group to compare
      g.samples<-grep(aG[g],Metadata$Group,value=FALSE,perl=TRUE)
      g.pos<-grep(aG[g],G,value=FALSE,perl=TRUE)
      ## NOTE: Below use of as.matrix is to avoid error when ncol = 1
      g.NSAFtable<-as.matrix(NSAF[,g.samples])     
      g.lnNSAFtable<-as.matrix(lnNSAF[,g.samples])
      g.log2.norm.PSM<-as.matrix(log(FM$PSMs.Norm[,g.samples],base=2))
      g.log2.norm.PSM.mean<-rowMeans(g.log2.norm.PSM)
        ### Abundances are taken accordingly to what is selected by group.analysis
        ### Group lnNSAF is always taken by the mean of the samples
        ### Mean of Summations are NOT performed in case a single sample per group
        if(length(c.samples)>1){
          c.NSAF<-rowMeans(c.NSAFtable)
          c.lnNSAF<-rowMeans(c.lnNSAFtable)
          if(group.analysis=="sum"){
            c.PSMNorm<-rowSums(FM$PSMs.Norm[,c.samples])
          }else{
            c.PSMNorm<-rowMeans(FM$PSMs.Norm[,c.samples])
          }
        }else{
          message(paste("NOTE: Samples in",control.name,"<=1"))
          c.PSMNorm<-FM$PSMs.Norm[,c.samples]
          c.NSAF<-c.NSAFtable          
          c.lnNSAF<-c.lnNSAFtable
        }
        if(length(g.samples)>1){
          g.NSAF<-rowMeans(g.NSAFtable)          
          g.lnNSAF<-rowMeans(g.lnNSAFtable)
          if(group.analysis=="sum"){
            g.PSMNorm<-rowSums(FM$PSMs.Norm[,g.samples])          
          }else{
            g.PSMNorm<-rowSums(FM$PSMs.Norm[,g.samples])
          }
        }else{
          message(paste("NOTE: Samples in",aG[g],"<=1"))
          g.PSMNorm<-FM$PSMs.Norm[,g.samples]
          g.NSAF<-g.NSAFtable          
          g.lnNSAF<-g.lnNSAFtable
        }

        #### Declaring subroutines for statistical tests ----
          ##### T-test 
          t.test.fun<-function(TB1,TB2,TB1means,TB2means){
            n1<-sapply(1:nrow(as.matrix(TB1)),function(x){sum(!is.na(TB1[x,]))}) # nr of cases for each protein on each group
            n2<-sapply(1:nrow(as.matrix(TB2)),function(x){sum(!is.na(TB2[x,]))}) # nr of cases for each protein on each group
            SSE1<-sapply(1:nrow(as.matrix(TB1)),function(x){sum((as.numeric(na.omit(TB1[x,]))-TB1means[x])^2)})
            SSE2<-sapply(1:nrow(as.matrix(TB2)),function(x){sum((as.numeric(na.omit(TB2[x,]))-TB2means[x])^2)})
            SSE<-SSE1+SSE2
            Sp2<-SSE/(n1+n2-2) # Pooled variance of the two groups (n1+n2-2 degrees of freedom)
            ESE<-sqrt((Sp2/n1)+(Sp2/n2))
            t.test<-(TB1means-TB2means)/ESE
              p.values<-sapply(1:length(t.test),function(x){
                if(t.test[x]<0){
                  2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=TRUE)                
                }else{
                  2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=FALSE)
                }
              })
              p.values.adj<-as.numeric(p.adjust(p.values),method="BH")
              results<-cbind(p.values,p.values.adj)
              colnames(results)<-c("t-test p-values","T-test adj p-values (FDR)")
              return(results)           
          }
          moderated.t.test.fun<-function(TB1,TB2){
            Dat<-cbind(TB1,TB2)
            Dat[is.na(Dat)]<-0
            group<-c(rep("A",ncol(TB1)),rep("B",ncol(TB2)))
            factor<-factor(group)
            des<-model.matrix(~factor)
            fit<-lmFit(Dat,des)
            fit2<-eBayes(fit)
            res<-topTable(fit2,number=nrow(Dat))
            moderated.p.values<-res[,4]
            moderated.p.values.adj<-res[,5]
            p.vals<-cbind(moderated.p.values,moderated.p.values.adj)
            colnames(p.vals)<-c("moderated t-test (p)","moderated t-test (FDR)")
            return(p.vals)
          }
          #### Complementary error function
          # erfc<-function(x) 2*pnorm(-sqrt(2)*x,lower.tail=FALSE)
          erfc.fun<-function(x,mean,sd) 2*pnorm(-sqrt(2)*x,mean=mean,sd=sd,lower.tail=TRUE)
      ## log2(Norm Adjusted PSMs) analysis -----------------------
        message("-- Proceeding with log2(Norm Adjusted PSMs) analysis")    
        ### Getting log2 fold change
          logFC.PSMs<-c.log2.norm.PSM.mean-g.log2.norm.PSM.mean
          ### Getting t-test
            if(length(c.samples)<=1 | length(g.samples)<=1){
              message("NOTE: t-test won't be calculated due to number of samples for group or control <=1")
              log2NormPSManalysis<-cbind(logFC.PSMs)
              colnames(log2NormPSManalysis)<-c(paste("log2 Norm PSMs FC(",control.name,"-",aG[g],sep=""))
            }else{
              #### NOTE: assuming equal variances and different groups size          
                t.test<-t.test.fun(c.log2.norm.PSM,g.log2.norm.PSM,c.log2.norm.PSM.mean,g.log2.norm.PSM.mean)
                mod.t.test<-moderated.t.test.fun(c.log2.norm.PSM,g.log2.norm.PSM)
              log2NormPSManalysis<-cbind(logFC.PSMs,t.test,mod.t.test)
              colnames(log2NormPSManalysis)<-c(paste("log2 Norm PSMs FC(",control.name,"-",aG[g],sep=""),paste(control.name,"vs",aG[g],"t-test (p) log 2 Norm PSM",sep=""),paste(control.name,"vs",aG[g],"t-test (FDR) log 2 Norm PSM",sep=""),paste(control.name,"vs",aG[g],"moderated t-test (p) log 2 Norm PSM",sep=""),paste(control.name,"vs",aG[g],"moderated t-test (FDR) log 2 Norm PSM",sep=""))
            }
            FM$log2.PSMs.Norm.analysis<-cbind(FM$log2.PSMs.Norm.analysis,log2NormPSManalysis)
      ## Abundance Ratio Analysis --------------------------------
        message("-- Proceeding with Abundance Ratio analysis")   
        Rcg<-c.PSMNorm/g.PSMNorm # Control over corresponding group
        Rgc<-g.PSMNorm/c.PSMNorm # Group over Control
        Gtest<-as.numeric(2*(Rcg*log(Rcg/((Rcg+Rgc)/2))+Rgc*log(Rgc/((Rcg+Rgc)/2))) ) # G-test
        pval<-as.numeric(pchisq(Gtest,df=1,lower.tail = FALSE)) # G-test p-value by chi-squared
        ### False Discovery Rate (FDR) of the G-test p-value by Benjamin-Hochberg method
        pvalFDR<-as.numeric(p.adjust(pval),method="BH")
        ### Complementary error function over adjusted normal distribution for log2 protein abundance ratio (R)
        #### Getting corrected log2(Abundance ratio Control:Group) over normal distribution mean
          l2Rcg<-log(Rcg,base=2)
          fitCG<-fitdistr(l2Rcg,densfun="normal")
          mCG<-fitCG$estimate[1] # Mean
          sdCG<-fitCG$estimate[2] # Standard Deviation
          Rcg.corrected<-(l2Rcg-mCG)/(sdCG/sqrt(2))
        #### Getting corrected log2(Abundance ratio Group:Control) over normal distribution mean
          l2Rgc<-log(Rgc,base=2)
          fitGC<-fitdistr(l2Rgc,densfun="normal")
          mGC<-fitGC$estimate[1] # Mean
          sdGC<-fitGC$estimate[2] # Standard Deviation
          Rgc.corrected<-(l2Rgc-mGC)/(sdGC/sqrt(2))

          RcgERFC<-erfc.fun(abs(Rcg.corrected),mCG,sdCG)
          erfcFDR<-as.numeric(p.adjust(RcgERFC),method="BH")
        #### Creating final Abundance Matrix
        ABUNDANCE<-matrix(c(Rcg.corrected,Rgc.corrected,Gtest,pval,pvalFDR,RcgERFC,erfcFDR),nc=7,dimnames = list(rownames(FM$PSMs.filtered),c(paste("log2(R) ",control.name,":",aG[g]," corrected",sep=""),paste("log2(R) ",aG[g],":",control.name," corrected",sep=""),"G-test","G-test p-value","adjusted G-test (FDR)","ERFC","adjusted ERFC (FDR)")))
        FM$Abundance.Ratio.analysis<-cbind(FM$Abundance.Ratio.analysis,ABUNDANCE)
      ## NSAF ANALYSIS -------------------------------------------
        message("-- Proceeding with NSAF analysis")   
        ### NSAF fold-change (FC) -------------------
          # fold change (on subtraction format because the use of logarithms)
          lnFC<-c.lnNSAF-g.lnNSAF 
          log2FC<-log(c.NSAF,base=2)-log(g.NSAF,base=2)
        #### Fitting FC to normal distribution -----
          fitlnFC<-fitdistr(lnFC,densfun="normal")  
          FCEstimates[g.pos,1]<-fitlnFC$estimate[1]
          FCEstimates[g.pos,2]<-fitlnFC$estimate[2]
        ### NSAF and lnNSAF t-test --------------------------
          if(length(c.samples)<=1 | length(g.samples)<=1){
            message("NOTE: t-test won't be calculated due to number of samples for group or control <=1")
            log2NormPSManalysis<-cbind(logFC.PSMs)
            colnames(log2NormPSManalysis)<-c(paste("log2 Norm PSMs FC(",control.name,"-",aG[g],sep=""))
            STATslnNSAF<-as.matrix(lnFC)
            STATsNSAF<-as.matrix(log2FC)            
            rownames(STATslnNSAF)<-rownames(FM$PSMs.filtered)
            colnames(STATslnNSAF)<-paste(control.name,"/",aG[g]," FC",sep="")
            rownames(STATsNSAF)<-rownames(FM$PSMs.filtered)
            colnames(STATsNSAF)<-paste(control.name,"/",aG[g]," FC",sep="")            
          }else{
            #### T-test for each protein among two different groups ##
            #### NOTE: assuming equal variances and different groups size
            t.test.NSAF<-t.test.fun(c.NSAFtable,g.NSAFtable,c.NSAF,g.NSAF)
            t.test.lnNSAF<-t.test.fun(c.lnNSAFtable,g.lnNSAFtable,c.lnNSAF,g.lnNSAF)
            ### lnNSAF Moderated T-test (Limma package)------------
            mod.t.test.NSAF<-moderated.t.test.fun(c.NSAFtable,g.NSAFtable)
            mod.t.test.lnNSAF<-moderated.t.test.fun(c.lnNSAFtable,g.lnNSAFtable)
            STATsNSAF<-matrix(c(log2FC,t.test.NSAF[,1],t.test.NSAF[,2],mod.t.test.NSAF[,1],mod.t.test.NSAF[,2]),nc=5,dimnames=list(rownames(FM$PSMs.filtered),c(paste(control.name,"/",aG[g]," NSAF log2(FC)",sep=""),"t-Test (p)","t-Test (FDR)","Moderated t-Test (p)","Moderated t-Test (FDR)")))            
            STATslnNSAF<-matrix(c(lnFC,t.test.lnNSAF[,1],t.test.lnNSAF[,2],mod.t.test.lnNSAF[,1],mod.t.test.lnNSAF[,2]),nc=5,dimnames=list(rownames(FM$PSMs.filtered),c(paste(control.name,"/",aG[g],"lnNSAF FC",sep=""),"t-Test (p)","t-Test (FDR)","Moderated t-Test (p)","Moderated t-Test (FDR)")))
          }
        FM$NSAF.analysis<-cbind(FM$NSAF.analysis,STATsNSAF)         
        FM$lnNSAF.analysis<-cbind(FM$lnNSAF.analysis,STATslnNSAF)
    }
        colnames(FCEstimates)<-c(paste("Mean fitted lnNSAF Fold Change (over ",control.name,")",sep=""),paste("SD fitted lnNSAF Fold Change (over ",control.name,")",sep=""))
        FM$NDE<-cbind(FM$NDE,FCEstimates)
  }else{
    message("- No comparative analysis")
  }

  message("------------------------------------------------------------------------")
  message("Both NSAF and lnNSAF have bee analysed. Remember that NSAF is used mostly")
  message("in volcano plot with the log2 of its FC, whilst lnNSAF is used more often")
  message("in QQplots with its ratio/sd(ratio) values.")   
  message("Methods for the analysis have been taking from: Swearingen et al. PLOS Neglected Tropical Diseases (2017)")
  message("https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0005791")
  message("------------------------------------------------------------------------") 
  return(FM)

}
