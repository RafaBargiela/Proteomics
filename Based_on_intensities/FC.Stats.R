# Author: Rafael Bargiela, PhD. Instituto de Catálisis y Petroleoquímica, ICP (CSIC). Last update: 4th of November, 2025

# Fold change calculation and t.test statistics
################################################
  # This function used the transformed, filtered, normalized and also imputed data returned from functions l2.filter.Norm and Imputation.
  # Here we get the fold change among two designated groups and calculate statistical significant differences among each or the proteins, based on 
  # t-test and moderate t-test based on limma package

    # M.Imputed: Previously filtered and normalized matrix (usually log-2 transformed). Usually, data was also previously imputed.
    # groups.vector: Vector assigning groups to each of the samples in matrix (columns). Only TWO different groups allowed.
    ## WARNING: Regard that groups.vector is not equal to the vector used on previous functions
    # LongAnalysisM : Data.frame containing Intensity, Time, Individuals and Group columns to considered in case on longitudinal analysis. Default NULL
    # method : Method for random effect measure and removal in case of longitudinal analysis. Default limma. WARNING: lme4 is not ready yet
FC.Stats<-function(M.filt.Norm.Imputed,Metadata,groups.vector,LongAnalysis=FALSE,method="limma"){
  message("- NOTE: Make sure cols in M are in the same order than rows in Metadata -")
    message("NOTE: Make sure groups vector to compare exists in Metadata$Group -")
  M<-M.filt.Norm.Imputed # Previously filtered and normalized matrix (usually log-2 transformed)
  FC.Stats<-{}
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
        block<-Metadata[,"Group"]
        dup_corrs <- duplicateCorrelation(M,design = design, block = block)
        mdl.fit<-Humanzee::ruv_mixed_model(M,ndups = 1,design = design, block = block,correlation = dup_corrs$cons)
        M.fit <- t( design %*% t(mdl.fit$coef) ) + mdl.fit$resid
        FC.Stats$Mixed.Model.Fit.limma<-M.fit
        message("Using Imputed data after random effect removal")
        M<-FC.Stats$Mixed.Model.Fit.limma
      }else{
        stop("Wrong method selection for Longitudinal Analysis")
      }
  }else{
    message("Using Imputed data without previous random effect removal")
  }

  G<-unique(groups.vector)
  g1sams<-grep(paste("^",G[1],"$",sep=""),Metadata[,"Group"])
  print(paste("Group 1 is ",G[1],": ",paste(colnames(M)[g1sams],collapse=" "),sep=""))
  G1<-M[,g1sams]
  g2sams<-grep(paste("^",G[2],"$",sep=""),Metadata[,"Group"])
  print(paste("Group 2 is ",G[2],": ",paste(colnames(M)[g2sams],collapse=" "),sep=""))   
  G2<-M[,g2sams]
  # G1<-M[,grep(G[1],groups.vector,value=FALSE)]
  # G2<-M[,grep(G[2],groups.vector,value=FALSE)]    
  G1av<-sapply(1:nrow(G1),function(x){mean(G1[x,],na.rm = TRUE)})
  G2av<-sapply(1:nrow(G2),function(x){mean(G2[x,],na.rm = TRUE)})
  names(G1av)<-rownames(G1)
  names(G2av)<-rownames(G2)
  FC<-G1av-G2av
  # T-test for each protein among two different groups ##
  # NOTE: assuming equal variances and different groups size
  n1<-sapply(1:nrow(G1),function(x){sum(!is.na(G1[x,]))}) # nr of cases for each protein on each group
  n2<-sapply(1:nrow(G2),function(x){sum(!is.na(G2[x,]))})  # If there aren't NA values n1 and n2 are equal for all proteins
  
  SSE1<-sapply(1:nrow(G1),function(x){sum((as.numeric(na.omit(G1[x,]))-G1av[x])^2)})
  SSE2<-sapply(1:nrow(G2),function(x){sum(as.numeric((na.omit(G2[x,]))-G2av[x])^2)})
  SSE<-SSE1+SSE2
  Sp2<-SSE/(n1+n2-2) # Pooled variance of the two groups (n1+n2-2 degrees of freedom)
  ESE<-sqrt((Sp2/n1)+(Sp2/n2))
  t.test<-FC/ESE
  p.values<-sapply(1:length(t.test),function(x){
    if(t.test[x]<0){
      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=TRUE)                
    }else{
      2*pt(q=t.test[x],df=(n1[x]+n2[x]-2),lower.tail=FALSE)
    }
  })

  p.values.adj<-p.adjust(p.values,method="BH")

  # Limma method
  require(limma)
  Dat<-data.frame(M[,c(g1sams,g2sams)])
  rownames(Dat)<-rownames(M)
  Dat[is.na(Dat)]<-0
  # groups<-factor(groups.vector)
  groups<-factor(Metadata[c(g1sams,g2sams),"Group"])
  des<-model.matrix(~groups)
  fit<-lmFit(Dat,des)
  fit2<-eBayes(fit)
  res<-topTable(fit2,number=nrow(Dat))
  res<-res[rownames(Dat),]
  moderated.p.values<-res[,4]
  moderated.p.values.adj<-res[,5]
  FC.Stats$DEA<-matrix(c(G1av,G2av,FC,p.values,p.values.adj,moderated.p.values,moderated.p.values.adj),nc=7,dimnames=list(rownames(M),c(paste(G[1],"log2 mean"),paste(G[2],"log2 mean"),"Fold change","t.test p-values","t-test adj. p-values","moderated-t.test p-values","moderated-t.test adj. p-values")))
  


  return(FC.Stats) 
}
