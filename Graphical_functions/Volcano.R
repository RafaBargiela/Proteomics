## FUNCTION TO CREATE VOLCANO PLOTS FROM PROTEOMICS DATA ANALYSIS #####

# FC : log2 FC matrix.
  # logPV : -log 10 of the adjusted p-values.
  # legend.cols: Logical. If TRUE, it shows the legend corresponding to each color. Usually related with FC and p-value threshold
  # x.legcol and y.legcol: Initial x and y point where to draw the legend.cols
  # x.legcollabsep: Separation between dots and text in legend.cols
  # y.legcolabbottom: y bottom point for the legend, been y.legcol the top.
  # legcol.cex: size for legend.cols labels
  # legcol.points.cex: size for legend.cols dots

volcano<-function(FC,logPv,FCthreshold=1,PVthreshold=2,labFCthreshold=NA,lablogPvthreshold=NA,vol.xlim=c(-5,5),vol.ylim=c(0,5),xlim.ticks=seq(vol.xlim[1],vol.xlim[2],2.5),volpch=21,points.labs=rownames(FC),highlight=NULL,points.cex=1.5,points.labs.cex=1,volcol=hcl.colors(sum(length(FCthreshold),length(PVthreshold))+1,"plasma"),legend.cols=TRUE,x.legcol=vol.xlim[1]+(abs(vol.xlim[1])*0.01),y.legcol=vol.ylim[2]-(abs(vol.ylim[2])*0.01),x.legcolabsep=((abs(vol.xlim[1])+abs(vol.xlim[2]))*0.02),y.legcolabbottom=(vol.ylim[2]-vol.ylim[1])/1.5,legcol.cex=1,legcol.points.cex=points.cex,volXlabCorr=rep(0.05,nrow(FC)),volYlabCorr=rep(0,nrow(FC)),legend.dots=FALSE,LegdotsLabs="",xlab=list("log\u2082(Fold change)",cex=1.2,font=2),ylab=list("-log\u2081\u2080(p-value)",cex=1.2,font=2),overRepLab="Overrepresented",x.overRepLab=vol.xlim[1],y.overRepLab=vol.ylim[1],overReplab.cex=1.2,underRepLab.cex=1.2,x.underRepLab=vol.xlim[2],y.underRepLab=vol.ylim[1],underRepLab="Underrepresented",xlegcol=vol.xlim[1]){
  

  # Chechings -----------------------------------------------------------------
    require(TeachingDemos)
    if(length(volcol)!=sum(length(FCthreshold),length(PVthreshold))+1){
      warning("- Number of colors for dots different than sum of thresholds divisions\n",paste("Colors:",length(volcol),"\n"),paste("Thresholds:",sum(length(FCthreshold),length(PVthreshold))+1,"\n"))
    }
    if(length(volXlabCorr)!=nrow(FC)){
      warning(paste("Nr of  X-axis corrections",length(volXlabCorr),"for labels must be equal to row names of FC",nrow(FC)))
    }
    if(length(volYlabCorr)!=nrow(FC)){
      warning(paste("Nr of  X-axis corrections",length(volYlabCorr),"for labels must be equal to row names of FC",nrow(FC)))
    }
    if(legend.dots==TRUE & length(volpch)==1){
      warning(paste("Only 1 dot type specified (volpch)\n","Samples to compare against the reference sample (first column): ",ncol(FC),"\n",sep=""))
    }       
  # End checkings -------------------------------------------------------------
    MostExpDF<-data.frame()
    MostExp<-vector("character")
    plot(NA,NA,xlim=vol.xlim,ylim=vol.ylim,ylab=ylab,xlab=xlab,axes=FALSE)
    axis(1,xlim.ticks,labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    axis(2,seq(vol.ylim[1],vol.ylim[2],1),labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
    arrows(c(-FCthreshold,FCthreshold),vol.ylim[1],c(-FCthreshold,FCthreshold),vol.ylim[2],code=0,lty=3,lwd=2,col="grey50")
    arrows(vol.xlim[1],PVthreshold,vol.xlim[2],PVthreshold,code=0,lty=3,lwd=2,col="grey50")
    rect(0,0,vol.xlim[1],vol.ylim[2],border=NA,col=rgb(t(col2rgb("thistle2")),alpha=100,maxColorValue=255),xpd=TRUE)
    text(x.overRepLab,y.overRepLab,labels=overRepLab,adj=c(0,0),cex=overReplab.cex,font=2,xpd=TRUE,col="mediumorchid3")
    rect(0,0,vol.xlim[2],vol.ylim[2],border=NA,col=rgb(t(col2rgb("grey90")),alpha=100,maxColorValue=255),xpd=TRUE)
    text(x.underRepLab,y.underRepLab,labels=underRepLab,adj=c(1,0),cex=overReplab.cex,font=2,xpd=TRUE,col="grey50")            
    for(r in 1:nrow(FC)){
      if(length(highlight)>0){
          if(sum(grepl(rownames(FC)[r],highlight))>0){
            message(paste(rownames(FC)[r],"will be highlighted"))
            next
          }else{

          }
      }
      #
      for(c in 1:ncol(FC)){
        n<-1
        topFC<-sum(abs(FC[r,c])>=FCthreshold,na.rm=TRUE)
        topPV<-sum(abs(logPv[r,c])>=PVthreshold,na.rm=TRUE)
        if(topFC==0 & topPV==0){
          n<-1
        }else{
          if((topFC>0 & topPV==0) | (topFC==0 & topPV>0)){
            n<-n+1
          }else{
            n<-n+topFC+topPV
          }
        } 

        points(FC[r,c],logPv[r,c],pch=volpch[c],col=volcol[n],bg=rgb(t(col2rgb(volcol[n])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)
        if(n==3){
          # In case we want to limit the labels to different thresholds
          if(is.na(labFCthreshold)==FALSE){
            if(labFCthreshold>abs(FC[r,c])){
              lab.plot<-FALSE
            }else{
              lab.plot<-TRUE
            }
          }else{
            lab.plot<-TRUE
          }
          if(is.na(lablogPvthreshold)==FALSE){
            if(lablogPvthreshold>logPv[r,c]){
              lab.plot<-FALSE
            }else{
              lab.plot<-TRUE
            }
          }else{
            lab.plot<-TRUE
          }       
          if(is.list(volXlabCorr)==TRUE){
              xc<-volXlabCorr[[c]][r]
            }else{
              xc<-volXlabCorr[r]
            }
              if(is.list(volYlabCorr)==TRUE){
              yc<-volYlabCorr[[c]][r]
            }else{
              yc<-volYlabCorr[r]
            } 
          if(FC[r,c]<0){
            adj<-c(1,0.5)
            xc<-xc*(-1)
          }else{
            adj<-c(0,0.5)
          }
          if(lab.plot==TRUE){
            text(FC[r,c]+xc,logPv[r,c]+yc,labels=points.labs[r],adj=adj,cex=points.labs.cex,font=4,xpd=TRUE)            
          }else{
            print(paste(points.labs[r],"label not printed"))
          }
          MostExp<-append(MostExp,rownames(FC)[r],length(MostExp))      
        }
      }
    }

    # Highlighting specific proteins
    if(length(highlight)>0){
      for(e in 1:length(highlight)){
        for(c in 1:ncol(FC)){
          n<-1
          topFC<-sum(abs(FC[highlight[e],c])>=FCthreshold)
          topPV<-sum(abs(logPv[highlight[e],c])>=PVthreshold)
          if(topFC==0 & topPV==0){
            n<-1
          }else{
            if((topFC>0 & topPV==0) | (topFC==0 & topPV>0)){
              n<-n+1
            }else{
              n<-n+topFC+topPV
            }
          } 
          points(FC[highlight[e],c],logPv[highlight[e],c],pch=volpch[c],col="white",bg="transparent",cex=points.cex+0.2,lwd=2)
          points(FC[highlight[e],c],logPv[highlight[e],c],pch=volpch[c],col=volcol[n],bg=rgb(t(col2rgb(volcol[n])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)         
          p<-grep(highlight[e],rownames(FC))
          if(is.list(volXlabCorr)==TRUE){
              xc<-volXlabCorr[[c]][p]
            }else{
              xc<-volXlabCorr[p]
            }
              if(is.list(volYlabCorr)==TRUE){
              yc<-volYlabCorr[[c]][p]
            }else{
              yc<-volYlabCorr[p]
            } 
          if(FC[highlight[e],c]<0){
            adj<-c(1,0.5)
            xc<-xc*(-1)
          }else{
            adj<-c(0,0.5)
          }
          if(length(points.labs)==length(highlight)){
            shadowtext(FC[highlight[e],c]+xc,logPv[highlight[e],c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=2,xpd=TRUE,col="black",bg="white",r=0.05) 

          }else{          
            shadowtext(FC[highlight[e],c]+xc,logPv[highlight[e],c]+yc,labels=points.labs[p],adj=adj,cex=1,font=2,xpd=TRUE,col="black",bg="white",r=0.05)              
          }
        }
      }
    }
    box(lwd=2)
    # Sample type Legend, when using more than one dot shape
      if(legend.dots==TRUE){
        xi<-vol.xlim[1]
        y<-vol.ylim[2]
        yw<-(y-(y-(y*0.17)))/length(LegdotsLabs)
        for(p in 1:length(volpch)){
          points(xi,y,pch=volpch[p],cex=2,lwd=2,bg="white",col="black",xpd=TRUE)
          text(xi+0.3,y,labels=LegdotsLabs[p],font=2,cex=1.2,xpd=TRUE,adj=c(0,0.5))
          y<-y-yw
        } 
      }
    # Thresholds Legend related to colors
      if(legend.cols==TRUE){
        legcolsLabs<-c(paste("log\u2082(FC)<",FCthreshold[1],", p>",10^(-PVthreshold[1]),sep=""),paste("log\u2082(FC)>",FCthreshold[1]," OR p<",10^(-PVthreshold[1]),sep=""),paste("log\u2082(FC)>",FCthreshold[1]," & p<",10^(-PVthreshold[1]),sep=""))
        FCleg<-vector("character",length(FCthreshold))
        PVleg<-vector("character",length(PVthreshold))
        FCleg[1]<-as.character(FCthreshold[1])
        PVleg[1]<-as.character(10^(-PVthreshold[1]))
        if(length(FCthreshold)>1){
          for(th in 2:length(FCthreshold)){
            FCleg[th]<-paste("log\u2082(FC)>",FCthreshold[th],sep="")
          }
        }
        if(length(PVthreshold)>1){
          for(th in 2:length(PVthreshold)){
            PVleg[th]<-paste("p<",10^(-PVthreshold[th]),sep="")
          }
        }
        FCPVleg<-paste(FCleg,"&",PVleg)
        if(length(FCPVleg)>1){
          legcolsLabs<-append(legcolsLabs,FCPVleg[2:length(FCPVleg)],length(legcolsLabs))
        }
        y.legcolabsep<-(y.legcol-y.legcolabbottom)/length(volcol)
        for(col in 1:length(volcol)){
          points(x.legcol,y.legcol,pch=volpch[1],cex=legcol.points.cex,col=volcol[col],bg=rgb(t(col2rgb(volcol[col])),alpha=100,maxColorValue=255),lwd=2)
          if(x.legcol<0){
            adj<-c(0,0.5)
            sep<-x.legcolabsep
          }else{
            adj<-c(1,0.5)
            sep<-x.legcolabsep*(-1)
          }
          text(x.legcol+sep,y.legcol,labels=legcolsLabs[col],adj=adj,font=2,cex=legcol.cex)
          y.legcol<-y.legcol-y.legcolabsep
        }
      }
    MostExM<-FC[unique(MostExp),]
    return(MostExM)
}
