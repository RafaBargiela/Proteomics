## FUNTION TO CREATE QUANTILE-QUANTILE plots from different expression data using FC and FC/sd ##################

# Alog: Abundances matrix, usually in logarithmns. NOTE: First column will correspond to top pink side of the diagram. Ensure you are using the same columns order than in the matrix returned by FC.Stats function.
  # FCsd:  Matrix with the log2(fold-change)/SD per element (protein,gene,transcript)
  # SD:  Standard deviation of the log2(FC) for all elements
  # MostExpLabs : Show labs for elements most expressed above a threshold
  # MostExpThreshold : Threshold to retrieved Most Expressed elements (genes, transcripts, proteins)
   # MostExpThresholdLabs : Threshold to show labs for Most Expressed elements (genes, transcripts, proteins). Usually, same than MostExpThreshold, but increase in case threshold shows to much labels on the graph.
  # y.legsep: Percentage of separataion on y-axis for Color legend

QQPlot<-function(Alog,FCsd,SD,sdrange=c(1,2,3),q.xlim=c(-5,5),q.ylim=c(-5,5),x.tick=seq(q.xlim[1],q.xlim[2],2.5),y.tick=seq(q.ylim[1],q.ylim[2],2.5),cex.axis=1.2,x.axis.labs=TRUE,y.axis.labs=TRUE,QXlabCorr=rep(0.05,nrow(FCsd)),QYlabCorr=rep(0,nrow(FCsd)),highlight=NULL,points.labs=rownames(Alog),points.cex=1.5,Qpch=21,Qcol=hcl.colors(length(sdrange)+1,"plasma"),legend=FALSE,leglabs=NULL,leglabs.cex=1.2,legpoints.cex=1.2,y.legsep=4,LegdotsLabs="",labcex=1,over.x.corr=0,over.y.corr=0,under.x.corr=0,under.y.corr=0,MostExpProt=FALSE,MostExpThreshold=2,MostExpThresholdLabs=MostExpThreshold,yl.corr=0,xl.corr=0,overRepLab="Overrepresented",underRepLab="Underrepresented",xlab="",ylab=""){
  
  # CHECKINGS ----------------------------------------------------------------------------
  if(ncol(Alog)!=(ncol(FCsd)+1)){
    stop(paste("ERROR: number of columns in Abundance Matrix must be equal to number of columns in Fold Change matrix +1: ",ncol(Alog),"!=",ncol(FCsd),sep=""))
   }
  # END CHECKINGS ------------------------------------------------------------------------

  require(TeachingDemos)
  plot(NA,NA,xlim=q.xlim,ylim=q.ylim,axes=FALSE,xlab=xlab,ylab=ylab)
    axis(1,x.tick,labels=x.axis.labs,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=cex.axis,las=1,lwd=2)
    axis(2,y.tick,labels=y.axis.labs,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=cex.axis,las=1,lwd=2)
    polygon(c(q.xlim[1],q.xlim[2],q.xlim[1]),c(q.ylim[1],q.ylim[2],q.ylim[2]),border=NA,col=rgb(t(col2rgb("thistle2")),alpha=100,maxColorValue=255),xpd=TRUE)
  polygon(c(q.xlim[1],q.xlim[2],q.xlim[2]),c(q.ylim[1],q.ylim[1],q.ylim[2]),border=NA,col=rgb(t(col2rgb("grey90")),alpha=100,maxColorValue=255),xpd=TRUE)
  text(q.xlim[2]+over.x.corr,q.ylim[2]+over.y.corr,labels=overRepLab,adj=c(1,1),cex=1.2,font=2,xpd=TRUE,col="mediumorchid3")
  text(q.xlim[2]+under.x.corr,q.ylim[1]-0.1+under.y.corr,labels=underRepLab,adj=c(1,0),cex=1.2,font=2,xpd=TRUE,col="grey50")                              
    A<-Alog[,1]
    if(ncol(Alog)>2){
      B<-sapply(2:ncol(Alog),function(x){
        v<-vector("numeric")
        v<-append(v,Alog[,x],length(v))
        v
      })
    }else{
      B<-Alog[,2]
    }

    # Calculating quantiles ----------------------------------------------------
    # Q1<-Alog2[round(seq(0.1,0.9,0.1)*(length(Alog2[,1])-1)),1]
    Q1<-quantile(A,seq(0.1,0.9,0.1))
    # Q2<-B[round(seq(0.1,0.9,0.1)*(length(B)-1))]
    Q2<-quantile(B,seq(0.1,0.9,0.1))
    arrows(q.xlim[1]-0.5,Q2,q.xlim[2]+0.5,Q2,lwd=1,lty=3,code=0,col="grey50")
    arrows(Q1,q.ylim[1]-0.5,Q1,q.ylim[2]+0.5,lwd=1,lty=3,code=0,col="grey50")
    text(quantile(A,.10),q.ylim[1]-0.25,labels="10%",font=2,cex=0.9,adj=c(1,0.5))
    text(quantile(A,.9),q.ylim[1]-0.25,labels="90%",font=2,cex=0.9,adj=c(0,0.5))
    text(q.xlim[2],quantile(B,.10),labels="10%",font=2,cex=0.9,adj=c(0.5,1))
    text(q.xlim[2],quantile(B,.9),labels="90%",font=2,cex=0.9,adj=c(0.5,0))
    # Printing dots-------- ----------------------------------------------------
    MostExpDF<-data.frame()
    MostExp<-vector("character")
    for(r in 1:nrow(Alog)){
      for(c in 2:ncol(Alog)){ # Regards that it starts from the second column
        n<-sum(abs(FCsd[r,c-1])>=sdrange)
        col<-n+1
        points(Alog[r,1],Alog[r,c],pch=Qpch[c-1],cex=points.cex,lwd=2,col=Qcol[col],bg=rgb(t(col2rgb(Qcol[col])),alpha=100,maxColorValue=255))
          if(is.list(QXlabCorr)==TRUE){
            xc<-QXlabCorr[[c-1]][r]
          }else{
            xc<-QXlabCorr[r]
          }
          if(is.list(QYlabCorr)==TRUE){
            yc<-QYlabCorr[[c-1]][r]
          }else{
            yc<-QYlabCorr[r]
          }
        if(MostExpProt==TRUE){
            if(abs(FCsd[r,c-1])>=MostExpThreshold){   # Printing labels for most expressed proteins
              if(abs(FCsd[r,c-1])>=MostExpThresholdLabs){
                adj<-c(0,0.5)
                if(Alog[r,1]<Alog[r,c]){      # Most expressed proteins on B are adjusted to the left
                  adj<-c(1,0.5)
                  xc<-xc*(-1)
                }
                text(Alog[r,1]+xc,Alog[r,c]+yc,labels=points.labs[r],font=4,cex=labcex,adj=adj,xpd=TRUE)
              }
            MostExp<-append(MostExp,rownames(Alog)[r],length(MostExp))
          } 
        }         
      }
    }
    # Highlighting specific proteins
    if(length(highlight)>0){
      for(e in 1:length(highlight)){
        for(c in 2:ncol(Alog)){
          if(abs(FCsd[highlight[e],c-1])<SD){
              col<-1
            }else{
              col<-(sum(sdrange<=abs(FCsd[highlight[e],c-1])))+1
            }
          points(Alog[highlight[e],1],Alog[highlight[e],c],pch=Qpch[c-1],col="white",bg="transparent",cex=points.cex+0.2,lwd=2)         
          points(Alog2[highlight[e],1],Alog2[highlight[e],c],pch=Qpch[c-1],col=Qcol[col],bg=rgb(t(col2rgb(Qcol[col])),alpha=100,maxColorValue=255),cex=points.cex,lwd=2)

          p<-grep(highlight[e],rownames(Alog))
          if(is.list(QXlabCorr)==TRUE){
              xc<-QXlabCorr[[c-1]][p]
            }else{
              xc<-QXlabCorr[p]
            }
              if(is.list(QYlabCorr)==TRUE){
              yc<-QYlabCorr[[c-1]][p]
            }else{
              yc<-QYlabCorr[p]
            } 
          if(FCsd[highlight[e],c-1]<0){
            adj<-c(1,0.5)
            xc<-xc*(-1)
          }else{
            adj<-c(0,0.5)
          }
          if(length(points.labs)==length(highlight)){
            # print(highlight[e])
            shadowtext(Alog2[highlight[e],1]+xc,Alog2[highlight[e],c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=4,xpd=TRUE,col="black",bg="white",r=0.05)      
          }else{
            # print(points.labs[highlight[e]])
            shadowtext(Alog[p,1]+xc,Alog[p,c]+yc,labels=points.labs[highlight[e]],adj=adj,cex=1,font=4,xpd=TRUE,col="black",bg="white",r=0.05)              
          }
        }
      }
    }


    # Legend
      if(legend==TRUE){
        yl<-q.ylim[2]+yl.corr
        yw<-(abs(q.ylim[1])+abs(q.ylim[2]))*y.legsep/100
        xw<-((abs(q.xlim[1])-0.2)+abs(q.xlim[2]))*2/100
        if(is.null(leglabs)==TRUE){
          leglabs<-vector("character",length(sdrange)+1)
          for(r in 1:(length(sdrange)+1)){
            if(r==1){
              leglabs[r]<-"FC < \u03c3"
            }else{
              if(r==2){
                leglabs[r]<-"FC > \u03c3"               
              }else{
                leglabs[r]<-paste("FC > ",r-1,"\u03c3",sep="")
              }
            }
          }
          # print(leglabs)
        }
        for(c in 1:length(Qcol)){
          points(q.xlim[1]-0.2,yl,pch=19,col=Qcol[c],cex=legpoints.cex)
          text(q.xlim[1]+xw-0.2,yl,labels=leglabs[c],font=2,cex=leglabs.cex,xpd=TRUE,adj=c(0,0.5))
          yl<-yl-yw
        }
        if(length(Qpch)>1){
        xi<-q.xlim[1]+xl.corr
        y<-yl-0.1
        yw<-(y-2)/length(LegdotsLabs)
        for(p in 1:length(Qpch)){
          points(xi,y,pch=Qpch[p],cex=2,lwd=2,bg="white",col="black",xpd=TRUE)
          text(xi+0.3,y,labels=LegdotsLabs[p],font=2,cex=1.2,xpd=TRUE,adj=c(0,0.5))
          y<-y-yw
        }   
        }     
      }
    box(lwd=2)
    MostExM<-FCsd[unique(MostExp),]
    return(MostExM)
}
