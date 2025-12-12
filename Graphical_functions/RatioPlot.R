## RATIO PLOT  function to create graphics based on differential expression of proteins using p-values and anbundace ratio between groups ###########33

# PSMs : Numeric vector with total PSMs per protein among both samples
  # R :    Numeric vector with Abundance Ratio A:B per protein
  # labels: Logical, if TRUE labels for proteins with both adjusted p-values <0.05 with be shown.

### Ratio plot using PSMs abundance ratio (R) and statistical test p-values
RatioPSMplot<-function(PSMs,R,xlim,ylim, x.tick=seq(xlim[1],xlim[2],2),y.tick=seq(ylim[1],ylim[2],2), cols=c("skyblue3","goldenrod3","grey50"),
                       labels=FALSE,labPSMlim=5, labXcorr=rep(0.15,length(R)), labYcorr=rep(0,length(R)),labAdj=matrix(rep(c(0,0.5),length(Rab)),nc=2,byrow=TRUE), ... ){
  
  R<-R[names(sort(abs(R)))]
  plot(NA,NA,xlim=xlim,ylim=ylim,axes=FALSE, ... )
  axis(1,x.tick,labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
  axis(2,y.tick,labels=TRUE,tick=TRUE,lwd.ticks=2,font.axis=2,cex.axis=1.2,las=1,lwd=2)
  arrows(x.tick,ylim[1]-0.1,x.tick,ylim[2]+0.1,code=0,lwd=1,lty=3,col="grey80") 
  arrows(xlim[1]-0.1,y.tick,xlim[2]+0.1,y.tick,code=0,lwd=1,lty=3,col="grey80") 
  box(lwd=2)
  
  if(is.null(names(labXcorr))) {names(labXcorr)<-names(R)}
  if(is.null(names(labYcorr))) {names(labYcorr)<-names(R)} 
  
  oe<-vector("character") # Vector to store OverExpressed proteins
  for(p in 1:length(R)){
    c<-3
    if(FDRs[names(R)[p],1]<=0.05 && FDRs[names(R)[p],2]<=0.05) {c<-1}
    if(FDRs[names(R)[p],1]<=0.05 && FDRs[names(R)[p],2]>0.05) {c<-2}
    xp<-PSMs[names(R)[p]]
    yp1<-R[p]
    points(xp,yp1,pch=19,cex=1.1,col=cols[c], ... )
    if(c==1 && labels==TRUE && xp>labPSMlim){
      xa<-labXcorr[names(R)[p]]
      if(labAdj[1]==1){
        xa<--xa
      }
      text(xp+xa,yp1+labYcorr[names(R)[p]],labels=names(R)[p],adj=labAdj[p], ... )
      oe<-append(oe,names(R)[p],length(oe))
    }
  }
  
  # Adding legend
  xl<-xlim[1]
  yl<-ylim[2]
  yw<-(abs(ylim[1])+abs(ylim[2]))*5/100
  xw<-(abs(xlim[1])+abs(xlim[2]))*2.5/100
  llabs<-c("FDR < 0.05","ERFC FDR > 0.05","G-test & ERFC FDR > 0.05")
  for(c in 1:length(cols)){
    points(xl,yl,pch=19,col=cols[c],cex=1.2, ... )
    text(xl+xw,yl,labels=llabs[c],adj=c(0,0.5),cex=1, ... )
    yl<-yl-yw
  }
  
  return(oe)
}
