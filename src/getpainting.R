source("FinestructureLibrary.R")
if(!require("fdrtool")){
    install.packages("fdrtool")
}
if(!require("poibin")){
    install.packages("poibin")
}
if(!require("gap")){
    install.packages("gap")
}

require("fdrtool")
require("poibin")
require("gap")
if(!exists("loadrds")){ ## load the results of a previous run to speed things up
  print("Recreating results. If you want to restore them instead, try \"loadrds<-TRUE\"")
  loadrds<-FALSE
}
if(!exists("ploidy")){ # Check ploidy
  print("Assuming ploidy=2, specify via \"ploidy<-1\" to change")
  ploidy<-2
}
if(capabilities()[["X11"]]){ ## allow creation of PNG without X11
  mypng<-png
}else{
  require("Cairo")
  mypng<-CairoPNG
}


source("paintingfns.R")

####################
# Read the data

## params only needed for GT analyses
#paramfile<-"wolfdog_hybrid.paramfile.txt"
#paramlist<-readParam(paramfile)

## Location of the data lists

########################################
## Read the data

if(!exists("alllist")){
    alllist<-readPaintings(recomblist,copyprobslist,mypops=mypops,
                           indsperfile=indsperfile,ploidy=ploidy,loadeach=loadrds)
}

# And this is how we reload from the save RDS objects
# alllist<-readPaintings(recomblist,copyprobslist,indsperfile=1,ploidy=2,loadeach=TRUE)

if(!exists("chromosomegap")) chromosomegap<-5e5 # how far we treat different chromosomes as being

## Get the mean painting; this is useful as it provides the reference info for the whole set of SNPs
if(!loadrds){
  meanpaintingall<-getMeanPainting(alllist,pops=mypops,chromosomegap = chromosomegap,naends=0)
  saveRDS(meanpaintingall,file=paste0(rootname,"_meanpainting.RDS"))
}else{
  meanpaintingall<-readRDS(file=paste0(rootname,"_meanpainting.RDS"))
}
cgaps<-which(diff(meanpaintingall[,1])>0)
cgappos<-(meanpaintingall[cgaps,2]+meanpaintingall[cgaps+1,2])/2

chromotextpos<-sapply(1:max(meanpaintingall[,"chromo"]),function(x){
  trange<-range(meanpaintingall[meanpaintingall[,"chromo"]==x,"pos"])
  mean(trange)
})

## This is how we would extract the probability that all hidden states were identical and 1
rowprod<-function(x){
  apply(x,1,prod)
}
## This is how we would extract the probability that all hidden states were identical and 0
rowprod0<-function(x){
  apply(1-x,1,prod)
}

if(loadrds){
  prodpaintingall<-readRDS(file=paste0(rootname,"_prodpaintingall.RDS"))
  prodpaintingall0<-readRDS(file=paste0(rootname,"_prodpaintingall0.RDS"))
}else{
  prodpaintingall<-getMeanPainting(alllist,pops=mypops,chromosomegap = chromosomegap,naends=0,fn=rowprod)
  prodpaintingall0<-getMeanPainting(alllist,pops=mypops,chromosomegap = chromosomegap,naends=0,fn=rowprod0)
  saveRDS(prodpaintingall,file=paste0(rootname,"_prodpaintingall.RDS"))
  saveRDS(prodpaintingall0,file=paste0(rootname,"_prodpaintingall0.RDS"))
}
pvalprod<-prodpaintingall
pvalprod0<-prodpaintingall0
pvalprod[,-(1:3)]<-1 -pvalprod[,-(1:3)]
pvalprod0[,-(1:3)]<-1 -pvalprod0[,-(1:3)]

#apply(prodpaintingall0[,-(1:3)],2,max)

pvalprodwindows<-getSnpRangesPval(pvalprod,rep(0.05,length(popnames)),snpwindow = 200)
pvalprod0windows<-getSnpRangesPval(pvalprod0,rep(0.05,length(popnames)),snpwindow = 200)

#####
if(!loadrds){
  writeSignifRegions(pvalprod0windows,rootname = paste0(rootname,"_exclusion_"))
  writeSignifRegions(pvalprodwindows,rootname = paste0(rootname,"_fixation_"))

#####
  for(adon in 1:2) {
    if(adon==1){
      td<-pvalprod
      tp<-pvalprodwindows
      trootname<-paste0(rootname,"_GenomeWideFixation")
    }else{
      td<-pvalprod0
      tp<-pvalprod0windows
      trootname<-paste0(rootname,"_GenomeWideExclusion")
    }
#    for(i in which(apply(td[,-(1:3)],2,min)<0.05)) {
    for(i in 1:(dim(td)[2]-3)) {
      jj<-i
      mypng(file=paste0(trootname,popnames[i],"_fdr.png"),height=500,width=1600)
      par(mfrow=c(1,1),cex=1.2)
      plot(td[,2],-log10(td[,jj+3]),type="n",frame.plot=F,xlab="Genome Position",ylab="-log10(Probability)",main=colnames(td)[jj+3])
      tmin<-1
      if(!is.null(tp[[jj]]$minpval)) tmin<-min(tp[[jj]]$minpval)
      mtext(paste("threshold p<0.05 retaining",length(tp[[jj]]$snps),"SNPs over",dim(tp[[jj]][["snpranges"]])[1],"regions, with minimum prob",format(tmin,digits=3)))
      abline(v=cgappos,col="grey",lwd=2)
      abline(h=-log10(0.05),lwd=2)
      lines(td[,2],-log10(td[,jj+3]))
      text(chromotextpos, rep(max(-log10(td[,jj+3])),max(td[,"chromo"])),
           labels=1:max(td[,"chromo"]),adj=c(0.5,1))
      dev.off()
    }
  }  
}

###############################################
########## poisson-binomial approach
source("paintingfns.R")
if(loadrds){
  signiflist<-readRDS(file=paste0(rootname,"_signiflist.RDS"))
}else{
  signiflist<-getPvalMatrix(alllist,names(mypops),naends=0) # 
  saveRDS(signiflist,file=paste0(rootname,"_signiflist.RDS"))
}
#signifmat<-signiflist$signifmat

####### qq plots
if(!loadrds){
  require("gap")
  mypng(file=paste0(rootname,"_qq.png"),height=800,width=400*dim(signiflist$signifmat)[2]/2)
  par(mfrow=c(2,dim(signiflist$signifmat)[2]/2),cex=1.5)
  for(a in 1:2){
    for(i in 1:(dim(signiflist$signifmat)[2]/2)) {
      ii<-(i-1)*2+a
      qqunif(signiflist$signifmat[,ii],main=dimnames(signiflist$signifmat)[[2]][ii])
    }
  }
  dev.off()
}

######### mean painting
if(!loadrds){
#######################
  for(i in 1:npops){
    mypng(file=paste0(rootname,"_MeanPainting",popnames[i],".png"),height=600,width=1600)
    par(mfrow=c(1,1),cex=1.2,mar=c(4,5,5,1))
      jj<-i
      plot(meanpaintingall[,2],meanpaintingall[,jj+3],type="n",frame.plot=F,xlab="Genome Position",ylab="Mean Painting",main=colnames(meanpaintingall)[jj+3])
      abline(v=cgappos,col="grey",lwd=2)
      lines(meanpaintingall[,2],meanpaintingall[,jj+3])
      text(chromotextpos, rep(max(-log(meanpaintingall[,jj+3])),max(meanpaintingall[,"chromo"])),
           labels=1:max(meanpaintingall[,"chromo"]),adj=c(0.5,1))
    dev.off()
  }


  for(i in 1:npops){
    mypng(file=paste0(rootname,"_Distribution",popnames[i],".png"),height=600,width=800)
    par(mfrow=c(1,1),cex=1.2,mar=c(4,5,5,1))
    jj<-i
    meanx<-signiflist$p0[,i]
    tx<-seq(0,1,length.out = length(meanx))
    tx<-(0:length(meanx))#/length(meanx)
    dx<-diff(tx)[1]/length(meanx)
    ty<-dpoibin(tx,meanx)/dx
    thist<-hist(meanpaintingall[,jj+3],breaks=tx/length(meanx),plot=FALSE)
    plot(thist$mids/dx,thist$density,xlab="Value",ylab="Density (log scale)",ylim=c(0.0001,max(c(ty,thist$density))),main=colnames(meanpaintingall)[jj+3],type="h",log="y")
    if(sum(tx*ty/length(meanx))/sum(ty)>0.5){tloc="topleft"
    }else{tloc="topright"     }
    legend(tloc,legend=c("Data","Poisson-Binomial Fit"),col=1:2,text.col=1:2,lty=1)
    lines(tx,ty,col=2)
    dev.off()
  }

  mypng(file=paste0(rootname,"_ecdf.png"),height=500,width=500)
  plot(ecdf(signiflist$p0[,1]),xlim=c(0,1),main="Empirical CDF")
  for(i in 2:dim(signiflist$p0)[2]) lines(ecdf(signiflist$p0[,i]),col=i)
  legend("left",legend=popnames,col=1:length(popnames),text.col=1:length(popnames),lty=1,pch=19)
  dev.off()
  
}

## Let x be the probability of the SNP in question
## and p0 be the genome wide probability
## then k~PoiBin(x)
## and k0+PoiBin(p0)
## We want p_high=prob(k>k0) for extreme HIGH values
## and p_low=prob(k<k0) = 1-prob(k>=k0) for extreme LOW values
# prob(k>k0) = prob(k=k' & k0<k') =prob(k=k') prob(k0<k')
#            = sum_{k'=0}^N ppoibin(k',x) dpoibin(k'-1,p0)
# prob(k>=k0) = prob(k=k' & k0<=k') =prob(k=k') prob(k0<=k')
#            = sum_{k'=0}^N ppoibin(k',x) dpoibin(k',p0)
library("fdrtool")
source("paintingfns.R")
if(loadrds){
  fdrlist<-readRDS(file=paste0(rootname,"_fdrlist.RDS"))
}else{
  fdrlist<-getlfdr(signiflist$signifmat,fdrmax=0.2)
  saveRDS(fdrlist,file=paste0(rootname,"_fdrlist.RDS"))
}
fdrvals<-fdrlist$fdrvals
# fdrlist$n50 # how many significant snps we find

########################
## Now we have a list of significant SNPs we might even believe in
# Make a nice version of the significance matrix
pvalsall<-cbind(meanpaintingall[,1:3],signiflist$signifmat)
source("paintingfns.R")
if(loadrds){
  pvalsnpwindows<-readRDS(file=paste0(rootname,"_pvalsnpwindows.RDS"))
}else{
  pvalsnpwindows<-getSnpRangesPval(pvalsall,fdrlist$thresh,snpwindow = 200)
  saveRDS(pvalsnpwindows,file=paste0(rootname,"_pvalsnpwindows.RDS"))
}

nsignifregions<-sapply(pvalsnpwindows,function(x){dim(x[["snpranges"]])[1]} )

gwsigthresh<-rep(0.05/dim(na.omit(signiflist$signifmat))[1],dim(signiflist$signifmat)[2])

gweffective<-mean(apply(signiflist$signifmat,2,function(x,myd=100){
  x<-na.omit(x)
#  tseq<-seq(1,ceiling(length(x)/myd-1)*myd+1,by=myd)
#  ttmat<-x[tseq]
#  tacf<-acf(ttmat,lag.max=10000)
  ttacf<-acf(x,lag.max = 60,plot=FALSE)
#  tmin<-min(which(tacf$acf<0.05))
#  teff<-tmin*myd
#  (gwsigthresh[1]*teff)
  gwsigthresh[1]*(1+ttacf$acf[2])/(1-ttacf$acf[2])
}))


source("paintingfns.R")
getgwregions<-function(gwpvalsnpwindows){
  lapply(gwpvalsnpwindows,function(x){x$minpval})
}

###################
## PUT THIS INSIDE AN IF *************88
if(!loadrds){
#if(TRUE) {
  gweffpvalsnpwindows<-getSnpRangesPval(pvalsall,rep(gweffective,length(gwsigthresh)),snpwindow = 200)
#  gwefflist<-getgwr(signifmat,gweffective)
  gweffregions<-getgwregions(gweffpvalsnpwindows) 
  saveRDS(gweffective,file=paste0(rootname,"_gweffective.RDS"))
  saveRDS(gweffpvalsnpwindows,file=paste0(rootname,"_gweffpvalsnpwindows.RDS"))
  saveRDS(gweffregions,file=paste0(rootname,"_gweffregions.RDS"))
}else{
  gweffective<-readRDS(file=paste0(rootname,"_gweffective.RDS"))
  gweffpvalsnpwindows<-readRDS(file=paste0(rootname,"_gweffpvalsnpwindows.RDS"))
  gweffregions<-readRDS(file=paste0(rootname,"_gweffregions.RDS"))
}
##

ngweffsignifregions<-sapply(gweffpvalsnpwindows,function(x){dim(x[["snpranges"]])[1]} )
ngweffsignifsnps<-sapply(gweffpvalsnpwindows,function(x){length(x[["snps"]])[1] })
ngweffpvals<-sapply(gweffpvalsnpwindows,function(x){
    if(is.null(x[["minpval"]])) return(1); 
  min(x[["minpval"]][1])} )


##########
if(!loadrds){
#if(TRUE) {
gwpvalsnpwindows<-getSnpRangesPval(pvalsall,gwsigthresh,snpwindow = 200)
#  gwlist<-getgwr(signifmat,gwective)
  gwregions<-getgwregions(gwpvalsnpwindows) 
  saveRDS(gwsigthresh,file=paste0(rootname,"_gwsigthresh.RDS"))
  saveRDS(gwpvalsnpwindows,file=paste0(rootname,"_gwpvalsnpwindows.RDS"))
  saveRDS(gwregions,file=paste0(rootname,"_gwregions.RDS"))
}else{
  gwsigthresh<-readRDS(file=paste0(rootname,"_gwsigthresh.RDS"))
  gwpvalsnpwindows<-readRDS(file=paste0(rootname,"_gwpvalsnpwindows.RDS"))
  gwregions<-readRDS(file=paste0(rootname,"_gwregions.RDS"))
}
##

ngwsignifregions<-sapply(gwpvalsnpwindows,function(x){dim(x[["snpranges"]])[1]} )
ngwsignifsnps<-sapply(gwpvalsnpwindows,function(x){length(x[["snps"]])[1] })
ngwpvals<-sapply(gwpvalsnpwindows,function(x){
    if(is.null(x[["minpval"]])) return(1); 
  min(x[["minpval"]][1])} )


if(!loadrds){
  fdrregions<-getgwregions(pvalsnpwindows)
  saveRDS(fdrregions,file=paste0(rootname,"_fdrregions.RDS"))
}else{
  fdrregions<-readRDS(file=paste0(rootname,"_fdrregions.RDS"))
}

if(!loadrds){
#if(TRUE) {
#######################
  for(i in 1:npops){
    mypng(file=paste0(rootname,"_GenomeWidePvals",popnames[i],"_fdr.png"),height=1000,width=1600)
    par(mfrow=c(2,1),cex=1.2,mar=c(4,5,5,1))
    for(j in 1:2){
      jj<-(i-1)*2+j
      plot(pvalsall[,2],-log10(pvalsall[,jj+3]),type="n",frame.plot=F,xlab="Genome Position",ylab="-log10(Probability)",main="")
      mtext(colnames(pvalsall)[jj+3],line=3)
      mtext(paste("FDR threshold p<",format(fdrlist$thresh[jj],digits=3),"retaining",fdrlist$n50[jj],"SNPs over",nsignifregions[jj],"regions, with lfdr in",format(fdrlist$minlfdr[jj],digits=3),":",format(fdrlist$maxlfdr[jj],digits=3)))
      mtext(paste("GW threshold p<",format(gwsigthresh[jj],digits=3),"retaining",ngwsignifsnps[jj],"SNPs over",ngwsignifregions[jj],"regions, with min pval ",format(ngwpvals[jj],digits=3)),line=1)
      mtext(paste("GWeff threshold p<",format(gweffective,digits=3),"retaining",ngweffsignifsnps[jj],"SNPs over",ngweffsignifregions[jj],"regions, with min pval ",format(ngweffpvals[jj],digits=3)),line=2)

      
      abline(v=cgappos,col="grey",lwd=2)
      abline(h=-log10(fdrlist$thresh[jj]),lwd=2)
      abline(h=-log10(gwsigthresh[jj]),lwd=2,col=2)
      abline(h=-log10(gweffective),lwd=2,col=3)
      lines(pvalsall[,2],-log10(pvalsall[,jj+3]))
      text(chromotextpos, rep(max(-log(pvalsall[,jj+3])),max(pvalsall[,"chromo"])),
           labels=1:max(pvalsall[,"chromo"]),adj=c(0.5,1))
    }
    dev.off()
  }
  
  writeSignifRegions(pvalsnpwindows,fdrregions = fdrregions ,rootname = paste0(rootname,"_poibinranges_"))
  writeSignifRegions(gwpvalsnpwindows,fdrregions = gwregions ,rootname = paste0(rootname,"_gwranges_"))
  writeSignifRegions(gweffpvalsnpwindows,fdrregions = gweffregions ,rootname = paste0(rootname,"_gweffranges_"))
}



# save.image(paste0(rootname,"_getpainting.RData"))


######################################
###
print("Trying to create painting figures - this requires FinestructureLibrary.R")

some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-c(some.colors,rgb(0,0,seq(1-0.025,0,-0.025)))

nhaps<-dim(alllist[[1]][[1]])[2]
ninds<-nhaps/ploidy
nchr<-length(alllist)
nsnpschr<-sapply(alllist,function(x){dim(x[[1]])[1]})

if(ploidy==1) {print("Haploid mode")
}else if(ploidy==2) {print("Diploid mode")
}else stop("Incorrect ploidy!")

for(mypop in 1:length(mypops)) {
  print(paste("Processing population",names(mypops)[mypop]))
  tcmmat<-sapply(alllist,function(x){colMeans(x[[mypop]])})
  tcm<-colMeans(tcmmat)
  indmat<-tcmmat[seq(1,ninds)*ploidy-(ploidy-1),]
  if(ploidy==2) indmat<-indmat + tcmmat[seq(1,ninds)*ploidy,]
  indmat<-t(indmat/ploidy)
  colnames(indmat)<-paste0("IND",1:ninds)
  rownames(indmat)<-paste0("CHR",1:nchr)
  if(ploidy==1) colnames(indmat)<-paste0("HAP",1:ninds)
  
  tind<-colMeans(indmat)
  tindorder<-order(tind)
  tchromorder<-rep(ploidy*tindorder,each=ploidy)
  if(ploidy==2)tchromorder<-tchromorder +c(-1,0)
  
  tpopfrac<-colMeans(indmat[,tindorder])

  pdf(paste0(rootname,"_ChromosomeMeanPainting",names(mypops)[mypop],".pdf"),height=10,width=15)
  plotFinestructure(indmat[,tindorder],labelsx=colnames(indmat)[tindorder],labelsy=rownames(indmat),labelsatx=1:ninds,labelsaty=1:nchr,cols=some.colorsEnd,ignorebelow = 0,main=paste("Mean",names(mypops)[mypop],"Painting averaged over each chromosome"),layoutd = 0.1,layoutf=0.1,cex.axis=1,colscale=c(0,1),scalenum=11)
  title(ylab="Chromosome",xlab=paste0("Individual (sorted by ",names(mypops)[mypop]," Fraction)"),line=7)
  axis(3,1:ninds,format(tpopfrac,digits=2),las=2)
  axis(4,1:nchr,format(apply(indmat[,tindorder],1,mean),digits=2),las=2)
  dev.off()


  pdf(paste0(rootname,"_ChromosomeECDF_",names(mypops)[mypop],".pdf"),height=6,width=10)
  for(i in 1:nchr){
#    tmp<-indmat[i,tindorder]/max(indmat[i,])
    tmp<-indmat[i,tindorder]
    if(i==1){
      plot(ecdf(tmp),xlab=paste(names(mypops)[mypop],"Fraction"),ylab=paste("CDF of",names(mypops)[mypop],"Fraction"),main=paste("ECDF of Chromosomes"),xlim=c(0,1.2))
    }
    else lines(ecdf(tmp),col=i,lty=i,pch=i)
  }
  legend("right",legend=paste("Chr",1:nchr),col=1:nchr,lty=1:nchr,pch=1:nchr,text.col=1:nchr,bty="n")
  abline(v=c(0,1))
  dev.off()

  sindmat<-indmat[,tindorder]
  tcor<-matrix(1,nrow=dim(sindmat)[2],ncol=dim(sindmat)[2])
  for(i in 1:(dim(sindmat)[2])-1) for(j in (1+1):dim(sindmat)[2]){
    tcor[i,j]<-tcor[j,i]<-cor(sindmat[,i],sindmat[,j])
  }

  pdf(paste0(rootname,"_ChromosomePaintingCorrelation_",names(mypops)[mypop],".pdf"),height=10,width=12)
  plotFinestructure(tcor,labelsx=colnames(sindmat),labelsatx=1:ninds,cols=some.colorsEnd,ignorebelow = -1,main=paste("Correlation between",names(mypops)[mypop]," chromosome paintings"),layoutd = 0.05,layoutf=0.1,cex.axis=1)
  title(xlab=paste0("Individual (sorted by ",names(mypops)[mypop]," Fraction)"),ylab=paste0("Individual (sorted by ",names(mypops)[mypop]," Fraction)"),line=6)
  dev.off()

#################
  ## Still needs bringing up to standards
  pdir<-paste0(rootname,"Paintings",names(mypops)[mypop])
  system(paste("mkdir -p",pdir),intern=TRUE) ## WON'T WORK ON WINDOWS

  tmatrange<-ceiling(max(c(
    -log10(min(signiflist$signifmat[,(mypop-1)*2+1:2])),
    -log10(gwsigthresh),
    -log10(fdrlist$thresh[fdrlist$thresh>0]) ))
    )
  
  for(i in 1:length(alllist)) {
    print(paste("Chromosome",i))
    tx<-round(seq(1,dim(alllist[[i]][[mypop]])[1],length.out=min(10000,dim(alllist[[i]][[1]])[1])))
    
    tmat<-signiflist$signifmat[meanpaintingall[,1]==i,(mypop-1)*2+1:2]

    mypng(paste0(pdir,"/",names(mypops)[mypop],"Painting_Chr",i,".png"),height=1200,width=2000)

    layout(matrix(c(2,1,3,1),nrow=2,ncol=2,byrow=TRUE),heights=c(0.75,0.25),widths = c(0.9,0.1))
    
    ## scale
    par(mar=c(5,5,5,3))
    colindex<-t(matrix(seq(0,1,length.out=100),ncol=1,nrow=100)) # colour scale
    image(1,1:100,colindex,xaxt="n",yaxt="n",xlab="",ylab="",col=some.colors,zlim=c(0,1))
    scalelocs <- seq(0,1,by=0.2)
    scalephysicalpos <- seq(1,100,length.out=6)
    axis(2,at=scalephysicalpos,labels=scalelocs,las=2,cex.axis=2)
    
    par(mar=c(3,5,5,1))
    image(tx,1:dim(alllist[[i]][[mypop]])[2], as.matrix(alllist[[i]][[mypop]][tx,tchromorder]),xlab="",ylab="",main=paste0("Painting for ",names(mypops)[mypop]," on Chromosome ",i,": Yellow=LOW, Blue=HIGH"),axes=F,col=some.colors,zlim=c(0,1),cex.main=1.5)
    axis(1,round(seq(1,dim(alllist[[i]][[mypop]])[1],length.out=11)),cex.axis=1.5)
    axis(2,2*(0:(ninds-1))+1.5,1:ninds,las=1,labels=colnames(indmat)[tindorder])
    if(ploidy==2) abline(h=2*(1:(ninds-1))+0.5)

    
    par(xaxs="i",mar=c(5,5,0,1))
    plot(range(tx),c(0,tmatrange),type="n",xlab=paste0("SNP number, Red=low, Blue=High"),ylab="-log10(Pval)",axes=F,cex.lab=1.5)
    axis(1,round(seq(1,dim(alllist[[i]][[mypop]])[1],length.out=11)),cex.axis=1.5)
#    axis(2,at=seq(1,20,length.out = 5),
    axis(2,at=seq(1,tmatrange,length.out = 5),
         labels=format(seq(0,tmatrange,length.out=5),digits=2),las=1)
#    for(j in 1:2) lines(tx, -log10(tmat[tx,j])/max(-log10(tmat[,1:2]))*20 ,col=c("red","blue")[j])
    abline(h=-log10(fdrlist$thresh[(mypop-1)*2+1:2]),lwd=2,lty=2,col=c("red","blue"))
    abline(h=-log10(gwsigthresh[(mypop-1)*2+1]),lwd=2,lty=2)
    abline(h=-log10(gweffective),lwd=2,lty=2,col="green")
    par(xpd=NA)
    text(x=dim(alllist[[i]][[mypop]])[1],y=-log10(gwsigthresh[(mypop-1)*2+1]),"GW",adj=0)
    text(x=dim(alllist[[i]][[mypop]])[1],y=-log10(gweffective),"GWeff",adj=0)
    text(x=dim(alllist[[i]][[mypop]])[1],y=-log10(fdrlist$thresh[(mypop-1)*2+1]),"FDR",adj=0,col="red")
    text(x=dim(alllist[[i]][[mypop]])[1],y=-log10(fdrlist$thresh[(mypop-1)*2+2]),"FDR",adj=0,col="blue")
    par(xpd=FALSE)
    for(j in 1:2) lines(tx, -log10(tmat[tx,j]),col=c("red","blue")[j])

    dev.off()
  }

}
#########################################



