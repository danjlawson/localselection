matColSums<-function(mat,poplist)
## Sums over columns in a matrix, by grouping all columns listed in poplist
## e.g. if mat is M*N matrix and poplist is length K, returns a M*K matrix
## the names of poplist are used to assign names to the returned matrix
{
        res<-matrix(0,nrow=dim(mat)[1],ncol=length(poplist))
        colindex<-lapply(poplist,function(x){which(colnames(mat)%in%x)})
        res<-t(apply(mat,1,function(x){
                sapply(colindex,function(y){sum(x[y])})
        }))
        colnames(res)<-names(poplist)
        rownames(res)<-rownames(mat)
        res
}

getPrA<-function(x){
  sum(x[-1]*pxAgivenX)
}

getSignifSNPs<-function(mysnps){
  mysnps<-c(mysnps,tail(mysnps,1))
  ret<-matrix(0,ncol=2,nrow=0)
  tdiff<-diff(mysnps)
  while(length(tdiff)>0){
    tnext<-(which(tdiff>1))[1]
    if(is.na(tnext)) tnext<-length(tdiff)
    ret<-rbind(ret,as.numeric(c(mysnps[1],mysnps[tnext])))
    tdiff<-tdiff[-(1:(tnext))]
    mysnps<-mysnps[-(1:(tnext))]
  }
  ret
}

which.is.max<-function(x){
  if(any(is.nan(x))) return(1)
  ret<-which.max(x)
  if(length(ret)>1) return(sample(ret,1))
  if(length(ret)<0) return(1)
  ret
}

getMonteCarloThresholdList<-function(alllist,nreps,element,cutoff=0.01,naends=0,verbose=TRUE,vverbose=FALSE,fn=rowMeans,...) {
counter<-1
myreps<-sapply(alllist,function(x) {
  if(verbose) print(paste("Performing monte carlo on chromosome",counter,"of",length(alllist)))
  counter<<-counter+1
  tx<<-x[[element]]
  if(naends>0) tx<-tx[-c(1:naends,dim(tx)[1]-(1:naends)+1),]
  tmp<-getMonteCarloThreshold(tx,nreps=nreps ,cutoff=cutoff,verbose=vverbose,fn=fn,...)
  return(c(tmp$min,tmp$max))
})
min=apply(myreps[1:nreps,,drop=FALSE],1,min)
max=apply(myreps[nreps+1:nreps,,drop=FALSE],1,max)
tthresh<-c(quantile(min,cutoff),quantile(max,cutoff))

return(list(thresh=tthresh,min=min,max=max))
}

getMonteCarloThreshold<-function(data,nreps,cutoff=0.01,aslist=TRUE,verbose=TRUE,fn=rowMeans,...){
  nsnpsrep<-dim(data)[1]
  mydistmax<-numeric(nreps)
  mydistmin<-numeric(nreps)
  for(rep in 1:nreps){
    if(verbose) print(paste("Running replicate",rep))
    trep<-data[1:nsnpsrep,]
    for(i in 1:dim(data)[2]) {
      tsample<-(1:nsnpsrep+sample(1:nsnpsrep,1))%%(nsnpsrep)+1
      trep[,i]<-trep[tsample,i]
    }
    trm<-fn(trep,...)
    mydistmax[rep]<-max(trm)
    mydistmin[rep]<-min(trm)
  }
  tthresh<-c(quantile(mydistmin,cutoff),quantile(mydistmax,cutoff))
  if(aslist){
    return(list(thresh=tthresh,min=mydistmin,max=mydistmax,sample=trm))
  }
  tthresh
}

spoibin<-function(kk, pp, method = "DFT-CF", wts = NULL){
  ## Discrete Survivor for the poisson binomial
  ## which is Pr(x >= kk|pp)
  1-ppoibin(kk,pp,method = method,wts=wts) + dpoibin(kk,pp,wts=wts)
}
  

getPvalMatrix<-function(alllist,popnames,naends=0,verbose=TRUE){
  require("poibin")

  tpp<-getP0(alllist,naends) # genome wide probabilities of painting from each pop
  tpppoibin_geq<-apply(tpp,2,function(x){spoibin(0:length(x),x)}) # prob that x>=K
  tpppoibin_leq<-apply(tpp,2,function(x){ppoibin(0:length(x),x)})  # prob that x<=K

  mypoibin<-function(x,arg) {
# probability of "arg"  weighted by the probability that each x occurs in the poisson binomial
    if(any(is.na(x))) return(NA)
    x[x<0]<-0
    x[x>1]<-1
    tmp<-dpoibin(0:length(x),x)* arg
    tmp[tmp<0]<-0
    tmp[tmp>1]<-1
    sum(tmp)
  }

  npops<-length(alllist[[1]])
  if(verbose) print("Obtaining p-values for low haplotype counts")
  poibin_low<-sapply(1:npops,function(x){
    getPvalList(alllist,element=x,
                fn=mypoibin,naends=naends,
                arg=tpppoibin_leq[,x])
  })
  if(verbose) print("Obtaining p-values for high haplotype counts")
  poibin_high<-sapply(1:npops,function(x){
    getPvalList(alllist,element=x,
                fn=mypoibin,naends=naends,
                arg=tpppoibin_geq[,x])
  })
  
## Extract as a matrix
  for(i in 1:npops){
    tx<-unlist(poibin_low[,i])
    tx_high<-unlist(poibin_high[,i])
    if(i==1) signifmat<-matrix(nrow=length(tx),ncol=2*npops)
    tx[tx<0]<-0
    tx[tx>1]<-1
    tx_high[tx_high<0]<-0
    tx_high[tx_high>1]<-1
    signifmat[,2*(i-1)+1]<-tx
    signifmat[,2*(i-1)+2]<-tx_high
  }
  colnames(signifmat)<-paste0(rep(popnames,each=2),rep(c("_low","_high"),times=npops))
  list(signifmat=signifmat,p0=tpp,tpppoibin_leq=tpppoibin_leq,
       tpppoibin_geq=tpppoibin_geq,fn=mypoibin)
}

getPvalList<-function(alllist,element,naends=0,verbose=TRUE,fn,arg) {
  counter<-1
  sapply(alllist,function(x) {
    if(verbose)print(paste("processing chromosome",counter,"(element",element,")"))
    counter<<-counter+1
    tx<-x[[element]]
#    if(naends>0) tx<-tx[-c(1:naends,dim(tx)[1]-(1:naends)+1),]
    if(naends>0) tx[c(1:naends,dim(tx)[1]-(1:naends)+1),]<-NA
    tmp<-getPval(tx,fn=fn,arg=arg)
    return(tmp)
  })
}

getPval<-function(tx,fn,...){
  apply(tx,1,fn,...)
}

NormaliseEntropy<-function(tentropy,tev=NULL){
  if(any(is.null(tev))) tev<-as.numeric(tentropy)
  tev<-na.omit(tev)
  (tentropy[]-  min(tev))/diff(range(tev))
}

readDonorTable<-function(donortabfile,recippop=NULL){
  ttab<-read.table(donortabfile,stringsAsFactors=FALSE)
  tnames<-unique(ttab[,2])
  if(!is.null(recippop)) tnames<-tnames[-which(tnames==recippop)]
  tlist<-lapply(tnames,function(x){
    ttab[ttab[,2]==x,1]
  })
  names(tlist)<-tnames
  tlist
}

readHaps<-function(tfiles,mypops,nhapstot,nsnps,ninds,ploidy=2,verbose=TRUE){
  nhapstot<-ploidy*ninds
  header<-read.table(tfiles[1],nrows=1,as.is=T)[1,]
  paintingmat<-matrix(0,nrow=nsnps,ncol=nhapstot)
  colnames(paintingmat)<-paste0("hap",1:nhapstot)
  paintinglist<-list()
  
  for(i in 1:length(mypops)) paintinglist[[i]]<-paintingmat

  hapon<-1
  for (tfon in 1:length(tfiles) ) {
    tf<-tfiles[tfon]
    if(verbose) print(paste("Computing population paintings in file number",tfon,"(",tf,")"))
    for(tindon in 1:indsperfile) {
      if(verbose&&indsperfile>1) print(paste("Computing population painting for individual ",tindon,"from that file"))
      hap1<-read.table(tf,skip=1+(nsnps+1)*(ploidy*tindon-ploidy)+1,nrows=nsnps,as.is=T)
      colnames(hap1)<-header
      if (ploidy==2) {
        hap2<-read.table(tf,skip=1+(nsnps+1)*(2*tindon-1)+1,nrows=nsnps,as.is=T)
        colnames(hap2)<-header
        painting2<-matColSums(hap2,mypops)
      }
      painting1<-matColSums(hap1,mypops)
      for(i in 1:dim(painting1)[2]){
        painting1[painting1[,i]>1,i]<-1
        painting1[painting1[,i]<0,i]<-0
        paintinglist[[i]][, hapon]<-painting1[,i]
        if (ploidy==2) {
          painting2[painting1[,i]>1,i]<-1
          painting2[painting1[,i]<0,i]<-0
          paintinglist[[i]][, hapon+1]<-painting2[,i]
        }
      }
      hapon<-hapon+ploidy
      if(tfon==1) for(i in 1:dim(painting1)[2]) rownames(paintinglist[[i]])<-hap1[,1]
    }
  }
  for(i in 1:length(paintinglist)){
    paintinglist[[i]]<-apply(paintinglist[[i]],2,rev)
  }
  paintinglist
}

getMeanPainting<-function(alllist,pops=NA,chromosomegap=100000,naends=0,fn=rowMeans){
  
  snpsperfile<-sapply(alllist,function(x){dim(x[[1]])[1]})
  meanpaintingfull<-matrix(nrow=0,ncol=length(alllist[[1]]))
  positionsfull<-numeric()
  positionsorig<-numeric()
  chromosomesfull<-numeric()
  
  tlastpos<-0
  for (recfileon in 1:length(alllist) ) {
    paintinglist1<-alllist[[recfileon]]
    tpositions<-as.numeric(rownames(paintinglist1[[1]]))
    meanpainting1<-lapply(paintinglist1,fn)
    meanpainting<-matrix(nrow=length(tpositions),ncol=length(meanpainting1))
    for(i in 1:length(meanpainting1)) meanpainting[,i]<-meanpainting1[[i]]
    if(naends>0) {
      meanpainting[1:naends,]<-NA
      meanpainting[dim(meanpainting)[1]-(1:naends)+1,]<-NA
    }
    meanpaintingfull<-rbind(meanpaintingfull,meanpainting)
    positionsfull<-c(positionsfull,tpositions+tlastpos)
    positionsorig<-c(positionsorig,tpositions)
    chromosomesfull<-c(chromosomesfull,rep(recfileon,length(tpositions)))
    tlastpos<-max(tpositions)+chromosomegap+tlastpos
  }
  ret<-cbind(chromosomesfull,positionsfull,positionsorig,meanpaintingfull)
  if(any(is.na(pops))) {
    colnames(ret)<-c("chromo","pos","pos0",paste0("pop",1:(dim(ret)[2]-1)))
  }else{
    colnames(ret)<-c("chromo","pos","pos0",names(pops))
  }
  ret
}



getSnpRangesPval<-function(pvalsall,tthresh,snpwindow=500){
  nppops<-dim(pvalsall)[2]-3
  ret<-list()
  if(length(tthresh)!=dim(pvalsall)[2]-3) stop("Must provide 1 p-value threshold per population column")
  for(popon in 1:nppops){
    ret[[popon]]<-list()
    tw<-pvalsall[,popon+3]<=tthresh[popon]
    if(any(na.omit(tw))) {
      ret[[popon]][["snps"]]<-which(tw)
      ret[[popon]][["snpranges"]]<-getSignifSNPs(ret[[popon]][["snps"]])
      ret[[popon]][["chrom"]]<-apply(ret[[popon]][["snpranges"]],1,function(x){
        pvalsall[x[1]]
      })
      ret[[popon]][["snpwindow"]]<-t(apply(ret[[popon]][["snpranges"]],1,function(x){
        tchrom<-pvalsall[x[1]]
        c(max(x[1]-snpwindow,min(which(pvalsall[,"chromo"]==tchrom))),min(x[2]+snpwindow,max(which(pvalsall[,"chromo"]==tchrom))))
      }))
      ret[[popon]][["psnpranges"]]<-matrix(apply(ret[[popon]][["snpranges"]],2,function(x){pvalsall[x,"pos"]}),ncol=2)
      ret[[popon]][["psnpwindow"]]<-matrix(apply(ret[[popon]][["snpwindow"]],2,function(x){pvalsall[x,"pos"]}),ncol=2)
      ret[[popon]][["p0snpranges"]]<-matrix(apply(ret[[popon]][["snpranges"]],2,function(x){pvalsall[x,"pos0"]}),ncol=2)
      ret[[popon]][["p0snpwindow"]]<-matrix(apply(ret[[popon]][["snpwindow"]],2,function(x){pvalsall[x,"pos0"]}),ncol=2)
      ret[[popon]][["minpval"]]<-apply(ret[[popon]][["snpranges"]],1,function(x){
        min(pvalsall[x[1]:x[2],popon+3])
      })
    }else{
      ret[[popon]][["p0snpwindow"]]<-ret[[popon]][["p0snpranges"]]<-ret[[popon]][["psnpwindow"]]<-ret[[popon]][["psnpranges"]]<-ret[[popon]][["snpwindow"]]<-ret[[popon]][["snpranges"]]<-matrix(nrow=0,ncol=2)
    }
  }
  names(ret)<-colnames(pvalsall)[-(1:3)]
  ret
}

getSnpRanges<-function(meanpaintingall,tthresh,snpwindow=500,type="max"){
  nppops<-dim(meanpaintingall)[2]-3
  ret<-list()
  if(class(tthresh)=="numeric") tthresh<-matrix(tthresh,ncol=2)
  for(popon in 1:nppops){
    ret[[popon]]<-list()
    if(type=="max") ret[[popon]][["snps"]]<-which(meanpaintingall[,popon+3]>tthresh[popon,2])
    else  ret[[popon]][["snps"]]<-which(meanpaintingall[,popon+3]<tthresh[popon,1])
    ret[[popon]][["snpranges"]]<-getSignifSNPs(ret[[popon]][["snps"]])
    ret[[popon]][["chrom"]]<-apply(ret[[popon]][["snpranges"]],1,function(x){
      meanpaintingall[x[1]]
    })
    ret[[popon]][["snpwindow"]]<-t(apply(ret[[popon]][["snpranges"]],1,function(x){
      tchrom<-meanpaintingall[x[1]]
      c(max(x[1]-snpwindow,min(which(meanpaintingall[,"chromo"]==tchrom))),min(x[2]+snpwindow,max(which(meanpaintingall[,"chromo"]==tchrom))))
    }))
    ret[[popon]][["psnpranges"]]<-apply(ret[[popon]][["snpranges"]],2,function(x){meanpaintingall[x,"pos"]})
    ret[[popon]][["psnpwindow"]]<-apply(ret[[popon]][["snpwindow"]],2,function(x){meanpaintingall[x,"pos"]})
    ret[[popon]][["p0snpranges"]]<-apply(ret[[popon]][["snpranges"]],2,function(x){meanpaintingall[x,"pos0"]})
    ret[[popon]][["p0snpwindow"]]<-apply(ret[[popon]][["snpwindow"]],2,function(x){meanpaintingall[x,"pos0"]})
  }
  ret
}

getHighFreq<-function(alllist,above=0.5){
  for(i in 1:length(alllist[[1]])) {
    tres<-unlist(lapply(alllist,function(x){
      rowMeans(x[[i]]>0.5)
    }))
    if(i==1) ret<-matrix(tres,ncol=1)
    else  ret<-cbind(ret,tres)
  }
  ret
}

myseq<-function(trange,length.out=10) {
  seq(trange[1],trange[2],length.out=length.out)
}

getBf<-function(x,alpha=0.05){ # Bonferoni
  return(alpha/length(x))
}

getBH<-function(x,alpha=0.05,ess=TRUE){ # returns a threshold below which we accept, above and equal to which we reject
  m<-length(x)
  ox<-order(x)
  if(ess) {
    require(coda)
    m<-effectiveSize(x)
  }
  
  ilist<-1:length(x)
  myw<-which(alpha*ilist/m-x[ox][ilist]>0)
  if(length(myw)==0) return(alpha/m)
  mythresh1<-max(myw)
  return((x[ox][mythresh1]+x[ox][mythresh1+1])/2)
}

getMeanProbabilites<-function(alllist,popnames=NULL) {
  ## Get the genome wide probabilities of being in each state, per haplotype
  snpsperfile<-sapply(alllist,function(x){dim(x[[1]])[1]})
  tcmslist<-lapply(1:length(alllist[[1]]),function(y){
    sapply(alllist,function(x){colMeans(x[[y]])})
  })
  ttcmslist<-lapply(tcmslist,function(x){
    t(snpsperfile*t(x)/sum(snpsperfile))
  })
  ret<-sapply(ttcmslist,rowSums)
  if(any(is.null(popnames))) popnames<-paste0("Pop",1:dim(ret)[2])
  colnames(ret)<-popnames
  ret
}

getpoibin<-function(tpp){
  tpppoibin_high<-apply(tpp,2,function(x){ppoibin((-1):(length(x)-1),x)})
  tpppoibin_low<-apply(tpp,2,function(x){ppoibin(0:length(x),x)})
  list(low=tpppoibin_low,high=tpppoibin_high)
}

mypoibin_testlow<-function(x,arg) { # prob that k>k0
  require("poibin")
  x[x>1]<-1
  x[x<0]<-0
  1-sum(dpoibin(0:length(x),x)* arg)
}
mypoibin_testhigh<-function(x,arg){
  require("poibin")
  x[x>1]<-1
  x[x<0]<-0
  (sum(dpoibin(0:length(x),x) * arg))
}

getpoibinPvals<-function(alllist,tpppoibin){
  poibin_low<-sapply(1:length(alllist[[1]]),function(x){
    getPvalList(alllist,element=x,
                fn=mypoibin_testlow,
                arg=tpppoibin[["low"]][,x]) 
  })
  poibin_high<-sapply(1:length(alllist[[1]]),function(x){
    getPvalList(alllist,element=x,
                fn=mypoibin_testhigh,
                arg=tpppoibin[["high"]][,x]) 
  })
  
  ## Extract as a matrix
  for(i in 1:dim(poibin_high)[2]){
    tx<-unlist(poibin_low[,i])
    tx_high<-unlist(poibin_high[,i])
    if(i==1) signifmat<-matrix(nrow=length(tx),ncol=2*dim(poibin_high)[2])
    signifmat[,2*(i-1)+1]<-tx
    signifmat[,2*(i-1)+2]<-tx_high
  }
  popnames<-colnames(tpppoibin[["high"]])
  colnames(signifmat)<-paste0(rep(popnames,each=2),c("_low","_high"))
  signifmat
}

##############################
readParam<-function(paramfile){
  ## Read a GT param file and its associated ID file
  paramdata<-read.table(paramfile,sep=":",as.is=T)
  paramdata[,2]<-sub(" ","",paramdata[,2])
  gtout<-paste0(paramdata[paramdata[,1]=="save.file.main",2],".txt")
  idfile<-paramdata[paramdata[,1]=="input.file.ids",2]
  copyvector_file<-paramdata[paramdata[,1]=="input.file.copyvectors",2]
  recippop<-paramdata[paramdata[,1]=="target.popname",2]
  donorpops<-strsplit(paramdata[paramdata[,1]=="copyvector.popnames",2]," ")[[1]]
  iddata<-read.table(idfile,as.is=T)
  list(iddata=iddata,paramdata=paramdata,gtout=gtout,idfile=idfile,copyvector_file=copyvector_file,recippop=recippop,donorpops=donorpops)
}

readPaintings<-function(recomblist,copyprobslist,mypops,
                        indsperfile=1,ploidy=2,
                        saveeach=TRUE,loadeach=FALSE,keeplist=TRUE,verbose=TRUE) {
  if(loadeach & saveeach) saveeach<-FALSE
  urecomb<-unique(recomblist)
  nfiles<-length(urecomb)
  
  ninds<-length(copyprobslist)/nfiles*indsperfile
  nhapstot<-ninds*ploidy
  
  if(keeplist) alllist<-list()
  for (recfileon in 1:length(urecomb) ) {
    nsnps<-dim(read.table(urecomb[recfileon],header=T,as.is=T))[1]
    procrecfiles<-which(recomblist==urecomb[recfileon])
    tfiles<-copyprobslist[procrecfiles]

    if(verbose & length(urecomb)>1) print(paste("Processing recombinination file",recfileon,"of",length(urecomb),"called",urecomb[recfileon],"containing",length(tfiles),"files"))

    if(loadeach){
      paintinglist<-readRDS(file=sub(".recombfile",".summedpainting.RDS",urecomb[recfileon]))
    }else{
      paintinglist<-readHaps(tfiles,mypops,nhapstot,nsnps,ninds,ploidy=ploidy,verbose=verbose)
      if(saveeach) saveRDS(paintinglist,file=sub(".recombfile",".summedpainting.RDS",urecomb[recfileon]))
    }
    if(keeplist) alllist[[recfileon]]<-paintinglist
  }
  alllist
}

getP0<-function(alllist,naends=0){
 # genome-wide probability of being breed_dog, wolf, village dog
  ## 
  snpsperfile<-sapply(alllist,function(x){dim(x[[1]])[1]})
  npops<-length(alllist[[1]])
  tcmslist<-lapply(1:npops,function(y){
    sapply(alllist,function(x){
      tx<-x[[y]]
      if(naends>0) tx<-tx[-c(1:naends,dim(tx)[1]-naends:1+1),]
      if(dim(tx)[1]==0) return(colMeans(x[[y]]))
      colMeans(tx)
    })
  })
  ttcmslist<-lapply(tcmslist,function(x){
    t(snpsperfile*t(x)/sum(snpsperfile))
  })
  tpp<-sapply(ttcmslist,rowSums)
  tpp
}

getgwr<-function(signifmat,tthresh){
  myna<-which(is.na(signifmat[,1]))
  fdrvals<-apply(signifmat[-myna,],2,function(x){
    x[x<0]<-0
    x[x>1]<-1
    myz<- qnorm(x)
    torder<-order(x,decreasing=F)

    tsignif=data.frame(index=numeric(),lfdr=numeric(),z=numeric(),p=numeric())
    default=list(signif=tsignif,lfdr=numeric(),p=numeric())
    myn<-length(na.omit(myz<tthresh))
    if(myn==0) return(default)
    
    torder<-order(myz,decreasing=F)
    tlfdr<-x[torder[1:myn]]

    list(signif=x,lfdr=x,p=x[torder[1:myn]])
  })
  list(fdrvals=fdrvals)
}

getlfdr<-function(signifmat,fdrmax=1-1e6){
  require("fdrtool")
  mykeep<-which(!is.na(signifmat[,1]))
  fdrvals<-apply(signifmat[mykeep,],2,function(x) {
    x[x<0]<-0
    x[x>1]<-1
    tsignif=data.frame(index=numeric(),lfdr=numeric(),z=numeric(),p=numeric())

    myz<- qnorm(x)
    tmp<-fdrtool(myz,plot=FALSE)
    default=list(signif=tsignif,lfdr=tmp$lfdr,lfdrmax=1,pmax=0,n50=0,fdrtool=tmp)
    tdull<-which(myz>0 & tmp$lfdr<1)
    tmine<-which(myz<0 & tmp$lfdr<1)
    if(length(tmine)==0) return(default)
    myfrac<-length(tmine)/(length(tdull)+length(tmine))
    
    n50<-floor(myfrac*as.numeric(floor((1-tmp$param[1,"eta0"])/2*length(tmp$pval))))

    if(n50==0) return(default)
    
    tmp$lfdr[myz>0]<-1
    torder<-order(myz,decreasing=F)
    tlfdr<-tmp$lfdr[torder[1:n50]]
    tsignif<-data.frame(index=torder[1:n50],lfdr=tmp$lfdr[torder[1:n50]],z=myz[torder[1:n50]],p=x[torder[1:n50]])
    tsignif<-tsignif[tsignif$lfdr<fdrmax,,drop=FALSE]
    n50<-dim(tsignif)[1]
    if(n50==0) return(default)

    list(signif=tsignif,lfdr=tmp$lfdr,lfdrmax=max(tsignif$lfdr),pmax=max(tsignif$p),n50=n50,fdrtool=tmp)
  })
  
  pthresh<-sapply(fdrvals,function(x){max(x$pmax)})
  n50<-sapply(fdrvals,function(x){x$n50})
  maxlfdr<-sapply(fdrvals,function(x){
    if(dim(x$signif)[1]==0)return(1);
    max(x$signif$lfdr) })
  minlfdr<-sapply(fdrvals,function(x){
    if(dim(x$signif)[1]==0)return(1);
    min(x$signif$lfdr)})
  list(fdrvals=fdrvals,thresh=pthresh,n50=n50,minlfdr=minlfdr,maxlfdr=maxlfdr)
}


writeSignifRegions<-function(pvalsnpwindows,rootname="poibinranges_",fdrregions=NULL) {
  for(i in 1:length(pvalsnpwindows)) {
    tdat<-cbind(pvalsnpwindows[[i]]$chrom,pvalsnpwindows[[i]]$p0snprange,pvalsnpwindows[[i]]$minpval,1+apply(pvalsnpwindows[[i]]$snpranges,1,diff))
    colnames(tdat)<-c("chromosome","snpstart","snpend","minp","nsnps")
    if(!any(is.null(fdrregions))) tdat<-cbind(tdat,minfdr=fdrregions[[i]])
    write.table(tdat,file=paste0(rootname,names(pvalsnpwindows)[i],".txt"),quote=F,row.names=F,col.names=T,sep=",")
  }
}
