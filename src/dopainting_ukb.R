
### READING THE NEW WAY
getsnppos<-function(files){
    ## Return the SNP positions (header) for each file in files
    ret=lapply(files,function(x)
        as.character(read.table(x,nrows=1,row.names=1)[1,]))
    names(ret)=names(files)
    ret
}
getnsnps<-function(filemat){
    ## Count the number of SNPs for every file in a matrix of files
    ## Where each row is a different population with the same assumed number of SNPs
    ## Also accepts a single file
    if((class(filemat)=="character") & (length(filemat)==1)) return(length(getsnppos(filemat)[[1]]))
    snppos=getsnppos(filemat[1,])
    sapply(snppos,length)
}
getnhaps<-function(file){
   length(count.fields(file, sep = "\n"))-1
}
processFile<-function(file,nhaps=NULL,nsnps=NULL,scale=1,verbose=1){
    ## Read the files with cpp code
    require("Rcpp")
    sourceCpp("rcppreadfile.cpp")
    if(all(is.null(nhaps))) {
        if(verbose>0) print(paste("Counting haplotypes..."))
        nhaps=getnhaps(file)
    }
    if(all(is.null(nsnps))) {
        if(verbose>0) print(paste("Counting SNPs..."))
        nsnps=getnsnps(file)
    }
    if(verbose>0) print(paste("Computing genome summaries..."))
    res=readfilecpp(file,nhaps,nsnps,scale,verbose)
    res$snpnames=as.numeric(res$snpnames)
    res$nsnps=nsnps
    res$nhaps=nhaps
    res
}
getmeanpainting<-function(filemat,nsnps=NULL,nhaps=NULL,scale=1/9,verbose=1,...){
    ## Run the data extraction on all files in the file matrix
    if(all(is.null(nhaps))) {
        if(verbose>0) print(paste("Counting haplotypes..."))
        nhaps=getnhaps(filemat[1,1])
    }
    if(all(is.null(nsnps))) {
        if(verbose>0) print(paste("Counting SNPs..."))
        nsnps=getnsnps(filemat)
    }
    res=lapply(1:dim(filemat)[2],function(chr){
        rres=lapply(1:dim(filemat)[1],function(pop){
            if(verbose>0) print(paste("Running chromosome",chr,
                                      "of",dim(filemat)[2],
                                      "with population",pop,"of",dim(filemat)[1]))
            processFile(filemat[pop,chr],
                        nhaps,nsnps[chr],
                        scale=scale,verbose=(verbose>1))
        })
        names(rres)=rownames(filemat)
    })
    return(res)
}
computesums = function(filepath,nsnps,nhaps,scale=1/9, sep=" ",verbose=T) {
    ## R approach: slow so this is for reference only
    sumsnps=rep(0,nsnps)
    sumhaps=rep(0,nhaps)
    con = file(filepath, "r")
    ## Skip the header
    myLine=scan(con,what="character",nlines=1,sep=sep,skip=0,quiet=TRUE)
    hapon=0
    while (length(
        myLine <- scan(con,what="character",nlines=1,sep=sep,skip=0,quiet=TRUE)
    ) > 0 ){
        hapon<-hapon + 1
        vals=as.numeric(myLine[-1])*scale
        sumhaps[hapon]=sum(vals)
        sumsnps=sumsnps + vals
    }
    close(con)
    return(list(sumsnps=sumsnps,
                sumhaps=sumhaps,
                nsnps=nsnps,
                nhaps=nhaps))
} 

##################
donorpops=c("Yamnaya","EHG","WHG","CHG","Farmer","African","EastAsian")
infiles=paste0("../yaoling_files/",donorpops,".50000inds.6.master_all_copyprobsperlocus.txt.gz")
infiles2=gsub(".50000inds",".tiny.1000inds",infiles)
infiles3=gsub(".50000inds",".tiny2.1000inds",infiles)

filemat=data.frame("chr6"=infiles)
##filemat=data.frame("chr6a"=infiles2,"chr6b"=infiles3)
rownames(filemat)=donorpops

##################
## Construction of smaller test data
for (i in 1:length(infiles)){
    f=infiles[i]
    ff=infiles2[i]
    cmd=paste0("gunzip -c ",f," | cut -f 1-9001 -d' ' | head -n 2001 | gzip -c > ",ff)
    system(cmd)
## Make a second file
    ff=infiles3[i]
    cmd=paste0("gunzip -c ",f," | cut -f 1,9002-18002 -d' ' | head -n 2001 | gzip -c > ",ff)
    system(cmd)
}
##################
filemat=data.frame("chr6a"=infiles2,"chr6b"=infiles3)
rownames(filemat)=donorpops

allres=getmeanpainting(filemat)
saveRDS(allres,file="chr6_meanpainting.RDS")
source("paintingfns.R")



### READING THE OLD WAY
## alllist=list()
## alllist[[1]]=
##     lapply(infiles2,function(x){
##         r=data.table::fread(x,header=T)
##         rn=r[,1]
##         r=as.matrix(r[,-1])
##         r=t(r)
##         colnames(r)=paste0(as.character(as.data.frame(rn)[,1]),c(".A",".B"))
##         r/9
##     })
## names(alllist[[1]])=donorpops
## loadrds<-FALSE
## system("mkdir -p ../testukb_out")
## rootname="../testukb_out"
## source("getpainting.R")


Rcpp::sourceCpp('rcppreadfile.cpp', verbose = T, rebuild = T)

resC=processFile("../yaoling_files/Yamnaya.tiny.1000inds.6.master_all_copyprobsperlocus.txt",scale=1/9,verbose=T)

resCgz=processFile("../yaoling_files/Yamnaya.tiny.1000inds.6.master_all_copyprobsperlocus.txt.gz",
                 scale=1/9,verbose=T)

resCgz=processFile("../yaoling_files/Yamnaya.50000inds.6.master_all_copyprobsperlocus.txt.gz",
                 nhaps,nsnps[fileon],scale=1/9)

resCgz=processFile(filemat[1,fileon],
                 nhaps,nsnps[fileon],scale=1/9)
