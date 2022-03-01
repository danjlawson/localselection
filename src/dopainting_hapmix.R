
#### Doing this from H
rootname<-"../hapmixtestout/test.ANCPROB"
rootnamein<-"../hapmixtest/test.ANCPROB"
system("mkdir -p ../hapmixtestout")

## Lastly, we need to know the populations to be contrasted
mypops<-list()
for(i in 1:2) mypops[[i]]<-paste0("IND",(i-1)*20+1:20)
names(mypops)<-c("EUR","AFR")

indsperfile<-1 # and the number of individuals in those samples files
npops<-length(mypops)
popnames<-names(mypops)
chromosomegap<-1e-3
    
## Manual reading of data
allin<-paste0(rootnamein,".",1:2,"")
mydata<-lapply(allin,read.table)

alllist<-list()
for(i in 1:length(mydata)) {
  alllist[[i]]<-list() # afr
  alllist[[i]][[1]]<-list()
  alllist[[i]][[2]]<-list()
  alllist[[i]][[1]]<-mydata[[i]]
  alllist[[i]][[2]]<-1-mydata[[i]]
}
rm(mydata)
gc()


loadrds<-FALSE
source("getpainting.R")
