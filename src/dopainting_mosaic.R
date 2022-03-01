
#### Doing this from H
rootname<-"../mosaictestout/test.mosaic"
rootnamein<-"~/Downloads/bees_code/wildcat_mosaic/localanc_scot_2way_1-36_1-19_130_60_0.99_58.RData"
system("mkdir -p ../mosaictestout")

## Lastly, we need to know the populations to be contrasted
mypops<-list()
for(i in 1:2) mypops[[i]]<-paste0("IND",i)
names(mypops)<-c("dom","wildcat")

npops<-length(mypops)
popnames<-names(mypops)

## Manual reading of data
tmp=load(rootnamein)

alllist<-list()
for(i in 1:length(localanc)) {
  alllist[[i]]<-list()
  alllist[[i]][[1]]<-t(localanc[[i]][1,,])
  alllist[[i]][[2]]<-t(localanc[[i]][2,,])
  rownames(alllist[[i]][[1]])=g.loc[[i]]
  rownames(alllist[[i]][[2]])=g.loc[[i]]
}
#rm(mydata)
gc()


loadrds<-FALSE
source("getpainting.R")
