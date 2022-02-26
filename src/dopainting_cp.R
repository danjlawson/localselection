
######################
## This is how to do it for standard finestructure
ninds=140
copyprobslist=c(paste0("../subsetdata/test/stage7/test_stage7_tmp_mainrun.linked_file1_ind",1:ninds,".copyprobsperlocus.out.gz"),
                paste0("../subsetdata/test/stage2/test_stage7_tmp_mainrun.linked_file7_ind",1:ninds,".copyprobsperlocus.out.gz"))

recomblist=c(rep("../subsetdata/test.chr1.recomb",times=ninds),
             rep("../subsetdata/test.chr2.recomb",times=ninds))
source("paintingfns.R")

## Lastly, we need to know the populations to be contrasted
allids=read.table("../subsetdata/test.ids",as.is=T)[,1]
mypops<-list()
for(i in 1:2) mypops[[i]]<-allids[(1:10) + (i-1)*10]
names(mypops)<-c("EUR","AFR")
indsperfile<-1 # and the number of individuals in those samples files
npops<-length(mypops)
popnames<-names(mypops)

alllist<-readPaintings(recomblist,
                       copyprobslist,
                       mypops=mypops,
                       indsperfile=indsperfile,
                       ploidy=2,
                       loadeach=FALSE)

loadrds<-FALSE
system("mkdir -p ../subsetdataout")
rootname="../subsetdataout/testcp"
source("getpainting.R")
