
######################
## This is how to do it for standard finestructure
ninds=140
copyprobslist=c(paste0("../subsetdata/test/stage7/test_stage7_tmp_mainrun.linked_file1_ind",
                       1:ninds,".copyprobsperlocus.out.gz"),
                paste0("../subsetdata/test/stage7/test_stage7_tmp_mainrun.linked_file2_ind",
                       1:ninds,".copyprobsperlocus.out.gz"))

recomblist=c(rep("../subsetdata/test.chr1.recomb",times=ninds),
             rep("../subsetdata/test.chr2.recomb",times=ninds))
source("paintingfns.R")

## Lastly, we need to know the populations to be contrasted
allids=read.table("../subsetdata/test.ids",as.is=T)[,1]
mypops=readIds("../subsetdata/test.ids")
#mypops=mypops[c("ada","afr")]

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
