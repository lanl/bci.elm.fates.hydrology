#spagetti plots
#get the statistics
#outdir<-"/lustre/scratch3/turquoise/cxu/cmft/fast_analysis/runs"
wd<-getwd()
outdir<-"/lustre/scratch3/turquoise/rutuja/ACME/cases"
filebase<-"BCI.ICLM45ED.badger.intel.Cabfd064-Fc47f968."
finaltag<-"clm2.h0.2018-12.nc" 

start.n <- 1
stop.n <- 5000
id.arr <- c(start.n:stop.n)
nsam <- length(id.arr)
pb <- txtProgressBar(min = 0, max = nsam, style = 3)

file.rm<-F

filter.arr<-rep(F,nsam)

for(i in 1:nsam){
  setTxtProgressBar(pb, i)
  casename <- paste(filebase,id.arr[i],sep="")
  filename<-paste(outdir,"/",casename,"/run/", casename,".",finaltag,sep="")
  if(file.exists(filename)){
    filter.arr[i] <-T
  }
  
}
miss.arr<-id.arr[filter.arr==F]
write.table(as.data.frame(cbind(filter.arr)),"Filter.txt",row.names=F,col.names=F)
write.table(as.data.frame(cbind(miss.arr)),"Missing.txt",row.names=F,col.names=F)
length(which(filter.arr==T))
