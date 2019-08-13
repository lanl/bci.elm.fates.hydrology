rm(list = ls())
library(ncdf4)
# wd<-"/turquoise/usr/projects/veg/rutuja/ACME_cases/BCI/FATES.Sen/OutputExtract/" # on server
wd <- getwd() #on desktop
# source(file.path(wd, "1.0_fun_extract.R")) # on server
source(file.path(wd,"R","1.0_fun_extract.R")) # on desktop
# outdir <- "/lustre/scratch3/turquoise/rutuja/ACME/cases/" # on server
outdir <- "data-raw" # on desktop
filebase <- "BCI.ICLM45ED.badger.intel.Ccd65684b0-F64b5affb."
start.year <- 2008
end.year <- 2016
sam.start <- 1
sam.end <- 4
filterFile <- "Filter.txt"
### for vars with both one and two dimensions
var.vec <- c('QVEGT','QVEGE','QSOIL','QRUNOFF','H2OSOI','GPP')
scale.vec <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

#extractres_h0(sam.start, sam.end, outdir, filebase, var.vec, scale.vec, filterFile, 
#                start.year, end.year)
#print("h0 complete")  
  
extractres_h1(sam.start, sam.end, outdir, filebase, var.vec, scale.vec, filterFile,  
                     start.year, end.year)
		     
print("h1 complete")  
# extractres_h2(sam.start, sam.end, outdir, filebase, var.vec, scale.vec, filterFile, 
  #               start.year, end.year)

# load(file.path(wd, "data-raw", "extract", out.file.name))
## load(file.path(wd, "extract", out.file.name)) #on server
