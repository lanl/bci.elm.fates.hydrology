rm(list = ls())
library(ncdf4)
wd <-"/turquoise/usr/projects/veg/rutuja/ACME_cases/BCI/FATES.Sen/OutputExtract/" # on server
source(file.path(wd, "1.0_fun_extract.R")) # on server
# source(file.path("R","1.0_fun_extract.R")) # on desktop
outdir <- "/lustre/scratch3/turquoise/rutuja/ACME/cases" # on server
# current.folder <- "2019-10-14_5000/trial_leafage_fdrought" # on desktop
# outdir <- paste0("data-raw/", current.folder, "/extract")  # on desktop
# if(!dir.exists(file.path("data-raw/extract", current.folder, "extract"))) {dir.create(file.path("data-raw/extract", current.folder, "extract"))} # on desktop

filebase <- "BCI.ICLM45ED.badger.intel.Cabfd064-Fc47f968.full."
start.year <- 1985
end.year <- 2018
sam.start <- 1
sam.end <- 100
filterFile <- "Filter.txt"
### for vars with both one and two dimensions
#var.vec <- c('QVEGT','QVEGE','QSOIL','QRUNOFF','H2OSOI','GPP', 'SOILPSI', 'BTRAN')
#var.vec <- c('TWS','TWS_MONTH_END','ZWT','ZWT_PERCH','RAIN','QINTR','QDRIP','QDRAI','QCHARGE','QSOIL','QVEGT','QVEGE','QRUNOFF')
#var.vec.h0 <- c('TLAI')
#scale.vec.h0 <- c(1.0)
#var.vec.h1 <- c('H2OSOI','QRUNOFF', 'QVEGE', 'QVEGT', 'QSOIL', 'GPP', 'TWS', 'BTRAN', 'SOILPSI')
var.vec.h1 <- c('H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'GPP', 'TWS', 'ZWT', 'BTRAN', 'SOILPSI')
#scale.vec.h1 <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
scale.vec.h1 <- c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)
## Monthly
#extractres_h0(sam.start, sam.end, outdir, filebase, var.vec.h0, scale.vec.h0, filterFile, 
#                start.year, end.year)
#print("h0 complete")  

## Daily  
extractres_h1(sam.start, sam.end, outdir, filebase, var.vec.h1, scale.vec.h1, filterFile,
              start.year, end.year, current.folder)
print("h1 complete")
## Hourly 
# extractres_h2(sam.start, sam.end, outdir, filebase, var.vec.h2, scale.vec.h2, filterFile,
#               start.year, end.year)
# print("h2 complete")

## load(file.path(wd, "data-raw", "extract", out.file.name))
## load(file.path(wd, "extract", out.file.name)) #on server