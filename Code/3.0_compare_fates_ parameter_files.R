##-------------------
## Author: Rutuja Chitra-Tarak
## Date: 08/15/2019
## Title: in fates paramter file reduce pfts to desired number and corresponding values in variables
##-------------------

rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, easyNCDF, tidyverse)

filename1 <- "data-raw/parameter_file_name.nc"
filename2 <- paste0("data-raw/parameter_file_name_2pft.nc")
nc1 <- nc_open(filename1, write = T)
nc2 <- nc_open(filename2, write = T)

old.pfts.n <- nc1$dim$fates_pft$len

## Choose pfts to use:
pft.names <- ncvar_get(nc, "fates_pftname")
# [1] "broadleaf_evergreen_tropical_tree          " "needleleaf_evergreen_extratrop_tree        "
# [3] "needleleaf_colddecid_extratrop_tree       "  "broadleaf_evergreen_extratrop_tree         "
# [5] "broadleaf_hydrodecid_tropical_tree         " "broadleaf_colddecid_extratrop_tree        " 
# [7] "broadleaf_evergreen_extratrop_shrub        " "broadleaf_hydrodecid_extratrop_shrub       "
# [9] "broadleaf_colddecid_extratrop_shrub       "  "arctic_c3_grass                            "
# [11] "cool_c3_grass                              " "c4_grass                                   "
# pfts.to.keep: 
pft.index <- 1
pft.to.keep <- rep(pft.names[pft.index])
## change all pft names to the chosen pft
ncvar_put(nc, "fates_pftname", rep(pft.to.keep, old.pfts.n))

all.vars <- NcReadVarNames(nc1)
var.length <- vector()
for (i in 1: length(all.vars)){
  var.length[i] <- length(ncvar_get(nc1, all.vars[i]))
}

var.df <- data.frame(var = all.vars, var.length = var.length) %>% 
  mutate(pft.n.based.vars = if_else(var.length %% old.pfts.n == 0, TRUE, FALSE),
         multiplier = if_else(pft.n.based.vars == TRUE, var.length/old.pfts.n, 0))
View(var.df)

var.ch <- var.df %>% subset(pft.n.based.vars == TRUE)

View(var.ch) ## We could keep parameters associated with this/these pfts
## Assuming parameter values are stored in the same sequence so as to correspond to sequence of pfts in fates_pftname:
for (i in 1: length(var.ch$var)){
  values.vec <- ncvar_get(nc, as.character(var.ch$var[i]))
  multiplier <- var.ch$multiplier[i]
  values.to.keep <- rep(values.vec[pft.index], old.pfts.n*multiplier)
  ncvar_put(nc, vars.to.change[i], values.to.keep)
}
nc_close(nc)
