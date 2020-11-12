##----------------
# Generate surface data files with varying parameters 
# Author: Rutuja Chitra-Tarak
# Date: Aug 15, 2019
##-----------------

rm(list=ls())
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, easyNCDF, tidyverse)
##-----------------
## Generate parameter ensembles for all parameters together: fates, elm and surface, 
## then substitute generated ensemble streams for fates, elm and surface into copies of a sample paramater file corresponding to fates, elm and surface
##-----------------
# https://docs.google.com/document/d/1cUtioJ99MTMBZWbVA7gMx9ETX5Vr1l4al1Q_WY29Qx8/edit#
# S.N. | Parameter | Description  | Range (Units) | Method of determination with reference
# 1 | fates_leaf_BB_slope | stomatal slope parameter, as per Ball-Berry | 4 - 16 Unitless | The range of Ball Berry parameter fitted across several studies as compiled in Table 2 of (Medlyn et al. 2011)
# 2 | fates_leaf_slatop | Specific Leaf Area (SLA) at the top of canopy, projected area basis | 0.0042 - 0.0400 m^2/gC | Based on the observed range of Leaf Mass Area (LMA)--24.97 to 235.8 g per m2-- for individual tree variation across 51 tree species in Barro Colorado Island (Kunert et al. unpublished data), and assumption of 50% carbon content in biomass. https://github.com/EcoClimLab/HydraulicTraits
# 3 | fates_leaf_vcmax25top | Maximum carboxylation rate of Rubisco at 25 deg C, canopy top | 22.7 - 92.5 umol CO2/m^2/s | Observed range of Vcmax at 25 deg C for tropical tree species under relative radiation of 50% or greater, thereby excluding highly shaded leaves. Values derived with a conversion specific for CLM model are used. (Ali et al. 2015)
# 4 | fates_roota_par | CLM rooting distribution parameter | 5.9 - 7.4 per m | Parameter a in Eq. 2 in (Zeng 2001) that regulates the shape of the rooting profile. Range corresponds to this parameter specified for BATS (or IGBP) land cover classified Deciduous Broadleaf Trees (5.9 m-1) and Evergreen Broadleaf Trees  (7.4 m-1) as given in Table 1 (or 2) of Zeng 2001.
# 5 | fates_rootb_par | CLM rooting distribution parameter | 0.217 - 1.956  per m  | Parameter b in Eq. 2 in (Zeng 2001) that regulates the depth of the rooting profile. Chosen range of b is derived using this equation so as to fit the observed range of rooting depth (dr) of 2 - 18 m for Tropical Deciduous Forest (mean ± Standard Error (S.E.); 3.7 ± 0.5, n = 5 trees; min = 2, max = 4.7)  and Tropical Evergreen Forest (mean ± S.E.; 7.3 ± 0.5, n = 3 trees and 3 communities; min = 2, max = 18 m) combined (Canadell et al. 1996). Besides the direct observation of roots at 18 m included by (Nepstad et al. 1994) in Paragominas, eastern Amazonia that is included in the above study; in Tapajos, eastern Amazonia water extraction by roots was also inferred up to 18 m. (Davidson et al. Forest Science 2011)
# 6 | fates_smpsc  | Soil water potential at full stomatal closure | 113000 - 242000 mm | Based on observed range of -1.13 to -2.42 MPa for leaf turgor loss point for individual tree variation across 49 BCI tree species. (Kunert et al. unpublished data) https://github.com/EcoClimLab/HydraulicTraits
# 7 | aveDTB | Distance to bedrock | 3 - 18 m | Ben Turner; pers. comm.
# 8 | fmax | The maximum fractional saturated area  | 0.01 - 0.80 unitless | Empirical 
# 9 | HKSAT | Soil hydraulic conductivity profile | 0.007 - 0.014 mm/s for = < 12.5 cm; 5.56e-05, 0.0003 for > = 60 cm | To account for high macroporosity and direct flow paths in tropical soils (Kinner and Stallard 2004; Broedel et al. 2017; Tomasella and Hodnett 1996) For Conrad catchment observed median Ksat (Lower95%CI, Upper95%CI) for 12.5 cm depth is 38.3 mm/hr (25.4, 51.2, n = 75), while for 60 cm depth 0.7 (0.2, 1.2, n = 40) mm/hr. For reference a storm with 12.5 mm/hr rainfall intensity has a 0.2 probability of occuring in any given rainfall event (Kinner and Stallard 2004).
# 10 | HKSAT_ADJ | Adjusting factor for soil hydraulic conductivity | 1 - 8 unitless | To account for high macroporosity and direct flow paths in tropical soils that is not accounted for by small soil core samples  (Kinner and Stallard 2004; Broedel et al. 2017; Tomasella and Hodnett 1996)


## setting which fates.params to vary and use in creating ensembles
para.OI.fates <- c(
  "fates_leaf_BB_slope", # 4-16 unitless
  "fates_leaf_slatop", # 0.0042 - 0.0400 m^2/gC
  "fates_leaf_vcmax25top", # 22.7 - 92.5 umol CO2/m^2/s
  "fates_roota_par",# 5.9 - 7.4 per m
  "fates_rootb_par", # 0.217 - 1.956  per m
  "fates_smpsc", # 113000 - 242000 mm
  "fates_allom_l2fr" # 0.001 - 0.94 unitless (actually Range same after plant selection as below)
)
# "fates_allom_l2fr", # previously called "fates_froot_leaf"; but does not play a role until hydro is turned on
# For posterity: Iversen et al_FRED2_20180518_rootshoot.csv:
# Root mass fraction (RMF); Root biomass divided by total plant biomass.F00846 # Root mass fraction (fraction of root dry mass per whole plant dry mass)
# After selecting Angiosperms(F01291)/Broadleaved(F00041)/trees(F00032)

fates.info <- list(par.names = para.OI.fates,
             ## to create this grid, need to provide min max for vars:
             min.param = c(4, 0.0042, 22.7, 5.9, 0.217, 113000, 0.001),
             max.param = c(16, 0.0400, 92.5, 7.4, 1.956, 242000, 0.94))
## for elm:
#1 | fpi_max | Maximum interception fraction of precipitation | 0.05-0.44 unitless | Based on throughfall data from (Zimmermann et al. 2010) and precipitation data for BCI from STRI Physical Monitoring program, defined for precipitation events greater than 10 mm.

elm.info <- list(par.names = "fpi_max",
               ## to create this grid, need to provide min max for vars:
               min.param = c(0.05),
               max.param = c(0.44))
## for surf:
### Surface data params to vary simultaneously
# aveDTB | Distance to bedrock | 3 - 18 m | Ben Turner; pers. comm. 
# fmax | The maximum fractional saturated area  | 0.1 - 0.8 unitless | Empirical 
# HKSAT_ADJ | Adjusting factor for hydraulic conductivity | 1 - 8 unitless | 
surf.info <- list(par.names = c("FMAX", "aveDTB","HKSAT_ADJ"), 
               min.param = c(0.01, 3, 1), 
               max.param = c(0.80, 18, 8)) 
all.info <- list(par.names = c(fates.info$par.names, elm.info$par.names, surf.info$par.names),
                 min.param = c(fates.info$min.param, elm.info$min.param, surf.info$min.param),
                 max.param = c(fates.info$max.param, elm.info$max.param, surf.info$max.param))

# Generated parameters should be saved and not changed; unless a new parameter set is required.
n.sam <- 10000
n.par <- length(all.info$par.names)
# grid <- lhs::randomLHS(n.sam, n.par)
# all.params <- matrix(0, n.sam, n.par)
# 
# ## generating ensembles
# for (i in 1: n.sam) {
#   for (j in 1: n.par) {
#     all.params[i, j] <- qunif(grid[i, j], min = all.info$min.param[j], max = all.info$max.param[j])
#   }
# }
# 
# all.params.df <- data.frame(all.params); colnames(all.params.df) <- all.info$par.names
# summary(all.params.df)
# # fates_leaf_BB_slope fates_leaf_slatop  fates_leaf_vcmax25top fates_roota_par fates_rootb_par   fates_smpsc     fates_allom_l2fr  
# # Min.   : 4          Min.   :0.004202   Min.   :22.71         Min.   :5.900   Min.   :0.2171   Min.   :113012   Min.   :0.001043  
# # 1st Qu.: 7          1st Qu.:0.013150   1st Qu.:40.15         1st Qu.:6.275   1st Qu.:0.6518   1st Qu.:145255   1st Qu.:0.235811  
# # Median :10          Median :0.022101   Median :57.60         Median :6.650   Median :1.0864   Median :177503   Median :0.470530  
# # Mean   :10          Mean   :0.022100   Mean   :57.60         Mean   :6.650   Mean   :1.0865   Mean   :177500   Mean   :0.470500  
# # 3rd Qu.:13          3rd Qu.:0.031049   3rd Qu.:75.05         3rd Qu.:7.025   3rd Qu.:1.5212   3rd Qu.:209746   3rd Qu.:0.705258  
# # Max.   :16          Max.   :0.039999   Max.   :92.50         Max.   :7.400   Max.   :1.9559   Max.   :241993   Max.   :0.939961  
# # fpi_max             FMAX             aveDTB         HKSAT_ADJ   
# # Min.   :0.05002   Min.   :0.01006   Min.   : 3.001   Min.   :1.00  
# # 1st Qu.:0.14750   1st Qu.:0.20751   1st Qu.: 6.750   1st Qu.:2.75  
# # Median :0.24499   Median :0.40502   Median :10.501   Median :4.50  
# # Mean   :0.24500   Mean   :0.40500   Mean   :10.500   Mean   :4.50  
# # 3rd Qu.:0.34249   3rd Qu.:0.60246   3rd Qu.:14.249   3rd Qu.:6.25  
# # Max.   :0.43999   Max.   :0.79997   Max.   :17.999   Max.   :8.00
# write.csv(all.params.df, "data-raw/all.params_but_ksat.csv", row.names = FALSE)
all.params.df <- read.csv("data-raw/all.params_but_ksat.csv", header = TRUE)

#-------------------------------------------------
# Generate fates param files with varying parameters 
#-------------------------------------------------
## to generate 1 pft file:
# Go to:
#   /turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/external_models/fates/tools
# 
# Run:
#   python FatesPFTIndexSwapper.py --pft-indices=1 --fin=/turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/external_models/fates/parameter_files/fates_params_default_newcopy.nc --fout=/turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/external_models/fates/parameter_files/fates_params_default_1pft.nc
#================================================
basefile.fates <- "data-raw/fates_params_default_2pft.nc"
nc <- nc_open(basefile.fates, write = T)
all.vars <- NcReadVarNames(nc)

#================================================
#See what are the actual values of these variables in the parameter file:
fates.n.par <- length(fates.info$par.names)
compare.1 <- data.frame(var = fates.info$par.names, gen.val = fates.info$min.param)   
compare.1$file.val <- NA
for (j in 1:fates.n.par) {
  var.name <- fates.info$par.names[j]
  compare.1$file.val[j] <- as.double(ncvar_get(nc, var.name)[1])
}
compare.1
## so we need to match sign and order as in the file:
compare.1 <- compare.1 %>% mutate(reverse.flag = rep(FALSE, fates.n.par), reciprocal.flag = rep(FALSE, fates.n.par), scale = rep(1, fates.n.par))
compare.1$reverse.flag[which(compare.1$var == "fates_smpsc")] <- TRUE 
# if sla were in biomass terms (m^2/g biomass): # carbon from biomass to carbon (0.5): 
# compare.1$scale[which(compare.1$var == "fates_leaf_slatop")] <- 0.5 
# so after all the conversions generated and file values should be comparable:
compare.1
compare.1 <- compare.1 %>% mutate(
  converted.new = if_else(reverse.flag == TRUE, -gen.val, gen.val),
  converted.new = if_else(reciprocal.flag == TRUE, 1/converted.new, converted.new),
  converted.new = scale*converted.new)
compare.1
# var  gen.val  file.val reverse.flag reciprocal.flag scale converted.new
# 1   fates_leaf_BB_slope 4.00e+00  8.00e+00        FALSE           FALSE     1      4.00e+00
# 2     fates_leaf_slatop 4.20e-03  1.20e-02        FALSE           FALSE     1      4.20e-03
# 3 fates_leaf_vcmax25top 2.27e+01  5.00e+01        FALSE           FALSE     1      2.27e+01
# 4       fates_roota_par 5.90e+00  7.00e+00        FALSE           FALSE     1      5.90e+00
# 5       fates_rootb_par 2.17e-01  1.00e+00        FALSE           FALSE     1      2.17e-01
# 6           fates_smpsc 1.13e+05 -2.55e+05         TRUE           FALSE     1     -1.13e+05
# 7      fates_allom_l2fr 1.00e-03  1.00e+00        FALSE           FALSE     1      1.00e-03

pft.dim.lower.1 <- rep(1, fates.n.par)
pft.dim.upper.1 <- rep(2, fates.n.par) ## also for dry deciduous
tissue.dim.lower.1 <- rep(0, fates.n.par)
tissue.dim.upper.1 <- rep(0, fates.n.par)

#-------------------------------------------------
# create the files
# n.sam <- nrow(fast.para)
# fates.n.par<-2

fates.dir <- "fates.param.sam"
if(!dir.exists(fates.dir)) {dir.create(fates.dir)}
f.to.rm <- list.files(fates.dir)
if(length(f.to.rm) > 0) {file.remove(file.path(fates.dir, f.to.rm), showWarnings =  FALSE)}

files.copied <- vector()
for (i in 1:n.sam) {
  basefile.fates <- "data-raw/fates_params_default_2pft.nc"
  filename2 <- paste0("fates.param.sam/parameter_file_name", i, ".nc")
  file.copy(from = basefile.fates, to = filename2, overwrite = TRUE)
  files.copied <- c(files.copied, file.exists(basefile.fates))
}
length(files.copied)
#-------------------------------------------------
#set the parameter value
pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
for (i in 1:n.sam) {
  setTxtProgressBar(pb, i)
  filename2 <- paste0("fates.param.sam/parameter_file_name", i, ".nc")
  fates_para <- nc_open(filename2, write = T)
  for (j in 1:fates.n.par) {
    var.name <- fates.info$par.names[j]
    if (var.name != "NA") {
      para.val <- all.params.df[i, var.name] * compare.1$scale[j]
      if (compare.1$reverse.flag[j])
        para.val = -para.val
      if (compare.1$reciprocal.flag[j])
        para.val = 1.0 / para.val
      # if (conver.flag1[j])
      #   para.val = conversion(para.val, 0.01, 0.99)
      # if (conver.flag2[j])
      #   para.val = conversion(para.val, 0.015, 0.995)
      fates.para.values <- ncvar_get(fates_para, var.name)
      if (tissue.dim.lower.1[j] == 0) {
        fates.para.values [pft.dim.lower.1[j]:pft.dim.upper.1[j]] <- para.val
      } else {
        fates.para.values [pft.dim.lower.1[j]:pft.dim.upper.1[j], tissue.dim.lower.1[j]:tissue.dim.upper.1[j]] <-
        para.val
      }
      ncvar_put(fates_para, var.name, fates.para.values)
    } # if (var.name!="NA"){
  } #j
  nc_close(fates_para)
} #i

## Check whether a file has been converted as desired:
filename.test.1 <- paste0("fates.param.sam/parameter_file_name", 1, ".nc")
para.test.1 <- nc_open(filename.test.1, write = T)
## value in the file
file.value.1 <- ncvar_get(para.test.1, fates.info$par.names[1])[1]
## value that was to be substituted
desired.value.1 <- all.params.df[1, fates.info$par.names[1]] * compare.1$scale[1]
## make sure the two does match
file.value.1 == desired.value.1 # should be true 

# now remove any older fates.param.sam.zip file, then zip fates.param.sam
file.remove("fates.param.sam.zip")
f.to.zip.1 <- list.files(fates.dir)
zip(zipfile = "fates.param.sam.zip", files = file.path(fates.dir, f.to.zip.1))
# if this does not work, zip manually

# unzip(zipfile = "fates.param.sam.zip")

# copy over to the server
# scp ~/Work_at_LANL/Projects/clm_fates_hydrology/fates.param.sam.zip rutuja@wtrw.lanl.gov:ba-fe1:/turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/external_models/fates/parameter_files/

### End of fates files section--------------

#-------------------------------------------------
# Generate elm param files with varying parameters 
#-------------------------------------------------

#================================================
basefile.elm <- "data-raw/clm_params_c180301.nc"
nc.2 <- nc_open(basefile.elm, write = T)
all.vars <- NcReadVarNames(nc.2)

#================================================
#See what are the actual values of these variables in the parameter file:
elm.n.par <- length(elm.info$par.names)
compare.2 <- data.frame(var = elm.info$par.names, gen.val = elm.info$min.param)   
compare.2$file.val <- NA
for (j in 1:elm.n.par) {
  var.name <- elm.info$par.names[j]
  compare.2$file.val[j] <- as.double(ncvar_get(nc.2, var.name)[1])
}
compare.2
## so we need to match sign and order as in the file:
compare.2 <- compare.2 %>% mutate(reverse.flag = rep(FALSE, elm.n.par), reciprocal.flag = rep(FALSE, elm.n.par), scale = rep(1, elm.n.par))
# so after all the conversions generated and file values should be comparable:
compare.2
compare.2 <- compare.2 %>% mutate(
  converted.new = if_else(reverse.flag == TRUE, -gen.val, gen.val),
  converted.new = if_else(reciprocal.flag == TRUE, 1/converted.new, converted.new),
  converted.new = scale*converted.new)
compare.2
# var gen.val file.val reverse.flag reciprocal.flag scale converted.new
# 1 fpi_max    0.05     0.05        FALSE           FALSE     1          0.05

pft.dim.lower.2 <- rep(1, elm.n.par)
pft.dim.upper.2 <- rep(1, elm.n.par)
tissue.dim.lower.2 <- rep(0, elm.n.par)
tissue.dim.upper.2 <- rep(0, elm.n.par)
#-------------------------------------------------
# create the files
elm.dir <- "elm.param.sam"
if(!dir.exists(elm.dir)) {dir.create(elm.dir)}
f.to.rm.2 <- list.files(elm.dir)
if(length(f.to.rm) > 0) {file.remove(file.path(elm.dir, f.to.rm.2), showWarnings =  FALSE)}

files.copied <- vector()
for (i in 1:n.sam) {
  basefile.fates <- "data-raw/clm_params_c180301.nc"
  filename2 <- paste0("elm.param.sam/parameter_file_name", i, ".nc")
  file.copy(from = basefile.fates, to = filename2, overwrite = TRUE)
  files.copied <- c(files.copied, file.exists(filename2))
}
length(files.copied)
#-------------------------------------------------
#set the parameter value
pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
for (i in 1:n.sam) {
  setTxtProgressBar(pb, i)
  filename2 <- paste0("elm.param.sam/parameter_file_name", i, ".nc")
  elm_para <- nc_open(filename2, write = T)
  for (j in 1:elm.n.par) {
    var.name <- elm.info$par.names[j]
    if (var.name != "NA") {
      para.val <- all.params.df[i, var.name] * compare.2$scale[j]
      if (compare.2$reverse.flag[j])
        para.val = -para.val
      if (compare.2$reciprocal.flag[j])
        para.val = 1.0 / para.val
      # if (conver.flag1[j])
      #   para.val = conversion(para.val, 0.01, 0.99)
      # if (conver.flag2[j])
      #   para.val = conversion(para.val, 0.015, 0.995)
      elm.para.values <- ncvar_get(elm_para, var.name)
      if (tissue.dim.lower.2[j] == 0) {
        elm.para.values [pft.dim.lower.2[j]:pft.dim.upper.2[j]] <- para.val
      } else {
        elm.para.values [pft.dim.lower.2[j]:pft.dim.upper.2[j], tissue.dim.lower.2[j]:tissue.dim.upper.2[j]] <-
          para.val
      }
      ncvar_put(elm_para, var.name, elm.para.values)
    } # if (var.name!="NA"){
  } #j
  nc_close(elm_para)
} #i


## Check whethe a file has been converted as desired:
filename.test.2 <- paste0("elm.param.sam/parameter_file_name", 1, ".nc")
para.test.2 <- nc_open(filename.test.2, write = T)
## value in the file
file.value.2 <- ncvar_get(para.test.2, elm.info$par.names[1])[1]
## value that was to be substituted
desired.value.2 <- all.params.df[1, elm.info$par.names[1]] * compare.2$scale[1]
## make sure the two does match
file.value.2 == desired.value.2 # should be true 

# now remove any older surf.sam.zip file, then zip surf.sam
if(file.exists("elm.param.sam.zip")) {file.remove("elm.param.sam.zip")}
f.to.zip.2 <- list.files(elm.dir)
zip(zipfile = "elm.param.sam.zip", files = file.path(elm.dir, f.to.zip.2))
# if this does not work, zip manually
# unzip(zipfile = "elm.param.sam.zip")

# copy over to the server
# scp ~/Work_at_LANL/Projects/clm_fates_hydrology/elm.param.sam.zip rutuja@wtrw.lanl.gov:ba-fe1:/turquoise/usr/projects/veg/rutuja/ACME_cases/param/
### End of elm files section--------------

#-------------------------------------------------
# Generate surface data files with varying parameters 
#-------------------------------------------------

## also adding HKSAT:
# From 3.2.5_Surfdata_Ksat_obs_and_bootstrapped_param.R
ksat <- read.csv("data-raw/Ten thousand bootstrapped Ksat profiles within 95CI of observed data.csv", header = TRUE)
## ksat n.sam (ncol - first column) should match n.sam for the create.files code to work;
ncol(ksat) - 1 == n.sam #should be TRUE
# if not TRUE, go to 3.2.5_Surfdata_Ksat_obs_and_bootstrapped_param.R,
# make n.sam in that file == n.sam in this file and regenerate the above ksat file
## storing with all params
## Only selecting ksat at 12.5 and 60 where maximum contrast is available
params <- all.params %>% mutate(HKSAT_12.5 = as.numeric(ksat[4, 2:ncol(ksat)]), # at 12.5 cm
                                HKSAT_60 = as.numeric(ksat[7, 2:ncol(ksat)])) # at 60 cm
write.csv(params, "data-raw/params.csv", row.names = FALSE)

surf.info$par.names[4] <- c("HKSAT")
surf.n.par <- length(surf.info$par.names)
#-------------------------------------------------
#create the files
# the base file "data-raw/surfdata_bci_panama_v1_c171113.nc" should have data for HKSAT for 10 depths; insert a sample if not present.
#HKSAT =  HKSAT =
# // HKSAT(0,0)
# 0.01,
# // HKSAT(1,0)
# 0.01,
# // HKSAT(2,0)
# 0.01,
# // HKSAT(3,0)
# 0.01,
# // HKSAT(4,0)
# 0.01,
# // HKSAT(5,0)
# 0.007,
# // HKSAT(6,0)
# 0.003,
# // HKSAT(7,0)
# 0.0001,
# // HKSAT(8,0)
# 0.0001,
# // HKSAT(9,0)
# 0.0001 ;

#-------------------------------------------------
#create the files

surf.dir <- "surf.sam"
if(!dir.exists(surf.dir)) {dir.create(surf.dir)}
f.to.rm.3 <- list.files(surf.dir)
if(length(f.to.rm.3) > 0) {file.remove(file.path(surf.dir, f.to.rm.3), showWarnings =  FALSE)}
f.to.rm.3 <- list.files(surf.dir)
if(length(f.to.rm.3) > 0) {file.remove(file.path(surf.dir, f.to.rm.3), showWarnings =  FALSE)}

for (i in 1:n.sam) {
  filename1 <- "data-raw/surfdata_bci_panama_v1_c171113.nc"
  filename2 <- paste0("surf.sam/surfdata_bci_panama_v1_c171113.", i, ".nc")
  file.copy(filename1, filename2, overwrite = T)
}
#-------------------------------------------------
#set the parameter value
pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
for (i in 1:n.sam) {
  setTxtProgressBar(pb, i)
  filename2 <- paste0("surf.sam/surfdata_bci_panama_v1_c171113.", i, ".nc")
  surf_para <- nc_open(filename2, write = T)
  for (j in 1:surf.n.par) {
    var.name <- surf.info$par.names[j]
    if (var.name == "HKSAT") {
      surf.para.values <- ncvar_get(surf_para, var.name)
      surf.para.values <- ksat[1:10, 1+i] ## need to substitute ksat for only first 10 depths, skipping the first column named depth
      ncvar_put(surf_para, var.name, surf.para.values)
    } else {
      if (var.name != "NA") {
        surf.para.values <- ncvar_get(surf_para, var.name)
        surf.para.values <- all.params.df[i, var.name]
        ncvar_put(surf_para, var.name, surf.para.values)
      } 
    }
  } #j
  nc_close(surf_para)
} #i
## Check whether a file has been converted as desired:
filename.test.3 <- paste0("surf.sam/surfdata_bci_panama_v1_c171113.", 1, ".nc")
para.test.3 <- nc_open(filename.test.3, write = T)
## value in the file
file.value.3 <- ncvar_get(para.test.3, surf.info$par.names[1])[1]
## value that was to be substituted
desired.value.3 <- all.params.df[1, surf.info$par.names[1]]
## make sure the two does match
file.value.3 == desired.value.3 # should be true 
## for ksat compare following
ncvar_get(para.test.3, surf.info$par.names[4])
## value that was to be substituted
ksat[1:10, 1+1]

## substituting AVA tower texture and organic matter

txr.4 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                    PCT_SAND = c(11.5, 11.5, 11.5, 7.0, 6.5, 6.0, 11.0, 13.5,	14.0, 13),
                    PCT_CLAY = c(60.0, 60.0, 60.0, 64.5, 67, 72.0, 69.0, 70, 13.0, 9), # earlier chosen 70 instead of 45 for 1.04, given Matteo's data at AVA TOWER
                    ORGANIC = c(12.1, 12.1, 12.1, 4.6, 3.2, 1.8, 0, 0, 0, 0))
txr.par.names <- colnames(txr.4)[-1]
txr.n.par <- length(txr.par.names)
write.csv(txr.4, file = file.path("data-raw", "AVA_soil_texture_in_CLM_10_depths.csv"), row.names = FALSE)

pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
for (i in 1:n.sam) {
  setTxtProgressBar(pb, i)
  filename2 <- paste0("surf.sam/surfdata_bci_panama_v1_c171113.", i, ".nc")
  surf_para <- nc_open(filename2, write = T)
  for (j in 1:txr.n.par) {
    var.name <- txr.par.names[j]
    if (var.name != "NA") {
      surf.para.values <- ncvar_get(surf_para, var.name)
      surf.para.values <- txr.4[, txr.par.names[j]]
      ncvar_put(surf_para, var.name, surf.para.values)
    } 
  } #j
  nc_close(surf_para)
} #i
## Check whether a file has been converted as desired:
filename.test.4 <- paste0("surf.sam/surfdata_bci_panama_v1_c171113.", 1, ".nc")
para.test.4 <- nc_open(filename.test.4, write = T)
## value in the file
file.value.4 <- ncvar_get(para.test.4, txr.par.names[1])
## value that was to be substituted
desired.value.4 <- txr.4[, txr.par.names[1]]
## make sure the two does match
file.value.4 == desired.value.4 # should be true 

# now remove any older surf.sam.zip file, then zip surf.sam
file.remove("surf.sam.zip")
f.to.zip.3 <- list.files(surf.dir)
zip(zipfile = "surf.sam.zip", files = file.path(surf.dir, f.to.zip.3))
# unzip(zipfile = "surf.sam.zip")

# scp ~/Work_at_LANL/Projects/clm_fates_hydrology/surf.sam.zip rutuja@wtrw.lanl.gov:ba-fe1:/turquoise/usr/projects/veg/rutuja/ACME_cases/BCI/bci_0.1x0.1_v4.0i/
### End of surface data files generation--------------



## For LAI & leaf age

# From Sobrado et al. 1994
# Stage | Collection/when stage began | Months spent in each stage | yr units | Amax mumol/m2/s | proportion to max Amax
# Expanding leaves | September | 1 | 0.08 | 1.6 | 0.21
# Young | NA | 2 | 0.17 | 7.7 | 1
# Mature leaves | April | 5 | 0.42 | 4.5 | 0.58
# Senescent | July | 5 | 0.42 | 1.8 |  0.23

# So far two pfts:
# fates_leaf_long =
#   0.08, 0.08,
#   0.17, 0.17,
#   0.42, 0.42,
#   0.42, 0.42;

# vcmax.base = 50
# prop <- c(0.21, 1, 0.58, 0.23)
# prop*50
# fates_leaf_vcmax25top = 
#   10.5, 10.5,
#   50.0, 50.0,
#   29.0, 29.0,
#   11.5, 11.5 ;
#   
  
# Evergreen species leaf longevity > 12 months Sobrado



#### R code to create soil_depth_arr & fmax_arr:
# depth <- rep(3:10, each = 5)
# fmax <- rep(c(0.2, 0.3, 0.4, 0.5, 0.6), times = 8)
# paste(depth, collapse = ",")
# paste(fmax, collapse = ",")


# #### On server:
# cd $CASE_ROOT/BCI/bci_0.1x0.1_v4.0i/
#   
#   set soil_depth_arr={3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10}
# set fmax_arr = {0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6,0.2,0.3,0.4,0.5,0.6}
# 
# foreach case_i (`seq 1 40`)
# cp xx.cdl new.cdl
# 
# set depth_i=$soil_depth_arr[$case_i]
# set fmax_i=$fmax_arr[$case_i]
# 
# set surf_file_CLONE = "surfdata_bci_panama_v1_c171113_aveDTB${depth_i}.fmax.${fmax_i}.nc" 
# 
# sed -i "s/FMAX = 0.2/FMAX = ${fmax_i}/g" new.cdl
# sed -i "s/aveDTB = 3/aveDTB = ${depth_i}/g" new.cdl
# 
# ncgen -o $surf_file_CLONE new.cdl
# end
# 
# ## new
# total.n <- 2000
# fates.n <- 100; surf.n <- 20
# surf.par_arr <- rep(c(1:surf.n), each = fates.n); paste(surf.par_arr, collapse=",")
# param_arr <- rep(c(1:fates.n), length.out = total.n); paste(param_arr, collapse=",")
# ## so that all surf.par are covered in first few iterations, this should be changed to:
# surf.par_arr <- rep(c(1:surf.n), length.out = total.n); paste(surf.par_arr, collapse=",")
# param_arr <- rep(c(1:fates.n), each = surf.n); paste(param_arr, collapse=",")
