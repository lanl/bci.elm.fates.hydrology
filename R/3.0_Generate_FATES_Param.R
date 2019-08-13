# setwd("C:/Users/242621/Documents/FAST")

## setting which params to vary and use in creating ensembles
info <- list(par.names = c("lma", # min max [30, 135] obtained from Knobby's data:
                           #https://github.com/rutujact/HydraulicTraits/blob/master/data/Panama/processed_trait_data/Panama_all_traits_table_species_level.csv
                           #Massoud et al. SLA range: [0.010,0.014], Knobby's data SLA range would be [0.007, 0.033]
                           "roota_par", # Massoud et al range: [5.95 8.05] 
                           "rootb_par", # Massoud et al range: [0.85 1.15] 
                           "rootshoot", # "fates_allom_l2fr", or "froot_leaf" Massoud et al range: [0.85, 1.15]
                           #"pitlp_node_leaf", # min max [1.2, 2.8 Mpa] obtained from Knobby's data
                           #https://github.com/rutujact/HydraulicTraits/blob/master/data/Panama/processed_trait_data/Panama_all_traits_table_species_level.csv"fates_smpsc", # same range as tlp
                           "fates_vcmax25top", # Ali et al 2015 range [10?, 100]# range used in Massoud et al for Para, Amazon was [51, 69]
                           "smpsc", #Soil water potential for stomata closure. # Using  min max [1.2, 2.8 Mpa] obtained from Knobby's data; # Massoud et al range [-2.93e4, -2.16e4]
                           # "smpso", # Soil water potential for stomata opening. Massoud et al range [-7.59e4 -5.61e4]
                           # "pitlp_node_stem",
                           # "pitlp_node_troot", 
                           # "pitlp_node_aroot",
                           "bbslope"), # Stomata conductance slope, # Massoud et al range [7.65 10.35]
             ## to create this grid, need to provide min max for vars:
             min.param = c(30, 5.95, 0.85, 0.85, 5, 1.2, 7.65),
             max.param = c(135, 8.05, 1.15, 1.15, 100, 2.8, 10.35)
             )

## could get min max for some of the ensembles from the past
fast.para.name = "HYDRAULICTRAITS.FAST.Sam"
fast.para = read.table(file.path("data-raw", fast.para.name) ,
                       header = T,
                       sep = "\t")
head(fast.para)
fast.para.available <- max.fast <- min.fast <- vector()
for (i in 1: 3) { # since bbslope is not present in FAST file; this will have to be changed if par.names is lengthend 
  col.name <- info$par.names[i]
  fast.para.available[i] <- col.name
  max.fast[i] <- ceiling(max(fast.para[, col.name]))
  min.fast[i] <- min(fast.para[, col.name])
}
fast.para.available; min.fast; max.fast
# [1] "lma"       "roota_par" "rootb_par"
# [1] 21.40954  6.00100  1.00100
# [1] 302   7   2
# Generated parameters should be saved and not changed; unless a new parameter set is required.
n.sam <- 100
n.par <- length(info$par.names)
# grid <- lhs::randomLHS(n.sam, n.par)
# params <- matrix(0, n.sam, n.par)

# ## generating ensembles
# for (i in 1: n.sam) {
#   for (j in 1: n.par) {
#     params[i, j] <- qunif(grid[i, j], min = info$min.param[j], max = info$max.param[j])
#   }
# }
# 
# params.df <- data.frame(params); colnames(params.df) <- info$par.names
# summary(params.df)
# lma           roota_par       rootb_par        rootshoot      fates_vcmax25top     smpsc          bbslope      
# Min.   : 30.35   Min.   :5.955   Min.   :0.8502   Min.   :0.8520   Min.   : 5.92    Min.   :1.209   Min.   : 7.670  
# 1st Qu.: 56.49   1st Qu.:6.472   1st Qu.:0.9264   1st Qu.:0.9261   1st Qu.:29.11    1st Qu.:1.606   1st Qu.: 8.330  
# Median : 82.09   Median :6.992   Median :0.9999   Median :0.9996   Median :52.73    Median :2.003   Median : 8.995  
# Mean   : 82.50   Mean   :6.999   Mean   :1.0000   Mean   :1.0000   Mean   :52.55    Mean   :1.999   Mean   : 9.001  
# 3rd Qu.:108.74   3rd Qu.:7.526   3rd Qu.:1.0748   3rd Qu.:1.0750   3rd Qu.:76.35    3rd Qu.:2.395   3rd Qu.: 9.664  
# Max.   :134.30   Max.   :8.035   Max.   :1.1486   Max.   :1.1493   Max.   :99.53    Max.   :2.785   Max.   :10.338 

# write.csv(params.df, "data-raw/params.csv", row.names = FALSE)
params.df <- read.csv("data-raw/params.csv", header = TRUE)
# ================================================
# test the parameter values
# setwd("./FASTPara")
# nfname = "parameter_file_name.nc"
# fates_para <- nc_open(file.path("data-raw", nfname), write = F)
# list_var <- names(fates_para$var)
# #write.table(file="varname.txt",list_var)
# var.name <- "fates_hydr_avuln_gs"
# fates_para$var[[var.name]]$unit
# fates_para$var[[var.name]]$longname
# para_values <- ncvar_get(fates_para, var.name)
# para_values
# nc_close(fates_para)
#================================================
reverse.flag <- rep(F, n.par)
reciprocal.flag <- rep(F, n.par)
reciprocal.flag[1] <- T # from lma to sla
scale <- rep(1, n.par)
scale[1] <- 0.5 # carbon from biomass to carbon (0.5) for lma
scale[6] <- 10000 # smpsc from Mpa to mm
# conversion <- function(y, l, u) {
#   x <- l + y * (u - l)
# }
# conver.flag1 <- rep(F, n.par) #,(x-lbound)/(ubound-lbound)	0.01	0.99
# conver.flag1[18:21] <- T
# conver.flag2 <- rep(F, n.par) #,(x-lbound)/(ubound-lbound)	0.015	0.995
# conver.flag2[30:33] <- T
# para.OI.fast <- colnames(fast.para)[1:n.par]
para.OI.fast <- info$par.names
para.OI.fates <- c(
  "fates_leaf_slatop",
  "fates_roota_par", 
  "fates_rootb_par",
  "fates_allom_l2fr", # previously called "fates_froot_leaf"
  "fates_leaf_vcmax25top",
  "fates_smpsc",
  # "fates_hydr_pitlp_node",
  # "fates_hydr_pitlp_node",
  # "fates_hydr_pitlp_node",
  "fates_leaf_BB_slope"
)

pft.dim.lower <- rep(1, n.par)
pft.dim.upper <- rep(12, n.par)
tissue.dim.lower <- rep(0, n.par)
tissue.dim.upper <- rep(0, n.par)
# 
# tissue.dim.lower[5] <- 1
# tissue.dim.upper[5] <- 2

#-------------------------------------------------
#create the files
# n.sam <- nrow(fast.para)
# n.par<-2
fates.dir <- "param.sam"
if(!dir.exists(fates.dir)) {dir.create(fates.dir)}
f.to.rm <- list.files(fates.dir)
if(length(f.to.rm) > 0) {file.remove(file.path(fates.dir, f.to.rm), showWarnings =  FALSE)}

for (i in 1:n.sam) {
  filename1 <- "data-raw/parameter_file_name.nc"
  filename2 <- paste0("param.sam/parameter_file_name", i, ".nc")
  file.copy(filename1, filename2, overwrite = T)
}
#-------------------------------------------------
#set the parameter value
pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
for (i in 1:n.sam) {
  setTxtProgressBar(pb, i)
  filename2 <- paste0("param.sam/parameter_file_name", i, ".nc")
  fates_para <- nc_open(filename2, write = T)
  for (j in 1:n.par) {
    var.name <- para.OI.fates[j]
    if (var.name != "NA") {
      fast.val <- params.df[i, j] * scale[j]
      if (reverse.flag[j])
        fast.val = -fast.val
      
      if (reciprocal.flag[j])
        fast.val = 1.0 / fast.val
      
      # if (conver.flag1[j])
      #   fast.val = conversion(fast.val, 0.01, 0.99)
      # if (conver.flag2[j])
      #   fast.val = conversion(fast.val, 0.015, 0.995)
      
      fates.para.values <- ncvar_get(fates_para, var.name)
      if (tissue.dim.lower[j] == 0) {
        fates.para.values [pft.dim.lower[j]:pft.dim.upper[j]] <- fast.val
      } else {
        fates.para.values [pft.dim.lower[j]:pft.dim.upper[j], tissue.dim.lower[j]:tissue.dim.upper[j]] <-
          fast.val
      }
      ncvar_put(fates_para, var.name, fates.para.values)
      #fates.para.values <- ncvar_get(fates_para,var.name)
      # if (var.name == "fates_hydr_p50_gs") {
      #   var.name = "fates_hydr_avuln_gs"
      #   fast.val.avuln = -4 * 60.15 * (-fast.val) ^ (-1.25) / 100 * fast.val
      #   fates.para.values <- ncvar_get(fates_para, var.name)
      #   if (tissue.dim.lower[j] == 0) {
      #     fates.para.values [pft.dim.lower[j]:pft.dim.upper[j]] <-
      #       fast.val.avuln
      #   } else {
      #     fates.para.values [pft.dim.lower[j]:pft.dim.upper[j], tissue.dim.lower[j]:tissue.dim.upper[j]] <-
      #       fast.val.avuln
      #   }
      #   ncvar_put(fates_para, var.name, fates.para.values)
      # }
      # if (var.name == "fates_hydr_pitlp_node") {
      #   var.name = "fates_hydr_epsil_node"
      #   fates.epsil.values <- ncvar_get(fates_para, var.name)
      #   epsil = mean(fates.epsil.values[pft.dim.lower[j]:pft.dim.upper[j], tissue.dim.lower[j]:tissue.dim.upper[j]])
      #   fast.val.pinot =  fast.val * epsil / (epsil - fast.val)
      #   var.name = "fates_hydr_pinot_node"
      #   fates.para.values <- ncvar_get(fates_para, var.name)
      #   if (tissue.dim.lower[j] == 0) {
      #     fates.para.values [pft.dim.lower[j]:pft.dim.upper[j]] <-
      #       fast.val.pinot
      #   } else {
      #     fates.para.values [pft.dim.lower[j]:pft.dim.upper[j], tissue.dim.lower[j]:tissue.dim.upper[j]] <-
      #       fast.val.pinot
      #   }
      #   ncvar_put(fates_para, var.name, fates.para.values)
      # }
    } # if (var.name!="NA"){
  } #j
  nc_close(fates_para)
} #i


# since I did not store the parameters ##-------
# extract params for the six pars from the files generated
# n.sam <- 100
# n.par <- length(info$par.names)
# para.OI.fates <- c(
#   "fates_leaf_slatop",
#   "fates_roota_par", 
#   "fates_rootb_par",
#   "fates_allom_l2fr", # previously called "fates_froot_leaf"
#   "fates_hydr_pitlp_node",
#   "fates_leaf_BB_slope"
# )
# 
# params.from.files <-  as.data.frame(matrix(0, n.sam, n.par))
# colnames(params.from.files) <- para.OI.fates 
# 
# pb <- txtProgressBar(min = 0, max = n.sam, style = 3)
# for (i in 1:n.sam) {
#   setTxtProgressBar(pb, i)
#   filename2 <- paste0("data-raw/param.sam/parameter_file_name", i, ".nc")
#   fates_para <- nc_open(filename2, write = T)
#   for (j in 1:n.par) {
#     var.name <- para.OI.fates[j]
#     fates.para.values <- ncvar_get(fates_para, var.name)
#     if (length(dim(fates.para.values) == 1)) {
#       params.from.files[i, j] <-  fates.para.values[1]
#     } else {
#       params.from.files[i, j] <-  fates.para.values[1, 1]
#     }
#     } #j
#   nc_close(fates_para)
#   } #i
# 
# head(params.from.files)
# write.csv(params.from.files, "data-raw/params.from.files.csv", row.names = FALSE)
# params <- params.from.files
# params$fates_leaf_slatop <- 1/params$fates_leaf_slatop
# params$fates_leaf_slatop <- 2*params$fates_leaf_slatop
# params$fates_hydr_pitlp_node <- -params$fates_hydr_pitlp_node
# colnames(params) <- info$par.names
#--------
# now remove any older surf.sam.zip file, then zip surf.sam
file.remove("param.sam.zip")
f.to.zip <- list.files(surf.dir)
zip(zipfile = "param.sam.zip", files = file.path(surf.dir, f.to.zip))
# unzip(zipfile = "param.sam.zip")