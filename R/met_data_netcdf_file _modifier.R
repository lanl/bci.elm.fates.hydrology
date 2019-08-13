
# rm(list = ls())
# graphics.off()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, tidyverse)
path.folder <- "data-raw/ConvertMetCSVtoCLM-expand-format/NCOut/"
new.folder <- paste0(path.folder, 'bci_0.1x0.1_met.v5.1.test')
dir.create(new.folder)
files2copy <- list.files(path = paste0(path.folder, 'bci_0.1x0.1_met.v5.1/CLM1PT_data/'), 
                         full.names = TRUE)
file.copy(files2copy, new.folder, 
          overwrite = TRUE, recursive = TRUE)

baseStartY <- 2015
baseEndY <- 2015
n.baseY <- baseEndY - baseStartY + 1
startY <- 2016
endY <- 2016
n.FutureY <- endY - startY + 1

# deltaT = 0 #temperature increase
# deltaT.ann <- deltaT / n.FutureY
# var.name <- 'TBOT'
########################
#copy the file
#-------------------------------------------------
#create the files
# pb <- txtProgressBar(min = 1, max = n.FutureY, style = 3)

for (i in 1:n.FutureY) {
   Sys.sleep(0.000001)
   # update progress bar
   # setTxtProgressBar(pb, i)
   Yi = startY + i - 1
   modi <- i %% n.baseY
   if (modi == 0)
      modi = n.baseY
   baseYi = baseStartY + modi - 1
   for (j in 1:12) {
      if (j < 10) {
         filename1 <- paste0(new.folder, "/", baseYi, "-0", j, ".nc")
         filename2 <- paste0(new.folder, "/", Yi, "-0", j, ".nc")
      } else {
         filename1 <- paste0(new.folder, "/", baseYi, "-", j, ".nc")
         filename2 <- paste0(new.folder, "/", Yi, "-", j, ".nc")
      }
      #if(!file.exists(filename2))
      file.copy(filename1,
                filename2,
                overwrite = T,
                copy.mode = TRUE)
      # cmd = paste("chmod g+w ", filename2)
      # t1 <- try(system(cmd, intern = TRUE))
      climate_data <- nc_open(filename2, write = T)
      # var.values <- ncvar_get(climate_data, var.name)
      # var.values <- var.values + deltaT.ann * i
      # ncvar_put(climate_data, var.name, var.values)
      ncatt_get(climate_data,"time","units")$value
      if (j < 10) {
         timeunit <- paste0('days since ', Yi, "-0", j, "-01 00:00:00")
      } else {
         timeunit <- paste0('days since ', Yi, "-", j, "-01 00:00:00")
      }
      ncatt_put(climate_data, 'time', "units", timeunit, prec = "text")
      nc_close(climate_data)
      #}
   }
}

# files2zip <- list.files(path = new.folder, full.names = TRUE)
# zip(zipfile = paste0(path.folder, 'bci_0.1x0.1_met.v5.1.test.zip'), files = files2zip)

