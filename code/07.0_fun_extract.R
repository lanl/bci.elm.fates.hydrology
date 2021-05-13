## works both for for vars like H2OSOI that have a vector of values (by depth) and also vars with one with  a single value

extractres_h0 <-
  function(sam.start,
           sam.end,
           outdir,
           filebase,
           var.vec,
           scale.vec,
           filterFile,
           start.year,
           end.year,
           current.folder) {
    nmonth <- 12
    nyears <- end.year - start.year + 1
    ncol <- nmonth * nyears
    # filter.arr <-
    #    read.table(file.path(wd, filterFile), header = F) # on server
    filter.arr <-
      read.table(file.path("data-raw/extract/", current.folder, filterFile), header = F) # on desktop
    # sam.vec <-
    #  c(sam.start:sam.end)[filter.arr$V1[c(sam.start:sam.end)]]
    sam.vec <-
      c(sam.start:sam.end)[filter.arr$V1]

    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h0.", start.year, "-")
    monstr <- sprintf("%02d", 1)
    filename <-
      paste0(outdir,
             "/",
             casename,
             "/run/",
             casename,
             ".",
             filetag,
             monstr,
             ".nc")
    nc <- nc_open(filename, write = F)
    res.arr <- vector("list", length = length(var.vec))
    var.length <- vector()
    for (v in 1:length(var.vec)) {
      var.name <- var.vec[v]
      val <- ncvar_get(nc, var.name)
      var.length[v] <- length(val)
      for (k in 1:var.length[v]) {
        res.arr[[v]][[k]] <- matrix(NA, nsam, ncol)
      }
    }
    
    ##----
    pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    cnames <- NULL
    for (yr in start.year:end.year) {
      for (j in 1:12) {
        cnames <- c(cnames, paste0(yr, "-", sprintf("%02d", j)))
      }
    }
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
      for (yr in start.year:end.year) {
        filetag <- paste("clm2.h0.", yr, "-", sep = "")
        for (j in 1:nmonth) {
          #months
          monstr <- sprintf("%02d", j)
          filename <-
            paste0(outdir,
                   "/",
                   casename,
                   "/run/",
                   casename,
                   ".",
                   filetag,
                   monstr,
                   ".nc")
          nc <- nc_open(filename, write = F)
          for (v in 1:length(var.vec)) {
            var.name <- var.vec[v]
            scale <- scale.vec[v]
            val <- ncvar_get(nc, var.name)
            yri <- yr - start.year
            indx <- j + nmonth * yri
            for (k in 1:var.length[v]) {
              res.arr[[v]][[k]][i, indx] <- val[k] * scale
            }
          }
          nc_close(nc)
        }
      }
      if(i %% 50 == 0) print(i)
    }
    for (v in 1:length(var.vec)) {
      var.name <- var.vec[v]
      for (k in 1:var.length[v]) {
        res.arr[[v]][[k]] <- as.data.frame(res.arr[[v]][[k]])
        colnames(res.arr[[v]][[k]]) <- cnames
        rownames(res.arr[[v]][[k]]) <- sam.vec
      }
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h0.extract.Rdata")
      # save(var.res.arr, file = file.path(wd, "extract", out.file.name)) # on server
      save(var.res.arr, file = file.path("data-raw", "extract", current.folder, "extract", out.file.name)) # on desktop
    }
  }

### FOR DAILY output -------
extractres_h1 <-
  function(sam.start,
           sam.end,
           outdir,
           filebase,
           var.vec,
           scale.vec,
           filterFile,
           start.year,
           end.year,
           current.folder) {
    cnames <- seq(from = as.Date(paste0(start.year, "-01-01")),
                       to = as.Date(paste0(end.year, "-12-31")),
                       by = 1)
    ## removing Feb 29 because CLM only produces output for 365 days not 366
    cnames <- cnames[format(cnames, "%m-%d") != "02-29"]
    ncol <- length(cnames)
    
    # filter.arr <- read.table(file.path(wd, filterFile), header = F) # on server
    filter.arr <-
      read.table(file.path("data-raw/extract/", current.folder, filterFile), header = F) # on desktop
    
    ### temporarily turning some cases off, since some cases from 1:100 are not well run:
    #filter.arr$V1[c(1:100, 740)] <- FALSE
    sam.vec <-
      c(sam.start:sam.end)#[filter.arr$V1[c(sam.start:sam.end)]] # turn this on after trial is finished
    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h1.", start.year, "-")
    filename <- paste0(outdir,
                       "/",
                       casename,
                       "/run/",
                       casename,
                       ".",
                       filetag,
                       "01-01-00000.nc")
    nc <- nc_open(filename, write = F)
    res.arr <- vector("list", length = length(var.vec))
    var.dim <- vector()
    for (v in 1:length(var.vec)) {
      var.name <- var.vec[v]
      val <- ncvar_get(nc, var.name)
      if (length(dim(val)) == 1) {
        var.dim[v] <- 1
      } else {
        var.dim[v] <- dim(val)[1]
      }
      for (k in 1:var.dim[v]) {
        res.arr[[v]][[k]] <- matrix(NA, nsam, ncol)
      }
    }
    ##----
    #pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    for (i in 1:nsam) {
      sample <- sam.vec[i]
      #setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
      for (yr in start.year:end.year) {
        filetag <- paste("clm2.h1.", yr, "-", sep = "")
        filename <-
          paste0(outdir,
                 "/",
                 casename,
                 "/run/",
                 casename,
                 ".",
                 filetag,
                 "01-01-00000.nc")
        nc <- nc_open(filename, write = F)
        for (v in 1:length(var.vec)) {
          var.name <- var.vec[v]
          scale <- scale.vec[v]
          val <- ncvar_get(nc, var.name)
          index.start <- (yr - start.year) * 365 + 1
          index.end <- index.start + 365 - 1
          if (var.dim[v] == 1) {
            res.arr[[v]][[1]][i, index.start:index.end] <- val * scale
          } else {
            for (k in 1:var.dim[v]) {
              res.arr[[v]][[k]][i, index.start:index.end] <- val[k,] * scale
            }
          } 	
        }
        nc_close(nc)
      }
      print(sample)
      # if (i %% 200 == 0)
      #  print(sample)
    }
    for (v in 1:length(var.vec)) {
      var.name <- var.vec[v]
      for (k in 1:var.dim[v]) {
        res.arr[[v]][[k]] <- as.data.frame(res.arr[[v]][[k]])
        colnames(res.arr[[v]][[k]]) <- cnames
        rownames(res.arr[[v]][[k]]) <- sam.vec
      }
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h1.extract.Rdata")
      # save(var.res.arr, file = file.path(wd, "extract", out.file.name)) # on server
      save(var.res.arr, file = file.path("data-raw", "extract", current.folder, "extract", out.file.name)) # on desktop
    }
  }
# load(file.path(wd, "data-raw", "extract", out.file.name))
### FOR HOURLY output -------
extractres_h2 <-
  function(sam.start,
           sam.end,
           outdir,
           filebase,
           var.vec,
           scale.vec,
           start.year,
           end.year,
           current.folder) {
    cnames <-
      seq(from = as.POSIXct(paste0(start.year, "-01-01 00:00:00", tz = "UTC")),
               to = as.POSIXct(paste0(end.year, "-12-31 23:00:00", tz = "UTC")),
               by = "hour")
    cnames <- cnames[format(cnames, "%m-%d") != "02-29"]
    ncol <- length(cnames)
    
    # filter.arr <-
    #   read.table(file.path(wd, filterFile), header = F) # on server
    filter.arr <-
      read.table(file.path("data-raw/extract/", current.folder, filterFile), header = F) # on desktop
    sam.vec <-
      sam.vector[filter.arr$V1[sam.vector]]
    nsam <- length(sam.vec)
    
    ## to get the length of value vector (for H2OSOI the number of depths)#-----
    casename <- paste(filebase, sam.vec[1], sep = "")
    filetag <- paste0("clm2.h2.", start.year, "-")
    filename <- paste0(outdir,
                       "/",
                       casename,
                       "/run/",
                       casename,
                       ".",
                       filetag,
                       "01-01-00000.nc")
    nc <- nc_open(filename, write = F)
    val <- ncvar_get(nc, var.name)
    res.arr <- list()
    if (length(dim(val)) == 1) {
      var.dim <- 1
    } else {
      var.dim <- dim(val)[1]
    }
    for (k in 1:var.dim) {
      res.arr[[k]] <- as.data.frame(matrix(NA, nsam, ncol))
    }
    
    ##----
    pb <- txtProgressBar(min = 0, max = nsam, style = 3)
    # days_in_year <- function(year) {
    #   365 + (year %% 4 == 0) - (year %% 100 == 0) + (year %% 400 == 0)
    # }
    for (i in 1:nsam) {
      sample <- sam.vec[sample]
      setTxtProgressBar(pb, i)
      casename <- paste(filebase, sample, sep = "")
        for (yr in start.year:end.year) {
          filetag <- paste("clm2.h2.", yr, "-", sep = "")
          filename <-
            paste0(outdir,
                   "/",
                   casename,
                   "/run/",
                   casename,
                   ".",
                   filetag,
                   "01-01-00000.nc")
          nc <- nc_open(filename, write = F)
          val <- ncvar_get(nc, var.name)
          index.start <- 24 * (yr - start.year) * 365 + 1
          index.end <- index.start + 365 * 24 - 1
          if (var.dim == 1) {
            res.arr[[1]][i, index.start:index.end] <- val * scale
          } else {
            for (k in 1:var.dim) {
              res.arr[[k]][i, index.start:index.end] <- val[k,] * scale
            }
          }
          nc_close(nc)
        }
      if (i %% 50 == 0) print(i)
    }
    for (k in 1:var.dim) {
      colnames(res.arr[[k]]) <- cnames
      rownames(res.arr[[k]]) <- sam.vec
      # data is stored as a list
      var.res.arr <- res.arr[[v]]
      out.file.name <- paste0(var.name, ".h2.extract.Rdata")
      # save(var.res.arr, file = file.path(wd, "extract", out.file.name)) # on server
      save(var.res.arr, file = file.path("data-raw", "extract", current.folder, "extract", out.file.name)) # on desktop
    }
  }
###-------
# load(file.path(wd, "data-raw", "extract", out.file.name))
## old code which only works for vars with single value and stores into tables by nsam
# extractres_h0 <- function(nsam, outdir, filebase, var.name, filterFile, ###----
#                           scale, start.year, end.year, out.file.name){
#   #outdir<-"/lustre/scratch3/tursampleuoise/rutuja/ACME/cases"
#   #filebase<-"BCI.ICLM45ED.badger.intel.Cabfd064cb-Fc47f968e."
#   #nsam<-8
#   nmonth <- 12
#   nyears <- end.year - start.year + 1
#   ncol <- nmonth*nyears
#   #scale = 30*24*3600 #rain from mm/s -> mm/month
#   #filter.arr<-rep(F,nsam)
#   res.arr <- as.data.frame(matrix(0, nsam, ncol))
#   pb <- txtProgressBar(min = 0, max = nsam, style = 3)
#   cnames <- NULL
#   for (yr in start.year: end.year){
#     for(j in 1:12){
#       cnames <- c(cnames, paste(var.name, "y", yr, "m", j, sep=""));
#     }
#   }
#   #filter.arr <- read.table(file.path(wd, filterFile), header = F)
#   filter.arr <- read.table(file.path(wd, "data-raw", filterFile), header = F)
#   filter.range <- filter.arr[1: nsam, 1] == 1
#
#   for(i in 1: nsam){
#     setTxtProgressBar(pb, i)
#     casename <- paste(filebase, i, sep ="")
#       for (yr in start.year:end.year){
#         filetag <- paste("clm2.h0.", yr, "-", sep = "")
#         for(j in 1:nmonth){ #months
#           monstr <- sprintf("%02d", j)
#           filename <- paste(outdir, "/", casename, "/run/", casename, ".",
#                             filetag, monstr,".nc", sep = "")
#           nc <- nc_open(filename, write = F)
#           val <- ncvar_get(nc, var.name)
#           yri <- yr - start.year
#           indx <- j + nmonth*yri
#           res.arr[i, indx] <- val*scale
#           nc_close(nc)
#         }
#       }
#   }
#   colnames(res.arr) <- cnames
#   #write.table(res.arr[filter.range,],file.path(wd,"extract",out.file.name),sep="\t",row.names=F,col.names=T)
#   write.table(res.arr[filter.range,], file.path(wd, "data-raw", "extract", out.file.name),
#               sep="\t", row.names = F, col.names = T)
# }
