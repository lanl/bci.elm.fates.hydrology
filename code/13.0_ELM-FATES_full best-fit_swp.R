#*******************************************
# Generating growth correlates from best-fit swp
# Author: Rutuja Chitra-Tarak
# First written: Oct 3, 2019
#*******************************************
# New model:
#   growthsp = a+ b * Btransp, 
# where Btransp=i=1i=zRootFracz, sp * f(swpz, tlpsp), 
# for z soil layers.

## creating table of Btran
library(groundhog)
groundhog.folder <- paste0("groundhog.library")
if(!dir.exists(file.path(groundhog.folder))) {dir.create(file.path(groundhog.folder))}
set.groundhog.folder(groundhog.folder)
groundhog.day = "2021-01-01"
pkgs=c("ncdf4", "easyNCDF", "lubridate", "tidyverse", "data.table", "doParallel", "foreach")
groundhog.library(pkgs, groundhog.day)

# graphics info
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

current.folder <- "2019-10-14_5000"
top.few <- 100
params.top.few.df <- read.csv(file.path("results", current.folder,
                                          paste0("params.top.few_", top.few, ".csv")), header = TRUE)
params.top.few <- params.top.few.df$x

nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
# depths
soil.depths <- as.numeric(nc$var[['H2OSOI']]$dim[[2]]$vals) # 15 depths in m

mini.date <- seq(from = as.Date("1990-01-01"), to = as.Date("2018-12-31"), by = "day")

#*******************************************
## get best-fit SOILPSI or swp ----
#*******************************************
load(file.path("data-raw/extract", current.folder, "best-fits/extract/SOILPSI.h1.extract.Rdata"))
# load(file.path("data-raw/extract", current.folder, "extract/SOILPSI.h1.extract.Rdata"))
psi.obs <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% gather(key = "n.best.par.sam", "psi", -date, -depth)
for(i in 2: length(var.res.arr)){
  psi.obs.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
    mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
    gather(key = "n.best.par.sam", "psi", -date, -depth)
  psi.obs <- psi.obs %>% bind_rows(psi.obs.i)
}
rm(psi.obs.i)
rm(var.res.arr) # large file
psi <- psi.obs %>% 
  mutate(n.best.par.sam = as.numeric(n.best.par.sam)) %>%
  left_join(data.frame(n.best.par.sam = 1:top.few, par.sam = params.top.few), 
            by = "n.best.par.sam") %>%
  select(-n.best.par.sam) %>%
  subset(par.sam %in% params.top.few) %>%
  mutate(date = as.Date(date), depth = signif(round(depth, 2), 2)) %>% 
  subset(date %in% mini.date) %>%   # As depths below 3 m are absent from some par.sam and has 0s associated with them; 
  # substituting those 0s with NAs, as there are no 0s otherwise 
  mutate_at(vars(psi), ~na_if(., -15)) %>% 
  subset(!is.na(psi)) %>% droplevels()
write.csv(psi, file = file.path("results", current.folder, "psi_model_daily_all_depths_params.top.few_full.csv"), row.names = FALSE)

usethis::use_data(psi, overwrite = TRUE)

psi.mean <- psi %>%
  group_by(date, depth) %>% summarise(psi = mean(psi, na.rm = TRUE)) %>%
  arrange(depth) %>% subset(!is.na(psi)) %>% droplevels()
usethis::use_data(psi.mean, overwrite = TRUE)

#*******************************************
#### BTRAN: Obs versus model #-------
#*******************************************
# BTRAN Daily #------
#*******************************************
load(file.path("data-raw/extract", current.folder, "best-fits/extract/BTRAN.h1.extract.Rdata"))

mod.btran.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date)) %>%
  as.data.frame()

rm(var.res.arr) # large file

btran.d <- mod.btran.d
head(btran.d[, 1:6]); summary(btran.d[, 1:6])
btran.d.long <- gather(btran.d, key = "n.best.par.sam", "value", -date) %>% 
  mutate(n.best.par.sam = as.numeric(n.best.par.sam)) %>%
  left_join(data.frame(n.best.par.sam = 1:top.few, par.sam = params.top.few), 
            by = "n.best.par.sam") %>%
  select(-n.best.par.sam)
btran <- btran.d.long %>% subset(date >= as.Date("1990-01-01") & 
                                              par.sam %in% params.top.few)
## mean + 95% CI
btran.stat <- btran %>% 
  group_by(date) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            upper.CI = quantile(value, probs = 0.95),
            lower.CI = quantile(value, probs = 0.05))
usethis::use_data(btran, overwrite = TRUE)
usethis::use_data(btran.stat, overwrite = TRUE)


## get growth factor = f(swp) ----

traits.indi <- read.csv("data-raw/hydraulic_traits_panama_kunert.csv") # Nobby's data
tlp <- traits.indi %>% group_by(sp) %>% select(-idividual, -ind_ID) %>%
  summarise(tlp = mean(mean_TLP_Mpa, na.rm = TRUE)) %>% 
  ungroup(sp) 
tlp.mean.positive <- -mean(tlp$tlp, na.rm = TRUE)

## IN Kuert's tlp data: swp ranges from -1.13 to -2.42 MPa 
# at -2.42 MPa f(swp) = 0
# from 0 to 0.5 f(swp) = 1, 
# f(swp) = linearly decreases from 0 to 1 f(swp) =  1/(2.42-0.5)*2.42 - 1/(2.42-0.5)*swp
# f(swp) =  1.260417 - 0.5208333*swp
## such that at -2.42 MPa f(swp) = 0

rel.swp.gfac <- data.frame(swp = seq(from = 0, to = -min(psi$psi), length.out = 1000)) %>%
  mutate(gfac = if_else(swp >=0 & swp < 0.5, 1, 
                        if_else(swp > tlp.mean.positive, 0,
                        if_else(swp == 0.5, 1, (tlp.mean.positive - swp)/(tlp.mean.positive - 0.5)))))
jpeg(file.path("figures", current.folder, "Growth_factor_swp_relationship.jpeg"), quality = 10, width = 360, height = 240)
ggplot(rel.swp.gfac, aes(y = gfac, x = swp)) + 
  geom_point() + scale_x_continuous(breaks = c(0, round(tlp.mean.positive, 1), 5, 10)) +
  ylab("Growth Factor") + xlab("Soil Water Potential (-MPa)") +
  ggtitle("SWP - Growth factor Relationship")
dev.off()

swp.gfac <- psi %>% mutate(psi.positive = -psi) %>%
  mutate(gfac = if_else(psi.positive >= 0 & psi.positive < 0.5, 1,
                         if_else(psi.positive == 0.5, 1, if_else(psi.positive > tlp.mean.positive, 0, 
                                                          (tlp.mean.positive - psi.positive)/(tlp.mean.positive - 0.5)))))
head(swp.gfac); range(swp.gfac$swp, na.rm =TRUE)
write.table(swp.gfac, file = file.path("results", current.folder, "swp.gfac_for_best-fit.txt"), row.names = FALSE)

usethis::use_data(swp.gfac, overwrite = TRUE)

#*******************************************
## soil water content ----
#*******************************************

load(file.path("data-raw/extract", current.folder, "best-fits/extract/H2OSOI.h1.extract.Rdata"))
# load(file.path("data-raw/extract", current.folder, "extract/SOILPSI.h1.extract.Rdata"))
swc.obs <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% 
  gather(key = "n.best.par.sam", "swc", -date, -depth)
for(i in 2: length(var.res.arr)){
  swc.obs.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
    mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
    gather(key = "n.best.par.sam", "swc", -date, -depth)
  swc.obs <- swc.obs %>% bind_rows(swc.obs.i)
}
rm(swc.obs.i)
rm(var.res.arr) # large file
swc <- swc.obs %>% 
  mutate(n.best.par.sam = as.numeric(n.best.par.sam)) %>%
  left_join(data.frame(n.best.par.sam = 1:top.few, par.sam = params.top.few), 
            by = "n.best.par.sam") %>%
  select(-n.best.par.sam) %>%
  subset(par.sam %in% params.top.few) %>%
  mutate(date = as.Date(date), depth = signif(round(depth, 2), 2)) %>% 
  subset(date %in% mini.date) %>%   # As depths below 3 m are absent from some par.sam and has 0s associated with them; 
  # substituting those 0s with NAs, as there are no 0s otherwise 
  mutate_at(vars(swc), ~na_if(., 0)) %>% 
  subset(!is.na(swc)) %>% droplevels()
write.csv(swc, file = file.path("results", current.folder, "swc_model_daily_all_depths_params.top.few_full.csv"), row.names = FALSE)

usethis::use_data(swc, overwrite = TRUE)


## Plant available water content ----

### swc - water content at wilting point
### water content at wilting point, using Clapp & Hornberger eqn 1978:

# In CLM:
#   /turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/biogeophysSoilStateType.F90

# this%watmin_col(c,lev) = &
#   this%watsat_col(c,lev)*(-min_liquid_pressure/this%sucsat_col(c,lev))**(-1._r8/this%bsw_col(c,lev))

# real(r8), pointer :: watsat_col           (:,:) ! col volumetric soil water at saturation (porosity)
# real(r8), pointer :: sucsat_col           (:,:) ! col minimum soil suction (mm) (nlevgrnd)
# real(r8), pointer :: bsw_col              (:,:) ! col Clapp and Hornberger "b" (nlevgrnd)

# real(r8), parameter :: min_liquid_pressure = -10132500._r8 ! Minimum soil liquid water pressure [mm]
# this%watsat_col(c,lev)    = 0.51_r8
# this%bsw_col(c,lev)       = 10_r8
# this%sucsat_col(c,lev)    = 200.0_r8

#thus, swc_wilt = watmin = 0.51*(10132500/200)^(-1/10) = 0.173

swc_wilt <- 0.17 # (volumetric water content)
load(file = file.path("data/swc.rda"))
swc.df <- swc %>% data.table()
paw <- swc.df[, ':='(paw = c(swc - swc_wilt))][, swc := NULL]
rm(swc, swc.df)
usethis::use_data(paw, overwrite = TRUE)


#*******************************************
## get best-fit GPP ----
#*******************************************
load(file.path("data-raw/extract", current.folder, "best-fits/extract/GPP.h1.extract.Rdata"))
# load(file.path("data-raw/extract", current.folder, "extract/SOILPSI.h1.extract.Rdata"))
mod.gpp.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date)) %>%
  as.data.frame()
mod.gpp.d[, -1] <- mod.gpp.d[, -1]*24*60*60 # converting from gC/m^2/s to gC/m^2/d

rm(var.res.arr) # large file

gpp.d <- mod.gpp.d
head(gpp.d[, 1:6]); summary(gpp.d[, 1:6])
gpp.d.long <- gather(gpp.d, key = "n.best.par.sam", "value", -date) %>% 
  mutate(n.best.par.sam = as.numeric(n.best.par.sam)) %>%
  left_join(data.frame(n.best.par.sam = 1:top.few, par.sam = params.top.few), 
            by = "n.best.par.sam") %>%
  select(-n.best.par.sam)
gpp <- gpp.d.long %>% subset(date >= as.Date("1990-01-01") & 
                                   par.sam %in% params.top.few)

write.csv(gpp, file = file.path("results", current.folder, "gpp_model_daily_all_depths_params.top.few_full.csv"), row.names = FALSE)

usethis::use_data(gpp, overwrite = TRUE)
devtools::document()

# install.packages("sinew")
# devtools::install_github("mdlincoln/docthis")
## To create documentaiton the first time
# file.create("R/data.R")
# sinew::makeOxygen(psi)
# sinew::makeOxygen(btran)
# sinew::makeOxygen(btran.stat)
# sinew::makeOxygen(psi.mean)
# sinew::makeOxygen(swc)
# sinew::makeOxygen(paw)
# sinew::makeOxygen(swp.gfac)
# sinew::makeOxygen(gpp)
