##----------------
# Choosing Best-fit Simulated ecosystem soil depth for hydrological water-balance
# Author: Rutuja Chitra-Tarak
# Date: May 8, 2019
##-----------------

#### Tasks:
## Load model outputs for all soil depths: 
# Convert .nc file from server to xxx
## Retrieve ET, Qrunoff and soil moisture comparable to depths 10, 40 and 100 cm to match against obervations
## Check fit for each soil depth
## Choose the best

rm(list=ls())
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(smooth, ncdf4, tidyverse, bci.hydromet, data.table, avasoilmoisture, hydroGOF, SemiPar)
## Somehow SemiPar does not get installed with Pacman
# install.packages("~/Downloads/SemiPar_1.0-4.2-2.tar", repos = NULL, type="source")
# library(SemiPar)
# pacman::p_load(tidyverse)
# H2OSOI	volumetric soil water (vegetated landunits only)	mm3/mm3

# QRUNOFF	total liquid runoff (does not include QSNWCPICE)	mm/s

## Total evapotransipration:
# QSOIL	Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)	mm/s
# QVEGE	canopy evaporation	mm/s
# QVEGT	canopy transpiration	mm/s
##--------------------------
## Setting ggplot theme
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             strip.background = element_blank()
)
##--------------------------

current.folder <- "2019-10-14_5000"
figures.folder <- paste0(current.folder, "/best-fits")
top.few <- 100
params.top.few.cond.df <- read.csv(file.path("results", current.folder, paste0("params.top.few.cond_", top.few, ".csv")), header = TRUE)
params.top.few <- params.top.few.cond.df$x
col1 <- c("Observed" = "#f04546", "Observed Gap-filled" = "purple", 
          "Observed Point Location" = "#f04546", "Observed Plot-wide" = "black",
          "Observed Point Location2" = "#3591d1", 
          "Simulated" = "#3591d1", # "Simulated" = "green", 
          "Best-fit RMSE" = "#3591d1", "Best-fit NSE" = "#3591d1", "Best-fit R-squared" = "purple", 
          "Observed at Ava-Tower" = "#f04546")

## for plotting RMSE
params.obj.top.few.cond <- read.csv(file.path("results", current.folder, paste0("params.obj.top.few.cond_", top.few, ".csv")), header = TRUE)

###
#### Runoff: Obs versus model #-------
####
# Daily #------
####
## QRUNOFF total liquid runoff (does not include QSNWCPICE)
## QUNOFF = QOVER + QDRAI 
# QOVER: surface runoff mm/s
# QDRAI: sub-surface drainage mm/s
load(file.path("data-raw/extract", current.folder, "extract/QRUNOFF.h1.extract.Rdata"))
mod.qrunoff.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date))
mod.qrunoff.d[, 2:ncol(mod.qrunoff.d)] <- mod.qrunoff.d[, 2:ncol(mod.qrunoff.d)]*24*60*60/3 # converting mm/s to mm/day
obs.qrunoff.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day
  rename(obs = flow_conrad) %>% subset(date > min(mod.qrunoff.d$date, na.rm = TRUE))

qrunoff.d <- obs.qrunoff.d %>% full_join(mod.qrunoff.d, by = "date")
#head(qrunoff.d); #summary(qrunoff.d)
rm(var.res.arr) # large file
# letting four years pass to allow for model initialisation before model-obs comparison
qrunoff.d.sub <- qrunoff.d %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4 + 1)) 

qrunoff.d.long <- gather(qrunoff.d.sub, key = "par.sam", "value", -date, -obs) %>% 
  subset(par.sam %in% params.top.few)

p.qrun.d.1 <- ggplot(qrunoff.d.long,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.05) +
  scale_colour_manual(name = "", values = col1) +
  geom_line(aes(y = obs, colour = "Observed"), size = 0.5) +
  ylab(expression("Stream Discharge (mm."*day^-1*")")) +
  xlab("Time") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.7, 0.95), 
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) #+
 # ggtitle("QRUNOFF: Observed vs Simulated_Daily")
ggsave("QRUNOFF_Obs_vs_model_daily.jpeg", plot = p.qrun.d.1, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units='in')

####
# Monthly #------
####

qrunoff.m.long <- qrunoff.d.long %>%
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m")),
         value.na = if_else(is.na(obs), obs, value)) %>%
  group_by(yrmo, par.sam) %>%
  ## when all in a month are NAS, na.rm =TRUE makes the sum zero
  # Correcting that
  summarise(obs = if_else(all(is.na(obs)), sum(obs), sum(obs, na.rm = TRUE)),
            value = if_else(all(is.na(obs)), sum(value.na), sum(value.na, na.rm = TRUE))) %>%
  as.data.frame()

## rsq for the ensemble with best rmse
rmse.qrun.m <- round(mean(1 - params.obj.top.few.cond$qrunoff_monthly), 2)

qrunoff.m.long <- qrunoff.m.long %>% 
  mutate(date = as.Date(paste0(yrmo, "-01")))
p.qrun.m1 <- ggplot(qrunoff.m.long,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  scale_colour_manual(name = "", values = col1) +
  ylab(expression("Stream Discharge (mm."*month^-1*")")) +
  xlab("Time") + 
  geom_text(aes(x = qrunoff.m.long$date[300], y = 500, label = rmse.qrun.m), 
            size = 6, vjust = "inward", hjust = "inward") +
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.3, 0.9), 
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) #+
  # ggtitle("QRUNOFF: Observed vs Simulated_Monthly")
ggsave("QRUNOFF_Obs_vs_model_monthly.jpeg", plot = p.qrun.m1, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units='in')

#--------------
####
#### Evapotranspiration: Obs versus model #-------
####
# Daily #------
####
load(file.path("data-raw/extract", current.folder, "extract/QVEGE.h1.extract.Rdata"))
mod.qvege.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") # in mm/s
load(file.path("data-raw/extract", current.folder, "extract/QVEGT.h1.extract.Rdata"))
mod.qvegt.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") # in mm/s
load(file.path("data-raw/extract", current.folder, "extract/QSOIL.h1.extract.Rdata"))
mod.qsoil.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") # in mm/s
rm(var.res.arr) # large file
## combining evaporation from soil and plants and transpiration from plants
flow.columns <- (mod.qvege.d[, -1] + mod.qvegt.d[, -1] + mod.qsoil.d[, -1])*24*60*60 # from mm/s to mm/day
mod.qet.d <-  cbind.data.frame(mod.qvege.d[, 1], flow.columns) %>% 
  mutate(date = as.Date(date))
## Observed ET from flux tower
# AET.flag.day has NAs substituted where insufficient (< 50%) actual data for the day
obs.qet.d <- bci.hydromet::forcings %>% select(date, AET, AET.flag.day) %>% # in mm/day
  rename(obs = AET.flag.day, 
         obs.2 = AET) %>% # this is gap filled AET 
  subset(date > "2012-07-02" & date < "2017-09-01") ## date after which ET data is present
head(obs.qet.d)
# joining obs with mod
qet.d <- obs.qet.d %>% full_join(mod.qet.d, by = "date") 
# does not have to be done for ET because observed ET is 3.5 years ahead of the run start date, which is 2008-01-1
qet.d.sub <- qet.d %>% subset(date >= c(min(mod.qet.d$date, na.rm = TRUE) + 365*4 + 1))
qet.d.long <- gather(qet.d.sub, key = "par.sam", "value", -date, -obs, -obs.2) %>% 
  subset(par.sam %in% params.top.few)

qet.d.long <- qet.d.long %>% subset(date >= as.Date(min(obs.qet.d$date)) & date < as.Date(max(obs.qet.d$date)) + 1) %>% droplevels() %>%
 mutate(obs.3 = if_else(is.na(obs), obs.2, obs))
qet.d.long.sub <- qet.d.long 

p.et1 <- ggplot(qet.d.long.sub,
                aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.1) +
  scale_colour_manual(name = "", values = col1) +
  ylab(expression("Evapotranspiration (mm."*day^-1*")")) +
  xlab("Time") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        plot.margin = margin(5.1, 4.1, 4.1, 2.1, "pt"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
        legend.position = c(0.75, 0.95), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) 
p.et1.a.d <- p.et1 + 
  geom_line(aes(x= date, y = obs.2, color = "Observed"), size = 0.5) #+ 
  # ggtitle("Evapotranspiration: Observed vs Simulated_Daily")
ggsave("ET_Obs_vs_model_daily_all_years.jpeg", plot = p.et1.a.d, 
       path = file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units='in')
p.et1.b.d <- p.et1 + 
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.5) #+ 
  # ggtitle("Evapotranspiration: Observed vs Simulated_Daily")
ggsave("ET_Obs_vs_model_daily_all_years_with_gaps.jpeg", plot = p.et1.b.d, 
       path = file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units='in')
######
##### lines and points format
######

col2 <- c("Observed" = "#f04546", "Observed Gap-filled" = "purple", "Simulated" = "#3591d1")
p.et1.p <- p.et1 +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs.2, color = "Observed"), size = 0.2, linetype = 5) +  
  geom_point(aes(x= date, y = obs.3, color = "Observed Gap-filled"), size = 1, shape = 1) +  
  geom_point(aes(x= date, y = obs, color = "Observed"), size = 1, shape = 1) +  
  theme(legend.position = "top") +
  scale_colour_manual(name = "", values = col2)
ggsave(paste0("ET_Obs_vs_model_daily_all_years_points_lines.jpeg"), plot = p.et1.p, 
       path = file.path("figures", figures.folder), device = "jpeg", height = 3, width = 15, units='in')

####
# Monthly #------
## using 

# qet.m.long <- qet.d.long %>%  
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
#   group_by(yrmo, par.sam) %>%
#   ## using gap-filled AET
#   summarise(obs = sum(obs.2), ## removing na.rm = TRUE so as to only get sum for months when all observations present
#             value = sum(value, na.rm = T), date = min(date)) %>% 
#   as.data.frame()
qet.m.long <- qet.d.long %>%
  subset(date >= as.Date(min(obs.qet.d$date)) & date < as.Date(max(obs.qet.d$date)) + 1) %>%
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m")),
         value.na = if_else(is.na(obs), obs, value)) %>%
  group_by(yrmo, par.sam) %>%
  ## when all in a month are NAS, na.rm =TRUE makes the sum zero
  # Correcting that
  summarise(obs = if_else(all(is.na(obs)), sum(obs), sum(obs, na.rm = TRUE)),
            value = if_else(all(is.na(obs)), sum(value.na), sum(value.na, na.rm = TRUE))) %>%
  as.data.frame()

# letting four years pass to allow for model initialisation before model-obs comparison
qet.m <- qet.m.long %>% 
  pivot_wider(names_from = par.sam, values_from = value) %>% as.data.frame()

qet.m.long.sub <- qet.m.long %>% mutate(date = as.Date(paste0(yrmo, "-01")))

rmse.qet.m <- round(mean(1 - params.obj.top.few.cond$qet_monthly), 2)

# rmse.label <- as.character(as.expression(italic(r)^2~"="~rmse.max))
p.et.m1 <- ggplot(qet.m.long.sub, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  scale_colour_manual(name = "", values = col1) +
  ylab(expression("Evapotranspiration (mm."*month^-1*")")) +
  xlab("Time") +
  geom_text(aes(x = qet.m.long.sub$date[300], y = 130, label = rmse.qet.m), 
            size = 6, vjust = "inward", hjust = "inward") +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.25, 0.95), 
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) #+
  # ggtitle("Evapotranspiration: Observed vs Simulated_Monthly")
ggsave("ET_Obs_vs_model_monthly.jpeg", plot = p.et.m1, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units='in')

###
#### GPP: Obs versus model #-------
####
# Daily #------
####
load(file.path("data-raw/extract", current.folder, "extract/GPP.h1.extract.Rdata"))

mod.gpp.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date)) %>%
  as.data.frame()
mod.gpp.d[, -1] <- mod.gpp.d[, -1]*24*60*60 # converting from gC/m^2/s to gC/m^2/d
summary(mod.gpp.d[, 1:3])

rm(var.res.arr) # large file
bci.tower <- read.csv("data-raw/BCI_v3.1.csv", header = TRUE)
bci.tower <- as.data.frame(bci.tower[-1, ])
bci.tower$datetime <- strptime(bci.tower$date, format = "%m/%d/%Y %H:%M")
bci.tower$datetime <- as.POSIXct(bci.tower$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "")
str(bci.tower$datetime)

obs.gpp.d <- bci.tower %>% select(datetime, gpp) %>% # in mumol/m2/2: units must be mumol/m2/s
  mutate(date = as.Date(format(datetime, "%Y-%m-%d")),  
         gpp.mumol = as.numeric(as.character(gpp))) %>%
  group_by(date) %>% summarise(obs = mean(gpp.mumol, na.rm = T)) %>%
  mutate(obs = obs*12*1e-06*24*60*60) %>% # gC/m2/d %>%
  subset(date != is.na(date))
# mutate(obs = as.numeric(bigleaf::umolCO2.to.gC(obs))) ## gives the same result #-------
# gC/m^2/d
# to convert mumol/m2/s to gC/m^2/d
# Cmol*mumol2mol*kg2g/days2seconds 
# constants	Cmol - molar mass of carbon (kg mol-1) 
# umol2mol - conversion micromole (umol) to mol (mol) 
# kg2g - conversion kilogram (kg) to gram (g) 
# days2seconds - seconds per day
# mean gpp.mumol = 10696 in mumol/m2/s #----
# so mean gpp in gC/m^2/d
# conversion.eq : 
summary(obs.gpp.d)
gpp.d <- obs.gpp.d %>% full_join(mod.gpp.d, by = "date") %>% as.data.frame()
head(gpp.d[, 1:6]); summary(gpp.d[, 1:6])
# letting four years pass to allow for model initialisation before model-obs comparison
gpp.d.sub <- gpp.d %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4 + 1))
gpp.d.long <- gather(gpp.d, key = "par.sam", "value", -date, -obs) %>% 
  subset(par.sam %in% params.top.few)
gpp.d.long.sub <- gpp.d.long %>% subset(date >= min(obs.gpp.d$date) & date < max(obs.gpp.d$date) + 1)

p.gpp.d1 <- ggplot(gpp.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.1) +
  scale_colour_manual(name = "", values = col1) +
  geom_line(aes(y = obs, colour = "Observed"), size = 0.5) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Observed vs Simulated_Daily")
ggsave( "GPP_Obs_vs_model_daily.jpeg", plot = p.gpp.d1, path = file.path("figures", figures.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

######
##### lines and points format
######

p.gpp1 <- ggplot(gpp.d.long.sub, aes(x = date, y = value)) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        plot.margin = margin(5.1, 4.1, 4.1, 2.1, "pt"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) 
p.gpp1.p <- p.gpp1 +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.2, linetype = 5) +  
  geom_point(aes(x= date, y = obs, color = "Observed"), size = 1, shape = 1) +  
  theme(legend.position = "top") +
  scale_colour_manual(name = "", values = col2)
ggsave(paste0("GPP_Obs_vs_model_daily_all_years_points_lines.jpeg"), plot = p.gpp1.p, path = file.path("figures", figures.folder), device = "jpeg", height = 3, width = 15, units='in')

####
# Monthly #------
####

gpp.m.long <- gpp.d.long %>%
  # letting four years pass to allow for model initialisation before model-obs comparison
  subset(date >= as.Date(min(obs.gpp.d$date)) & date < as.Date(max(obs.gpp.d$date) + 1)) %>%
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m")),
         value.na = if_else(is.na(obs), obs, value)) %>%
  group_by(yrmo, par.sam) %>%
  ## when all in a month are NAS, na.rm =TRUE makes the sum zero
  # Correcting that
  summarise(obs = if_else(all(is.na(obs)), sum(obs), sum(obs, na.rm = TRUE)),
            value = if_else(all(is.na(obs)), sum(value.na), sum(value.na, na.rm = TRUE))) %>%
  as.data.frame()

gpp.m.long <- gpp.m.long %>% mutate(date = as.Date(paste0(yrmo, "-01")))
p.gpp.m.1 <- ggplot(gpp.m.long,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  scale_colour_manual(name = "", values = col1) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Observed vs Simulated_Monthly")
ggsave("GPP_Obs_vs_model_monthly.jpeg", plot = p.gpp.m.1, path = file.path("figures", figures.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

##
#--------------
####
#### Soil Moisture: Obs versus model #-------
####
# Daily #------
# bci.hydromet::TDR
# bci.hydromet::sm.historic
####
nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
# depths
soil.depths <- nc$var[['H2OSOI']]$dim[[2]]$vals 
rm(nc)
# 15 depths
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934
# so depths 4 [ 0.11 m], 6[0.36 m] and 8[1.03 m] roughly correspond to 10, 40 and 100 cm
# in units m3/m3

# On Thu, Sep 14, 2017 at 4:55 AM, Rutuja Chitra-Tarak <arutuj@gmail.com> wrote:
#   Hi Matteo, 
# A quick question: what are the depths at which you have soil moisture sensors deployed? Both vertically and horizontally?
#   Thanks,
# r
# The vertical (probe 1-2-3 are representative of 0-15 cm, the horizontal 3-4-5 are 10, 40 and 100 cm)
# 
# MTT

load(file.path("data-raw/extract", current.folder, "extract/H2OSOI.h1.extract.Rdata"))
mod.swc <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% 
  gather(key = "par.sam", "value", -date, -depth) %>% 
  subset(par.sam %in% params.top.few) %>% droplevels()
for(i in 1: length(var.res.arr)){
  mod.swc.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
    mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
    gather(key = "par.sam", "value", -date, -depth) %>% 
    subset(par.sam %in% params.top.few) %>% droplevels()
  mod.swc <- mod.swc %>% bind_rows(mod.swc.i)
}
rm(mod.swc.i, var.res.arr) # large file

### adding stephan Kupers data: #-------
swp.df <- read.table(
  file.path("data-raw/Kupers_et_al/Input/BCI_Soil_moisture_mapping.txt"),
  na.strings = c("NA", ""),
  header = T,
  row.names = NULL,
  check.names = F
)

summary(swp.df)
unique(swp.df$depth)
range(swp.df$swc, na.rm = T)
swp.df$depth.old <- swp.df$depth
swp.df$depth <- cut(swp.df$depth.old, breaks = c(0, 30, 70, 120), labels = c(10, 40 , 100)) # It should really be c(15, 40 , 100
summary(swp.df)
swp.df[is.na(swp.df$depth),]
swp.df <- swp.df[!is.na(swp.df$depth),]
  
steph.sm <- swp.df %>%
  select(date, depth, swc) %>% 
  mutate(date = as.Date(date, format = "%m/%d/%Y"),
         data = "Observed Plot-wide",
         depth = as.numeric(as.character(depth)),
         ## swc data is gravimetric. Need to be volumetric to comapre with model
         ## with 0.8 g/cm^3 bulk dry density at 15 cm depth. Data by Matteo: data/rawHydrological model Barro Colorado Island.docx
         swc = swc*0.8/100, # converting from % to v/v fraction 
         yrmo = format(date, format = "%Y-%m")) 
steph.sm <- steph.sm %>% 
  group_by(depth, yrmo) %>% 
  mutate(mean.date = mean(date, na.rm = TRUE)) %>%
  ungroup(depth, yrmo)

str(steph.sm); summary(steph.sm)
steph.quant <- steph.sm %>% 
  group_by(depth, mean.date) %>% 
  summarise(mean = mean(swc, na.rm = TRUE), 
            lower = quantile(swc, na.rm = TRUE, probs = c(0.05)),
            upper = quantile(swc, na.rm = TRUE, probs = c(0.95)),
            q.25 = quantile(swc, na.rm = TRUE, probs = c(0.25)),
            q.75 = quantile(swc, na.rm = TRUE, probs = c(0.75))) %>%
  ungroup(depth, mean.date)

####-----
swc.vertical <- avasoilmoisture::vertical %>% rename(obs = swc) %>% 
  group_by(date, depth) %>% summarise(obs = mean(obs, na.rm = TRUE)) %>% 
  mutate(date.depth = paste0(date, ".", depth))
swc.vertical.mean <-  swc.vertical %>% group_by(date) %>% summarise(obs = mean(obs, na.rm = TRUE)) 

mod.swc.d.vert <- mod.swc %>%  
  # subset(date >= min(as.Date(swc.vertical$date))) %>% 
  subset(date >= min(as.Date(swc.vertical$date)) & depth %in% soil.depths[1:4]) %>% ## date filter since otherwise vector memory exhausted
  mutate(date.depth = paste0(date, ".", depth)) %>% 
  group_by(date, par.sam) %>%
  summarise_at("value", list(~ mean(., na.rm = T))) %>% as.data.frame()

## joining vertical obs with mod
swc.d.long.vert <- mod.swc.d.vert %>%
  right_join(swc.vertical.mean, by = "date") %>% 
  as.data.frame() 
head(swc.d.long.vert); summary(swc.d.long.vert)

# with depths 10, 40 & 100 cm
mod.swc.d <- mod.swc %>% subset(depth %in% soil.depths[c(4, 6, 8)]) %>%
  mutate(depth = as.character(depth))
mod.swc.d$depth <- recode(mod.swc.d$depth, "0.118865065276623" = 10, 
                          "0.366065800189972" = 40, "1.03802704811096" = 100)

mod.swc.d <- mod.swc.d %>% mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()

obs.swc.d <- avasoilmoisture::horizontal %>% 
  rename(obs = swc) %>% group_by(date, depth) %>% summarise_at(vars(obs), ~mean(., na.rm = TRUE)) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() #%>%
# bind_rows(swc.single) %>% mutate(depth = as.numeric(depth))
## To compare with plot-wide data points
obsplot.swc.d <- steph.quant %>% 
  rename(date = mean.date, obs = mean) %>%
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() 
swcplot.d.long <- mod.swc.d %>% 
  right_join(obsplot.swc.d %>%
              select(date.depth, obs), by = "date.depth") %>% 
  as.data.frame() 
head(swcplot.d.long); summary(swcplot.d.long)

swc.d.long <- mod.swc.d %>% 
  right_join(obs.swc.d %>% select(-date, -depth), by = "date.depth") %>% 
  as.data.frame()
head(swc.d.long); summary(swc.d.long)

swc.range.best <- mod.swc.d %>% subset(date >= min(obs.swc.d$date)) %>% 
  group_by(depth) %>% 
  summarise(lower = min(value, na.rm = TRUE), upper = max(value, na.rm = TRUE))
swc.d.long <- swc.d.long %>% left_join(swc.range.best, by = "depth")
swc.d.long <- swc.d.long %>% subset(!is.na(depth))
swc.d.long <- swc.d.long %>% 
  group_by(depth) %>% 
  mutate(sat = scales::rescale(obs, to = c(lower[1], upper[1]))) %>% as.data.frame() 

swc.d.long.sub <- subset(swc.d.long, !depth %in% c(5, 300)) %>% droplevels()
swc.d.vert.long.sub <- subset(swc.d.long.vert, date >= as.Date("2012-07-01") & date < as.Date("2019-01-01")) 

data.backto.steph <- mod.swc.d %>% subset(date >= as.Date("2014-01-01"))

swc.d.long.sub <- swc.d.long.sub %>% 
  mutate(depth.plot = recode_factor(as.factor(depth), `10` = "0.1 m", `40` = "0.4 m",
                                    `100` = "1 m"))  %>% 
  transform(depth.plot = factor(depth.plot, 
                                levels = c("0.1 m", "0.4 m", "1 m"), ordered = TRUE))

steph.quant <- steph.quant %>%  
  mutate(depth.plot = recode_factor(as.factor(depth), `10` = "0.1 m", `40` = "0.4 m",
                                    `100` = "1 m"))  %>% 
  transform(depth.plot = factor(depth.plot, 
                                levels = c("0.1 m", "0.4 m", "1 m"), ordered = TRUE))


label.depths.hor <- round(1 - c(mean(params.obj.top.few.cond$sat10_daily), mean(params.obj.top.few.cond$sat40_daily),
                                            mean(params.obj.top.few.cond$sat100_daily)), 2)
dat_text.hor <- data.frame(label = label.depths.hor, depth.plot = c("0.1 m", "0.4 m", "1 m"))
 
label.depths.plot <- round(1 - c(mean(params.obj.top.few.cond$swcplot10_daily), mean(params.obj.top.few.cond$swcplot40_daily),
                                                 mean(params.obj.top.few.cond$swcplot100_daily)), 2)
dat_text.plot <- data.frame(label = label.depths.plot, depth.plot = c("0.1 m", "0.4 m", "1 m"))

g1 <- ggplot(swc.d.long.sub %>% 
               subset(par.sam == params.top.few[1] & depth %in% c(10, 40, 100)), aes(x = date)) +
  geom_line(data = data.backto.steph, aes(y = value, group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  ## adding Observed normalised
  geom_line(aes(y = sat, color = "Observed Point Location"), size = 0.5) +
  scale_colour_manual(name = "", values = col1) + ylim(0.18, 0.57) +
  ylab(expression('Soil Water Content ('*cm^3~cm^-3*')')) +
  xlab("Time") +
  geom_text(data = dat_text.plot, mapping = aes(x = data.backto.steph$date[200], 
                                                y = 0.57, label = label),
            size = 6, vjust = "inward", hjust = "inward") +
  geom_text(data = dat_text.hor, mapping = aes(x = data.backto.steph$date[200], 
                                           y = 0.2, label = label), 
            size = 6, vjust = "inward", hjust = "inward") +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.background = element_rect(fill = "transparent")) #+
  #ggtitle("Soil Water Content: Observed vs Ensemble simulations by depth")
p.h.swc.steph <- g1 + facet_grid(depth.plot ~ .) + 
  theme(legend.position = "top")
p.h.swc.steph.mean <- p.h.swc.steph +
  geom_point(data = steph.quant, aes(x = mean.date, y = mean, color = "Observed Plot-wide"),
             shape = 21, color = "white", fill = "black", alpha = 0.9, size = 1) +
  geom_errorbar(data = steph.quant, aes(x = mean.date, ymin = lower, ymax = upper, color = "Observed Plot-wide"), width = 0.1) +
  theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y")) 
ggsave("swc_Obs_vs_model_daily_horizontal_with_steph_mean.jpeg", plot = p.h.swc.steph.mean, 
       path = file.path("figures", figures.folder), device = "jpeg", height = 4.5, width = 6, units ='in')

### 
## With vertical:
rmse.v.swc.d <- round(mean(1 - params.obj.top.few.cond$swc.vert_daily), 2)
  
p.v.swc <- ggplot(swc.d.vert.long.sub, aes(x = date)) +
  geom_line(aes(y = value, group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  ## adding Observed normalised
  geom_line(aes(y = obs, color = "Observed Point Location"), size = 0.4) +
  scale_colour_manual(name = "", values = col1) + 
  ylab(expression('Soil Water Content ('*cm^3*cm^-3*')')) +
  xlab("Time") + 
  geom_text(aes(x = swc.d.vert.long.sub$date[200], y = 0.57, label = rmse.v.swc.d), 
            size = 6, vjust = "inward", hjust = "inward") +
  # ggtitle("Soil Water Content: Observed vs Ensemble simulations by depth") +
  theme(legend.justification = c(1, 1),  legend.position = c(0.85, 0.95), legend.direction = "horizontal",
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"))
ggsave("swc_Obs_vs_model_daily_vertical.jpeg", plot = p.v.swc, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 3, width = 6, units ='in')

## Plotting jointly

p.h.swc.steph.mean.joint <- ggplot(swc.d.long.sub %>% 
               subset(par.sam == params.top.few[1] & depth %in% c(10, 40, 100)), aes(x = date)) +
  geom_line(data = data.backto.steph, aes(y = value, group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(y = sat, color = "Observed"), size = 0.5) +
  scale_colour_manual(name = "", values = col1) + ylim(0.18, 0.57) +
  scale_fill_manual(name = "", values = col1) +
  ylab(expression('Soil Water Content ('*cm^3~cm^-3*')')) +
  xlab("Time") +
  geom_text(data = dat_text.plot, mapping = aes(x = data.backto.steph$date[50], 
                                                y = 0.57, label = label),
            size = 4, vjust = "inward", hjust = "inward") +
  geom_text(data = dat_text.hor, mapping = aes(x = data.backto.steph$date[50], 
                                               y = 0.2, label = label), 
            size = 4, vjust = "inward", hjust = "inward") +
  theme(axis.text = element_text(size = 14, face = "plain"),
        strip.text.y = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.background = element_rect(fill = "transparent")) +
  facet_grid(depth.plot ~ .) + 
  theme(legend.position = "top") +
  geom_errorbar(data = steph.quant, aes(x = mean.date, ymin = lower, ymax = upper), 
                color = "black", width = 0.2) +
  geom_point(data = steph.quant, aes(x = mean.date, y = mean, 
                                     fill = "Observed"),
             shape = 21, color = "black", alpha = 1, size = 2) +
  # theme(axis.text.x = element_text(face = "plain", angle = 90, vjust = 1, hjust = 1)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(legend.position = "none")

swc.d.vert.long.sub <- swc.d.vert.long.sub %>% mutate(depth.plot = "0-0.15 m")
p.v.swc.joint <- ggplot(swc.d.vert.long.sub, aes(x = date)) +
  geom_line(aes(y = value, group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  ## adding Observed normalised
  geom_line(aes(y = obs, color = "Observed"), size = 0.4) +
  scale_colour_manual(name = "", values = col1) + 
  ylab(expression('VWC ('*cm^3/cm^3*')')) +
  xlab("Time") + 
  geom_text(aes(x = swc.d.vert.long.sub$date[200], y = 0.57, label = rmse.v.swc.d), 
            size = 4, vjust = "inward", hjust = "inward") +
  theme(legend.position = "none") +
  facet_grid(depth.plot ~ .) + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        strip.text.y = element_text(size = 14, face = "plain")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y"))

p.et.m1.joint <- ggplot(qet.m.long.sub, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  scale_colour_manual(name = "", values = col1) +
  ylab("Evapotranspiration (mm/month)") +
  xlab("Time") + 
  geom_text(aes(x = qet.m.long.sub$date[200], y = 130, label = rmse.qet.m), 
            size = 4, vjust = "inward", hjust = "inward") +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.25, 0.95), 
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(legend.position = "none")

p.qrun.m1.joint <- ggplot(qrunoff.m.long %>% subset(date < as.Date("2018-07-01")),
                          aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  scale_colour_manual(name = "", values = col1) +
  ylab("Stream Discharge (mm/month)") +
  xlab("Time") + 
  geom_text(aes(x = qrunoff.m.long$date[200], y = 500, label = rmse.qrun.m), 
            size = 4, vjust = "inward", hjust = "inward") +
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.5, 0.8), 
        legend.key = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y"))

right_col <- cowplot::plot_grid(p.qrun.m1.joint, p.et.m1.joint, labels = c('B', 'C'), 
                                label_size = 14, nrow = 2, #label_x = 0, label_y = 0,
                                hjust = 1) #vjust = -0.5
joint.plot <- cowplot::plot_grid(p.h.swc.steph.mean.joint, right_col, labels = c('A', ''), 
                                 label_size = 14, ncol = 2, rel_widths = c(1, 1.2))
ggsave("hydro-calib.jpeg", plot = joint.plot, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 6, width = 9, units ='in')

p.et.m1.joint.shortlab <-  p.et.m1.joint + ylab(expression("ET (mm month"^-1*")"))
p.qrun.m1.joint.shortlab <-  p.qrun.m1.joint + ylab(expression("Discharge (mm month"^-1*")"))
p.h.swc.steph.mean.joint.shortlab <- p.h.swc.steph.mean.joint +  ylab(expression('VWC ('*cm^3~cm^-3*')'))

right_col2 <- cowplot::plot_grid(p.v.swc.joint, p.qrun.m1.joint.shortlab, 
                                 p.et.m1.joint.shortlab, 
                                                labels = c('B', 'C', 'D'), 
                                label_size = 14, nrow = 3, #label_x = 0, label_y = 0,
                                hjust = 1) #vjust = -0.5
joint.plot2 <- cowplot::plot_grid(p.h.swc.steph.mean.joint.shortlab, right_col2, labels = c('A', ''), 
                                 label_size = 14, ncol = 2, rel_widths = c(1, 1))
ggsave("hydro-calib2.jpeg", plot = joint.plot2, path = 
         file.path("figures", figures.folder), device = "jpeg", height = 6, width = 9, units ='in')
