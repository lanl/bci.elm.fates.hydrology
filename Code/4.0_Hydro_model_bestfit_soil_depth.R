##----------------
# Choosing best-fit ecosystem soil depth for hydrological water-balance
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
graphics.off()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, tidyverse, bci.hydromet, data.table, avasoilmoisture)

nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
#ncprint(paste("The file has", nc$nvars, "variables"))
#nc_close(nc)
# depths
nc$var[['H2OSOI']]$dim[[2]]$vals # 15 depths
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934
# so depths 4 [ 0.11 m], 6[0.36 m] and 8[1.03 m] roughly correspond to 15, 40 and 100 cm
# in units
nc$var[['H2OSOI']]$dim[[2]]$units # meters
## the time axis:  
nc$var[['H2OSOI']]$dim[[3]]$vals
## it's the time series of 
nc$var[['H2OSOI']]$dim[[3]]$units
# "days since 2003-01-01 00:00:00"

## data itself can be obtained thus
sm <- ncvar_get(nc, "H2OSOI")	# by default, reads ALL the data
print(data1)
dim(data1)
qvegt <- ncvar_get(nc, "QVEGT")
qvege <- ncvar_get(nc, "QVEGE")
qsoil <- ncvar_get(nc, "QSOIL")
qrunoff <- ncvar_get(nc, "QRUNOFF")

# H2OSOI	volumetric soil water (vegetated landunits only)	mm3/mm3

# QRUNOFF	total liquid runoff (does not include QSNWCPICE)	mm/s

# QSOIL	Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)	mm/s
# QVEGE	canopy evaporation	mm/s
# QVEGT	canopy transpiration	mm/s

bci.hydromet::forcings$AET.flag.day
bci.hydromet::forcings$flow_conrad
bci.hydromet::BCI.volumetric_soil_moisture_by_depth
##--------------------------
##--------------------------
nsim = 2 # for aveDTB 3:10 m
####
#### Runoff: Obs versus model #-------
####
# Daily #------
####
load("data-raw/extract/QRUNOFF.h1.extract.Rdata")

mod.qrunoff.d <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date))
mod.qrunoff.d[, -1] <- mod.qrunoff.d[, -1]*24*60*60 # converting mm/s to mm/day
obs.qrunoff.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day
  rename(obs = flow_conrad)
qrunoff.d <- obs.qrunoff.d %>% full_join(mod.qrunoff.d, by = "date")
head(qrunoff.d); summary(qrunoff.d)

rsq.table <- data.frame("soil_depth" = c(1: nsim) + 2, "qrunoff_daily" = rep(NA, nsim))
for (i in 1: nsim) {
  test <- cor.test(qrunoff.d$obs, qrunoff.d[, i + 2], 
                   method = "spearman")
  rsq.table$qrunoff_daily[i] <- round(as.numeric(test$estimate)^2, 3)
}
qrunoff.d.long <- gather(qrunoff.d, key = "nsim", "value", -date, -obs) %>%
  mutate(nsim = as.numeric(nsim), soil_depth = nsim + 2) 
ggplot(qrunoff.d.long %>% subset(date >= as.Date("2003-01-01")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth, group = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  # #facet_grid(soil_depth ~ .) +
  ylab("QRUNOFF [mm/day]") +
  xlab("Date") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain", angle = 90)) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Obs vs Model. Daily")
ggsave(file.path("figures/QRUNOFF_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
####
## from monthly files #-----
# load("data-raw/extract/QRUNOFF.h0.extract.Rdata") ##---
# mod.qrunoff.m <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "yrmo") %>%# in mm/s
#   as.data.frame()
# mod.qrunoff.m[, -1] <- mod.qrunoff.m[, -1]*24*60*60*30 # converting mm/s to mm/month
# obs.qrunoff.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day   
#   rename(obs = flow_conrad)
# obs.qrunoff.m <- obs.qrunoff.d %>% 
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
#   group_by(yrmo) %>% summarise(obs = sum(obs, na.rm = T)) %>% as.data.frame() # in mm/month
# qrunoff.m <- obs.qrunoff.m %>%
#   full_join(mod.qrunoff.m, by = "yrmo")
# head(qrunoff.m); summary(qrunoff.m) ###----
## from daily files #-------
qrunoff.m.long <- qrunoff.d.long %>%  
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, soil_depth) %>% 
  summarise(obs = sum(obs, na.rm = T), value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
qrunoff.m <- spread(qrunoff.m.long %>% select(-date), key = soil_depth, value = value)
summary(qrunoff.m)
rsq.table$qrunoff_monthly <- NA
for (i in 1: nsim) {
  test <- cor.test(qrunoff.m$obs, qrunoff.m[, i + 2], 
                   method = "spearman")
  rsq.table$qrunoff_monthly[i] <- round(as.numeric(test$estimate)^2, 3)
}
rsq.table
ggplot(qrunoff.m.long %>% subset(date >= as.Date("2003-01-01")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth, group = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  #facet_grid(soil_depth ~ .) +
  ylab("QRUNOFF [mm/month]") +
  xlab("Date") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain", angle = 90)) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Obs vs Model. Monthly")
ggsave(file.path("figures/QRUNOFF_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')

#--------------------
####
#### Evapotranspiration: Obs versus model #-------
####
# Daily #------
####
load("data-raw/extract/QVEGE.h1.extract.Rdata")
mod.qvege.d <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "date") # in mm/s
load("data-raw/extract/QVEGT.h1.extract.Rdata")
mod.qvegt.d <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "date") # in mm/s
load("data-raw/extract/QSOIL.h1.extract.Rdata")
mod.qsoil.d <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "date") # in mm/s
## combining evaporation from soil and plants and transpiration from plants
flow.columns <- (mod.qvege.d[, -1] + mod.qvegt.d[, -1] + mod.qsoil.d[, -1])*24*60*60 # from mm/s to mm/day
mod.qet.d <-  cbind.data.frame(mod.qvege.d[, 1], flow.columns) %>% 
  mutate(date = as.Date(date))
## observed ET from flux tower
obs.qet.d <- bci.hydromet::forcings %>% select(date, AET) %>% # in mm/day
  rename(obs = AET)
qet.d <- obs.qet.d %>% full_join(mod.qet.d, by = "date") 
head(qet.d); summary(qet.d)
qet.d.cor <- qet.d %>% subset(!is.na(obs))
rsq.table$qet_daily <- NA
for (i in 1: nsim) {
  test <- cor.test(qet.d.cor$obs, qet.d.cor[, i + 2], 
                   method = "spearman")
  rsq.table$qet_daily[i] <- round(as.numeric(test$estimate)^2, 3)
}
rsq.table
qet.d.long <- gather(qet.d, key = "nsim", "value", -date, -obs) %>%
  mutate(nsim = as.numeric(nsim), soil_depth = nsim + 2) 
ggplot(qet.d.long %>% subset(date >= as.Date("2012-06-01") & date < as.Date("2017-06-01")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth, group = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  #facet_grid(soil_depth ~ .) +
  ylab("ET [mm/day]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain")) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Obs vs Model. Daily")
ggsave(file.path("figures/ET_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
qet.m.long <- qet.d.long %>%  
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, soil_depth) %>% 
  summarise(obs = sum(obs, na.rm = T), value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
qet.m <- spread(qet.m.long %>% select(-date), key = soil_depth, value = value)
summary(qet.m)
rsq.table$qet_monthly <- NA
for (i in 1: nsim) {
  test <- cor.test(qet.m$obs, qrunoff.m[, i + 2], 
                   method = "spearman")
  rsq.table$qet_monthly[i] <- round(as.numeric(test$estimate)^2, 3)
}
rsq.table
ggplot(qet.m.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2016-12-31")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth, group = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  #facet_grid(soil_depth ~ .) +
  ylab("ET [mm/month]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain")) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Obs vs Model. Monthly")
ggsave(file.path("figures/ET_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')

#--------------
####
#### Soil Moisture: Obs versus model #-------
####
# Daily #------
avasoilmoisture::tdr.daily
bci.hydromet::sm.historic
####
nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
# depths
nc$var[['H2OSOI']]$dim[[2]]$vals # 15 depths
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934
# so depths 4 [ 0.11 m], 6[0.36 m] and 8[1.03 m] roughly correspond to 10, 40 and 100 cm
# in units m3/m3

load("data-raw/extract/H2OSOI.h1.extract.Rdata")
mod.swc.d.10 <- setDT(as.data.frame(t(res.arr[[4]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 10) %>% gather(key = "nsim", "value", -date, -depth)
mod.swc.d.40 <- setDT(as.data.frame(t(res.arr[[6]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 40) %>% gather(key = "nsim", "value", -date, -depth)
mod.swc.d.100 <- setDT(as.data.frame(t(res.arr[[8]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 100) %>% gather(key = "nsim", "value", -date, -depth)

mod.swc.d <- mod.swc.d.10 %>% rbind(mod.swc.d.40) %>% rbind(mod.swc.d.100) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()
obs.swc.d <- avasoilmoisture::horizontal %>% 
  rename(obs = swc) %>% mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()
swc.d.long <- obs.swc.d %>% select(-date, -depth) %>%
  full_join(mod.swc.d, by = "date.depth") %>% select(-date.depth) %>% 
  as.data.frame() %>% subset(!is.na(nsim))
head(swc.d.long); summary(swc.d.long)
swc.d.10 <- swc.d.long %>% subset(depth == "10" & !is.na(nsim)) %>% select(-depth) %>%
  spread(key = nsim, value = value) 
swc.d.40 <- swc.d.long %>% subset(depth == "40" & !is.na(nsim)) %>% select(-depth) %>%
  spread(key = nsim, value = value) 
swc.d.100 <- swc.d.long %>% subset(depth == "100" & !is.na(nsim)) %>% select(-depth) %>%
  spread(key = nsim, value = value) 

rsq.table$swc100_daily <- rsq.table$swc40_daily <- rsq.table$swc10_daily <- NA
for (i in 1: nsim) {
  test <- cor.test(swc.d.10$obs, swc.d.10[, i + 2], 
                   method = "spearman")
  rsq.table$swc10_daily[i] <- round(as.numeric(test$estimate)^2, 3)
  test <- cor.test(swc.d.40$obs, swc.d.40[, i + 2], 
                   method = "spearman")
  rsq.table$swc40_daily[i] <- round(as.numeric(test$estimate)^2, 3)
  
  test <- cor.test(swc.d.100$obs, swc.d.100[, i + 2], 
                   method = "spearman")
  rsq.table$swc100_daily[i] <- round(as.numeric(test$estimate)^2, 3)
  
}
rsq.table
swc.d.long <- swc.d.long %>%
  mutate(nsim = as.numeric(nsim), soil_depth = nsim + 2) %>% subset(!is.na(soil_depth))
ggplot(swc.d.long %>% subset(soil_depth == 3) %>% subset(date >= as.Date("2008-01-01") & date < as.Date("2017-01-01")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  facet_grid(depth ~ .) +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Date") + theme_bw() +
  theme(axis.text.x = element_text(size = 16, face = "plain", angle = 90)) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 month", labels = function(x) format(x, "%b%y")) +
  ggtitle("Soil Water Content: Obs vs Model. Daily")
ggsave(file.path("figures/swc_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
####
## from monthly files #-----
# load("data-raw/extract/QRUNOFF.h0.extract.Rdata") ##---
# mod.swc.m <- setDT(as.data.frame(t(res.arr[[1]])), keep.rownames = "yrmo") %>%# in mm/s
#   as.data.frame()
# mod.swc.m[, -1] <- mod.swc.m[, -1]*24*60*60*30 # converting mm/s to mm/month
# obs.swc.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day   
#   rename(obs = flow_conrad)
# obs.swc.m <- obs.swc.d %>% 
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
#   group_by(yrmo) %>% summarise(obs = sum(obs, na.rm = T)) %>% as.data.frame() # in mm/month
# swc.m <- obs.swc.m %>%
#   full_join(mod.swc.m, by = "yrmo")
# head(swc.m); summary(swc.m) ###----
## from daily files #-------
swc.m.long <- swc.d.long %>%  subset(!is.na(obs)) %>%
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, soil_depth, depth) %>% 
  summarise(obs = sum(obs, na.rm = T), value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
swc.m.10 <- swc.m.long %>% subset(depth == "10" & !is.na(nsim)) %>% select(-depth, -date) %>%
  spread(key = soil_depth, value = value) 
swc.m.40 <- swc.m.long %>% subset(depth == "40" & !is.na(nsim)) %>% select(-depth, -date) %>%
  spread(key = soil_depth, value = value) 
swc.m.100 <- swc.m.long %>% subset(depth == "100" & !is.na(nsim)) %>% select(-depth, -date) %>%
  spread(key = soil_depth, value = value) 

rsq.table$swc100_monthly <- rsq.table$swc40_monthly <- rsq.table$swc10_monthly <- NA
for (i in 1: nsim) {
  test <- cor.test(swc.m.10$obs, swc.m.10[, i + 2], 
                   method = "spearman")
  rsq.table$swc10_monthly[i] <- round(as.numeric(test$estimate)^2, 3)
  test <- cor.test(swc.m.40$obs, swc.m.40[, i + 2], 
                   method = "spearman")
  rsq.table$swc40_monthly[i] <- round(as.numeric(test$estimate)^2, 3)
  
  test <- cor.test(swc.m.100$obs, swc.m.100[, i + 2], 
                   method = "spearman")
  rsq.table$swc100_monthly[i] <- round(as.numeric(test$estimate)^2, 3)
  
}
rsq.table
ggplot(swc.m.long %>% subset(date >= as.Date("2012-06-01") & date < as.Date("2017-06-01")),
       aes(x = date, y = value)) +
  geom_line(aes(color = soil_depth, group = soil_depth), show.legend = F, size = 1) +
  scale_colour_gradient(guide = guide_colourbar(reverse = TRUE)) +
  geom_line(aes( x= date, y = obs), colour = "red", size = 1) +
  #facet_grid(soil_depth ~ .) +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Date") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain", angle = 90)) +
  # theme(axis.text = tex, strip.text.y = tex, title = tex) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Soil Water Content: Obs vs Model. Monthly")
ggsave(file.path("figures/swc_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')


## Plotting rsq table
rsq.long <- rsq.table %>% gather(key = "variable", value = "rsq", -soil_depth)
ggplot(rsq.long, aes(x = soil_depth, y = rsq)) +
  geom_bar(aes(fill = variable), stat = "identity", show.legend = F) +
  # scale_fill_gradient(guide = guide_fill(reverse = TRUE)) +
  facet_grid(variable ~ ., scales = "free_y") +
  ylab("Rsq") +
  xlab("Soil Depth [m]") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain")) +
  ggtitle("Rsq between Observation and Model By Soil Depth")
ggsave(file.path("figures/Rsq between Observation and Model By Soil Depth.jpeg"), height = 8.94, width = 12.9, units='in')
