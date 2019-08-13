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
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, tidyverse, bci.hydromet, data.table, avasoilmoisture, hydroGOF)

# H2OSOI	volumetric soil water (vegetated landunits only)	mm3/mm3

# QRUNOFF	total liquid runoff (does not include QSNWCPICE)	mm/s

## Total evapotransipration:
# QSOIL	Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew)	mm/s
# QVEGE	canopy evaporation	mm/s
# QVEGT	canopy transpiration	mm/s
##--------------------------
##--------------------------
current.folder <- "4000"
params.from.files <- read.csv("data-raw/params.from.files.csv", header = TRUE)
## but these are only hundred
params <- params.from.files
# params <- read.csv("data-raw/params.csv", header = TRUE)

for (i in 1: (40 - 1)) {
  params <- bind_rows(params, params.from.files)
}
params$par.sam <- 1:nrow(params)
params$soil_depth <- rep(c(3:10), each = 5*100)
#soil depth array # in R # 
params$fmax <-rep(c(0.2, 0.3, 0.4, 0.5, 0.6), times = 8*100)

head(params)
filterFile <- read.table("data-raw/Filter.txt")
nrow(filterFile)
filterFile$V1[c(1:100, 740)] <- FALSE 

col1 <- c("observed" = "#f04546", "simulated" = "#3591d1", "best-fit" = "green")

# 
# ### test ----
# params <- params.from.files[c(3, 3, 3, 3, 3), ]
# params$nsam <- 1
# nsam = 4
# filterFile <- read.table("data-raw/Filter.txt")
# filterFile <- data.frame(V1 = filterFile$V1[1])
### ---

par.on <- params[filterFile$V1, ]
nsam <- length(par.on$par.sam) #length(which(filterFile$V1 == TRUE))

###
#### Runoff: Obs versus model #-------
####
# Daily #------
####
load(file.path("data-raw/extract", current.folder, "extract/QRUNOFF.h1.extract.Rdata"))
mod.qrunoff.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date))
mod.qrunoff.d[, -1] <- mod.qrunoff.d[, -1]*24*60*60 # converting mm/s to mm/day
obs.qrunoff.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day
  rename(obs = flow_conrad)
qrunoff.d <- obs.qrunoff.d %>% full_join(mod.qrunoff.d, by = "date")
#head(qrunoff.d); #summary(qrunoff.d)

rmse.table <- data.frame(par.sam = par.on$par.sam, "qrunoff_daily" = rep(NA, nsam))
nse.table <- data.frame(par.sam = par.on$par.sam, "qrunoff_daily" = rep(NA, nsam))
for (i in 1: nsam) {
  obs <- as.numeric(qrunoff.d$obs);
  sim <- as.numeric(qrunoff.d[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$qrunoff_daily[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$qrunoff_daily[i] <- NSE(obs, sim, na.rm = TRUE)
}
summary(rmse.table); summary(nse.table)
qrunoff.d.long <- gather(qrunoff.d, key = "par.sam", "value", -date, -obs) 
nse.max.run.d <- max(nse.table$qrunoff_daily)
max.par.n.run.row.d <- which(nse.table$qrunoff_daily == nse.max.run.d)
max.par.n.run.d <- nse.table$par.sam[max.par.n.run.row.d]
qrunoff.d.long.sub <- qrunoff.d.long %>% subset(date >= as.Date("2008-01-01"))

ggplot(qrunoff.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.1) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam %in% max.par.n.run.d), 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_text(data = qrunoff.d.long.sub %>% subset(par.sam %in% max.par.n.run.d),
            aes(x = as.Date("2010-01-01"), y = 100, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.run.d, 2))))) +
  scale_colour_manual(name = "", values = col1) +
  geom_line(aes(y = obs, colour = "observed"), size = 0.5) +
  ylab("QRUNOFF [mm/day]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Obs vs Model. Daily")
ggsave(file.path("figures", current.folder, "QRUNOFF_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
####
## from monthly files #-----
# load(file.path("data-raw/extract", current.folder, "extract/QRUNOFF.h0.extract.Rdata")) ##---
# mod.qrunoff.m <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%# in mm/s
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
  group_by(yrmo, par.sam) %>% 
  summarise(obs = sum(obs, na.rm = T), value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
qrunoff.m <- spread(qrunoff.m.long %>% select(-date), key = par.sam, value = value)
summary(qrunoff.m[, 1:3])
rmse.table$qrunoff_monthly <- NA; nse.table$qrunoff_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(qrunoff.d$obs); sim <- as.numeric(qrunoff.d[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$qrunoff_monthly[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$qrunoff_monthly[i] <- NSE(obs, sim, na.rm = TRUE)
}

summary(rmse.table);  summary(nse.table)
nse.max.run.m <- max(nse.table$qrunoff_monthly)
max.par.n.run.row.m <- which(nse.table$qrunoff_monthly == nse.max.run.m)
max.par.n.run.m <- nse.table$par.sam[max.par.n.run.row.m]
qrunoff.m.long.sub <- qrunoff.m.long %>% subset(date >= as.Date("2008-01-01"))
nse.text <- as.character(as.expression(italic(NSE[MAX])~" = "~nse.max.run))
ggplot(qrunoff.m.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam %in% max.par.n.run.m),
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 0.7) +
  geom_text(data = qrunoff.m.long.sub %>% subset(par.sam %in% max.par.n.run.m),
            aes(x = as.Date("2010-01-01"), y = 600, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.run.m, 2))))) +
  scale_colour_manual(name = "", values = col1) +
  ylab("QRUNOFF [mm/month]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Obs vs Model. Monthly")
ggsave(file.path("figures", current.folder, "QRUNOFF_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')

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
## combining evaporation from soil and plants and transpiration from plants
flow.columns <- (mod.qvege.d[, -1] + mod.qvegt.d[, -1] + mod.qsoil.d[, -1])*24*60*60 # from mm/s to mm/day
mod.qet.d <-  cbind.data.frame(mod.qvege.d[, 1], flow.columns) %>% 
  mutate(date = as.Date(date))
## observed ET from flux tower
obs.qet.d <- bci.hydromet::forcings %>% select(date, AET) %>% # in mm/day
  rename(obs = AET)
qet.d <- obs.qet.d %>% full_join(mod.qet.d, by = "date") 
head(qet.d); summary(qet.d)
rmse.table$qet_daily <- NA; nse.table$qet_daily <- NA
for (i in 1: nsam) {
  obs <- as.numeric(qet.d$obs); sim <- as.numeric(qet.d[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$qet_daily[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$qet_daily[i] <- NSE(obs, sim, na.rm = TRUE)
}
summary(rmse.table); summary(nse.table)

qet.d.long <- gather(qet.d, key = "par.sam", "value", -date, -obs)
nse.max.et.d <- max(nse.table$qet_daily)
max.par.n.et.row.d <- which(nse.table$qet_daily == nse.max.et.d)
max.par.n.et.d <- nse.table$par.sam[max.par.n.et.row.d]
qet.d.long.sub <- qet.d.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31")) %>% droplevels()
qet.d.long.sub.best.fit <- qet.d.long.sub %>% subset(par.sam %in% max.par.n.et.d)

p.et1 <- ggplot(qet.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.1) +
  geom_line(data = qet.d.long.sub.best.fit, 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 0.5) +
  geom_text(data = qet.d.long.sub %>% subset(par.sam %in% max.par.n.et.d),
            aes(x = as.Date("2013-01-01"), y = 5, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.et.d, 2))))) +
  scale_colour_manual(name = "", values = col1) +
  ylab("ET [mm/day]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Obs vs Model. Daily")
p.et1 +   geom_line(data = qet.d.long.sub.best.fit, 
                    aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 0.5)
ggsave(file.path("figures", current.folder, "ET_Obs_vs_model_daily_all_years.jpeg"), height = 8.94, width = 12.9, units='in')

qet.d.long.sub.2 <- qet.d.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2015-02-15"))
qet.d.long.sub.best.fit.2 <- qet.d.long.sub.2 %>% subset(par.sam %in% max.par.n.et.d) %>% droplevels()
p.et1 %+% qet.d.long.sub.2 +
  geom_line(data = qet.d.long.sub.best.fit.2, 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 0.5) +
  scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y"))
ggsave(file.path("figures", current.folder, "ET_Obs_vs_model_daily_2012-2014.jpeg"), height = 8.94, width = 12.9, units='in')

qet.d.long.sub.3 <- qet.d.long %>% subset(date >= as.Date("2015-05-01") & date < as.Date("2017-08-31")) 
qet.d.long.sub.best.fit.3 <- qet.d.long.sub.3 %>% subset(par.sam %in% max.par.n.et.d)
p.et1 %+% qet.d.long.sub.3 +
  geom_line(data = qet.d.long.sub.best.fit.3, 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 0.5) +
  scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y"))
ggsave(file.path("figures", current.folder, "ET_Obs_vs_model_daily_2015-2017.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
qet.m.long <- qet.d.long %>%  
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, par.sam) %>% 
  summarise(obs = sum(obs, na.rm = T), value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
qet.m <- spread(qet.m.long %>% select(-date), key = par.sam, value = value)
summary(qet.m[, 1:5])
nse.table$qet_monthly <- NA; nse.table$qet_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(qet.d$obs); sim <- as.numeric(qet.d[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$qet_monthly[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$qet_monthly[i] <- NSE(obs, sim, na.rm = TRUE)
}
summary(rmse.table);  summary(nse.table)
nse.max.et.m <- max(nse.table$qet_monthly, na.rm = T)
max.par.n.et.row.m <- which(nse.table$qet_monthly == nse.max.et.m)
max.par.n.et.m <- nse.table$par.sam[max.par.n.et.row.m]

qet.m.long.sub <- qet.m.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31"))
qet.m.long.sub.best.fit <- qet.m.long.sub %>% subset(par.sam %in% max.par.n.et.m)
qet.m.long.sub.1 <- qet.m.long.sub %>% subset(par.sam %in% max.par.n.et.m[1])
# rmse.label <- as.character(as.expression(italic(r)^2~"="~rmse.max))
p.et.m1 <- ggplot(qet.m.long.sub, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
  geom_line(data = qet.m.long.sub.best.fit, aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes( x= date, y = obs, color = "observed"), size = 1) +
  geom_text(data = qet.m.long.sub %>% subset(par.sam %in% max.par.n.et.m),
            aes(x = as.Date("2013-01-01"), y = 60, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.et.m, 2))))) +
  scale_colour_manual(name = "", values = col1) +
  ylab("ET [mm/month]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Obs vs Model. Monthly")
p.et.m1
ggsave(file.path("figures", current.folder, "ET_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')

# 
# m1 <- lm(value ~ obs, data = qet.m.long.sub.best.fit %>% subset(par.sam  == 627 ))
# summary(m1)
# ggplot(qet.m.long.sub,
#        aes(x = obs, y = value)) +
#   geom_point(aes(color = "simulated"), show.legend = F) +
#   geom_point(data = qet.m.long.sub.best.fit, 
#             aes(color = "best-fit"), show.legend = F) +
#   geom_smooth(data = qet.m.long.sub.best.fit, method = "lm") +
#   ylab("Simulated Best-Fit [mm/month]") +
#   xlab("Observed [mm/month]") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain")) +
#   # scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Evapotranspiration: Obs vs Model. Best-FitMonthly_scatter")


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

bci.tower <- read.csv("data-raw/BCI_v3.1.csv", header = TRUE)
bci.tower <- as.data.frame(bci.tower[-1, ])
bci.tower$datetime <- strptime(bci.tower$date, format = "%m/%d/%Y %H:%M")
bci.tower$datetime <- as.POSIXct(bci.tower$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "")
str(bci.tower$datetime)

obs.gpp.d <- bci.tower %>% select(datetime, gpp) %>% # in mumol/m2/2: units must be mumol/m2/s
  mutate(date = as.Date(format(datetime, "%Y-%m-%d")),  
         gpp.mumol = as.numeric(as.character(gpp))) %>%
  group_by(date) %>% summarise(obs = mean(gpp.mumol, na.rm = T)) %>%
  mutate(obs = obs*12*1e-06*24*60*60) # gC/m2/d
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

rmse.table$gpp_daily <- NA; nse.table$gpp_daily <- NA
for (i in 1: nsam) {
  obs <- as.numeric(gpp.d$obs); sim <- as.numeric(gpp.d[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$gpp_daily[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$gpp_daily[i] <- NSE(obs, sim, na.rm = TRUE)
}
summary(rmse.table); summary(nse.table)
gpp.d.long <- gather(gpp.d, key = "par.sam", "value", -date, -obs) 
nse.max.gpp.d <- max(nse.table$gpp_daily)
max.par.n.gpp.row.d <- which(nse.table$gpp_daily == nse.max.gpp.d)
max.par.n.gpp.d <- nse.table$par.sam[max.par.n.gpp.row.d]

gpp.d.long.sub <- gpp.d.long %>% subset(date >= as.Date("2008-01-01"))

ggplot(gpp.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.1) +
  geom_line(data = gpp.d.long.sub %>% subset(par.sam %in% max.par.n.gpp.d), 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  scale_colour_manual(name = "", values = col1) +
  geom_line(aes(y = obs, colour = "observed"), size = 0.5) +
  geom_text(data = gpp.d.long.sub %>% subset(par.sam %in% max.par.n.gpp.d),
            aes(x = as.Date("2013-01-01"), y = 10, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.gpp.d, 2))))) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Obs vs Model. Daily")
ggsave(file.path("figures", current.folder, "GPP_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')

####
# Monthly #------
####
## from monthly files #-----
# load(file.path("data-raw/extract", current.folder, "extract/GPP.h0.extract.Rdata")) ##---
# mod.gpp.m <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%# in mm/s
#   as.data.frame()
# mod.gpp.m[, -1] <- mod.gpp.m[, -1]*24*60*60 # converting rom gC/m^2/s to gC/m^2/d
# obs.gpp.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day   
#   rename(obs = flow_conrad)
# obs.gpp.m <- obs.gpp.d %>% 
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
#   group_by(yrmo) %>% summarise(obs = sum(obs, na.rm = T)) %>% as.data.frame() # in mm/month
# gpp.m <- obs.gpp.m %>%
#   full_join(mod.gpp.m, by = "yrmo")
# head(gpp.m); summary(gpp.m) ###----
## from daily files #-------
gpp.m.long <- gpp.d.long %>%  
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, par.sam) %>% 
  summarise(obs = mean(obs, na.rm = T), value = mean(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
gpp.m <- spread(gpp.m.long %>% select(-date), key = par.sam, value = value)
summary(gpp.m[, 1:3])
nse.table$gpp_monthly <- NA; nse.table$gpp_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(gpp.m$obs); sim <- as.numeric(gpp.m[, i + 2])
  variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
  rmse.table$gpp_monthly[i] <- round(as.numeric(variance)^0.5, 3)
  nse.table$gpp_monthly[i] <- NSE(obs, sim, na.rm = TRUE)
}

summary(rmse.table);  summary(nse.table)
nse.max.gpp.m <- max(nse.table$gpp_monthly)
max.par.n.gpp.row.m <- which(nse.table$gpp_monthly == nse.max.gpp.m)
max.par.n.gpp.m <- nse.table$par.sam[max.par.n.gpp.row.m]
gpp.m.long.sub <- gpp.m.long %>% subset(date >= as.Date("2008-01-01"))
ggplot(gpp.m.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
  geom_line(data = gpp.m.long.sub %>% subset(par.sam %in% max.par.n.gpp.m), 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes( x= date, y = obs, color = "observed"), size = 0.7) +
  geom_text(data = gpp.m.long.sub %>% subset(par.sam %in% max.par.n.gpp.m),
            aes(x = as.Date("2009-01-01"), y = 11, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.gpp.m, 2))))) +
  scale_colour_manual(name = "", values = col1) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.9, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Obs vs Model. Monthly")
ggsave(file.path("figures", current.folder, "GPP_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')

# 
# ##
# ## Monthly from daily files #-------
# ## from monthly files 
# load(file.path("data-raw/extract", current.folder, "extract/GPP.h0.extract.Rdata")) ##---
# mod.gpp.m <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%
#   as.data.frame() 
# head(mod.gpp.m[,1:3])
# mod.gpp.m[, -1] <- mod.gpp.m[, -1]*24*60*60 # converting from gC/m^2/s to gC/m^2/d
# 
# obs.gpp.d <- obs.gpp.d %>%  
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) 
# obs.gpp.m <- obs.gpp.d %>%
#   group_by(yrmo) %>% summarise(obs = mean(obs, na.rm = T))
# gpp.m <- obs.gpp.m %>% full_join(mod.gpp.m, by = "yrmo") %>% as.data.frame()
# head(gpp.m[,1:3])
# 
# gpp.m.long <-  gather(gpp.m, key = "par.sam", "value", -yrmo, -obs) %>%
#   as.data.frame()
# 
# str(gpp.m.long)
# gpp.m.long$date <- as.Date(strptime(paste0(gpp.m.long$yrmo, "-01"), "%Y-%m-%d"))
# gpp.m <- spread(gpp.m.long, key = par.sam, value = value)
# summary(gpp.m[, 1:3])
# nse.table$gpp_monthly <- NA; nse.table$gpp_monthly <- NA
# for (i in 1: nsam) {
#   obs <- as.numeric(gpp.m$obs); sim <- as.numeric(gpp.m[, i + 2])
#   variance <- var(obs, sim, na.rm = TRUE, use = "complete.obs")
#   rmse.table$gpp_monthly[i] <- round(as.numeric(variance)^0.5, 3)
#   nse.table$gpp_monthly[i] <- NSE(obs, sim, na.rm = TRUE)
# }
# 
# summary(rmse.table);  summary(nse.table)
# nse.max.gpp.m <- max(nse.table$gpp_monthly)
# max.par.n.gpp.row.m <- which(nse.table$gpp_monthly == nse.max.gpp.m)
# max.par.n.gpp.m <- nse.table$par.sam[max.par.n.gpp.row.m]
# gpp.m.long.sub <- gpp.m.long %>% subset(date >= as.Date("2008-01-01"))
# 
# # gpp.m.long.sub <- gpp.m.long %>% subset(date >= as.Date("2008-01-01"))
# ## add best-fits for ET
# ggplot(gpp.m.long,
#               aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
#   geom_line(data = gpp.m.long %>% subset(par.sam %in% max.par.n.gpp.m), 
#             aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
#   geom_line(aes( x= date, y = obs, color = "observed"), size = 0.7) +
#   scale_colour_manual(name = "", values = col1) +
#   ylab("GPP [gC/m^2/d]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = c(0.8, 0.5), legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
#   ggtitle("GPP: Modelled Monthly")
# ggsave(file.path("figures", current.folder, "GPP_model_monthly_bestfit.jpeg"), height = 8.94, width = 12.9, units='in')
###--

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
nc$var[['H2OSOI']]$dim[[2]]$vals # 15 depths
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
mod.swc.d.1 <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 1) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.3 <- setDT(as.data.frame(t(var.res.arr[[2]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 3) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.6 <- setDT(as.data.frame(t(var.res.arr[[3]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 6) %>% gather(key = "par.sam", "value", -date, -depth)

mod.swc.d.5 <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>% 
  bind_rows(setDT(as.data.frame(t(var.res.arr[[2]])), keep.rownames = "date")) %>% 
  bind_rows(setDT(as.data.frame(t(var.res.arr[[3]])), keep.rownames = "date")) %>% 
  bind_rows(setDT(as.data.frame(t(var.res.arr[[4]])), keep.rownames = "date")) %>% 
  group_by(date) %>% summarise_all(mean, na.rm = TRUE) %>%
  mutate(date = as.Date(date), depth = 5) %>% gather(key = "par.sam", "value", -date, -depth) %>% as.data.frame()

mod.swc.d.10 <- setDT(as.data.frame(t(var.res.arr[[4]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 10) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.40 <- setDT(as.data.frame(t(var.res.arr[[6]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 40) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.100 <- setDT(as.data.frame(t(var.res.arr[[8]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 100) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.300 <- setDT(as.data.frame(t(var.res.arr[[10]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 300) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.500 <- setDT(as.data.frame(t(var.res.arr[[11]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 500) %>% gather(key = "par.sam", "value", -date, -depth)
mod.swc.d.800 <- setDT(as.data.frame(t(var.res.arr[[12]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = 800) %>% gather(key = "par.sam", "value", -date, -depth)

swc.vertical <- avasoilmoisture::vertical %>% rename(obs = swc) %>% 
  group_by(date, depth) %>% summarise(obs = mean(obs, na.rm = TRUE)) %>% 
  mutate(date.depth = paste0(date, ".", depth))

str(swc.vertical)
mod.swc.d.shallow <- mod.swc.d.1 %>% bind_rows(mod.swc.d.3) %>% bind_rows(mod.swc.d.6) %>% 
  bind_rows(mod.swc.d.10) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()
str(mod.swc.d.shallow)
swc.d.shallow.sub <- mod.swc.d.shallow %>% subset(date >= as.Date("2016-01-01") & date < as.Date("2016-12-31"))
swc.vertical.sub <- swc.vertical %>% subset(date >= as.Date("2016-01-01") & date < as.Date("2016-12-31"))
# 
# nsam.on = 4 ## see txr.x in fsurdata_parameters.R for corresponding texture used
# ggplot(swc.d.shallow.sub %>% subset(par.sam  == nsam.on), aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = as.factor(depth)), show.legend = F, size = 1) +
#   geom_line(data = swc.vertical.sub, aes(x= date, y = obs, color = as.factor(depth)), linetype = "longdash", size = 1) +
#   # scale_colour_manual(name = "", values = col1) +
#   # facet_grid(depth ~ .) +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = "top", legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Obs vs Model. Daily")
# ggsave(file.path(paste0(paste0("figures/", current.folder, "/txr/swc_Obs_vs_model_daily_0-15cm_txr.", nsam.on,".jpeg")), height = 8.94, width = 12.9, units='in')
# 
mod.swc.d <- mod.swc.d.5 %>% bind_rows(mod.swc.d.10) %>% bind_rows(mod.swc.d.40) %>% 
  bind_rows(mod.swc.d.100) %>% bind_rows(mod.swc.d.300) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()
obs.swc.d <- avasoilmoisture::horizontal %>% 
  rename(obs = swc) %>% group_by(date, depth) %>% summarise(obs = mean(obs, na.rm = TRUE)) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() #%>%
  #bind_rows(swc.single) %>% mutate(depth = as.numeric(depth))

swc.d.long <- obs.swc.d %>% select(-date) %>%
  full_join(mod.swc.d %>% select(-depth), by = "date.depth") %>% 
  # separate(date.depth, c("date", "depth"), sep = ".", extra = "merge") %>% 
  as.data.frame() #%>% subset(!is.na(nsam))
head(swc.d.long); summary(swc.d.long)

### swc normalised:
norm.index <- function(x, ...){x/max(x, ...)}
mod.sat.d <- mod.swc.d %>% group_by(depth, par.sam) %>%
  mutate_at(c("value"), norm.index, na.rm = TRUE) %>% as.data.frame() 
## mod soil moisture also picks at Nov 23/24 2016
### tdr normalised or saturation.index(uses swc)
## because maximas are overestimated due to macropore structure, choose value after two days
norm.index2 <- function(x, na.rm = TRUE) (x/max(x, na.rm = TRUE))
obs.tdr.daily <- avasoilmoisture::horizontal %>% 
  rename(obs = tdr) %>% group_by(date, depth) %>% # this is directly normalised tdr
  summarise(obs = mean(obs, na.rm = TRUE)) %>% as.data.frame() 
obs.tdr.daily.spread <- obs.tdr.daily %>% spread(key = depth, value = obs)
obs.sat.daily.spread <- obs.tdr.daily.spread
col.max <- as.numeric(apply(obs.tdr.daily.spread[,-1],2, max, na.rm = T)) 
obs.sat.daily.spread[,2] <-  obs.tdr.daily.spread[, 2]/col.max[1]
obs.sat.daily.spread[,3] <-  obs.tdr.daily.spread[, 3]/col.max[2]
obs.sat.daily.spread[,4] <-  obs.tdr.daily.spread[, 4]/col.max[3]
  # mutate_at(c("obs"), norm.index2) %>%
obs.sat.d <- obs.sat.daily.spread %>%
  gather(key = depth, value = obs, -date) %>%
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() %>%
  mutate(depth = as.numeric(depth))
#bind_rows(swc.single) %>% mutate(depth = as.numeric(depth))

sat.d.long <- obs.sat.d %>% select(-date) %>%
  full_join(mod.sat.d %>% select(-depth), by = "date.depth") %>% 
  # separate(date.depth, c("date", "depth"), sep = ".", extra = "merge") %>% 
  as.data.frame() #%>% subset(!is.na(nsam))
head(sat.d.long); summary(sat.d.long)

sat.d.long <- sat.d.long %>% subset(!is.na(depth))
# sat.d.long <- sat.d.long %>% subset(depth != 5)
sat.d.long.sub <- sat.d.long %>% subset(date >= as.Date("2012-01-01") & date < as.Date("2017-05-31"))

sat.d.10 <- sat.d.long %>% subset(depth == "10" & !is.na(par.sam)) %>% select(-depth) %>%
  spread(key = par.sam, value = value) 
sat.d.40 <- sat.d.long %>% subset(depth == "40" & !is.na(par.sam)) %>% select(-depth) %>%
  spread(key = par.sam, value = value) 
sat.d.100 <- sat.d.long %>% subset(depth == "100" & !is.na(par.sam)) %>% select(-depth) %>%
  spread(key = par.sam, value = value) 

rmse.table$sat100_daily <- rmse.table$sat40_daily <- rmse.table$sat10_daily <- NA
nse.table$sat100_daily <- nse.table$sat40_daily <- nse.table$sat10_daily <- NA
for (i in 1: nsam) {
  variance1 <- var(sat.d.10$obs, sat.d.10[, i + 3], use = "complete.obs")
  rmse.table$sat10_daily[i] <- round(as.numeric(variance1)^0.5, 3)
  variance2 <- var(sat.d.40$obs, sat.d.40[, i + 3], use = "complete.obs")
  rmse.table$sat40_daily[i] <- round(as.numeric(variance2)^0.5, 3)
  variance3 <- var(sat.d.100$obs, sat.d.100[, i + 3], use = "complete.obs")
  rmse.table$sat100_daily[i] <- round(as.numeric(variance3)^0.5, 3)
  nse.table$sat10_daily[i] <- NSE(sat.d.10$obs, sat.d.10[, i + 3], na.rm = TRUE)
  nse.table$sat40_daily[i] <- NSE(sat.d.40$obs, sat.d.40[, i + 3], na.rm = TRUE)
  nse.table$sat100_daily[i] <- NSE(sat.d.100$obs, sat.d.100[, i + 3], na.rm = TRUE)
}
summary(rmse.table); summary(nse.table)
swc.d.long <- swc.d.long %>% subset(!is.na(depth))
# swc.d.long <- swc.d.long %>% subset(depth != 5)
swc.d.long.sub <- swc.d.long %>% subset(date >= as.Date("2012-01-01") & date < as.Date("2017-05-31"))

nse.max.10.d <- max(nse.table$sat10_daily, na.rm = TRUE)
nse.max.40.d <- max(nse.table$sat40_daily, na.rm = TRUE)
nse.max.100.d <- max(nse.table$sat100_daily, na.rm = TRUE)
max.par.n.swc.row.d <- c(which(nse.table$sat10_daily == nse.max.10.d),  
               which(nse.table$sat40_daily == nse.max.40.d), 
               which(nse.table$sat100_daily == nse.max.100.d)) 
max.par.n.swc.d <- nse.table$par.sam[max.par.n.swc.row.d]
# #nse.text <- as.character(as.expression(italic(NSE[10, MAX])~"="~nse.max.10))
# p1  <- ggplot(swc.d.long.sub, aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 1) +
#   geom_line(data = swc.d.long.sub %>% subset(par.sam %in% max.par.n.swc.d),
#             aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
#   geom_line(aes(x= date, y = obs, color = "observed"), size = 1) +
#   # geom_text(x = 0.2, y = 0.8, label = nse.text) +
#   scale_colour_manual(name = "", values = col1) +
#   facet_grid(depth ~ .) +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = "top", legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Obs vs Model. Daily")
# p1
# ggsave(file.path("figures", current.folder, "swc_Obs_vs_model_daily_with_extra_sm.jpeg"), height = 8.94, width = 12.9, units='in')
# #p1 %+% sat.d.long.sub

swc.d.long.sub.mini <- subset(swc.d.long.sub, date >= as.Date("2016-01-01") & 
                                date < as.Date("2017-05-31") & depth !=5)
p1 %+% swc.d.long.sub.mini 
ggplot(swc.d.long.sub.mini, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
  geom_line(data = swc.d.long.sub.mini %>% subset(par.sam %in% max.par.n.swc.d), 
            aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "observed"), size = 1) +
  geom_text(data = swc.d.long.sub %>% subset(par.sam %in% max.par.n.swc.d[1]),
            aes(x = as.Date("2016-03-01"), y = 0.5, label = 
                  as.character(paste0("par.sam = ", par.sam)))) +
  scale_colour_manual(name = "", values = col1) +
  facet_grid(depth ~ .) +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = "top", legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
  ggtitle("Soil Water Content: Obs vs Model. Daily")
ggsave(file.path("figures", current.folder, "swc_Obs_vs_model_daily.jpeg"), height = 8.94, width = 12.9, units='in')
# # only looking at 2016 right now---
# mod.swc <- bind_rows(mod.swc.d.1) %>% bind_rows(mod.swc.d.3) %>% 
#   bind_rows(mod.swc.d.6) %>% bind_rows(mod.swc.d.10) %>% bind_rows(mod.swc.d.40) %>% bind_rows(mod.swc.d.100) %>%
#   bind_rows(mod.swc.d.300) %>% bind_rows(mod.swc.d.500) %>% bind_rows(mod.swc.d.800) %>%
#   mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()
# head(mod.swc)
# mod.swc <- mod.swc %>% subset(!is.na(depth))
# mod.swc.sub <- mod.swc %>% subset(date >= as.Date("2012-01-01") & date < as.Date("2017-05-31"))
# 
# m1 <- ggplot(mod.swc, aes(x = date, y = value)) +
#   geom_line(aes(group = depth, color = as.factor(depth)), show.legend = TRUE, size = 1) +
#   geom_line(data = obs.swc.d , aes(y = obs, group = depth, color = as.factor(depth)), show.legend = F, size = 0.2) +
#   scale_colour_discrete("Depth (cm)") +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = "top", legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Obs vs Model. Daily")
# m1
# m2 <- ggplot(mod.swc, aes(x = date, y = value)) +
#   geom_line(aes(group = depth, color = as.factor(depth)), show.legend = FALSE, size = 1) +
#   # scale_colour_gradient(high = "#132B43", low = "#f7fbff") + 
#   facet_grid(depth ~ .) +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = "top", legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Obs vs Model. Daily")
# m2

####
# Monthly #------
####
## from monthly files #-----
# load(file.path("data-raw/extract", current.folder, "extract/QRUNOFF.h0.extract.Rdata")) ##---
# mod.swc.m <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%# in mm/s
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
# swc.m.long <- swc.d.long %>%  subset(!is.na(obs)) %>%
#   mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
#   group_by(yrmo, par.sam, depth) %>% 
#   summarise(obs = mean(obs, na.rm = T), value = mean(value, na.rm = T), date = min(date)) %>% 
#   as.data.frame()
# swc.m.10 <- swc.m.long %>% subset(depth == "10" & !is.na(par.sam)) %>% select(-depth, -date) %>%
#   spread(key = par.sam, value = value) 
# swc.m.40 <- swc.m.long %>% subset(depth == "40" & !is.na(par.sam)) %>% select(-depth, -date) %>%
#   spread(key = par.sam, value = value) 
# swc.m.100 <- swc.m.long %>% subset(depth == "100" & !is.na(par.sam)) %>% select(-depth, -date) %>%
#   spread(key = par.sam, value = value) 
# ###
# rmse.table$swc100_monthly <- rmse.table$swc40_monthly <- rmse.table$swc10_monthly <- NA
# nse.table$swc100_monthly <- nse.table$swc40_monthly <- nse.table$swc10_monthly <- NA
# for (i in 1: nsam) {
#   variance <- var(swc.m.10$obs, swc.m.10[, i + 2], use = "complete.obs")
#   rmse.table$swc10_monthly[i] <- round(as.numeric(variance)^0.5, 3)
#   variance <- var(swc.m.40$obs, swc.m.40[, i + 2], use = "complete.obs")
#   rmse.table$swc40_monthly[i] <- round(as.numeric(variance)^0.5, 3)
#   variance <- var(swc.m.100$obs, swc.m.100[, i + 2], use = "complete.obs")
#   rmse.table$swc100_monthly[i] <- round(as.numeric(variance)^0.5, 3)
#   nse.table$swc10_monthly[i] <- NSE(swc.m.10$obs, swc.m.10[, i + 2], na.rm = TRUE)
#   nse.table$swc40_monthly[i] <- NSE(swc.m.40$obs, swc.m.40[, i + 2], na.rm = TRUE)
#   nse.table$swc100_monthly[i] <- NSE(swc.m.100$obs, swc.m.100[, i + 2], na.rm = TRUE)
# }
# summary(rmse.table); summary(nse.table)
# swc.m.long.sub <- swc.m.long %>% subset(date >= as.Date("2012-01-01") & date < as.Date("2017-05-31"))
# 
# ggplot(swc.m.long.sub, aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
#   geom_line(aes(x= date, y = obs, color = "observed"), size = 1) +
#   scale_colour_manual(name = "", values = col1) +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#   legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Obs vs Model. Monthly")
# ggsave(file.path("figures", current.folder, "swc_Obs_vs_model_monthly.jpeg"), height = 8.94, width = 12.9, units='in')


###
#### GPP: Obs versus model #-------
####
# Daily #------
####
# load(file.path("data-raw/extract", current.folder, "extract/GPP.h1.extract.Rdata"))
# 
# mod.btran.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
#   mutate(date = as.Date(date)) %>%
#   as.data.frame()
# btran.d <- mod.btran.d
# head(btran.d[, 1:6]); summary(btran.d[, 1:6])
# btran.d.long <- gather(btran.d, key = "par.sam", "value", -date) 
# btran.d.long.sub <- btran.d.long %>% subset(date >= as.Date("2008-01-01"))
# 
# b.d <- ggplot(btran.d.long.sub,
#        aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.1) +
#   scale_colour_manual(name = "", values = col1) +
#   ylab("GPP [unitless]") +
#   xlab("Date") + theme_bw() +
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
#   ggtitle("GPP: Model. Daily")
# b.d + geom_line(data = btran.d.long.sub %>% subset(par.sam %in% max.par.n.et.d), 
#                 aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5)
# ggsave(file.path("figures", current.folder, "GPP_model_daily_bestfit_et.d.jpeg"), height = 8.94, width = 12.9, units='in')
# b.d + geom_line(data = btran.d.long.sub %>% subset(par.sam %in% max.par.n.swc.d), 
#                 aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5)
# ggsave(file.path("figures", current.folder, "GPP_model_daily_bestfit_swc.d.jpeg"), height = 8.94, width = 12.9, units='in')

## Monthly from daily files #-------
## from monthly files #-----
load(file.path("data-raw/extract", current.folder, "extract/GPP.h0.extract.Rdata")) ##---
mod.btran.m <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%
  as.data.frame()
btran.m.long <- gather(mod.btran.m, key = "par.sam", "value", -yrmo) 
btran.m.long <- btran.m.long %>%  
  # mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, par.sam) %>% 
  summarise(value = mean(value, na.rm = T)) %>% 
  as.data.frame()
str(btran.m.long)
btran.m.long$date <- as.Date(strptime(paste0(btran.m.long$yrmo, "-01"), "%Y-%m-%d"))
btran.m <- spread(btran.m.long, key = par.sam, value = value)
summary(btran.m[, 1:3])
# btran.m.long.sub <- btran.m.long %>% subset(date >= as.Date("2008-01-01"))
## add best-fits for ET
b.m <- ggplot(btran.m.long,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "simulated"), show.legend = F, size = 0.2) +
  scale_colour_manual(name = "", values = col1) +
  ylab("GPP [unitless]") +
  xlab("Date") + theme_bw() +
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.5), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Modelled Monthly")
b.m + geom_line(data = btran.m.long %>% subset(par.sam %in% max.par.n.et.m), 
                  aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5)
ggsave(file.path("figures", current.folder, "GPP_model_monthly_bestfit_et.m.jpeg"), height = 8.94, width = 12.9, units='in')
b.m + geom_line(data = btran.m.long %>% subset(par.sam %in% max.par.n.swc.d), 
                aes(group = par.sam, color = "best-fit"), show.legend = F, size = 0.5)
ggsave(file.path("figures", current.folder, "GPP_model_monthly_bestfit_swc.d.jpeg"), height = 8.94, width = 12.9, units='in')



###------------------
## Plotting nse table
###------------------

nse.long <- nse.table %>% gather(key = "variable", value = "nse", -par.sam)
ggplot(nse.long, aes(x = soil_depth, y = nse)) +
  geom_bar(aes(fill = variable), stat = "identity", show.legend = F) +
  # scale_fill_gradient(guide = guide_fill(reverse = TRUE)) +
  facet_grid(variable ~ ., scales = "free_y") +
  ylab("nse") +
  xlab("Soil Depth [m]") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain")) +
  ggtitle("Nash-Sutcliffe Efficiency between Observation and Model By Soil Depth")
ggsave(file.path("figures", current.folder, "nse between Observation and Model By Soil Depth.jpeg"), height = 8.94, width = 12.9, units='in')

head(nse.table)
summary(nse.table)
params.nse <- par.on %>% full_join(nse.table, by = "par.sam")
params.nse[max.par.n.swc,]
params.nse[max.par.n.et,]
params.nse[max.par.n.run,]
params.nse[max.par.n.swc,]
## Parameters sensitivity to obervations
par.on.long <- par.on %>% gather(key = "variable", value = "var.value", -par.sam) 

par.gpp.m.long <- par.on.long %>% full_join(gpp.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
yrmo.vec <- c("2015-02", "2015-07", "2016-04", "2016-07")
ggplot(par.gpp.m.long, aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
  #geom_smooth(method = "auto") +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Evapotranspiration [mm/month]") +
  xlab("Variable") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain"),
        strip.text = element_text(size = 14, face = "plain")) +
  ggtitle("Sensitivity of monthly ET to parameters")
ggsave(file.path("figures", current.folder, "Sensitivity of monthly ET to parameters_",
                        yrmo.vec[1], ".jpeg"), height = 8.94, width = 12.9, units='in')

par.qet.m.long <- par.on.long %>% full_join(qet.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
# yrmo.vec <- c("2015-02", "2015-07", "2016-04", "2016-07")
ggplot(par.qet.m.long, aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
  #geom_smooth(method = "auto") +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Evapotranspiration [mm/month]") +
  xlab("Variable") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain"),
        strip.text = element_text(size = 14, face = "plain")) +
  ggtitle("Sensitivity of monthly ET to parameters")
ggsave(file.path("figures", current.folder, "Sensitivity of monthly ET to parameters_",
                        yrmo.vec[1], ".jpeg"), height = 8.94, width = 12.9, units='in')

par.qrunoff.m.long <- par.on.long %>% full_join(qrunoff.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
# par.qrunoff.m.long %>% subset(yrmo == yrmo.vec[2])
ggplot(par.qrunoff.m.long, aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("QRUNOFF [mm/month]") +
  xlab("Variable") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain"),
        strip.text = element_text(size = 14, face = "plain")) +
  ggtitle("Sensitivity of monthly runoff to parameters")
ggsave(file.path("figures", current.folder, "Sensitivity of monthly runoff to parameters_",
                        yrmo.vec[2], ".jpeg"), height = 8.94, width = 12.9, units='in')

par.swc.m.long <- par.on.long %>% full_join(swc.m.long %>% mutate(par.sam = as.integer(par.sam)) %>%
                                              subset(), by = "par.sam")

date.vec <- as.Date(c("2016-03-01", "2016-04-01", "2016-07-01", "2016-10-01"))
sen.swc.100 <- ggplot(par.swc.m.long %>% 
                        subset(depth == 100 & !is.na(variable)), # & date == date.vec[3] 
                      aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Variable") + theme_bw() +
  theme(axis.text.x = element_text(size = 14, face = "plain"),
        strip.text = element_text(size = 14, face = "plain")) +
  ggtitle("Sensitivity of monthly soil water content at 1 m to parameters")
sen.swc.100
ggsave(file.path("figures", current.folder, "Sensitivity of monthly soil water content at 1 m to parameters_",
                        date.vec[3], ".jpeg"), height = 8.94, width = 12.9, units='in')

sen.swc.100 %+% subset(par.swc.m.long, depth == 10 & !is.na(variable)) + #& date == date.vec[3]
  ggtitle("Sensitivity of monthly soil water content at 10 cm to parameters")
ggsave(file.path("figures", current.folder, "Sensitivity of monthly soil water content at 10 cm to parameters_",
                        date.vec[3], ".jpeg"), height = 8.94, width = 12.9, units='in')

write.csv(nse.table, file = "data/nse.table.csv", row.names = FALSE)
