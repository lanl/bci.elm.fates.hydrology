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


current.folder <- "2019-10-14_5000/trial"
if(!dir.exists(file.path("data-raw", current.folder))) {dir.create(file.path("data-raw", current.folder))}
if(!dir.exists(file.path("results", current.folder))) {dir.create(file.path("results", current.folder))}
if(!dir.exists(file.path("figures", current.folder))) {dir.create(file.path("figures", current.folder))}


params <- read.csv(file.path("data-raw/params.csv"), header = TRUE)
dim(params); str(params)

head(params)
total.n <- nrow(params)

# until folder 500 #------
# params <- fates.params
# for (i in 1: (surf.n - 1)) {
#   params <- bind_rows(params, fates.params)
# }
# head(surf.params)
# params <- params %>% mutate(FMAX = rep(surf.params$FMAX, each = fates.n),
#                             aveDTB = rep(surf.params$aveDTB, each = fates.n),
#                             HKSAT_ADJ = rep(surf.params$HKSAT_ADJ, each = fates.n),
#                             HKSAT_12.5 = rep(as.numeric(ksat[4, 1 + (1:surf.n)]), each = fates.n), # at 12.5 cm
#                             HKSAT_60 = rep(as.numeric(ksat[7, 1 + (1:surf.n)]), each = fates.n),
#                             fpi_max = rep(elm.params$fpi_max, length.out = total.n)) # at 60 cm
# ----
## Make sure to get matching Filter.txt file from server
filterFile <- read.table(file.path("data-raw/extract", current.folder, "Filter.txt"))
nrow(filterFile)

col1 <- c("Observed" = "#f04546", "Simulated" = "#3591d1", "Best-fit Simulated" = "green", 
          "Best-fit NSE" = "green", "Best-fit R-squared" = "purple", 
          "Observed at Ava-Tower" = "#f04546", "Observed Plot-wide" = "black")

# 
# ### test ----
# params <- params.from.files[c(3, 3, 3, 3, 3), ]
# params$nsam <- 1
# nsam = 4
# filterFile <- read.table("data-raw/Filter.txt")
# filterFile <- data.frame(V1 = filterFile$V1[1])
### ---
params <- params %>% mutate(par.sam = 1:total.n)
# par.sam.on <- params$par.sam[1000:5000][filterFile$V1[1000:5000] == TRUE]
# par.on <- params %>% subset(par.sam %in% par.sam.on) 
par.on.1 <- params[params$par.sam == 3638,]
par.on <- par.on.1[rep(row.names(par.on.1), 5),]
nsam <- length(par.on$par.sam) 
rmse.table <- data.frame(par.sam = par.on$par.sam, "qrunoff_daily" = rep(NA, nsam), n = c(7001, 7002, 7003, 7004, 7005))
rsq.table <- data.frame(par.sam = par.on$par.sam, "qrunoff_daily" = rep(NA, nsam), n = c(7001, 7002, 7003, 7004, 7005))
nse.table <- data.frame(par.sam = par.on$par.sam, "qrunoff_daily" = rep(NA, nsam), n = c(7001, 7002, 7003, 7004, 7005))
###
#### Runoff: Obs versus model #-------
####
# Daily #------
####
## QRUNOFF total liquid runoff (does not include QSNWCPICE)
## Since total it probably includes QOVER (surface runoff) + qh2osfc is runoff from surface water storage
load(file.path("data-raw/extract", current.folder, "extract/QRUNOFF.h1.extract.Rdata"))
mod.qrunoff.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date))
mod.qrunoff.d[, -1] <- mod.qrunoff.d[, -1]*24*60*60/3 # converting mm/s to mm/day
obs.qrunoff.d <- bci.hydromet::forcings %>% select(date, flow_conrad) %>% # in mm/day
  rename(obs = flow_conrad) %>% subset(date > min(mod.qrunoff.d$date, na.rm = TRUE))

qrunoff.d <- obs.qrunoff.d %>% full_join(mod.qrunoff.d, by = "date")
#head(qrunoff.d); #summary(qrunoff.d)
rm(var.res.arr) # large file
# letting four years pass to allow for model initialisation before model-obs comparison
qrunoff.d.sub <- qrunoff.d %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4))
qrunoff.d.long <- gather(qrunoff.d, key = "par.sam", "value", -date, -obs) 

for (i in 1: nsam) {
  obs <- as.numeric(qrunoff.d.sub$obs)
  sim <- as.numeric(qrunoff.d.sub[, i + 2])
  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(qrunoff.d[i + 2])) # same for all other tables
  # i <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  
  rmse.table$qrunoff_daily[i] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE)^0.5), 3)
  rsq.table$qrunoff_daily[i] <- round(as.numeric(correlation)^2, 3)
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$qrunoff_daily[i] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
}

summary(rsq.table); summary(nse.table)
nse.max.run.d <- max(nse.table$qrunoff_daily, na.rm = TRUE) 
rmse.min.run.d <- min(rmse.table$qrunoff_daily, na.rm = TRUE)
rsq.max.run.d <- max(rsq.table$qrunoff_daily, na.rm = TRUE)
max.par.n.run.row.d <- c(which(nse.table$qrunoff_daily == nse.max.run.d, which(rsq.table$qrunoff_daily == rsq.max.run.d))) 
max.par.n.run.d <- nse.table$par.sam[max.par.n.run.row.d]
## rsq for the ensemble with max nse
rsq.for.max.nse.run.d <- rsq.table$qrunoff_daily[which(nse.table$qrunoff_daily == nse.max.run.d)]

qrunoff.d.long.sub <- qrunoff.d.long %>% subset(date >= as.Date("2008-01-01"))
col2 <- c("Observed" = "#f04546",  "Half HKsatADJ" = "blue", "Full BBslope" = "purple", 
          "Half BBslope" = "green", "Half BBslope+ Half Vcmax" = "pink", "Half Vcmax" = "magenta")

p.qrun.d <- ggplot(qrunoff.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam == 7005), aes(group = par.sam, color = "Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam == 7001), aes(group = par.sam, color = "Half HKsatADJ"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam == 7004), aes(group = par.sam, color = "Half BBslope+ Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam == 7003), aes(group = par.sam, color = "Full BBslope"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.d.long.sub %>% subset(par.sam == 7002), aes(group = par.sam, color = "Half BBslope"), show.legend = F, size = 0.5) +
  geom_line(aes(y = obs, colour = "Observed"), size = 0.5) +
  scale_colour_manual(name = "", values = col2) +
  # geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.1) +
  # geom_line(data = qrunoff.d.long.sub %>% subset(par.sam %in% max.par.n.run.d),
  #           aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
  # geom_text(data = qrunoff.d.long.sub %>% subset(par.sam %in% max.par.n.run.d),
  #           aes(x = min(date) + 700, y = Inf, label =
  #                 as.character(paste0("NSE.max = ", round(nse.max.run.d, 2), "; R-squared = ", round(max(rsq.for.max.nse.run.d, na.rm = TRUE), 2),
  #                                     "\n RMSE.max = ", round(rmse.min.run.d, 2),
  #                                     "\n R-squared.max = ", round(rsq.max.run.d, 2)))), vjust = 2) +


  ylab("QRUNOFF [mm/day]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Observed vs Simulated_Daily")
ggsave("QRUNOFF_Obs_vs_model_daily.jpeg", plot = p.qrun.d, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')
rm(p.qrun.d, qrunoff.d.long.sub)
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
  summarise(obs = sum(obs), ## removing na.rm = TRUE so as to only get sum for months when all observations present
            value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
rm(qrunoff.d.long)
# letting four years pass to allow for model initialisation before model-obs comparison
qrunoff.m <- qrunoff.m.long %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4)) %>% 
  select(-date) %>% pivot_wider(names_from = par.sam, values_from = value) %>% as.data.frame()
summary(qrunoff.m[, 1:3])

rmse.table$qrunoff_monthly <- NA; nse.table$qrunoff_monthly <- NA; rsq.table$qrunoff_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(qrunoff.m$obs); sim <- as.numeric(qrunoff.m[, i + 2])
  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(qrunoff.m[i + 2])) # same for all other tables
  # i <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  
  rmse.table$qrunoff_monthly[i] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
  rsq.table$qrunoff_monthly[i] <- round(as.numeric(correlation)^2, 3)
  
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$qrunoff_monthly[i] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
}

summary(rmse.table);  summary(nse.table)
nse.max.run.m <- max(nse.table$qrunoff_monthly, na.rm = TRUE) 
rmse.min.run.m <- min(rmse.table$qrunoff_monthly, na.rm = TRUE)
rsq.max.run.m <- max(rsq.table$qrunoff_monthly, na.rm = TRUE)
nse.max.par.n.run.m <- nse.table$par.sam[which(nse.table$qrunoff_monthly == nse.max.run.m)]
rsq.max.par.n.run.m <- rsq.table$par.sam[which(rsq.table$qrunoff_monthly == rsq.max.run.m)]

## rsq for the ensemble with max nse
rsq.for.max.nse.run.m <- rsq.table$qrunoff_monthly[which(nse.table$qrunoff_monthly == nse.max.run.m)]

qrunoff.m.long.sub <- qrunoff.m.long %>% subset(date >= as.Date("2008-01-01"))
# nse.text <- as.character(as.expression(italic(NSE[MAX])~" = "~nse.max.run))
p.qrun.m <- ggplot(qrunoff.m.long.sub,
       aes(x = date, y = value)) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam == 7005), aes(group = par.sam, color = "Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam == 7001), aes(group = par.sam, color = "Half HKsatADJ"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam == 7004), aes(group = par.sam, color = "Half BBslope+ Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam == 7003), aes(group = par.sam, color = "Full BBslope"), show.legend = F, size = 0.5) +
  geom_line(data = qrunoff.m.long.sub %>% subset(par.sam == 7002), aes(group = par.sam, color = "Half BBslope"), show.legend = F, size = 0.5) +
  geom_line(aes(y = obs, colour = "Observed"), size = 0.5) +
  # geom_text(data = qrunoff.m.long.sub %>% subset(par.sam %in% nse.max.par.n.run.m),
  #           aes(x = min(date) + 700, y = Inf, label = 
  #                 as.character(paste0("NSE.max = ", round(nse.max.run.m, 2), "; R-squared = ", round(max(rsq.for.max.nse.run.m, na.rm = TRUE), 2), 
  #                                     "\n RMSE.max = ", round(rmse.min.run.m, 2),
  #                                     "\n R-squared.max = ", round(rsq.max.run.m, 2)))), vjust = 2) +
  scale_colour_manual(name = "", values = col2) +
  ylab("QRUNOFF [mm/month]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("QRUNOFF: Observed vs Simulated_Monthly")
ggsave("QRUNOFF_Obs_vs_model_monthly.jpeg", plot = p.qrun.m, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')
rm(p.qrun.m)

write.csv(qrunoff.d.sub, file = file.path("results", current.folder, "qrunoff.d.sub.csv"), row.names = FALSE)
write.csv(qrunoff.m, file = file.path("results", current.folder, "qrunoff.m.csv"), row.names = FALSE)
# do not remove qrunoff.m.long
rm(qrunoff.d, qrunoff.d.sub, qrunoff.d.long, qrunoff.d.long.sub, qrunoff.m, qrunoff.m.long.sub)
rm(obs.qrunoff.d, mod.qrunoff.d)
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
n.month <- length(unique(format(obs.qet.d$date[!is.na(obs.qet.d$obs)], format = "%Y-%m")))
df.factor <- 1
set.order <- 21 ## how many days to average over
sma.sp <- smooth::sma(obs.qet.d$obs.2, order = set.order)
# sp <- smooth.spline(obs.qet.d$obs.2, df = n.month*df.factor) # obs.qet.d$obs[!is.na(obs.qet.d$obs)]
# obs.qet.d$sp <- sp$y
obs.qet.d$sp <- as.numeric(sma.sp$fitted) 
plot(obs.2 ~ date, type = "l", data = obs.qet.d, col = "pink", ylim = range(obs.qet.d$obs.2, na.rm = TRUE))
lines(sma ~ date, type = "l", data = obs.qet.d, col = "red")
# joining obs with mod
qet.d <- obs.qet.d %>% full_join(mod.qet.d, by = "date") 
# letting four years pass to allow for model initialisation before model-obs comparison; 
# does not have to be done for ET because observed ET is 3.5 years ahead of the run start date, which is 2008-01-1
qet.d.sub <- qet.d %>% subset(date >= c(min(mod.qet.d$date, na.rm = TRUE) + 365*4))

# points(obs.qet.d$date, obs.qet.d$obs, col = "red")
dates.valid.to.compare <- seq.Date(from = c(min(qet.d.sub$date, na.rm = TRUE) + 30), to =  c(max(qet.d.sub$date, na.rm = TRUE) - 30), by = "day")
# qet.sp <- qet.d.sub
# for (i in 1:nsam){
#   rows.present <- !is.na(qet.d.sub[, i + 4])
#   # only works when no gaps
#   sp.i <- smooth::sma(obs.qet.d$obs.2, order = set.order)
#   # sp.i <- smooth.spline(qet.d.sub[rows.present, i + 4], df = n.month*df.factor)
#   sp.date <- data.frame(date = qet.d.sub$date[rows.present], y = sp.i$y)
#   # getting back dates when for which data are missing
#   sp.date.full <- qet.sp %>% select(date) %>% full_join(sp.date, by = "date") %>% dplyr::arrange(date)
#   qet.sp[, i + 4] <- sp.date.full$y
# }
 
qet.sma <- qet.d.sub 
for (i in 1:nsam){
  sma.i <- smooth::sma(qet.d.sub[, i + 4], order = set.order)
  qet.sma[, i + 4] <- as.numeric(sma.i$fitted)
}
#head(qet.d.sub); summary(qet.d.sub)
rsq.table$qet_daily <- NA; rmse.table$qet_daily <- NA; nse.table$qet_daily <- NA
# rsq.table$qet_spline <- NA; rmse.table$qet_spline <- NA; nse.table$qet_spline <- NA
qet.sma.sub <- qet.sma %>% subset(date %in% dates.valid.to.compare)
# qet.sp.sub <- qet.sp %>% subset(date %in% dates.valid.to.compare)
for (i in 1: nsam) {
  # obs <- as.numeric(qet.d.sub$obs); sim <- as.numeric(qet.d.sub[, i + 4])
  # obs <- as.numeric(qet.sp.sub$sp); sim <- as.numeric(qet.sp.sub[, i + 4])
  obs <- as.numeric(qet.sma.sub$sp); sim <- as.numeric(qet.sma.sub[, i + 4])
  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(qet.d.sub[i + 4])) # same for all other tables
  # i <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  
  rmse.table$qet_daily[i] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
  rsq.table$qet_daily[i] <- round(as.numeric(correlation)^2, 3)
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$qet_daily[i] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
}
summary(rmse.table); summary(nse.table); summary(rsq.table)
rmse.table; nse.table; rsq.table
nse.max.et.d <- max(nse.table$qet_daily, na.rm = TRUE)
rmse.min.et.d <- min(rmse.table$qet_daily, na.rm = TRUE)
rsq.max.et.d <- max(rsq.table$qet_daily, na.rm = TRUE)
# max.par.n.et.row.d <- c(which(nse.table$qet_daily == nse.max.et.d), which(rsq.table$qet_daily == rsq.max.et.d))
nse.max.par.n.et.d <- nse.table$par.sam[which(nse.table$qet_daily == nse.max.et.d)]
rsq.max.par.n.et.d <- nse.table$par.sam[which(rsq.table$qet_daily == rsq.max.et.d)]
## rsq for the ensemble with max nse
rsq.for.max.nse.et.d <- rsq.table$qet_daily[which(nse.table$qet_daily == nse.max.et.d)]
top.quant <- 0.95 # 
top.par.n.et.d <- unique(rsq.table$par.sam[rsq.table$qet_daily > as.numeric(quantile(rsq.table$qet_daily, prob = top.quant, na.rm = TRUE))],
                         nse.table$par.sam[nse.table$qet_daily > as.numeric(quantile(nse.table$qet_daily, prob = top.quant, na.rm = TRUE))])
                         
# qet.d.long <- gather(qet.d.sub, key = "par.sam", "value", -date, -obs, -obs.2)
# qet.d.long <- gather(qet.sp, key = "par.sam", "value", -date, -obs, -obs.2, -sp)
qet.d.long <- gather(qet.d.sub, key = "par.sam", "value", -date, -obs, -obs.2, -sp)
qet.d.long.sub <- qet.d.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31")) %>% droplevels() %>% 
  subset(par.sam %in% c(top.par.n.et.d, nse.max.par.n.et.d, rsq.max.par.n.et.d))
qet.d.long.sub.best.fit.nse <- qet.d.long.sub %>% subset(par.sam %in% nse.max.par.n.et.d)
qet.d.long.sub.best.fit.rsq <- qet.d.long.sub %>% subset(par.sam %in% rsq.max.par.n.et.d)

## for sma
qet.d.long.sma <- gather(qet.sma, key = "par.sam", "value", -date, -obs, -obs.2, -sp) %>%
  subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31")) %>% droplevels() 
qet.d.long.sub.sma <- qet.d.long.sma %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31")) %>% droplevels() %>% 
  subset(par.sam %in% c(top.par.n.et.d, nse.max.par.n.et.d, rsq.max.par.n.et.d))
qet.d.long.sub.best.fit.sma.nse <- qet.d.long.sub.sma %>% subset(par.sam %in% nse.max.par.n.et.d)
qet.d.long.sub.best.fit.sma.rsq <- qet.d.long.sub.sma %>% subset(par.sam %in% rsq.max.par.n.et.d)

p.et1.label.text <- geom_text(aes(x = min(date) + 700, y = Inf, label = 
                                    as.character(paste0("NSE.max = ", round(nse.max.et.d, 2), "; R-squared = ", round(max(rsq.for.max.nse.et.d, na.rm = TRUE), 2),
                                                        "\n RMSE.max = ", round(rmse.min.et.d, 2),
                                                        "\n R-squared.max = ", round(rsq.max.et.d, 2)))), vjust = 2)

col2 <- c("Observed" = "#f04546", "Simulated" = "#3591d1", "Best-fit Simulated" = "green")

p.et1 <- ggplot(qet.d.long.sub,
                aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.1) +
  scale_colour_manual(name = "", values = col1) +
  ylab("ET [mm/day]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        plot.margin = margin(5.1, 4.1, 4.1, 2.1, "pt"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Observed vs Simulated_Daily")

p.et1.a <- p.et1 + 
  geom_line(data = qet.d.long.sub.best.fit.nse, aes(group = par.sam, color = "Best-fit NSE"), show.legend = F, size = 0.5) +
  # geom_line(data = qet.d.long.sub.best.fit.rsq, aes(group = par.sam, color = "Best-fit R-squared"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs.2, color = "Observed"), size = 0.5) + 
  p.et1.label.text
ggsave("ET_Obs_vs_model_daily_all_years.jpeg", plot = p.et1.a, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')
p.et1.b <- p.et1 + 
  geom_line(data = qet.d.long.sub.best.fit.nse, aes(group = par.sam, color = "Best-fit NSE"), show.legend = F, size = 0.5) +
  # geom_line(data = qet.d.long.sub.best.fit.rsq, aes(group = par.sam, color = "Best-fit R-squared"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.5) + 
  p.et1.label.text
ggsave("ET_Obs_vs_model_daily_all_years_with_gaps.jpeg", plot = p.et1.b, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')
# p.et1.c <- p.et1 + geom_line(data = qet.d.long.sub.best.fit, aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
#   geom_line(aes(x= date, y = sp, color = "Observed"), size = 0.5) + 
#   p.et1.label.text
# ggsave("ET_Obs_vs_model_daily_all_years_splines.jpeg", plot = p.et1.c, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

col1 <- c("Observed" = "#f04546", "Simulated" = "#3591d1", "Best-fit Simulated" = "green", 
          "Best-fit NSE" = "green", "Best-fit R-squared" = "purple", 
          "Observed at Ava-Tower" = "#f04546", "Observed Plot-wide" = "black")

p.et1.d <- ggplot(qet.d.long.sma, aes(x = date, y = value)) +
  geom_line(data = qet.d.long.sma %>% subset(par.sam == 7005), aes(group = par.sam, color = "Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qet.d.long.sma %>% subset(par.sam == 7001), aes(group = par.sam, color = "Half HKsatADJ"), show.legend = F, size = 0.5) +
  geom_line(data = qet.d.long.sma %>% subset(par.sam == 7004), aes(group = par.sam, color = "Half BBslope+ Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = qet.d.long.sma %>% subset(par.sam == 7003), aes(group = par.sam, color = "Full BBslope"), show.legend = F, size = 0.5) +
  geom_line(data = qet.d.long.sma %>% subset(par.sam == 7002), aes(group = par.sam, color = "Half BBslope"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = sp, color = "Observed"), size = 0.5) +  
  theme(legend.position = "top") +
  scale_colour_manual(name = "", values = col2) +
  ylab("ET [mm/day]") +
  xlab("Date") + 
  # theme(axis.text = element_text(size = 14, face = "plain"),
  #       legend.text = element_text(size = 16, face = "plain"),
  #       plot.margin = margin(5.1, 4.1, 4.1, 2.1, "pt"),
  #       plot.title = element_text(hjust = 0.5, size = 14, face = "plain"),
  #       legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y"))
  # geom_line(data = qet.d.long.sub.best.fit.sma.nse, aes(group = par.sam, color = "Best-fit NSE"), show.legend = F, size = 0.5) +
  # geom_line(data = qet.d.long.sub.best.fit.sma.rsq, aes(group = par.sam, color = "Best-fit R-squared"), show.legend = F, size = 0.5) +

  # p.et1.label.text
ggsave(paste0("ET_Obs_vs_model_daily_all_years_sma_", set.order, ".jpeg"), plot = p.et1.d, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

# qet.d.long.sub.2 <- qet.d.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2015-02-15"))
# qet.d.long.sub.best.fit.2 <- qet.d.long.sub.2 %>% subset(par.sam %in% max.par.n.et.d) %>% droplevels()
# p.et1 %+% qet.d.long.sub.2 +
#   geom_line(data = qet.d.long.sub.best.fit.2, 
#             aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
#   geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.5) +
#   scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y")) +
#   p.et1.label.text
# ggsave(plot = file.path("figures", current.folder, "ET_Obs_vs_model_daily_2012-2014.jpeg"), height = 6.25, width = 8.94, units='in')
# 
# qet.d.long.sub.3 <- qet.d.long %>% subset(date >= as.Date("2015-05-01") & date < as.Date("2017-08-31")) 
# qet.d.long.sub.best.fit.3 <- qet.d.long.sub.3 %>% subset(par.sam %in% max.par.n.et.d)
# p.et1 %+% qet.d.long.sub.3 +
#   geom_line(data = qet.d.long.sub.best.fit.3, 
#             aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
#   geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.5) +
#   scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y")) +
#   p.et1.label.text
# ggsave(plot = file.path("figures", current.folder, "ET_Obs_vs_model_daily_2015-2017.jpeg"), height = 6.25, width = 8.94, units='in')
# 
rm(p.et1, p.et1.d) # large.file
####
# Monthly #------
## using 

qet.m.long <- qet.d.long %>%  
  mutate("yrmo" = paste0(format(date, "%Y"), "-", format(date, "%m"))) %>%
  group_by(yrmo, par.sam) %>%
  ## using gap-filled AET
  summarise(obs = sum(obs.2), ## removing na.rm = TRUE so as to only get sum for months when all observations present
            value = sum(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
# letting four years pass to allow for model initialisation before model-obs comparison
qet.m <- qet.m.long %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4)) %>% 
  select(-date) %>% 
  pivot_wider(names_from = par.sam, values_from = value) %>% as.data.frame()
summary(qet.m[, 1:5])
rmse.table$qet_monthly <- NA; nse.table$qet_monthly <- NA; rsq.table$qet_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(qet.m$obs); sim <- as.numeric(qet.m[, i + 2])

  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(qet.m[i + 2])) # same for all other tables
  # row.sam <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  rmse.table$qet_monthly[i] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
  rsq.table$qet_monthly[i] <- round(as.numeric(correlation)^2, 3)
  
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$qet_monthly[i] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
} 
summary(rmse.table);  summary(nse.table)
nse.max.et.m <- max(nse.table$qet_monthly, na.rm = TRUE) 
rmse.min.et.m <- min(rmse.table$qet_monthly, na.rm = TRUE) 
rsq.max.et.m <- max(rsq.table$qet_monthly, na.rm = TRUE)
nse.max.par.n.et.m <- nse.table$par.sam[which(nse.table$qet_monthly == nse.max.et.m)]
rsq.max.par.n.et.m <- rsq.table$par.sam[which(rsq.table$qet_monthly == rsq.max.et.m)] 
## nse for the ensemble with max rsq
nse.table$qet_monthly[which(rsq.table$qet_monthly == rsq.max.et.m)]
## rsq for the ensemble with max nse
rsq.for.max.nse.et.m <- rsq.table$qet_monthly[which(nse.table$qet_monthly == nse.max.et.m)]

qet.m.long.sub <- qet.m.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31"))
qet.m.long.sub.best.fit.nse <- qet.m.long.sub %>% subset(par.sam %in% nse.max.par.n.et.m)
qet.m.long.sub.best.fit.rsq <- qet.m.long.sub %>% subset(par.sam %in% rsq.max.par.n.et.m)
qet.m.long.sub.1 <- qet.m.long.sub %>% subset(par.sam %in% rsq.max.par.n.et.m[1])
# rmse.label <- as.character(as.expression(italic(r)^2~"="~rmse.max))
p.et.m1 <- ggplot(qet.m.long.sub, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(data = qet.m.long.sub.best.fit.nse, aes(group = par.sam, color = "Best-fit NSE"), show.legend = F, size = 0.5) +
  geom_line(data = qet.m.long.sub.best.fit.rsq, aes(group = par.sam, color = "Best-fit R-squared"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 1) +
  geom_text(data = qet.m.long.sub %>% subset(par.sam %in% rsq.max.par.n.et.m),
            aes(x = min(date) + 700, y = Inf, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.et.m, 2), "; R-squared = ", round(max(rsq.for.max.nse.et.m, na.rm = TRUE), 2), 
                                      "\n RMSE.max = ", round(rmse.min.et.m, 2),
                                      "\n R-squared.max = ", round(rsq.max.et.m, 2)))), vjust = 2) +
  scale_colour_manual(name = "", values = col1) +
  ylab("ET [mm/month]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("Evapotranspiration: Observed vs Simulated_Monthly")
ggsave("ET_Obs_vs_model_monthly.jpeg", plot = p.et.m1, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

write.csv(qet.sma, file = file.path("results", current.folder, "qet.sma.csv"), row.names = FALSE)
write.csv(qet.sp, file = file.path("results", current.folder, "qet.sp.csv"), row.names = FALSE)
write.csv(qet.d, file = file.path("results", current.folder, "qet.d.csv"), row.names = FALSE)
write.csv(qet.m, file = file.path("results", current.folder, "qet.m.csv"), row.names = FALSE)
# do not remove qet.m.long
rm(qet.d, qet.d.long, qet.d.long.sub, qet.d.long.sub.best.fit, qet.m, qet.m.long.sub, qet.m.long.sub.1, qet.m.long.sub.best.fit)
rm(obs.qet.d, mod.qvege.d, mod.qvegt.d, mod.qet.d)
rm(p.et.m1)
# 
# m1 <- lm(value ~ obs, data = qet.m.long.sub.best.fit %>% subset(par.sam  == 627 ))
# summary(m1)
# ggplot(qet.m.long.sub,
#        aes(x = obs, y = value)) +
#   geom_point(aes(color = "Simulated"), show.legend = F) +
#   geom_point(data = qet.m.long.sub.best.fit, 
#             aes(color = "Best-fit Simulated"), show.legend = F) +
#   geom_smooth(data = qet.m.long.sub.best.fit, method = "lm") +
#   ylab("Simulated Best-fit Simulated [mm/month]") +
#   xlab("Observed [mm/month]") + 
#   theme(axis.text = element_text(size = 14, face = "plain")) +
#   # scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Evapotranspiration: Observed vs Simulated_Best-fit SimulatedMonthly_scatter")


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
# letting four years pass to allow for model initialisation before model-obs comparison
gpp.d.sub <- gpp.d %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4))
gpp.d.long <- gather(gpp.d, key = "par.sam", "value", -date, -obs) 

rmse.table$gpp_daily <- NA; nse.table$gpp_daily <- NA
for (i in 1: nsam) {
  obs <- as.numeric(gpp.d.sub$obs); sim <- as.numeric(gpp.d.sub[, i + 2])
  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(gpp.d[i + 2])) # same for all other tables
  row.sam <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  rmse.table$gpp_daily[row.sam] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
  rsq.table$gpp_daily[row.sam] <- round(as.numeric(correlation)^2, 3)
  
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$gpp_daily[row.sam] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
}
summary(rmse.table); summary(nse.table)

nse.max.gpp.d <- max(nse.table$gpp_daily, na.rm = TRUE); 
rmse.min.gpp.d <- min(rmse.table$gpp_daily, na.rm = TRUE); 
rsq.max.gpp.d <- max(rsq.table$gpp_daily, na.rm = TRUE)
max.par.n.gpp.row.d <- c(which(nse.table$gpp_daily == nse.max.gpp.d), which(rsq.table$gpp_daily == rsq.max.gpp.d))
max.par.n.gpp.d <- nse.table$par.sam[max.par.n.gpp.row.d]
## rsq for the ensemble with max nse
rsq.for.max.nse.gpp.d <- rsq.table$gpp_daily[which(nse.table$gpp_daily == nse.max.gpp.d)]

gpp.d.long.sub <- gpp.d.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31"))

p.gpp.d <- ggplot(gpp.d.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.1) +
  geom_line(data = gpp.d.long.sub %>% subset(par.sam %in% max.par.n.gpp.d), 
            aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
  scale_colour_manual(name = "", values = col1) +
  geom_line(aes(y = obs, colour = "Observed"), size = 0.5) +
  geom_text(data = gpp.d.long.sub %>% subset(par.sam %in% max.par.n.gpp.d),
            aes(x = min(date) + 700, y = Inf, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.gpp.d, 2), "; R-squared = ", round(max(rsq.for.max.nse.gpp.d, na.rm = TRUE), 2),
                                      "\n RMSE.max = ", round(rmse.min.gpp.d, 2),
                                      "\n R-squared.max = ", round(rsq.max.gpp.d, 2)))), vjust = 2) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Observed vs Simulated_Daily")
ggsave( "GPP_Obs_vs_model_daily.jpeg", plot = p.gpp.d, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

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
  summarise(obs = mean(obs, na.rm = T), ## removing na.rm = TRUE so as to only get sum for months when all observations present
            value = mean(value, na.rm = T), date = min(date)) %>% 
  as.data.frame()
# letting four years pass to allow for model initialisation before model-obs comparison
gpp.m <- gpp.m.long %>% subset(date >= c(min(date, na.rm = TRUE) + 365*4)) %>% 
  select(-date) %>% 
  pivot_wider(names_from = par.sam, values_from = value) %>% as.data.frame()
summary(gpp.m[, 1:3])
nse.table$gpp_monthly <- NA; rmse.table$gpp_monthly <- NA; rsq.table$gpp_monthly <- NA
for (i in 1: nsam) {
  obs <- as.numeric(gpp.m$obs); sim <- as.numeric(gpp.m[, i + 2])

  correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
  
  par.sam.mod <- as.numeric(colnames(gpp.m[i + 2])) # same for all other tables
  row.sam <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  rmse.table$gpp_monthly[row.sam] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
  rsq.table$gpp_monthly[row.sam] <- round(as.numeric(correlation)^2, 3)
  
  ## not computing efficiency when observed variable was zero
  obs.drop0 <- obs[!obs == 0]; sim.drop0 <- sim[!obs == 0]
  # NSE = 1 - ( sum( (obs - sim)^2 ) / sum( (obs - mean(obs))^2 )
  nse.table$gpp_monthly[row.sam] <- NSE(obs.drop0, sim.drop0, na.rm = TRUE) 
}

summary(rmse.table);  summary(nse.table)
nse.max.gpp.m <- max(nse.table$gpp_monthly, na.rm = TRUE) 
rmse.min.gpp.m <- min(rmse.table$gpp_monthly, na.rm = TRUE)
rsq.max.gpp.m <- max(rsq.table$gpp_monthly, na.rm = TRUE)
max.par.n.gpp.row.m <- c(which(nse.table$gpp_monthly == nse.max.gpp.m), which(rsq.table$gpp_monthly == rsq.max.gpp.m))
max.par.n.gpp.m <- nse.table$par.sam[max.par.n.gpp.row.m]
## rsq for the ensemble with max nse
rsq.for.max.nse.gpp.m <- rsq.table$gpp_monthly[which(nse.table$gpp_monthly == nse.max.gpp.m)]

gpp.m.long.sub <- gpp.m.long %>% subset(date >= as.Date("2012-07-01") & date < as.Date("2017-08-31"))
p.gpp.m <- ggplot(gpp.m.long.sub,
       aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(data = gpp.m.long.sub %>% subset(par.sam %in% max.par.n.gpp.m), 
            aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.7) +
  geom_text(data = gpp.m.long.sub %>% subset(par.sam %in% max.par.n.gpp.m),
            aes(x = min(date) + 700, y = Inf, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.gpp.m, 2), "; R-squared = ", round(max(rsq.for.max.nse.gpp.m, na.rm = TRUE), 2), 
                                      "\n RMSE.max = ", round(rmse.min.gpp.m, 2),
                                      "\n R-squared.max = ", round(rsq.max.gpp.m, 2)))), vjust = 2) +
  scale_colour_manual(name = "", values = col1) +
  ylab("GPP [gC/m^2/d]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.9, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("GPP: Observed vs Simulated_Monthly")
ggsave("GPP_Obs_vs_model_monthly.jpeg", plot = p.gpp.m, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

write.csv(gpp.d.sub, file = file.path("results", current.folder, "gpp.d.sub.csv"), row.names = FALSE)
write.csv(gpp.m, file = file.path("results", current.folder, "gpp.m.csv"), row.names = FALSE)

rm(gpp.d.long, gpp.d.sub, gpp.m.long.sub, obs.gpp.d, mod.gpp.d)
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
# 
# correlation <- cor(obs, sim, use = "complete.obs", method = "pearson")
# 
#   rmse.table$gpp_monthly[i] <- round(as.numeric(mean((sim - obs)^2, na.rm = TRUE))^0.5, 3)
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
# ## add Best-fit Simulateds for ET
# ggplot(gpp.m.long,
#               aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
#   geom_line(data = gpp.m.long %>% subset(par.sam %in% max.par.n.gpp.m), 
#             aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
#   geom_line(aes( x= date, y = obs, color = "Observed"), size = 0.7) +
#   scale_colour_manual(name = "", values = col1) +
#   ylab("GPP [gC/m^2/d]") +
#   xlab("Date") + 
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = c(0.8, 0.5), legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
#   ggtitle("GPP: Modelled Monthly")
# ggsave(plot = file.path("figures", current.folder, "GPP_model_monthly_bestfit.jpeg"), height = 6.25, width = 8.94, units='in')
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
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% gather(key = "par.sam", "value", -date, -depth)
for(i in 1: length(var.res.arr)){
  mod.swc.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
    mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
    gather(key = "par.sam", "value", -date, -depth)
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
         swc = swc/100) # converting from % to v/v fraction 
str(steph.sm); summary(steph.sm)
####-----
swc.vertical <- avasoilmoisture::vertical %>% rename(obs = swc) %>% 
  group_by(date, depth) %>% summarise(obs = mean(obs, na.rm = TRUE)) %>% 
  mutate(date.depth = paste0(date, ".", depth))
swc.vertical.mean <-  swc.vertical %>% group_by(date) %>% summarise(obs = mean(obs, na.rm = TRUE)) 

mod.swc.d.vert <- mod.swc %>%  
  # subset(date >= min(as.Date(swc.vertical$date))) %>% 
  subset(date >=min(as.Date(swc.vertical$date)) & depth %in% soil.depths[1:4]) %>% ## date filter since otherwise vector memory exhausted
  mutate(date.depth = paste0(date, ".", depth)) %>% 
  group_by(date, par.sam) %>%
  summarise_at("value", list(~mean(., na.rm = T))) %>% as.data.frame()

## joining vertical obs with mod
swc.d.long.vert <- mod.swc.d.vert %>%
  full_join(swc.vertical.mean, by = "date") %>% 
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
#bind_rows(swc.single) %>% mutate(depth = as.numeric(depth))

swc.d.long <- mod.swc.d %>% 
  full_join(obs.swc.d %>% select(-date, -depth), by = "date.depth") %>% 
  # separate(date.depth, c("date", "depth"), sep = ".", extra = "merge") %>% 
  as.data.frame() #%>% subset(!is.na(nsam))
head(swc.d.long); summary(swc.d.long)

### swc normalised from 0 - 1
norm.index <- function(x, ...){(x - min(x, ...))/c(max(x, ...) - min(x, ...))}
mod.sat.d <- mod.swc.d %>% 
  subset(date %in% unique(obs.swc.d$date)) %>%
  group_by(depth) %>% 
  mutate_at(c("value"), norm.index, na.rm = TRUE) %>% as.data.frame() 
# because maximas are overestimated due to macropore structure, especially for the 100 cm depth normalised swc looks skewed by one value

max.cut <- as.numeric(quantile(obs.swc.d$obs[which(obs.swc.d$depth == 100)], prob = 0.99)) # qunaitle at 0.99 = 0.490
max(obs.swc.d$obs[which(obs.swc.d$depth == 100)]) # 0.490
obs.sat.d <- obs.swc.d %>% 
  mutate(obs = if_else(obs >=max.cut & depth == 100, max.cut, obs)) %>%
  group_by(depth) %>%
  mutate_at(c("obs"), norm.index, na.rm = TRUE) %>% 
  mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame()

# ## mod soil moisture also peaks at Nov 23/24 2016
# ### tdr normalised or saturation.index(uses swc)
# ## because maximas are overestimated due to macropore structure, choose value after two days
# # norm.index2 <- function(x) {(x/max(x, na.rm = TRUE))}
# obs.tdr.daily <- avasoilmoisture::horizontal %>% 
#   rename(obs = sat.index) %>% group_by(date, depth) %>% # this is directly normalised tdr
#   summarise(obs = mean(obs, na.rm = TRUE)) %>% as.data.frame() 
# obs.tdr.daily.spread <- obs.tdr.daily %>% pivot_wider(names_from = depth, values_from = obs) %>% as.data.frame()
# obs.sat.daily.spread <- obs.tdr.daily.spread
# # col.max <- as.numeric(apply(obs.tdr.daily.spread[,-1],2, max, na.rm = T)) 
# # obs.sat.daily.spread[,2] <-  obs.tdr.daily.spread[, 2]/col.max[1]
# # obs.sat.daily.spread[,3] <-  obs.tdr.daily.spread[, 3]/col.max[2]
# # obs.sat.daily.spread[,4] <-  obs.tdr.daily.spread[, 4]/col.max[3]
# # mutate_at(c("obs"), norm.index2) %>%
# # obs.sat.d <- obs.sat.daily.spread %>%
# #   gather(key = depth, value = obs, -date) %>%
# #   mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() %>%
# #   mutate(depth = as.numeric(depth))
# obs.sat.d <- obs.tdr.daily %>%
#   mutate(date.depth = paste0(date, ".", depth)) %>% as.data.frame() %>%
#   mutate(depth = as.numeric(depth))
# #bind_rows(swc.single) %>% mutate(depth = as.numeric(depth))

sat.d.long <- obs.sat.d %>% select(-date) %>%
  full_join(mod.sat.d %>% select(-depth), by = "date.depth") %>% 
  # separate(date.depth, c("date", "depth"), sep = ".", extra = "merge") %>% 
  as.data.frame() #%>% subset(!is.na(nsam))
head(sat.d.long); summary(sat.d.long)

sat.d.long <- sat.d.long %>% subset(!is.na(depth))
# sat.d.long <- sat.d.long %>% subset(depth != 5)
sat.d.long.sub <- sat.d.long %>% subset(date >= as.Date("2012-01-01"))

sat.d.10 <- sat.d.long %>% subset(depth == "10" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 
sat.d.40 <- sat.d.long %>% subset(depth == "40" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 
sat.d.100 <- sat.d.long %>% subset(depth == "100" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 

swc.d.10 <- swc.d.long %>% subset(depth == "10" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 
swc.d.40 <- swc.d.long %>% subset(depth == "40" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 
swc.d.100 <- swc.d.long %>% subset(depth == "100" & !is.na(par.sam)) %>% select(-depth) %>%
  pivot_wider(names_from = par.sam, values_from = value, -date.depth) %>% as.data.frame() 
swc.d.vert <- swc.d.long.vert %>% subset(!is.na(par.sam)) %>% 
  pivot_wider(names_from = par.sam, values_from = value) %>% as.data.frame() 

rmse.table$sat100_daily <- rmse.table$sat40_daily <- rmse.table$sat10_daily <- NA
nse.table$sat100_daily <- nse.table$sat40_daily <- nse.table$sat10_daily <- NA
rsq.table$sat100_daily <- rsq.table$sat40_daily <- rsq.table$sat10_daily <- NA

rmse.table$swc.vert_daily <- rmse.table$swc100_daily <- rmse.table$swc40_daily <- rmse.table$swc10_daily <- NA
nse.table$swc.vert_daily <- nse.table$swc100_daily <- nse.table$swc40_daily <- nse.table$swc10_daily <- NA
rsq.table$swc.vert_daily <- rsq.table$swc100_daily <- rsq.table$swc40_daily <- rsq.table$swc10_daily <- NA

sat.obs.10 <- sat.d.10$obs; sat.obs.40 <- sat.d.40$obs; sat.obs.100 <- sat.d.100$obs; 
swc.obs.10 <- swc.d.10$obs; swc.obs.40 <- swc.d.40$obs; swc.obs.100 <- swc.d.100$obs; swc.obs.vert <- swc.d.vert$obs;
for (i in 1: nsam) {
  par.sam.mod <- as.numeric(colnames(swc.d.10[i + 2])) # same for all other tables
  row.sam <- which(nse.table$par.sam == par.sam.mod) ## same for rmse
  swc.sim.10 <- swc.d.10[, i + 2]
  swc.sim.40 <- swc.d.40[, i + 2]
  swc.sim.100 <- swc.d.100[, i + 2]
  swc.sim.vert <- swc.d.vert[, i + 2]
  ## normalised
  sat.sim.10 <- sat.d.10[, i + 2]
  sat.sim.40 <- sat.d.40[, i + 2]
  sat.sim.100 <- sat.d.100[, i + 2]
  ## for rmse
  # rmse.table$swc10_daily[row.sam] <- round(as.numeric(mean((swc.sim.10 - swc.obs.10)^2, na.rm = TRUE))^0.5, 3)
  # rmse.table$swc40_daily[row.sam] <- round(as.numeric(mean((swc.sim.40 - swc.obs.40)^2, na.rm = TRUE))^0.5, 3)
  # rmse.table$swc100_daily[row.sam] <- round(as.numeric(mean((swc.sim.100 - swc.obs.100)^2, na.rm = TRUE))^0.5, 3)
  # rmse.table$swc.vert_daily[row.sam] <- round(as.numeric(mean((swc.sim.vert - swc.obs.vert)^2, na.rm = TRUE))^0.5, 3)
  # 
  # rmse.table$sat10_daily[row.sam] <- round(as.numeric(mean((sat.sim.10 - sat.obs.10)^2, na.rm = TRUE))^0.5, 3)
  # rmse.table$sat40_daily[row.sam] <- round(as.numeric(mean((sat.sim.40 - sat.obs.40)^2, na.rm = TRUE))^0.5, 3)
  # rmse.table$sat100_daily[row.sam] <- round(as.numeric(mean((sat.sim.100 - sat.obs.100)^2, na.rm = TRUE))^0.5, 3)
  
  # for rsq
  rsq.table$swc10_daily[row.sam] <- round(as.numeric(cor(swc.obs.10, swc.sim.10, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$swc40_daily[row.sam] <- round(as.numeric(cor(swc.obs.40, swc.sim.40, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$swc100_daily[row.sam] <- round(as.numeric(cor(swc.obs.100, swc.sim.100, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$swc.vert_daily[row.sam] <- round(as.numeric(cor(swc.obs.vert, swc.sim.vert, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$sat10_daily[row.sam] <- round(as.numeric(cor(sat.obs.10, sat.sim.10, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$sat40_daily[row.sam] <- round(as.numeric(cor(sat.obs.40, sat.sim.40, method = "pearson", use = "complete.obs"))^2, 3)
  rsq.table$sat100_daily[row.sam] <- round(as.numeric(cor(sat.obs.100, sat.sim.100, method = "pearson", use = "complete.obs"))^2, 3)
  
  ## for nse
  nse.table$swc10_daily[row.sam] <- NSE(swc.obs.10, swc.sim.10, na.rm = TRUE)
  nse.table$swc40_daily[row.sam] <- NSE(swc.obs.40, swc.sim.40, na.rm = TRUE)
  nse.table$swc100_daily[row.sam] <- NSE(swc.obs.100, swc.sim.100, na.rm = TRUE)
  nse.table$swc.vert_daily[row.sam] <- NSE(swc.obs.vert, swc.sim.vert, na.rm = TRUE)
  nse.table$sat10_daily[row.sam] <- NSE(sat.obs.10, sat.sim.10, na.rm = TRUE)
  nse.table$sat40_daily[row.sam] <- NSE(sat.obs.40, sat.sim.40, na.rm = TRUE)
  nse.table$sat100_daily[row.sam] <- NSE(sat.obs.100, sat.sim.100, na.rm = TRUE)
}
summary(rsq.table); summary(nse.table)

write.csv(nse.table, file = file.path("results", current.folder, "nse.table.csv"), row.names = FALSE)
write.csv(rsq.table, file = file.path("results", current.folder, "rsq.table.csv"), row.names = FALSE)

params.nse <- par.on %>% full_join(nse.table, by = "par.sam")
write.csv(params.nse, file = file.path("results", current.folder, "params.nse.csv"), row.names = FALSE)
params.rsq <- par.on %>% full_join(rsq.table, by = "par.sam")
write.csv(params.rsq, file = file.path("results", current.folder, "params.rsq.csv"), row.names = FALSE)
best.params.nse <- params.nse[which(params.nse$par.sam %in% c(max.par.n.run.d, max.par.n.et.d, max.par.n.swc.d)),]
write.csv(best.params.nse, file = file.path("results", current.folder, "best.params.nse.csv"), row.names = FALSE)


rsq.max.vert.d <- max(rsq.table$swc.vert_daily, na.rm = TRUE)
rsq.max.10.d <- max(rsq.table$sat10_daily, na.rm = TRUE)
rsq.max.40.d <- max(rsq.table$sat40_daily, na.rm = TRUE)
rsq.max.100.d <- max(rsq.table$sat100_daily, na.rm = TRUE)
nse.max.vert.d <- max(nse.table$swc.vert_daily, na.rm = TRUE)
nse.max.10.d <- max(nse.table$sat10_daily, na.rm = TRUE)
nse.max.40.d <- max(nse.table$sat40_daily, na.rm = TRUE)
nse.max.100.d <- max(nse.table$sat100_daily, na.rm = TRUE)
max.par.n.swc.row.d <- c(which(nse.table$swc.vert_daily == nse.max.vert.d),
                         which(nse.table$sat10_daily == nse.max.10.d),  
                         which(nse.table$sat40_daily == nse.max.40.d), 
                         which(nse.table$sat100_daily == nse.max.100.d)) 
max.par.n.swc.d <- nse.table$par.sam[max.par.n.swc.row.d]

## rsq for the ensemble with max nse
rsq.for.max.nse.10.d <- rsq.table$swc10_daily[which(nse.table$sat10_daily == nse.max.10.d)]
rsq.for.max.nse.40.d <- rsq.table$sat40_daily[which(nse.table$sat40_daily == nse.max.40.d)]
rsq.for.max.nse.100.d <- rsq.table$sat100_daily[which(nse.table$sat100_daily == nse.max.100.d)]
rsq.for.max.nse.vert.d <- rsq.table$swc.vert_daily[which(nse.table$swc.vert_daily == nse.max.vert.d)]
# > nse.max.10.d
# [1] -0.2555529
# > nse.max.40.d
# [1] -2.271846
# > nse.max.100.d
# [1] -7.493096
# #nse.text <- as.character(as.expression(italic(NSE[10, MAX])~"="~nse.max.10))
# p1  <- ggplot(swc.d.long.sub, aes(x = date, y = value)) +
#   geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 1) +
#   geom_line(data = swc.d.long.sub %>% subset(par.sam %in% max.par.n.swc.d),
#             aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
#   geom_line(aes(x= date, y = obs, color = "Observed"), size = 1) +
#   # geom_text(x = 0.2, y = 0.8, label = nse.text) +
#   scale_colour_manual(name = "", values = col1) +
#   facet_grid(depth ~ .) +
#   ylab("Soil Water Content [m3/m3]") +
#   xlab("Date") + 
#   theme(axis.text = element_text(size = 14, face = "plain"),
#         legend.text = element_text(size = 16, face = "plain"),
#         legend.position = "top", legend.background = element_rect(fill = "transparent")) +
#   scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
#   ggtitle("Soil Water Content: Observed vs Simulated_Daily")
# p1
# ggsave(plot = file.path("figures", current.folder, "swc_Obs_vs_model_daily_with_extra_sm.jpeg"), height = 6.25, width = 8.94, units='in')
# #p1 %+% sat.d.long.sub
swc.d.long <- swc.d.long %>% subset(!is.na(depth))
swc.d.long.sub <- swc.d.long %>% subset(date >= as.Date("2012-01-01"))
swc.d.long.sub.mini <- subset(swc.d.long.sub, date >= as.Date("2016-01-01") & !depth %in% c(5, 300)) %>% droplevels()
swc.d.vert.long.sub.mini <- subset(swc.d.long.vert, date >= as.Date("2012-07-01") & date < as.Date("2019-01-01") ) 

sat.d.long <- sat.d.long %>% subset(!is.na(depth))
sat.d.long.sub <- sat.d.long %>% subset(date >= as.Date("2012-01-01"))
sat.d.long.sub.mini <- subset(sat.d.long.sub, date >= as.Date("2016-01-01") & date < as.Date("2018-01-01")& 
                                !depth %in% c(5, 300)) %>% droplevels()

col2 <- c("Observed" = "#f04546",  "Half HKsatADJ" = "blue", "Full BBslope" = "purple", "Half BBslope" = "green", "Half BBslope+ Half Vcmax" = "pink")

label.depths <- paste0(paste0("NSE.max = ", round(c(nse.max.10.d, nse.max.40.d, nse.max.100.d), 2)), 
                       "; ", paste0("R-squared = ", round(c(max(rsq.for.max.nse.10.d, na.rm = TRUE), max(rsq.for.max.nse.40.d, na.rm = TRUE), max(rsq.for.max.nse.100.d, na.rm = TRUE)), 2),
                                    "; ", paste0("R-squared.max = ", round(c(rsq.max.10.d, rsq.max.40.d, rsq.max.100.d), 2))))
dat_text <- data.frame(label = label.depths, depth = c(10, 40, 100))
g1 <- ggplot(swc.d.long.sub.mini, aes(x = date, y = value)) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 0.5) +
  geom_line(data = swc.d.long.sub.mini %>% subset(par.sam == 7001), aes(group = par.sam, color = "Half HKsatADJ"), show.legend = F, size = 0.5) +
  geom_line(data = swc.d.long.sub.mini %>% subset(par.sam == 7004), aes(group = par.sam, color = "Half BBslope+ Half Vcmax"), show.legend = F, size = 0.5) +
  geom_line(data = swc.d.long.sub.mini %>% subset(par.sam == 7003), aes(group = par.sam, color = "Full BBslope"), show.legend = F, size = 0.5) +
  geom_line(data = swc.d.long.sub.mini %>% subset(par.sam == 7002), aes(group = par.sam, color = "Half BBslope"), show.legend = F, size = 0.5) +
  theme(legend.position = "top") +
  scale_colour_manual(name = "", values = col2) + ylim(0.18, 0.57) +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y")) +
  ggtitle("Soil Water Content: Observed vs Ensemble simulations by depth")
p.h.swc <- g1 +  facet_grid(depth ~ .)
  # geom_text(data  = dat_text, mapping = aes(x = min(swc.d.long.sub.mini$date) + 100, 
            #                                 y = 0.57, label = label),
            # hjust   = -0.1, vjust = 2) +
  # geom_point(data = steph.sm, aes(x = date, y = swc, color = "Observed Plot-wide"), 
  #                                            show.legend = F, size = 0.5, alpha = 0.3, fill = "gray") +
  # theme(legend.position = "top") +
  # geom_line(data = swc.d.long.sub.mini %>% subset(par.sam %in% max.par.n.swc.d), 
  #           aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5)
ggsave("swc_Obs_vs_model_daily_horizontal.jpeg", plot = p.h.swc, path = file.path("figures", current.folder), device = "jpeg", height = 7, width = 9, units ='in')
## With vertical:
p.v.swc <- g1 %+% swc.d.vert.long.sub.mini + ylim(0.2, 0.52) +
  theme(legend.justification = c(1, 1), legend.position=c(0.97, 0.95)) +
  # geom_point(data = steph.sm, aes(x = date, y = swc, color = "Observed Plot-wide"), 
  #            show.legend = F, size = 0.5, alpha = 0.3, fill = "gray") +
  geom_text(aes(x = min(date) + 700, y = 0.51, label = 
                  as.character(paste0("NSE.max = ", round(nse.max.vert.d, 2), "; R-squared = ", round(max(rsq.for.max.nse.vert.d, na.rm = TRUE), 2), 
                                      "\n R-squared.max = ", round(rsq.max.vert.d, 2)))), vjust = 2) +
  geom_line(data = swc.d.vert.long.sub.mini %>% subset(par.sam %in% max.par.n.swc.d[1]), 
            aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5)
ggsave("swc_Obs_vs_model_daily_vertical.jpeg", plot = p.v.swc, path = file.path("figures", current.folder), device = "jpeg", height = 5, width = 9, units ='in')
    # Horizontal normalised
p.h.sat <- ggplot(sat.d.long.sub.mini, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = "Simulated"), show.legend = F, size = 0.2) +
  geom_line(aes(x= date, y = obs, color = "Observed"), size = 1) +
  scale_colour_manual(name = "", values = col1) + 
  ylab("Normalised Soil Water Content [0-1, unitless]") +
  xlab("Date") + 
  theme(axis.text = element_text(size = 14, face = "plain"),
        legend.text = element_text(size = 16, face = "plain"),
        legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "6 months", labels = function(x) format(x, "%b%y")) +
  facet_grid(depth ~ .) + 
  geom_text(data  = dat_text, mapping = aes(x = min(sat.d.long.sub.mini$date) + 100, 
                                            y = 1.2, label = label),
            hjust   = -0.1, vjust = 2) +
  theme(legend.position = "top") +
  geom_line(data = sat.d.long.sub.mini %>% subset(par.sam %in% max.par.n.swc.d), 
            aes(group = par.sam, color = "Best-fit Simulated"), show.legend = F, size = 0.5) +
  ggtitle("Normalised Soil Water Content: Observed vs Ensemble simulations by depth")
ggsave("sat_Obs_vs_model_daily_horizontal.jpeg", plot = p.h.sat, path = file.path("figures", current.folder), device = "jpeg", height = 7, width = 9, units ='in')

rm(g1)
swc.d.long.all <- obs.swc.d %>% select(-date) %>%
  full_join(mod.swc.d.all %>% select(-depth), by = "date.depth") %>% 
  as.data.frame() 
head(swc.d.long.all); summary(swc.d.long.all)

## do not remove swc.d.long

write.csv(swc.d.long, file = file.path("results", current.folder, "swc.d.long.csv"), row.names = FALSE)
write.csv(swc.d.vert.long, file = file.path("results", current.folder, "swc.d.vert.long.csv"), row.names = FALSE)

rm(mod.swc.d, mod.swc.d.1, mod.swc.d.10, mod.swc.d.100, mod.swc.d.3, mod.swc.d.300, mod.swc.d.40, mod.swc.d.5, mod.swc.d.500, mod.swc.d.6, mod.swc.d.800, mod.swc.d.all)
rm(mod.sat.d)
rm(g1)
rm(swc.d.10, swc.d.100, swc.d.100, swc.d.40, swc.d.long.sub, swc.d.long.sub.mini, swc.d.shallow.sub)
rm(swc.d.vert, swc.d.vert.long.sub.mini, sat.d.long, sat.d.long.sub)
rm(obs.sat.d, obs.sat.daily.spread, obs.tdr.daily, obs.tdr.daily.spread, swc.vertical, swc.vertical.mean, obs.swc.d)
rm(bci.tower)


###------------------
## Plotting nse table
###------------------

# nse.long <- nse.table %>% gather(key = "variable", value = "nse", -par.sam)
# ggplot(nse.long, aes(x = soil_depth, y = nse)) +
#   geom_bar(aes(fill = variable), stat = "identity", show.legend = F) +
#   # scale_fill_gradient(guide = guide_fill(reverse = TRUE)) +
#   facet_grid(variable ~ ., scales = "free_y") +
#   ylab("nse") +
#   xlab("Soil Depth [m]") + 
#   theme(axis.text.x = element_text(size = 14, face = "plain")) +
#   ggtitle("Nash-Sutcliffe Efficiency between Observation and Model By Soil Depth")
# ggsave(plot = file.path("figures", current.folder, "nse between Observation and Model By Soil Depth.jpeg"), height = 6.25, width = 8.94, units='in')

qrunoff.m <- read.csv(file = file.path("results", current.folder, "qrunoff.m.csv"), header = TRUE)
qet.m <- read.csv(file = file.path("results", current.folder, "qet.m.csv"), header = TRUE)
swc.d.long  <- read.csv(file = file.path("results", current.folder, "swc.d.long.csv"), header = TRUE)
swc.d.long.vert <- read.csv(file = file.path("results", current.folder, "swc.d.vert.long.csv"))

qrunoff.m.long <- qrunoff.m %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "par.sam", names_prefix = "X", 
               values_to = "value") %>% mutate_at("par.sam", as.integer) %>% as.data.frame()
qet.m.long <- qet.m %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "par.sam", names_prefix = "X", 
               values_to = "value") %>% mutate_at("par.sam", as.integer) %>% as.data.frame()

## Parameters sensitivity to obervations
par.on.long <- par.on %>% gather(key = "variable", value = "var.value", -par.sam) 

yrmo.vec <- c("2015-02", "2015-07", "2015-10", "2016-04", "2016-07")
yrmo.on <- yrmo.vec[5]

# par.gpp.m.long <- par.on.long %>% full_join(gpp.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
# par.gpp.m.long.sub <- par.gpp.m.long %>% subset(yrmo == yrmo.on)
# # 
# sen.gpp <- ggplot(par.gpp.m.long.sub, aes(x = var.value, y = value)) +
#   geom_point(aes(colour = variable), show.legend = F) +
#   geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
#   #geom_smooth(method = "auto") +
#   facet_wrap(~ variable, scales = "free_x") +
#   ylab("GPP [gC/m^2/d]") +
#   xlab("Variable") + 
#   theme(
#         strip.text = element_text(size = 14, face = "plain")) +
#   ggtitle(paste0("Sensitivity of monthly GPP to parameters for ", yrmo.on))
# ggsave(paste0("Sensitivity of monthly GPP to parameters_",
#               yrmo.on, ".jpeg"), plot = sen.gpp, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

par.qet.m.long <- par.on.long %>% full_join(qet.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
par.qet.m.long.sub <- par.qet.m.long %>% subset(yrmo == yrmo.on)
sen.et <- ggplot(par.qet.m.long.sub, aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
  #geom_smooth(method = "auto") +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Evapotranspiration [mm/month]") +
  xlab("Variable") + 
  theme(
        strip.text = element_text(size = 14, face = "plain")) +
  ggtitle(paste0("Sensitivity of monthly ET to parameters for ", yrmo.on))
ggsave(paste0("Sensitivity of monthly ET to parameters_",
              yrmo.on, ".jpeg"), plot = sen.et, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

par.qrunoff.m.long <- par.on.long %>% full_join(qrunoff.m.long %>% mutate(par.sam = as.integer(par.sam)), by = "par.sam")
par.qrunoff.m.long.sub <- par.qrunoff.m.long %>% subset(yrmo == yrmo.on)
sen.qrun <- ggplot(par.qrunoff.m.long.sub, aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3), se = TRUE) +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("QRUNOFF [mm/month]") +
  xlab("Variable") + 
  theme(strip.text = element_text(size = 14, face = "plain")) +
  ggtitle(paste0("Sensitivity of monthly runoff to parameters for ", yrmo.on))
ggsave(paste0("Sensitivity of monthly runoff to parameters_",
              yrmo.on, ".jpeg"), plot = sen.qrun, path = file.path("figures", current.folder), device = "jpeg", height = 6.25, width = 8.94, units='in')

par.swc.d.long.sub <- par.on.long %>% full_join(swc.d.long %>% 
                                                  mutate(date = as.Date(date),
                                                         par.sam = as.integer(par.sam)), by = "par.sam")

date.vec <- as.Date(c("2016-03-01", "2016-04-01", "2016-07-01", "2016-10-01"))
date.on <- date.vec[1] 
par.swc.d.long.sub.date.100 <- par.swc.d.long.sub %>% 
  subset(depth == 100 & date == date.on & !is.na(variable)) %>%
  mutate(var.value = as.numeric(var.value))
sen.swc.100 <- ggplot(par.swc.d.long.sub.date.100, 
                      aes(x = var.value, y = value)) +
  geom_point(aes(colour = variable), show.legend = F) +
  geom_smooth(method = glm, formula = y ~ splines::bs(x, 3)) +
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Soil Water Content [m3/m3]") +
  xlab("Variable") + 
  theme(
        strip.text = element_text(size = 14, face = "plain")) 
p.sen.swc.100 <- sen.swc.100 +
  ggtitle(paste0("Sensitivity of daily soil water content at 1 m to parameters for ", date.on))
ggsave(paste0("Sensitivity of daily soil water content at 1 m to parameters_",
              date.on, ".jpeg"), plot = p.sen.swc.100, path = file.path("figures", current.folder), height = 6.25, width = 8.94, units='in')
par.swc.d.long.sub.date.10 <- par.swc.d.long.sub %>%
  subset(depth == 10 & date == date.on & !is.na(variable)) %>% 
  mutate(var.value = as.numeric(var.value))
p.sen.swc.10 <- sen.swc.100 %+% par.swc.d.long.sub.date.10 + 
  ggtitle(paste0("Sensitivity of daily soil water content at 10 cm to parameters for ", date.on))
ggsave(paste0("Sensitivity of daily soil water content at 10 cm to parameters_",
              date.on, ".jpeg"), plot = p.sen.swc.10, path = file.path("figures", current.folder), height = 6.25, width = 8.94, units='in')


max.nse <- round(apply(nse.table, 2, max, na.rm = T), 2); max.nse
write.csv(max.nse, file = file.path("results", current.folder, "max.nse.csv"), row.names = FALSE)

# par.sam   qrunoff_daily qrunoff_monthly       qet_daily     qet_monthly       gpp_daily     gpp_monthly     swc10_daily     swc40_daily 
# 500.00            0.30            0.43           -0.95           -0.50           -1.54           -1.47            0.37            0.85 
# swc100_daily  swc.vert_daily 
# -0.21            0.69 
max.rsq <- round(apply(rsq.table, 2, max, na.rm = T), 2); max.rsq
write.csv(max.rsq, file = file.path("results", current.folder, "max.rsq.csv"), row.names = FALSE)
# par.sam   qrunoff_daily qrunoff_monthly       qet_daily     qet_monthly       gpp_daily     gpp_monthly     swc10_daily     swc40_daily 
# 500.00            0.79            0.80            0.35            0.78            0.14            0.39            0.80            0.86 
# swc100_daily  swc.vert_daily 
# 0.75            0.80


round(best.params.nse, 2)
# fates_leaf_BB_slope fates_leaf_slatop fates_leaf_vcmax25top fates_roota_par fates_rootb_par fates_smpsc fates_allom_l2fr FMAX aveDTB HKSAT_ADJ
# 138               12.57              0.02                 70.73            6.64            0.82    150607.4             0.51 0.28   5.95      7.53
# 142               11.84              0.00                 24.52            7.24            1.27    233706.3             0.72 0.28   5.95      7.53
# 218               12.57              0.02                 70.73            6.64            0.82    150607.4             0.51 0.36  12.97      3.48
# 400               11.84              0.00                 24.52            7.24            1.27    233706.3             0.72 0.77   6.70      3.27
# 403                9.43              0.01                 39.76            6.12            0.69    126788.4             0.17 0.77   6.70      3.27
# 418               13.24              0.02                 55.45            7.04            0.58    157748.8             0.45 0.77   6.70      3.27
# HKSAT_12.5 HKSAT_60 fpi_max par.sam qrunoff_daily qrunoff_monthly qet_daily qet_monthly gpp_daily gpp_monthly swc10_daily swc40_daily swc100_daily
# 138       0.01        0    0.06     167          0.24            0.41     -1.21       -0.80     -1.67       -1.66        0.30        0.85        -0.96
# 142       0.01        0    0.44     172          0.26            0.42     -0.95       -0.50     -1.59       -1.53        0.33        0.85        -1.01
# 218       0.01        0    0.06     267          0.26            0.41     -1.25       -0.96     -1.68       -1.70        0.29        0.83        -0.21
# 400       0.01        0    0.44     472          0.28            0.43     -0.98       -0.66     -1.54       -1.47        0.37        0.83        -0.37
# 403       0.01        0    0.05     475          0.28            0.42     -1.30       -1.09     -1.62       -1.59        0.34        0.83        -0.33
# 418       0.01        0    0.05     490          0.28            0.42     -1.29       -1.08     -1.62       -1.58        0.35        0.83        -0.33
# swc.vert_daily
# 138           0.63
# 142           0.66
# 218           0.62
# 400           0.69
# 403           0.67
# 418           0.67
