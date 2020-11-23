##---------------------------
# Choosing Best-fit ELM-FATES model params for hydrological water-balance
# Author: Rutuja Chitra-Tarak
# First written: Oct 3, 2019
##---------------------------

rm(list=ls())
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, easyNCDF, lubridate, tidyverse, data.table, hydroGOF)
theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

###***************************
#### Load Observed ET -------
###***************************
load(file.path("data-raw/bci.hydromet/forcings.rda"))

###***************************
current.folder <- "2019-10-14_5000"
params.rmse <- read.csv(file = file.path("results", current.folder, "params.rmse.csv"), header = TRUE)
params.rsq <- read.csv(file = file.path("results", current.folder, "params.rsq.csv"), header = TRUE)

## see Root/individual-flux-best-fits.html to lead selection:
col.rmse.1 <-  c("qrunoff_daily", "qrunoff_monthly", "qet_yearly", "qet_monthly",  
                 "swc.vert_daily", "swcplot10_daily", "swcplot40_daily", "swcplot100_daily")
col.rmse.2 <-  c("qrunoff_monthly", "qet_monthly",  
               "swc.vert_daily", "swcplot10_daily", "swcplot40_daily", "swcplot100_daily")
col.rsq <-  c("sat10_daily", "sat40_daily", "sat100_daily")
col.select.1 <- c(col.rmse.1, col.rsq)
col.select.2 <- c(col.rmse.2, col.rsq)
## using normalised rmse for qet
### swc normalised from 0 - 1
norm.index <- function(x, ...){(x - min(x, ...))/c(max(x, ...) - min(x, ...))}
comple.norm.rmse <- params.rmse[, c("par.sam", col.rmse.1)] %>%
  mutate_at(col.rmse.1, norm.index, na.rm = TRUE) %>% 
  ## since best rmse is minimum, while rsq is maximum, taking rmse's complement, so highest of that would be best
  mutate_at(col.rmse.1, function(x) {1 - x}) %>% as.data.frame() 
# creating new objective function table
params.obj <- params.rsq %>%  
  select(-col.rmse.1, -qrunoff_yearly, -gpp_daily, -gpp_yearly, -qet_daily, -gpp_monthly, -swc10_daily, -swc40_daily, -swc100_daily) %>%
  left_join(comple.norm.rmse, by = "par.sam")

params.obj.wt.1 <- params.obj %>% mutate_at(col.select.1, function(x) {x/length(col.select.1)})
params.obj$obj.all.1 <- as.numeric(apply(params.obj.wt.1[, col.select.1], 1, sum, na.rm = TRUE))
params.obj <- params.obj %>% arrange(desc(obj.all.1))
head(params.obj)
# par.sam sat10_daily sat40_daily sat100_daily qrunoff_daily qrunoff_monthly qet_yearly qet_monthly
# 1    3473       0.837       0.879        0.762     0.9371069       0.8055378  0.9582425   0.9492230
# 2    2800       0.831       0.883        0.791     0.8574423       0.8698784  0.9663997   0.9447287
# swc.vert_daily swcplot10_daily swcplot40_daily swcplot100_daily   obj.all obj.all.2 obj.all.1
# 1      0.8461538       0.9876543       0.9523810           0.7000 0.8740272 0.8576611 0.8740272
# 2      0.8205128       0.9753086       0.9166667           0.6625 0.8653125 0.8549550 0.8653125

## but it's redundant to use qrunoff_daily as well as monthly, and qet_yearly as well as monthly.
# prioritising monthly
params.obj.wt.2 <- params.obj %>% mutate_at(col.select.2, function(x) {x/length(col.select.2)})
params.obj$obj.all.2 <- as.numeric(apply(params.obj.wt.2[, col.select.2], 1, sum, na.rm = TRUE))
params.obj <- params.obj %>% arrange(desc(obj.all.2))
head(params.obj)
# par.sam sat10_daily sat40_daily sat100_daily qrunoff_daily qrunoff_monthly qet_yearly qet_monthly
# 1    3130       0.844       0.889        0.798     0.6415094       0.8634763  0.8757850   0.8859670
# 2    2217       0.846       0.889        0.784     0.6624738       0.8895647  0.7701102   0.7702128
# swc.vert_daily swcplot10_daily swcplot40_daily swcplot100_daily   obj.all obj.all.2 obj.all.1
# 1      0.8205128       0.9753086       0.9404762           0.7250 0.8417305 0.8601934 0.8417305
# 2      0.8974359       0.9135802       1.0000000           0.7500 0.8338525 0.8599771 0.8338525
# obj.all.2 finds best fits with better match for swcplot100_daily & sat100_daily
## number of top params to select
# top.few <- as.integer(round(nrow(params.nse)*0.001, 0))
# params.top.all <- setNames(data.frame(matrix(ncol = length(col.select), nrow = top.few)), col.select)
# for (i in 1: length(col.select)){
#   chosen.col <- col.select[i]
#   params.nse <- params.nse[order(params.nse[, chosen.col], decreasing = TRUE),]
#   params.top.all[, chosen.col] <- params.nse$par.sam[1: top.few]
# }
# ## unique params across all
# # params.common <- intersect(as.vector(as.matrix(params.top.all)))
# params.top.few <- sort(unique(as.vector(as.matrix(params.top.all)))) #unique(params.top.all$qrunoff_monthly, params.top.all$qrunoff_daily)
top.few <- 100
params.obj <- params.obj %>%  
  rename(obj.func = obj.all.2) %>% 
  select(-obj.all.1) %>% 
  arrange(desc(obj.func))

params.top.few <- params.obj$par.sam[1:100]
length(params.top.few)
params.obj.top.few <- params.obj %>%  
  subset(par.sam %in% params.top.few)
summary(params.obj.top.few)
write.csv(params.obj, file = file.path("results", current.folder, paste0("params.obj.csv")), row.names = FALSE)
write.csv(params.top.few, file = file.path("results", current.folder, paste0("params.top.few_", top.few, ".csv")), row.names = FALSE)
write.csv(params.obj.top.few, file = file.path("results", current.folder, paste0("params.obj.top.few_", top.few, ".csv")), row.names = FALSE)
params.obj <- read.csv(file.path("results", current.folder, paste0("params.obj.csv")), header = TRUE)

params.obj.top.few <- read.csv(file.path("results", current.folder, paste0("params.obj.top.few_", top.few, ".csv")), header = TRUE)
## the above selection does not pass through even one of the plot-wide absolute soil moisture observations at 100 cm Depth 
## there are about 300 par.sam above 0.8 swcplot100_daily. So conditionally selecting on that:
params.obj.cond <- params.obj %>% arrange(desc(obj.func)) %>% subset(swcplot100_daily > 0.8)
params.top.few.cond <- params.obj.cond$par.sam[1:100]
length(params.top.few.cond)
params.obj.top.few.cond <- params.obj.cond %>%  
  subset(par.sam %in% params.top.few.cond)
write.csv(params.top.few.cond, file = file.path("results", current.folder, paste0("params.top.few.cond_", top.few, ".csv")), row.names = FALSE)
write.csv(params.obj.top.few.cond, file = file.path("results", current.folder, paste0("params.obj.top.few.cond_", top.few, ".csv")), row.names = FALSE)

equi <- params.obj %>% 
  select(-col.select.1) %>%
  subset(par.sam %in% params.obj$par.sam[1:1000]) %>%
  pivot_longer(
    cols = c(-par.sam, -obj.func),
    names_to = "parameter",
    values_to = "value",
    values_drop_na = TRUE
    ) %>%
  mutate(top.few = ifelse(par.sam %in% params.obj$par.sam[1:20], 'Top Few', 'The Rest'),
         top.few = fct_relevel(top.few, 'Top Few', 'The Rest'))
# p_load(wesanderson, ggsci)
# pal <- wes_palette("Zissou1", 100, type = "continuous")
# scale_color_gradientn(colours = pal) + 

equi.plot <- ggplot(equi, aes(x = value, y = obj.func)) +
  geom_point(aes(colour = top.few), size = 0.3) +
  scale_colour_manual(name = "", values = c("red", "grey")) +
  theme(legend.position = c(0.8, 0.1)) +
  facet_wrap(~ parameter, scales = "free_x") +
  ylab("Objective Function [0-1]") +
  xlab("Parameter") + 
  theme(strip.text = element_text(size = 14, face = "plain")) +
  ggtitle(paste0("Parameter Equifinality"))
ggsave(paste0("Parameter Equifinality.jpeg"), plot = equi.plot, path = file.path("figures", current.folder, "Sensitivity"), device = "jpeg", height = 6.25, width = 8.94, units='in')

# tot.nse <- colSums(params.nse, na.rm = TRUE)
# likelihood <- params.nse %>% 
#   mutate(qrunoff_monthly = ifelse(par.sam %in% params.top.few, qrunoff_monthly/as.numeric(tot.nse["qrunoff_monthly"]), 0))
# summary(likelihood$qrunoff_monthly[likelihood$par.sam %in% params.top.few])

nc <- nc_open(file.path("data-raw/extract", current.folder, "BCI.ICLM45ED.badger.intel.Cabfd064-Fc47f968.1050.clm2.h0.2018-01.nc"), 
               readunlim = FALSE)
all.vars <- NcReadVarNames(nc)
var.info <- setNames(data.frame(matrix(ncol = 3, nrow = length(all.vars))), 
                     c("var", "longname", "units"))

for (i in 1: length(all.vars)){
  var.extract <- nc$var[[all.vars[i]]]
  if(!is.null(var.extract)){
    var.info$var[i] <- all.vars[i]
    var.info$longname[i] <- var.extract$longname
    var.info$units[i] <- var.extract$units
  }
}
View(var.info)
write.csv(var.info, file = file.path("results", current.folder, "sample_case_var.info_for_monthly_file.csv"), row.names = FALSE)

# ERRH2O total water conservation error mm
# H2OCAN "intercepted water" units mm
# WA "water in the unconfined aquifer (vegetated landunits only)" # default value 0, so not tracked?
# ZWT_PERCH "perched water table depth (vegetated landunits only)" m | but values extra-ordinarily large
# ZWT "water table depth (vegetated landunits only)" m | but values extra-ordinarily large
# H2OSFC "surface water depth" mm: # default value negligible 0.1 - 0001 mm/month
# TWS "total water storage" mm
# TWS_MONTH_END : "total water storage at the end of a month"
# RAIN atmospheric rain mm/s
# QINTR interception mm/s
# QSOIL Ground evaporation (soil/snow evaporation + soil/snow sublimation - dew) mm/s
# QVEGE canopy evaporation mm/s
# QVEGT canopy transpiration mm/s
# QCHARGE: aquifer recharge rate (vegetated landunits only) mm/s
# QDRAI: sub-surface drainage mm/s
# WT total water storage (unsaturated soil water + groundwater, veg landunits) mm
# (NOt Present; TWS looks equivalent)
# QRGWL: surface runoff at glaciers (liquid only), wetlands, lakes mm/s

# QRUNOFF total liquid runoff (does not include QSNWCPICE) mm/s
## So probably includes QOVER + QDRAI ## confirmed from QOVER
# QOVER: surface runoff mm/s

#delta TWS = RAIN - QINT -QDRIP - QSOIL - QVEGE - QVEGT - QRUNOFF - QCHARGE - QDRAI
wb <- list()
wb.vec <- c("TWS", "RAIN", "QINTR", "QDRIP", "TWS_MONTH_END", "QSOIL", #"ERRH2O", "WA", "H2OSFC", 
"QVEGE", "QVEGT", "QRUNOFF", "QCHARGE", "QDRAI", "ZWT", "ZWT_PERCH")

for(i in 1:length(wb.vec)){
  load(file.path("data-raw/extract", current.folder, "extract", paste0(wb.vec[i], ".h0.extract.Rdata")))
  wb[[i]] <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "yrmo") %>%
    as.data.frame()  
  if(!wb.vec[i] %in% c("ERRH2O", "WA", "H2OSFC", "TWS", "TWS_MONTH_END")){
      wb[[i]][, -1] <- wb[[i]][, -1]*24*60*60*30 # converting mm/s to mm/month
  }
  wb[[i]] <- wb[[i]] %>% 
    gather(key = "par.sam", "value", -yrmo) %>% 
    subset(par.sam %in% params.top.few) %>%
    mutate(date = as.Date(paste0(yrmo,"-01")),
           bioyear = as.numeric(if_else(as.numeric(format(date, format = "%m")) < 7, 
                             # if yes, call it the previous year; 
                             #using lubridate %m+% Add and subtract months to a date without exceeding the last day of the new month
                             format(date %m+% years(-1), "%Y"), format(date, "%Y")))
    ) %>% # removing first three bioyears allowing for model initialisation
    subset(!bioyear %in% seq(from = min(bioyear, na.rm = TRUE), length.out = 3))
}
names(wb) <- wb.vec

### by month
wb.table.empty.mo <- 
  setNames(data.frame(matrix(ncol = 1 + length(wb.vec), 
                             nrow = length(unique(wb[["QRUNOFF"]]$yrmo)))), 
           c("yrmo", wb.vec)) %>%
  mutate(yrmo = unique(wb[["QRUNOFF"]]$yrmo))

## wb.table by par.sam by month
wb.par.sam.mo <- list()
par.sam <- as.character(params.top.few)[!is.na(params.top.few)] # unique(wb.yr[["RAIN"]]$par.sam) #same as params.top.few
for(i in 1:length(par.sam)) {
  par.sam.i <- par.sam[i]
  wb.par.sam.mo[[par.sam.i]] <- wb.table.empty.mo
  for(j in 1: length(wb.vec)){
    web.vec.j <- wb.vec[j]
    var.yrmo <- wb[[web.vec.j]] %>% 
      subset(par.sam == par.sam.i)
    wb.par.sam.mo[[par.sam.i]][, web.vec.j] <- var.yrmo$value
  }
  wb.par.sam.mo[[par.sam.i]] <- wb.par.sam.mo[[par.sam.i]] %>% 
    mutate(QET = QSOIL + QVEGE + QVEGT,
           # Source = RAIN - QINTR + QDRIP,
           QOVER = QRUNOFF - QDRAI,
           Source = QDRIP, ## as QDRIP is given as throughfall
           Sink = QET + QRUNOFF, #+ QCHARGE,
           W.Balance = Source - Sink,
           del.TWS = TWS - c(NA, TWS[-nrow(wb.table.empty.mo)]),
           del.TWS_MONTH_END = c(NA, TWS_MONTH_END[-nrow(wb.table.empty.mo)]),
           del.ZWT = ZWT - c(NA, ZWT[-nrow(wb.table.empty.mo)]),
           del.ZWT_PERCH = c(NA, ZWT_PERCH[-1]- ZWT_PERCH[-nrow(wb.table.empty.mo)])
    )
}
View(head(wb.par.sam.mo[["3130"]]))


## wb table by year:
head(wb[["QRUNOFF"]])
wb.yr <- list()
value.lastmo <- function(vCol, dCol) {vCol[dCol == max(dCol, na.rm = TRUE)]}
value.fistmo <- function(vCol, dCol) {vCol[dCol == min(dCol, na.rm = TRUE)]}
for(i in 1:length(wb.vec)){
  if(wb.vec[i] %in% c("WA", "TWS", "TWS_MONTH_END", "ZWT", "ZWT_PERCH")){
    wb.yr[[i]] <- wb[[i]] %>% group_by(bioyear, par.sam) %>%
      summarise(value = value.lastmo(value, date)) # Storage at the end of year
  } else {
    wb.yr[[i]] <- wb[[i]] %>% group_by(bioyear, par.sam) %>%
      summarise(value = sum(value, na.rm = TRUE)) # Sum for the year
  }
}
names(wb.yr) <- wb.vec

wb.table.empty <- 
  setNames(data.frame(matrix(ncol = 1 + length(wb.vec), 
                             nrow = length(unique(wb.yr[["QRUNOFF"]]$bioyear)))), 
           c("bioyear", wb.vec)) %>%
  mutate(bioyear = unique(wb.yr[["QRUNOFF"]]$bioyear))

## wb.table by par.sam by year
wb.par.sam <- list()
par.sam <- as.character(params.top.few) # unique(wb.yr[["RAIN"]]$par.sam) #same as params.top.few
for(i in 1:length(par.sam)) {
  par.sam.i <- par.sam[i]
  wb.par.sam[[par.sam.i]] <- wb.table.empty
  for(j in 1: length(wb.vec)){
    web.vec.j <- wb.vec[j]
    var.yr <- wb.yr[[web.vec.j]] %>% 
                 subset(par.sam == par.sam.i)
    wb.par.sam[[par.sam.i]][, web.vec.j] <- var.yr$value
  }
  wb.par.sam[[par.sam.i]] <- wb.par.sam[[par.sam.i]] %>% 
    mutate(QET = QSOIL + QVEGE + QVEGT,
           # Source = RAIN - QINTR + QDRIP,
           QOVER = QRUNOFF - QDRAI,
           Source = QDRIP, ## as QDRIP is given as throughfall
           Sink = QET + QRUNOFF, #+ QCHARGE,
           W.Balance = Source - Sink,
           del.TWS = TWS - c(NA, TWS[-nrow(wb.table.empty)]),
           del.TWS_MONTH_END = c(NA, TWS_MONTH_END[-nrow(wb.table.empty)]),
           del.ZWT = ZWT - c(NA, ZWT[-nrow(wb.table.empty)]),
           del.ZWT_PERCH = c(NA, ZWT_PERCH[-1]- ZWT_PERCH[-nrow(wb.table.empty)])
    )
}
View(wb.par.sam[["3130"]])
## Adding Obs ET & Runoff
## Observed ET from flux tower
# AET.flag.day has NAs substituted where insufficient (< 50%) actual data for the day
# AET.flag.month has NAs substituted where insufficient (< 50%) actual data for the month

obs.df <- forcings %>% select(date, Precip, AET, AET.flag.day, AET.flag.month, flow_conrad) %>%
  rename(obs.ET = AET, obs.runoff = flow_conrad) %>% # in mm/day
  # check if month is before July
  mutate(date = as.Date(date),
         yrmo = format(date, format = "%Y-%m"),
         bioyear = as.numeric(if_else(as.numeric(format(date, format = "%m")) < 7, 
                           # if yes, call it the previous year; 
                           # using lubridate %m+% Add and subtract months to a date without exceeding the last day of the new month
                           format(date %m+% years(-1), "%Y"), format(date, "%Y")))) %>% 
  ## selecting data when obs.ET present, but note that QRUNOFF is available earlier
  subset(date > as.Date("2001-01-01"), na.rm = TRUE)

obs.yr <- obs.df %>% select(-date, -Precip) %>% group_by(bioyear) %>%
  ## selecting data when obs.ET present, but note that QRUNOFF is available earlier
  #subset(date > min(date[!is.na(obs.ET)], na.rm = TRUE)) %>%
  summarise_at(vars(obs.ET, obs.runoff), sum, na.rm = TRUE) %>%
  subset(bioyear %in% unique(wb.yr[["QRUNOFF"]]$bioyear)) #%>%
  # mutate(obs.ET = replace(obs.ET, bioyear %in% c(2007:2011), NA),
  #        AET.flag.day = replace(AET.flag.day, bioyear %in% c(2007:2011), NA), 
  #        AET.flag.month = replace(AET.flag.month, bioyear %in% c(2007:2011), NA))
obs.yrmo <- obs.df %>% select(-date) %>% group_by(yrmo) %>%
  summarise_at(vars(obs.ET, AET.flag.month, obs.runoff), sum) %>%
  subset(yrmo %in% unique(wb[["QRUNOFF"]]$yrmo))
## adding observations to wb.table by par.sam
# don't repeat this 
for(i in 1:length(par.sam)) {
  par.sam.i <- par.sam[i]
  wb.par.sam[[par.sam.i]] <- wb.par.sam[[par.sam.i]] %>%
    left_join(obs.yr, by = "bioyear") %>%
    mutate(ET.sim.obs = QET/obs.ET,
           runoff.sim.obs = QRUNOFF/obs.runoff) %>%
    subset(bioyear >= 2012) 
  wb.par.sam.mo[[par.sam.i]] <- wb.par.sam.mo[[par.sam.i]] %>%
    left_join(obs.yrmo, by = "yrmo") %>%
    mutate(ET.sim.obs = QET/obs.ET,
           runoff.sim.obs = QRUNOFF/obs.runoff) 
}

wb.table.all <- do.call(rbind, wb.par.sam) %>% group_by(bioyear) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  mutate(QDRAI.QOVER = round(QDRAI/QOVER, 0),
         QDRAI.obs.runoff = round(QDRAI/obs.runoff, 0)) %>%
  mutate_at(vars(-c(bioyear, ET.sim.obs, runoff.sim.obs, QDRAI.QOVER, QDRAI.obs.runoff)), list(~round(., 0))) %>%
  mutate_at(vars(ET.sim.obs, runoff.sim.obs, QDRAI.QOVER, QDRAI.obs.runoff), list(~round(., 1))) %>%
  mutate_at(vars(ET.sim.obs), ~replace(., bioyear %in% c(2017, 2018), NA))  %>%# since currently I have not calculated matching model flux totals when obs present
  subset(bioyear != 2018)
wb.table <- wb.table.all %>% select(-c(QSOIL, QVEGE, QVEGT, ZWT, ZWT_PERCH, del.ZWT, del.ZWT_PERCH))

View(wb.table)
mean.wb.table <- apply(wb.table, 2, FUN = function(x) {as.numeric(mean(x, na.rm = TRUE))}) %>%
  t() %>%
  as.data.frame() %>% 
  mutate_at(vars(-c(bioyear, ET.sim.obs, runoff.sim.obs, QDRAI.QOVER, QDRAI.obs.runoff)), list(~round(., 0))) %>%
  mutate_at(vars(ET.sim.obs, runoff.sim.obs, QDRAI.QOVER, QDRAI.obs.runoff), list(~round(., 1)))
  
wb.table <- wb.table %>% bind_rows(mean.wb.table)
wb.table$bioyear[nrow(wb.table)] <- "Average"
write.csv(wb.table, file = file.path("results", current.folder, "wb.table.csv"), row.names = FALSE)
write.csv(wb.table.all, file = file.path("results", current.folder, "wb.table.all.csv"), row.names = FALSE)

## for yrmo:
wb.table.yrmo.all <- do.call(rbind, wb.par.sam.mo) %>% group_by(yrmo) %>%
  summarize_all(mean, na.rm = TRUE) %>%
  mutate_at(vars(-c(yrmo, ET.sim.obs, runoff.sim.obs)), list(~round(., 0))) %>%
  mutate_at(vars(ET.sim.obs, runoff.sim.obs), list(~round(., 2))) %>%
  subset(!is.na(obs.ET)) #%>%
  # mutate(yr = as.numeric(format(as.Date(paste0(yrmo, "-1")), format = "%Y"))) %>%
  # subset(yr >= 2012) %>% select(-yr)
wb.table.yrmo <- wb.table.yrmo.all %>% select(-c(QSOIL, QVEGE, QVEGT, ZWT, ZWT_PERCH, del.ZWT, del.ZWT_PERCH, AET.flag.month))

# View(wb.table.yrmo)
# mean.wb.table.yrmo <- apply(wb.table.yrmo, 2, FUN = function(x) {as.numeric(mean(x, na.rm = TRUE))}) %>%
#   t() %>%
#   as.data.frame() %>% 
#   mutate_at(vars(-c(yrmo, ET.sim.obs, runoff.sim.obs)), list(~round(., 0))) %>%
#   mutate_at(vars(ET.sim.obs, runoff.sim.obs), list(~round(., 2)))

write.csv(wb.table.yrmo, file = file.path("results", current.folder, "wb.table.yrmo.csv"), row.names = FALSE)
write.csv(wb.table.yrmo.all, file = file.path("results", current.folder, "wb.table.yrmo.all.csv"), row.names = FALSE)
### HEATMAP for water balance data
wb.table.long <- gather(wb.table, key = variable, value = value, -bioyear)
p.wb.table <- ggplot(wb.table.long, aes(bioyear, variable, fill = value)) + 
  geom_tile() + #ylab("Depth [cm]") + xlab("Date") +   
  #scale_fill_continuous("SWC [% v/v]", direction = -1)
  scale_fill_viridis_c("", direction = -1, option = "plasma") +
  ggtitle("Water balance fluxes and source")

run <- wb.table.yrmo %>% select(yrmo, RAIN, QRUNOFF, QDRAI, QOVER, obs.runoff, runoff.sim.obs) %>%
  mutate(QDRAI.QOVER = round(QDRAI/QOVER, 1),
         QDRAI.obs.runoff = round(QDRAI/obs.runoff, 1)) %>%
  mutate(mo = as.numeric(format(as.Date(paste0(yrmo, "-1")), format = "%m"))) 
run.sub <- run %>% subset(mo %in% c(1:4)) %>% select(-mo)

mean.run.sub <- apply(run.sub %>% mutate(yrmo = NA), 2, FUN = function(x) {as.numeric(mean(x[!is.infinite(x)][!is.nan(x)], na.rm = TRUE))}) %>%
  t() %>%
  as.data.frame() %>%
  mutate_all(list(~round(., 0))) %>%
  mutate(yrmo = "Average")
run.sub <- run.sub %>% bind_rows(mean.run.sub)
write.csv(run.sub, file = file.path("results", current.folder, "run.table.dryssn.csv"), row.names = FALSE)

run.long <- run %>%
  gather(key = variable, value = value, -yrmo)
ggplot(run.long, aes(yrmo, variable, fill = value)) + 
  geom_tile() + #ylab("Depth [cm]") + xlab("Date") +   
  #scale_fill_continuous("SWC [% v/v]", direction = -1)
  scale_fill_viridis_c("", direction = -1, option = "plasma") +
  ggtitle("Runoff terms")

#*********-----------------------------
### HEATMAP for best-fit soil moisture
#*********-----------------------------

## also see Seaborn: https://towardsdatascience.com/python-seaborn-plots-in-r-using-reticulate-fb59cebf61a7
## Or plot_ly: https://plot.ly/r/heatmaps/ : ‘plot_ly’ is not available (for R version 3.6.1) 

nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
# depths
soil.depths <- as.numeric(nc$var[['H2OSOI']]$dim[[2]]$vals) # 15 depths in m
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934
## Plotting soil moisture
load(file.path("data-raw/extract", current.folder, "extract/H2OSOI.h1.extract.Rdata"))
swc.obs <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% gather(key = "par.sam", "value", -date, -depth)
for(i in 1: length(var.res.arr)){
    swc.obs.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
      mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
      gather(key = "par.sam", "value", -date, -depth)
    swc.obs <- swc.obs %>% bind_rows(swc.obs.i)
  }
rm(var.res.arr); rm(swc.obs.i)
mini.date <- seq(from = as.Date("2012-01-01"), to = as.Date("2018-12-31"), by = "day")
swc.d.long.sub <- swc.obs %>% subset(par.sam %in% params.top.few) %>%
  mutate(date = as.Date(date), depth = signif(round(depth, 2), 2)) %>% 
  subset(date %in% mini.date) %>%   # As depths below 3 m are absent from some par.sam and has 0s associated with them; 
  # substituting those 0s with NAs, as there are no 0s otherwise 
  mutate_at(vars(value), ~na_if(., 0)) %>% 
  subset(!is.na(value)) %>% droplevels()
swc.d.long.sub.mean <- swc.d.long.sub %>%
  group_by(date, depth) %>% 
  summarise_at(vars(value), mean, na.rm = TRUE) %>%
  arrange(depth) %>% subset(!is.na(value)) %>% droplevels()

# Heatmap 
swc.top.few <- ggplot(swc.d.long.sub.mean, aes(date, as.factor(-depth), fill = value*100)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Date") +   
  #scale_fill_continuous("SWC [% v/v]", direction = -1)
  scale_fill_viridis_c("SWC [% v/v]", direction = -1, option = "plasma") +
  ggtitle("Average SWC by depth across best-fit ensembles")
ggsave("swc_mean_across_params.top.few.jpeg", plot = swc.top.few, path = file.path("figures", current.folder), height = 8.94, width = 8.94, units='in')
rm(swc.top.few)
## by depth panels
theme_set(theme_bw())
theme_update(panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank()
)
p.mod.swc <- ggplot(swc.d.long.sub, aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = par.sam), show.legend = F, size = 0.2) +
  facet_grid(depth ~ .) +
  ylab("Soil Water Content [m3/m3]") + xlab("Date") + 
  scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Content for best-fit ensembles")
ggsave("swc_model_daily_all_depths_params.top.few.jpeg", plot = p.mod.swc, path = file.path("figures", current.folder), height = 12, width = 8.94, units='in')
rm(p.mod.swc)

###
#### SOILPSI: Obs versus model #-------
#### 
# Daily #------
####
# SOILPSI	soil water potential in each soil layer	MPa 
## var.info file says its MPa but is given in KPa since the range is -500 to -1500
load(file.path("data-raw/extract", current.folder, "extract/SOILPSI.h1.extract.Rdata"))
psi.obs <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date), depth = soil.depths[1]) %>% gather(key = "par.sam", "value", -date, -depth)
for(i in 1: length(var.res.arr)){
  psi.obs.i <- setDT(as.data.frame(t(var.res.arr[[i]])), keep.rownames = "date") %>%
    mutate(date = as.Date(date), depth = soil.depths[i]) %>% 
    gather(key = "par.sam", "value", -date, -depth)
  psi.obs <- psi.obs %>% bind_rows(psi.obs.i)
}
rm(psi.obs.i)
rm(var.res.arr) # large file
psi.d.long.sub <- psi.obs %>% subset(par.sam %in% params.top.few) %>%
  mutate(date = as.Date(date), depth = signif(round(depth, 2), 2)) %>% 
  subset(date %in% mini.date) %>%   # As depths below 3 m are absent from some par.sam and has 0s associated with them; 
  # substituting those 0s with NAs, as there are no 0s otherwise 
  mutate_at(vars(value), ~na_if(., -15)) %>% 
  subset(!is.na(value)) %>% droplevels()
write.csv(psi.d.long.sub, file = file.path("results", current.folder, "psi_model_daily_all_depths_params.top.few.csv"), row.names = FALSE)

psi.d.long.sub.mean <- psi.d.long.sub %>%
  group_by(date, depth) %>% summarise(value = mean(value, na.rm = TRUE)) %>%
  arrange(depth) %>% subset(!is.na(value)) %>% droplevels()

# Heatmap 
psi.top.few <- ggplot(psi.d.long.sub.mean, aes(date, as.factor(-depth), fill = value)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Date") +   
  #scale_fill_continuous("psi [% v/v]", direction = -1)
  scale_fill_viridis_c("PSI [MPa]", direction = -1, option = "plasma") +
  ggtitle("Average SOILPSI by depth across best-fit ensembles")
ggsave("psi_mean_across_params.top.few.jpeg", plot = psi.top.few, path = file.path("figures", current.folder), height = 8.94, width = 8.94, units='in')
rm(psi.top.few)
## by depth panels
## subsetting since too big to plot
require(scales);
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt", 
    transform = function(x) -sqrt(abs(x)), 
    inverse = function(x) x^2);
}
p.mod.psi <- ggplot(psi.d.long.sub, aes(x = date, y = -value)) +
  scale_y_continuous(trans="rev_sqrt", breaks = c(0, 0.5, 2, 5, 10, 15)) +
  geom_line(aes(group = par.sam, color = par.sam), show.legend = F, size = 0.2) +
  facet_grid(depth ~ .) +
  ylab("-Soil Water Potential [MPa]") + xlab("Date") + 
  scale_x_date(date_breaks = "3 months", labels = function(x) format(x, "%b%y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Potential for best-fit ensembles")
ggsave("psi_model_daily_all_depths_params.top.few.jpeg", plot = p.mod.psi, path = file.path("figures", current.folder), height = 12, width = 8.94, units='in')
rm(p.mod.psi)

###
#### BTRAN: Obs versus model #-------
####
# Daily #------
####
load(file.path("data-raw/extract", current.folder, "extract/BTRAN.h1.extract.Rdata"))

mod.btran.d <- setDT(as.data.frame(t(var.res.arr[[1]])), keep.rownames = "date") %>%
  mutate(date = as.Date(date)) %>%
  as.data.frame()

rm(var.res.arr) # large file

btran.d <- mod.btran.d
head(btran.d[, 1:6]); summary(btran.d[, 1:6])
btran.d.long <- gather(btran.d, key = "par.sam", "value", -date) %>% mutate(par.sam = as.numeric(par.sam))
btran.d.long.sub <- btran.d.long %>% subset(date %in% mini.date & par.sam %in% params.top.few)

b.btran <- ggplot(btran.d.long.sub,
                  aes(x = date, y = value)) +
  geom_violin(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.1) +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("BTRAN: Simulated Best-Fits")
ggsave("BTRAN_model_daily_bestfit_params.top.few.jpeg", plot = b.btran, file.path("figures", current.folder), device = "jpeg", height = 3, width = 8.94, units='in')

btran.d.long.sub.stat <- btran.d.long.sub %>% 
  group_by(date) %>%
  summarise(mean = mean(value, na.rm = TRUE),
         upper.CI = quantile(value, probs = 0.95),
         lower.CI = quantile(value, probs = 0.05))

b.btran.ci <- ggplot(btran.d.long.sub.stat,
                  aes(x = date, y = mean)) +
  geom_line(show.legend = F, size = 0.1) +
  geom_ribbon(aes(ymin = lower.CI, ymax = upper.CI), alpha=0.3) +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y")) +
  ggtitle("BTRAN: Simulated best-fits: mean with 95% Confidence Interval")
ggsave("BTRAN_model_daily_bestfit_params.top.few_CI.jpeg", plot = b.btran.ci, file.path("figures", current.folder), device = "jpeg", height = 3, width = 8.94, units='in')

## Plot Hydrograph
pdf(file.path("figures", current.folder, "Hydrograph.pdf"), width = 9, height = 5)
maxSF   <- max(c(obs.df$obs.runoff), na.rm = T)
maxPR   <- max(obs.df$Precip, na.rm = T)
par(mar = c(4, 4, 3, 4) + 0.1)
plot(obs.df$date, obs.df$obs.runoff,
     type = 'l', col = "red",
     ylim = c(0, 1.3 * maxSF),
     xaxs = "i", yaxs = "i",
     xlab = "Time", ylab = "Streamflow (mm/day)",
     main = "BCI Conrad Catchment Hydrograph", cex.main = 0.9)
# lines(obsDF0$date, obsDF0$streamflow, col = "black")
par(new = TRUE)
plot(x = obs.df$date, y = obs.df$Precip,
     type = "n", ylim = c(5 * maxPR, 0),
     xaxs = "i", yaxs = "i",
     axes = FALSE, xlab = "", ylab = "")
segments(x0 = obs.df$date, y0 = rep(0, nrow(obs.df)),
         x1 = obs.df$date, y1 = obs.df$Precip,
         lend = 2, lwd =1)

yrAxis  <- seq(0, ceiling(maxPR), length.out = 5)
axis(4, at = yrAxis, labels = paste0(yrAxis))
mtext("Precip (mm/day)", side = 4, line = 2, adj = 1)
dev.off()

