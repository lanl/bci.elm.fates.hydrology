## checking/cleaning climate data

## sourcing all .nc files
rm(list = ls())
graphics.off()
if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, tidyverse, bci.hydromet, data.table, chron)

path = "data-raw/ConvertMetCSVtoCLM-expand-format/NCOut/bci_0.1x0.1_met.v5.1/CLM1PT_data"

all.files = list.files(
  path = path,
  pattern = ".nc",
  recursive = FALSE,
  full.names = TRUE
)
file.names <- basename(all.files)
months <- sub(".nc", "", file.names)
total <- length(all.files)

time.df <- data.frame(filename = file.names, time.stamp = NA)
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1: length(file.names)){
  setTxtProgressBar(pb, i)
  filename <- all.files[i]
  nc <- nc_open(filename, write = F)
  time.df$time.stamp[i] <- ncatt_get(nc,"time","units")$value
}
View(time.df)
start.month = 1
## datetime column is created based on the start and end months of files available and filling up at hourly freq
# rnames <- seq(from = as.POSIXct(paste0(months[start.month], "-01 0:00"), tz = "America/Cayman"),
#               to = as.POSIXct(paste0(months[length(months)], "-31 23:00"), tz = "America/Cayman"), by = "1 hour")
# head(rnames)
#
# start.month = 384
ncin <- nc_open("data-raw/ConvertMetCSVtoCLM-expand-format/NCOut/bci_0.1x0.1_met.v5.1/CLM1PT_data/2016-01.nc", write = F)
val1 <- ncvar_get(ncin, "PRECTmms")
length(val1)
tunits1 <- ncatt_get(ncin,"time","units")
tunits1
# # get time
time1 <- ncvar_get(ncin,"time")
tail(time1)
nt1 <- dim(time1)
nt1

# convert time -- split the time units string into fields
tustr1 <- strsplit(tunits1$value, " ")
tdstr1 <- strsplit(unlist(tustr1)[3], "-")
tmonth1 <- as.integer(unlist(tdstr1)[2])
tday1 <- as.integer(unlist(tdstr1)[3])
tyear1 <- as.integer(unlist(tdstr1)[1])
chron(time1,origin=c(tmonth1, tday1, tyear1))
time.vec1 <- as.POSIXct(chron(time1, origin = c(tmonth1, tday1, tyear1)), tz = "GMT")
# the above takes in the given data as in GMT format, but converts it to Sys timezone
# so to convert back to GMT:
attr(time.vec1, "tzone") <- "GMT"
head(time.vec1)
tail(time.vec1)
## this converts to system time zone
#----

nc.nvars <- ncin$nvars
var.name <- vector()
for (j in 1: nc.nvars) {
  var.name[j] <- ncin$var[[j]]$name
}
var.name
# [1] "PSRF"     "PRECTmms" "FSDS"     "QBOT"     "RH"       "TBOT"     "WIND"     "ZBOT"     "LONGXY"   "LATIXY"   "EDGEE"
# [12] "EDGEW"    "EDGES"    "EDGEN"
## selected variables
vars <- var.name[c(1, 2, 3, 5, 7)]
vars
# "PSRF"     "PRECTmms" "FSDS"     "RH"       "WIND"
nvars <- length(vars)
## matrix to store individual file var data
pb <- txtProgressBar(min = 0, max = total, style = 3)
for (i in start.month: length(all.files)){ #l
  setTxtProgressBar(pb, i)
  filename <- all.files[i]
  nc <- nc_open(filename, write = F)
  ## time vector
  time <- ncvar_get(nc,"time")
  tunits <- ncatt_get(nc, "time", "units")
  tustr <- strsplit(tunits$value, " ")
  tdstr <- strsplit(unlist(tustr)[3], "-")
  tmonth <- as.integer(unlist(tdstr)[2])
  tday <- as.integer(unlist(tdstr)[3])
  tyear <- as.integer(unlist(tdstr)[1])
  time.vec <- as.POSIXct(chron(time, origin = c(tmonth, tday, tyear)), tz = "GMT")
  attr(time.vec, "tzone") <- "GMT"
  # form a dataframe to store extraction from the current nc file
  current.df <- setNames(data.frame(matrix(ncol = nvars + 1, nrow = length(time.vec))),
                         c("datetime", vars))
  current.df$datetime  <- time.vec
  for (j in 1: nvars) {
  var.name.on <- vars[j]
  val <- ncvar_get(nc, var.name.on)
  # extract and store a variable's data
  current.df[, var.name.on] <- val
  }
  # current.df <- current.df[1:length(val), ]
  # append to past extractions
  if (i == 1) {
    nc.dat <- current.df
  } else {
    nc.dat <- bind_rows(nc.dat, current.df)
  }
  nc_close(nc)
}

str(nc.dat)

full.datetime <- seq(from = as.POSIXct(nc.dat$datetime[1]),
              to = as.POSIXct(nc.dat$datetime[nrow(nc.dat)]), by = "1 hour")
length(full.datetime); nrow(nc.dat)

## because some hours in ncdat are not rounded to the hour, correcting that:
nc.dat$datetime[which(format(nc.dat$datetime, "%M") == "59")] <- nc.dat$datetime[which(format(nc.dat$datetime, "%M") == "59")]  + 01:01

str(full.datetime); str(nc.dat$datetime)
xx <- full.datetime[!as.character(full.datetime) %in% as.character(nc.dat$datetime)]
length(xx)
## nc.dat does not have leap year Feb 29; and those are the only days missing from nc.dat
## raw data has this date, but hardly ever rains that day, so should give similar monthly sum of rains
nc.dat$year <- format(nc.dat$datetime, "%Y")
nc.dat.yr <- nc.dat %>% subset(!is.na(year) & year != 2019) %>%
  group_by(year) %>%
  summarise(rain = sum(PRECTmms, na.rm = TRUE)*60*60) # mms*60*60 = mm/hr. Since data is given every hour summed over a year, it is in mm.
ggplot(nc.dat.yr, aes(x = year, y = rain)) +
  geom_bar(stat = "identity") + ggtitle("Annual rainfall at BCI from 1985 through 2018 from metv5.1 .nc files")
ggsave(file.path(paste0("figures/met/Annual rainfall at BCI from 1985 through 2018 from metv5.1 .nc file.jpeg")), height = 10, width = 20, units='in')

xx <- nc.dat %>% subset(year == 2017)
View(xx)
str(xx)
summary(xx)
xx.long <- gather(xx, key = "variable", value = "value", -datetime, -year)
ggplot(xx.long, aes(x = datetime, y = value)) +
  geom_point() +
  facet_wrap( ~ variable, scales = "free_y")
###-------------------
### Boris's data - 2018 substituted by rutuja
###-------------------
raw.dat <- read.csv("data-raw/BCI_1985_2018c_mod_2018substituted.csv")
str(raw.dat)
# converting to system time zone

raw.dat$datetime <-  as.POSIXct(raw.dat$DateTime, format = "%m/%d/%y %H:%M", tz = "America/Panama")
head(raw.dat$datetime)
# the above takes in the given data as in Panama format, but converts it to Sys timezone
# so to convert to GMT:
attr(raw.dat$datetime, "tzone") <- "GMT"
head(raw.dat$datetime)
# so somehow it doesnt work, so removing one hour
raw.dat$datetime <- raw.dat$datetime - 60*60
head(raw.dat)
raw.dat$year <- format(raw.dat$datetime, "%Y")
raw.dat.yr<- raw.dat %>% subset(!is.na(year) & year != 2019) %>%
  group_by(year) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = TRUE)) # m3/m2*1000 == L/m2 == mm
ggplot(raw.dat.yr, aes(x = year, y = rain)) +
  geom_bar(stat = "identity")
## 2018 seems really low: checking which months
## That was before substitution not after
raw.dat$moyr <- format(raw.dat$datetime, "%y-%m")
raw.dat.moyr <- raw.dat %>% subset(!is.na(moyr) & year %in% c(2016:2018)) %>%
  group_by(moyr) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = TRUE)) # m3/m2*1000 == L/m2 == mm
ggplot(raw.dat.moyr, aes(x = moyr, y = rain)) +
  geom_bar(stat = "identity" )
View(raw.dat.moyr)

# ###-------------------
# ### Ryan's data
# ###-------------------
# unix.dat <- read.csv("data-raw/BCI_1985_2018_unix2.csv", header = TRUE)
# # unix.dat <- read.csv("data-raw/BCI_1985_2018c_mod.csv", header = TRUE)
# head(unix.dat)
# unix.dat$datetime <-  as.POSIXct(unix.dat$DateTime, format = "%m/%d/%y %H:%M", tz = "GMT")
# dif <- unix.dat$datetime %in% raw.dat$datetime
#
# full <- seq(from = as.POSIXct(unix.dat$datetime[1], tz = "GMT"),
#                       to = as.POSIXct(unix.dat$datetime[nrow(unix.dat)], tz = "GMT"), by = "1 hour")
# head(full)
# # what is in full that is not in unix.dat
# setdiff(full, unix.dat$datetime)
# numeric(0)
## so unix.dat is continuous and has all dates and hours without gaps
#
# # checking whether raw.dat has the same date-time range after removing hours in may as in unix.dat
# extra <- seq(from = as.POSIXct("2018-05-31 20:00:00", tz = "GMT"),
#              to = as.POSIXct("2018-06-04 09:00:00", tz = "GMT"), by = "1 hour")
# raw.dat.trunc <- raw.dat %>% subset(!datetime %in% extra)
#
# length(raw.dat.trunc$datetime)
# 292892
# length(unix.dat$datetime)
# 292896
# because unix.dat has four hours appended to the beginning

raw.dat.1 <- raw.dat %>% select(datetime, Rainfall_mm_hr) %>%
  mutate(month = as.Date(paste0("01-", format(datetime, "%m-%y")), format = "%d-%m-%y")) %>%
  group_by(month) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = T)) %>%
  mutate(source = "raw from Boris")
# unix.dat.1 <- unix.dat %>% select(datetime, Rainfall_mm_hr) %>%
#   mutate(month = as.Date(paste0("01-", format(datetime, "%m-%y")), format = "%d-%m-%y")) %>%
#   group_by(month) %>%
#   summarise(rain = sum(Rainfall_mm_hr, na.rm = T)) %>%
#   mutate(source = "unix csv")
# both <- bind_rows(raw.dat.1, unix.dat.1) %>% as.data.frame()
# ggplot(both, aes(x = month, y = rain, colour = source)) +
#   geom_line(aes(group = source))
# head(raw.dat.1); head(unix.dat.1)
# tail(raw.dat.1); tail(unix.dat.1)
# ggplot(raw.dat, aes(x = datetime, y = Rainfall_mm_hr)) +
#   geom_line()

raw.dat.2 <- raw.dat %>% select(datetime, Rainfall_mm_hr) %>%
  mutate(month = as.Date(paste0("01-", format(datetime, "%m-%y")), format = "%d-%m-%y")) %>%
  group_by(month) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = T)) %>%
  mutate(source = "raw csv")

nc.dat.1 <- nc.dat %>% select(datetime, PRECTmms) %>%
  mutate(month = as.Date(paste0("01-", format(datetime, "%m-%y")), format = "%d-%m-%y")) %>%
  group_by(month) %>%
  summarise(rain = sum(PRECTmms, na.rm = T)*60*60) %>% # mms*60*60 = mm/hr.
  mutate(source = "netcdf")
summary(nc.dat.1)
both2 <- bind_rows(nc.dat.1, raw.dat.2) %>% as.data.frame()
ggplot(both2, aes(x = month, y = rain, colour = source)) +
  geom_line(aes(group = source)) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%b%y"))

both2 <- both2 %>% mutate(year = format(month, "%Y"))
ggplot(both2 %>% subset(year %in% c(2017)), aes(x = month, y = rain, colour = source)) +
  geom_line(aes(group = source)) +
  scale_x_date(date_breaks = "1 month", labels = function(x) format(x, "%b%y"))

raw.dat.day <- raw.dat %>% select(datetime, Rainfall_mm_hr) %>%
  mutate(date = as.Date(datetime)) %>% select(-datetime) %>%
  group_by(date) %>%
  summarise(rain = sum(Rainfall_mm_hr, na.rm = T)) %>%
  mutate(source = "raw csv")
head(raw.dat.day)
nc.dat.day <- nc.dat %>% select(datetime, PRECTmms) %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(rain = sum(PRECTmms, na.rm = T)*60*60) %>% # mms*60*60 = mm/hr.
  mutate(source = "netcdf")
head(nc.dat.day)
summary(nc.dat.day)

ts.sub <- seq(from = as.Date("2017-12-27"),
             to = as.Date("2017-12-31"), by = "1 day")

both2.day <- bind_rows(nc.dat.day, raw.dat.day) %>% as.data.frame()
both2.day.sub <- both2.day %>% subset(date %in% ts.sub)
ggplot(both2.day.sub, aes(x = date, y = rain, colour = source)) +
  geom_line(aes(group = source))

## some hours in nc.df are not rounded to hours in the tz conversion so hourly comparison not exact until; that is fixed
nc.dat.1.sub <- nc.dat %>%
  #subset(year %in% c(2015, 2016, 2017, 2018)) %>%
  select(datetime, PRECTmms) %>%
  mutate(month = as.Date(paste0("01-", format(datetime, "%m-%y")), format = "%d-%m-%y")) %>%
  group_by(month) %>%
  summarise(rain = sum(PRECTmms, na.rm = T)*60*60) %>% # mms*60*60 = mm/hr, and sum of it would be total mm/month
  mutate(source = "netcdf")
g1 <- ggplot(nc.dat.1.sub, aes(x = month, y = rain)) +
  geom_line(aes(group = source)) +
  geom_point() +
  scale_x_date(date_breaks = "12 months", labels = function(x) format(x, "%b%y")) +
  ylab(" Rainfall [mm]")
g1 +   ggtitle("Monthly rainfall at BCI from 1985 through 2018")
ggsave(file.path(paste0("figures/met/Monthly rainfall at BCI from 1985 through 2018.jpeg")), height = 10, width = 20, units='in')

nc.dat.1.sub$year <- format(nc.dat.1.sub$month, "%Y")
g1 %+% subset(nc.dat.1.sub,  year %in% c(2015, 2016, 2017, 2018)) +
  scale_x_date(date_breaks = "2 month", labels = function(x) format(x, "%b%y")) +
  ggtitle("Monthly rainfall at BCI from 2015 through 2018")
ggsave(file.path(paste0("figures/met/Monthly rainfall at BCI from 2015 through 2018.jpeg")), height = 8.94, width = 12.9, units='in')



