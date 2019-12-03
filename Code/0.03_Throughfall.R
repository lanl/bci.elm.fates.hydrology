# library(tidyverse)
# library(bci.hydromet)
# readr::read_table2("data-raw/Zimmermann et al_2009_WRR_TableA1", skip = 5)
# thr <- read.table("data-raw/Zimmermann et al_2009_WRR_TableA1", skip = 5, header =  FALSE)
# str(thr)
# colnames(thr) <- c("event", "date", "median", "mad", "mad/median", "skewness", "octile")
# thr <- thr %>% mutate(date = as.Date(date, format = "%d-%m-%y", tz = "America/Panama"))
# ## dates surrounding the throughfall dates
# 
# surr.dates <- vector()
# for( i in 1:nrow(thr)){
#   i.dates <- seq(from = thr$date[i] - 2, 
#                  to = thr$date[i] + 1, by = "1 day")
#   surr.dates <- c(surr.dates, i.dates) 
# }
# surround.dates <- unique(surr.dates); length(surround.dates)
# 
# ###-------------------
# ### Boris's met data file - 2018 substituted by rutuja
# ###-------------------
# # met <- read.csv("data-raw/BCI_1985_2018c_mod_2018substituted.csv")
# met <- read.csv("data-raw/Precip_hourly_2007_2008.csv")
# str(met)
# # converting to system time zone
# 
# # met$datetime <-  as.POSIXct(met$DateTime, format = "%m/%d/%y %H:%M", tz = "America/Panama") 
# met$datetime <-  as.POSIXct(met$Hr, format = "%m/%d/%y %H:%M", tz = "America/Panama") 
# 
# head(met$datetime)
# # so somehow it doesnt work, so removing one hour
# met$datetime <- met$datetime - 60*60
# head(met)
# met$date <- as.Date(met$datetime, tz = "America/Panama")
# met.sub <- met %>% mutate(rain = Precip) %>% 
#   select(datetime, date, rain) %>% subset(date %in% surround.dates)
# View(met.sub)
# ## adding throughfall at the end of the day
# head(thr)
# thr <- thr %>% mutate(datetime = as.character(paste0(date, " 23:00:00"), format = "%Y-%m-%d %H:%M:%S"))
# head(thr)
# str(thr)
# str(met.sub); nrow(met.sub)
# met.thr <- met.sub %>% mutate(datetime = as.character(datetime)) %>%
#   left_join(thr, by = "datetime") 
# nrow(met.thr)
# write.csv(met.thr, "data-raw/throughfall_rain_hourly_for_inspection_new.csv", row.names = FALSE)
# 
# rain <- met %>% 
#   group_by(date) %>%
#   summarise(rain = sum(Precip, na.rm = TRUE)) # m3/m2*1000 == L/m2 == mm
# head(rain)
# 
# rain.sub <- rain %>% subset(date %in% surround.dates)
# View(rain.sub)
# # bci.hydromet::forcings, 
# df <- left_join(rain.sub, thr, by = "date") %>% mutate(pc.thr = median/rain * 100)
# View(df)
# str(df)
# summary(df)
