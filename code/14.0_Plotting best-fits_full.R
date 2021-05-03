##---------------------------
# Generating growth correlates from best-fit swp
# Author: Rutuja Chitra-Tarak
# First written: Oct 3, 2019
##---------------------------
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
pkgs=c("lubridate", "tidyverse", "data.table")
groundhog.library(pkgs, groundhog.day)

theme_set(theme_bw())
theme_update(text = element_text(size=14),
             panel.grid.major.x = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank()
)

census.meds <- readr::read_rds("data-raw/census.mediandates.rds")
census.beg <- census.meds[3: length(census.meds)]
cut.breaks <- census.beg
cut.labels.2 <- paste0(seq(1990, 2010, by = 5), "-", seq(1995, 2015, by = 5))
require(scales);
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt", 
    transform = function(x) -sqrt(abs(x)), 
    inverse = function(x) x^2);
}

reverselog_trans <- function(base = exp(1)) {
  scales::trans_new(name = paste0("reverselog-", format(base)),
                    log_breaks(base = base),
                    domain = c(1e-100, Inf),
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x))
}
current.folder <- "2019-10-14_5000"
sub.folder <- "best-fits.full"
if(!dir.exists(file.path("figures", current.folder, "best-fits.full"))) {dir.create(file.path("figures", current.folder, "best-fits.full"))}

top.few <- 100
params.top.few.df <- read.csv(file.path("results", current.folder,
                                          paste0("params.top.few_", top.few, ".csv")), header = TRUE)
params.top.few <- params.top.few.df$x

##-----------------------------
## get best-fit SOILPSI or swp
##-----------------------------
load(file.path("data/psi.rda"))
load(file.path("data/psi.mean.rda"))
## by depth panels
plot.psi <- ggplot(psi, aes(x = date, y = -psi)) +
  scale_y_continuous(trans=reverselog_trans(10),
                     breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
                     labels = c(0.001, 0.01, 0.1, 1, 10, 100)) +
  annotation_logticks() +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = c(as.Date("1991-01-01"), census.beg[-1]), color = "gray") +
  facet_grid(depth ~ .) +
  ylab("-Soil Water Potential [MPa]") + xlab("Date") + 
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Potential for best-fit ensembles")
ggsave("psi_model_daily_all_depths_params.top.few_full.jpeg", plot = plot.psi, path = 
         file.path("figures", current.folder, sub.folder), height = 12, width = 8.94, units='in')

plot.psi.for.ERD <- plot.psi %+% subset(psi, date >= as.Date("1991-01-01")) 
ggsave("psi_model_daily_all_depths_params.top.few_full_beg_1991.jpeg", plot = plot.psi.for.ERD, path = 
         file.path("figures", current.folder, sub.folder), height = 12, width = 8.94, units='in')

# Heatmap 
plot.psi.mean <- ggplot(psi.mean, aes(date, as.factor(-depth), fill = psi)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Date") +
  scale_fill_viridis_c("PSI [MPa]", trans = "reverse", option = "plasma") +
  ggtitle("Average SOILPSI by depth across best-fit ensembles")
ggsave("psi_mean_across_params.top.few_full.jpeg", plot = plot.psi.mean, path = 
         file.path("figures", current.folder, sub.folder), height = 8.94, width = 8.94, units='in')

## interval summary
psi <- psi %>% 
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks, 
                            labels = cut.labels.2, right = TRUE))
psi.stat.1 <- psi %>% 
  group_by(date, depth) %>%
  summarise(mean = -mean(psi, na.rm = TRUE),
            upper.CI = -quantile(psi, probs = 0.95),
            lower.CI = -quantile(psi, probs = 0.05))

plot.psi.stat.1 <- ggplot(psi.stat.1,
                        aes(x = date, y = mean)) +
  geom_line(show.legend = F, size = 0.1, color = "blue") +
  geom_ribbon(aes(ymin = lower.CI, ymax = upper.CI), alpha = 0.3) +
  geom_vline(xintercept = census.beg) +
  scale_y_reverse() +
  facet_grid(depth ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  xlab("Date") + 
  theme(text = element_text(size = 12)) +
  ggtitle("PSI:Violins of interval means of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.psi.stat.1, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 12, width = 7, units='in')

psi.stat.2 <- psi %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, par.sam, depth) %>%
  summarise(mean = -mean(psi, na.rm = TRUE))

plot.psi.stat.2 <- ggplot(psi.stat.2,
                             aes(x = interval.yrs, y = mean)) +
  # geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_violin(fill = "grey") +
  scale_y_reverse() +
  facet_grid(depth ~ .) +
  ylab("-Soil Water Potential [MPa]") +
  xlab("Date") + 
  theme(text = element_text(size = 12)) +
  ggtitle("PSI:Violins of interval means of best-fit ensembles")
ggsave("psi_model_daily_bestfit_params.top.few_Violins_full.jpeg", plot = plot.psi.stat.2, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 12, width = 7, units='in')
psi.stat.3 <- psi %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, depth) %>%
  summarise(mean = -mean(psi, na.rm = TRUE))
# Heatmap 
plot.psi.mean.stat.3 <- ggplot(psi.stat.3, aes(interval.yrs, as.factor(-depth), fill = -mean)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Census Interval") +   
  scale_fill_viridis_c("PSI [MPa]", trans = "reverse", option = "plasma") +
  ggtitle("Average SOILPSI by interval and depth across best-fit ensembles")
ggsave("psi_mean_across_params.top.few_full_interval.jpeg", plot = plot.psi.mean.stat.3, path = 
         file.path("figures", current.folder, sub.folder), height = 8.94, width = 8.94, units='in')

###
#### BTRAN #-------
####
# Daily #------
####
load("data/btran.rda")
load("data/btran.stat.rda")

plot.btran <- ggplot(btran,
                  aes(x = date, y = value)) +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = census.beg, color = "gray") +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("BTRAN: for best-fit ensembles Daily")
ggsave("BTRAN_model_daily_bestfit_params.top.few_full.jpeg", plot = plot.btran, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 3, width = 8.94, units='in')
## mean + 95% CI
plot.btran.stat <- ggplot(btran.stat,
                     aes(x = date, y = mean)) +
  geom_line(show.legend = F, size = 0.1, color = "blue") +
  geom_ribbon(aes(ymin = lower.CI, ymax = upper.CI), alpha=0.3) +
  geom_vline(xintercept = census.beg) +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  theme(legend.text = element_text(size = 16, face = "plain"),
        legend.position = c(0.8, 0.9), legend.background = element_rect(fill = "transparent")) +
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("BTRAN: Simulated best-fits: mean with 95% Confidence Interval ")
ggsave("BTRAN_model_daily_bestfit_params.top.few_CI_full.jpeg", plot = plot.btran.stat, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 3, width = 8.94, units='in')

btran <- btran %>% 
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks, 
                            labels = cut.labels.2, right = TRUE))
btran.stat.2 <- btran %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, par.sam) %>%
  summarise(mean = mean(value, na.rm = TRUE))

plot.btran.stat <- ggplot(btran.stat.2,
                        aes(x = interval.yrs, y = mean)) +
  # geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_violin(fill = "grey") +
  theme(text = element_text(size = 12)) +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  ggtitle("BTRAN: Violins of interval means of best-fit ensembles")
ggsave("btran_model_daily_bestfit_params.top.few_Violins_full.jpeg", plot = plot.btran.stat, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 3, width = 6, units='in')

## also checking annual variation


btran <- btran %>% 
  mutate(year = format(date, "%Y"))
btran.stat.3 <- btran %>% 
  subset(!is.na(year)) %>%
  group_by(year, par.sam) %>%
  summarise(mean = mean(value, na.rm = TRUE))

plot.btran.stat.yr <- ggplot(btran.stat.3,
                          aes(x = year, y = mean)) +
  geom_violin(fill = "grey") +
  theme(text = element_text(size = 12)) +
  ylab("BTRAN [unitless]") +
  xlab("Date") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("BTRAN: Violins of annual means of best-fit ensembles")
ggsave("btran_model_daily_bestfit_params.top.few_Violins_full_yrly.jpeg", plot = plot.btran.stat.yr, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 3, width = 8.94, units='in')

###
#### Growth Factor #-------
####
# Daily #------
####
load("data/swp.gfac.rda")
plot.swp.gfac <- ggplot(swp.gfac, aes(x = date, y = gfac)) +
  # scale_y_continuous(trans="rev_sqrt", breaks = c(0, 0.5, 2, 5, 10, 15)) +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = census.beg, color = "gray") +
  facet_grid(depth ~ .) +
  ylab("Growth Factor [unitless, 0-1]") + xlab("Date") + 
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Growth Factor for best-fit ensembles")
ggsave("swp.gfac_model_daily_all_depths_params.top.few_full.jpeg", plot = plot.swp.gfac, path = 
         file.path("figures", current.folder, sub.folder), height = 12, width = 8.94, units='in')

## interval summary
swp.gfac <- swp.gfac %>% 
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks, 
                      labels = cut.labels.2, right = TRUE))

swp.gfac.stat <- swp.gfac %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, par.sam, depth) %>%
  summarise(mean = mean(gfac, na.rm = TRUE)) %>%
  droplevels()

plot.swp.gfac.stat <- ggplot(swp.gfac.stat,
                          aes(x = interval.yrs, y = mean)) +
  # geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_violin(fill = "grey") +
  facet_grid(depth ~ .) +
  ylab("Growth Factor [unitless, 0-1]") +
  xlab("Date") + 
  theme(text = element_text(size = 12)) +
  ggtitle("GrowthFactor:Violins of interval means of best-fit ensembles")
ggsave("swp.gfac_model_daily_bestfit_params.top.few_Violins_full.jpeg", plot = plot.swp.gfac.stat, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 12, width = 7, units='in')

swp.gfac.stat.2 <- swp.gfac %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, depth) %>%
  summarise(mean = mean(gfac, na.rm = TRUE)) %>%
  droplevels()

plot.swp.gfac.stat.2 <- ggplot(swp.gfac.stat.2, aes(interval.yrs, as.factor(-depth), fill = mean)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Census Interval") +   
  scale_fill_viridis_c("Growth\nFactor\n[unitless, 0-1]", trans = "reverse", option = "plasma") +
  ggtitle("Average Growth Factor by interval and depth across best-fit ensembles")
ggsave("swp.gfac_mean_across_params.top.few_full_interval.jpeg", plot = plot.swp.gfac.stat.2, path = 
         file.path("figures", current.folder, sub.folder), height = 8.94, width = 8.94, units='in')

##-----------------------------
## get soil water content
##-----------------------------
load(file.path("data/swc.rda"))

swc.mean <- swc %>%
  group_by(date, depth) %>% 
  summarise_at(vars(swc), mean, na.rm = TRUE) %>%
  arrange(depth) %>% subset(!is.na(swc)) %>% droplevels()


## by depth panels
plot.swc <- ggplot(swc, aes(x = date, y = swc)) +
  geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_vline(xintercept = census.beg, color = "gray") +
  facet_grid(depth ~ .) +
  ylab("Soil Water Content [mm3/mm3]") + xlab("Date") + 
  scale_x_date(date_breaks = "1 year", labels = function(x) format(x, "%Y")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Soil Water Content for best-fit ensembles")
ggsave("swc_model_daily_all_depths_params.top.few_full.jpeg", plot = plot.swc, path = 
         file.path("figures", current.folder, sub.folder), height = 12, width = 8.94, units='in')

# Heatmap 
plot.swc.mean <- ggplot(swc.mean, aes(date, as.factor(-depth), fill = swc)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Date") +   
  scale_fill_viridis_c("Soil\nWater\nContent\n[mm3/mm3]", trans = "reverse", option = "plasma") +
  ggtitle("Average SWC by depth across best-fit ensembles")
ggsave("swc_mean_across_params.top.few_full.jpeg", plot = plot.swc.mean, path = 
         file.path("figures", current.folder, sub.folder), height = 8.94, width = 8.94, units='in')

## interval summary

swc <- swc %>% 
  mutate(interval.yrs = cut(date, include.lowest = TRUE, breaks = cut.breaks, 
                            labels = cut.labels.2, right = TRUE))
swc.stat.2 <- swc %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, par.sam, depth) %>%
  summarise(mean = mean(swc, na.rm = TRUE))

plot.swc.stat.2 <- ggplot(swc.stat.2,
                          aes(x = interval.yrs, y = mean)) +
  # geom_line(aes(group = par.sam, color = as.factor(par.sam)), show.legend = F, size = 0.2) +
  geom_violin(fill = "grey") +
  facet_grid(depth ~ .) +
  ylab("Soil Water Content [mm3/mm3]") +
  xlab("Date") + 
  theme(text = element_text(size = 12)) +
  ggtitle("SWC:Violins of interval means of best-fit ensembles")
ggsave("swc_model_daily_bestfit_params.top.few_Violins_full.jpeg", plot = plot.swc.stat.2, 
       file.path("figures", current.folder, sub.folder), device = "jpeg", height = 12, width = 7, units='in')
# Heatmap 
swc.stat.3 <- swc %>% 
  subset(!is.na(interval.yrs)) %>%
  group_by(interval.yrs, depth) %>%
  summarise(mean = mean(swc, na.rm = TRUE))

plot.swc.mean.stat.2 <- ggplot(swc.stat.3, aes(interval.yrs, as.factor(-depth), fill = mean)) + 
  geom_tile() + ylab("Depth [m]") + xlab("Census Interval") +   
  scale_fill_viridis_c("Soil\nWater\nContent\n[mm3/mm3]", trans = "reverse", option = "plasma") +
  ggtitle("Average SWC by interval and depth across best-fit ensembles")
ggsave("swc_mean_across_params.top.few_full_interval.jpeg", plot = plot.swc.mean.stat.2, path = 
         file.path("figures", current.folder, sub.folder), height = 8.94, width = 8.94, units='in')
