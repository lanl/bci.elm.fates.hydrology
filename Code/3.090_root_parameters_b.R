## Estimating parameter b (fates_rootb_par) range for rooting profiles 


rm(list=ls())

if (!require("pacman")) install.packages("pacman"); library(pacman)
pacman::p_load(ncdf4, easyNCDF, tidyverse)

# FATES uses equation (2) in Zeng et al.
## y = 1 - 0.5(exp(-a*d) + exp(-b*d)) ....(eq 1)
# where d = depth and Y = root fraction, a = fates_roota_par, b = fates_rootb_par


# fates_roota_par: Parameter a in Eq. 2 in (Zeng 2001) that regulates the shape of the rooting profile. 
# Range 5.9 - 7.4 per m 
# Range corresponds to this parameter specified for BATS (or IGBP) land cover classified 
# Deciduous Broadleaf Trees (5.9 m-1) and Evergreen Broadleaf Trees  (7.4 m-1) as given in 
# Table 1 (or 2) of Zeng 2001.

# fates_rootb_par : Parameter b in Eq. 2 in (Zeng 2001) that regulates the depth of the rooting profile. 
# Chosen range of b is derived using this equation so as to fit the observed range of rooting depth (dr) 
# of 2 - 18 m for Tropical Deciduous Forest (mean ± Standard Error (S.E.); 
# 3.7 ± 0.5, n = 5 trees; min = 2, max = 4.7) and 
# Tropical Evergreen Forest (mean ± S.E.; 7.3 ± 0.5, n = 3 trees and 3 communities; min = 2, max = 18 m) 
# combined (Canadell et al. 1996). Besides the direct observation of roots at 18 m included by (Nepstad et al. 1994) 
# in Paragominas, eastern Amazonia that is included in the above study; 
# in Tapajos, eastern Amazonia water extraction by roots was also inferred up to 18 m. 
# (Davidson et al. Forest Science 2011)

# Rearranging eqn (1)
# exp(-a*d) + exp(-b*d) = 2*(1 - Y)
# exp(-b*d) = 2*(1 - Y) - exp(-a*d)
# b = -log(2*(1 - Y) - exp(-a*d))/d
# evaluate for Y = 0.99
# For a = 5.9 and d = 18 m
# exp(-b*18) = (1-0.99)*2 - exp(-5.9*18)
a <- 5.9; Y <- 0.99
d <- 18
b <- -log((1 - Y)*2 - exp(-a*d))/d
b
# 0.2173346
d <- 2
b <- -log((1 - Y)*2 - exp(-a*d))/d
# 1.956199

# For a = 7.4 and d = 18 m

# Therefore chosen range of b: 0.2173346 - 1.956199  per m

### Expanding these limits:
## for 90% ob roots to be in 20 cm, this would:
d <- 0.2; Y <- 0.9; 
# but for this combination to work a can't be 5.9 or 7.4 but at least 8.1 so that log is to be found for a positive number
a <- 8.1
b <- -log((1 - Y)*2 - exp(-a*d))/d
# 30.82599


# This depth is going to be the same for both `a`s, because as define the rooting profile, while bs define the depth of the rooting profile


####----------------------------------------------
#### Q1 : If all the roots of the forest were below 1 m: 
#### Will that reduce soil moisture variation at 1 m 
#### depth to match with the observations?
####----------------------------------------------

# To have most of the rooting depth by half a meter
d <- 0.5; Y <- 0.9999999
 # a needs to be > 12.... such that  (1 - Y)*2 - exp(-a*d) > 0 OR -log((1 - Y)*2)/d > a
-log((1 - Y)*2)/d
#30.8499
a <- -log((1 - Y)*2)/d
b <- -log((1 - Y)*2 - exp(-a*d))/d
#57.63159
d <- 1; Y <- 0.9999999
# a needs to be > 12.... such that  (1 - Y)*2 - exp(-a*d) > 0
-log((1 - Y)*2)/d
a <- -log((1 - Y)*2)/d
b <- -log((1 - Y)*2 - exp(-a*d))/d
# 50.59974
Y <- 0.9999999

rf.par <- data.frame(rf.sam = c(1,2,3)) %>%
  mutate(max.depth = c(0.5, 1, 2),
         a = -log((1 - Y)*2)/max.depth,
         b = -log((1 - Y)*2 - exp(-a*max.depth))/max.depth)

cum.rf.par <- data.frame(depth = rep(seq(0.1, 3, by = 0.1), each = nrow(rf.par))) 
rf.par.big <- rf.par
for (i in 2:c(nrow(cum.rf.par)/nrow(rf.par))) {
  rf.par.big <- rbind(rf.par.big, rf.par)
}
nrow(rf.par.big)
cum.rf.par <- cum.rf.par %>% bind_cols(rf.par.big) %>%
  mutate(cum.root.frac = 1 - 0.5*(exp(-a*depth) + exp(-b*depth))) %>%
  group_by(rf.sam) %>%
  mutate(root.frac = cum.root.frac - lag(cum.root.frac, default = 0)) %>%
  ungroup(rf.sam)

require(scales);
rev_sqrt_trans <- function() {
  scales::trans_new(
    name = "rev_sqrt", 
    transform = function(x) -sqrt(abs(x)), 
    inverse = function(x) x^2);
}
theme_set(theme_bw())
ggplot(data = cum.rf.par, aes(y = depth, x = root.frac)) +
  xlab("Root fraction") +
  geom_line(aes(group = rf.sam, color = as.factor(max.depth))) + 
  scale_colour_discrete(name = "Depth at \nCumulative Root Fraction = 0.9999999") + 
  theme(legend.position = c(0.5, 0.5)) +
  scale_y_continuous(trans="rev_sqrt", breaks = unique(cum.rf.par$depth)) +
  coord_cartesian(ylim = c(0.1, 2))
ggsave(file.path(paste0("figures/Root Fraction when Max Rooting depth around 1 m.jpeg")), height = 5, width = 5, units ='in')

       
ggplot(data = cum.rf.par, aes(y = depth, x = cum.root.frac)) +
  xlab("Cumulative Root Fraction") +
  geom_line(aes(group = rf.sam, color = as.factor(max.depth))) +
  scale_colour_discrete(name = "Depth at \nCumulative Root Fraction = 0.9999999") + 
  scale_y_continuous(trans="rev_sqrt", breaks = unique(cum.rf.par$depth)) +
  theme(legend.position = c(0.4, 0.5)) +
  scale_x_continuous(breaks = seq(0.1, 1, by = 0.1)) +
  coord_cartesian(ylim = c(0.1, 2))
ggsave(file.path(paste0("figures/Cumulative Root Fraction when Max Rooting depth around 1 m.jpeg")), height = 5, width = 5, units ='in')

####------------------------ End Q1 --------------------
  
## for the case that root fraction linearly increases with depth: 
l.18 <- data.frame(depth = seq(0.2, 18, by = 0.2)) %>%
  mutate(root.frac = rep(1/length(depth), length.out = length(depth))) %>%
  mutate(cum.root.frac = cumsum(root.frac)) %>%
  mutate(a = -log((1 - cum.root.frac)*2 - exp(-1.956199*depth))/depth)
l.18 <- l.18 %>%
  mutate(Y = 1 - 0.5*(exp(1.3178170*depth) + exp(-1.956199*depth)))
head(l.18)

View(l.18)
l.02 <- data.frame(depth = seq(0.2, 0.2, by = 0.2)) %>%
  mutate(root.frac = rep(1/length(depth), length.out = length(depth))) %>%
  mutate(cum.root.frac = cumsum(root.frac)) %>%
  mutate(a = -log((1 - cum.root.frac)*2 - exp(-1.956199*depth))/depth)

new <- data.frame(depth = seq(0.2, 18, by = 0.2),
                  root.frac)

# depths
####--------------------------------------------
#### Rooting profiles for inversion
####--------------------------------------------

nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)

soil.depths <- as.numeric(nc$var[['H2OSOI']]$dim[[2]]$vals)
# max.rds.df <- data.frame(max.rds = soil.depths[-length(soil.depths)]) %>% 
#   mutate(lag.max.rds = lag(max.rds, default = 0), 
#          diffby2 = (max.rds - lag.max.rds)/2,
#          half.depths = max.rds - diffby2)
# max.rds <- c(max.rds.df$max.rds, max.rds.df$half.depths)
exponents1plus <- data.frame(x = seq(-6, 6, by = 0.1)) %>% mutate(y = exp(x))
plot(data = exponents1plus, y~x)
# exponents.full <- rbind(exponents01, exponents1plus)
# plot(data = exponents.full, y~x)
max.rds.df <- data.frame(x = seq(-5.1, 3.1, by = 0.1)) %>% mutate(y = exp(x))
plot(data = max.rds.df, y~x)
max.rds <- c(max.rds.df$y)


profiles.list <- vector("list", length(max.rds))
profiles.dfs <- vector("list", length(max.rds))

exponents <- exponents1plus$y #seq(0.01, 30, by = 0.1) #c(seq(0.1, 1, length.out = 10), seq(1.2, 2.6, by = 0.2), seq(3, 15, length.out = 20))
nsam <- length(max.rds)*length(exponents)
## each rf.sam is associated with one maxD & power combination
for (i in 1: length(max.rds)) {
  maxD <- max.rds[i]
  profiles.list[[i]] <- vector("list", length(exponents))
  for (j in 1:length(exponents)) {
    rf.sam <- (i-1)*length(exponents) + j
    power <- exponents[j]
    profiles.list[[i]][[j]] <- data.frame(rf.sam = rf.sam, maxD = maxD, power = power,
                                     base = seq(0, 1, length.out = 100)) %>%
      mutate(abs.depth = base*maxD, 
             cum.root.frac = base^power,
             soil.depth.1 = soil.depths[1][match(signif(round(abs.depth, 3), 3), 
                                             signif(round(soil.depths[1], 3), 3))],
             soil.depth.rest = soil.depths[-1][match(signif(round(abs.depth, 2), 2), 
                                                    signif(round(soil.depths[-1], 2), 2))],
             soil.depth = if_else(is.na(soil.depth.1), soil.depth.rest, soil.depth.1),
             depth = ifelse(duplicated(soil.depth), NA, soil.depth))
  }
  profiles.dfs[[i]] <- do.call(rbind, profiles.list[[i]])
}

profiles <- do.call(rbind, profiles.dfs)
n.profiles <- unique(profiles$rf.sam)
length(n.profiles)
ggplot(profiles %>% subset(maxD == max.rds[length(max.rds)]), aes(y = abs.depth, x = cum.root.frac, color = maxD)) + 
  xlab("Depths [m]") +
  geom_line(aes(group = rf.sam)) + scale_y_continuous(trans="rev_sqrt", breaks = max.rds)
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion_one depth.jpeg")), height = 5, width = 8, units ='in')

ggplot(profiles, aes(y = abs.depth, x = cum.root.frac, color = maxD)) + 
  xlab("Depths [m]") +
  geom_line(aes(group = rf.sam)) + scale_y_continuous(trans="rev_sqrt", breaks = max.rds)
ggsave(file.path(paste0("figures/rooting_profiles_for_inversion.jpeg")), height = 5, width = 8, units ='in')

write.csv(profiles, "data-raw/root_profiles.csv", row.names = FALSE)

from <- c(NA, soil.depths); to <- c(NA, 1:length(soil.depths))
library(plyr)
profiles <- profiles %>% 
  mutate(soil.levels = plyr::mapvalues(depth, from, to),
  rf.sam.soil.levels = paste(rf.sam, soil.levels, sep = "."))
  # need to find depths that are equivalent to soil depths
  ## find which depths match which soil.depth and get those soil.depths in front of those depths
soil.depths[-length(soil.depths)]
unique(profiles$soil.depth)[-1]
unique(profiles$soil.levels)
View(profiles)
detach("package:plyr")
## but only first match would be enough
# so creating subgroup maxD.power.soildepths so that within those only first soil depth can be retained

profiles.sub1 <- profiles %>%
  subset(!is.na(soil.levels)) 
View(profiles.sub1)
## creating a table to work with ELM-FATES output, so all the soil depths present there are needed
# and a value for root fraction against it, 
## this is to be made for each nsam above.
## zeroes for root_frac when depth absent

pro.df <- data.frame(rf.sam = rep(1:nsam, each = length(soil.depths)),
                     depth = rep(soil.depths, times = nsam),
                     soil.levels = rep(1:length(soil.depths), times = nsam)) %>%
  mutate(rf.sam.soil.levels = paste(rf.sam, soil.levels, sep = ".")) %>%
  left_join(profiles.sub1 %>% select(rf.sam.soil.levels, cum.root.frac, maxD, base, power, abs.depth), 
            by = "rf.sam.soil.levels") 
## remove columns that vary within a rf.sam
cum.root.pro <- pro.df %>% select(-depth, -rf.sam.soil.levels, -abs.depth, -base) %>%
  ## maxD ab=nd opwer do not vary, but let's remove them too, because for some reason duplicate rows of rf.sam with NAs retained
  select(-maxD, -power) %>%
  pivot_wider(names_from = soil.levels, names_prefix = "level.", values_from = cum.root.frac)

cum.root.mat <- cum.root.pro %>% select(paste0("level.", 1:length(soil.depths)))
cum.root.mat.filled <- apply(cum.root.mat, 1, function(x){
  x <- as.numeric(x)
  n.na <- length(x) - length(x[!is.na(x)])  
  if ( n.na >= 1 ) {
    c(x[!is.na(x)], 1, rep(0, times = n.na - 1))
  } else{
    x
  }
}
)
cum.pro.df <- cbind(cum.root.pro$rf.sam, t(cum.root.mat.filled)) 
colnames(cum.pro.df) <- colnames(cum.root.pro)

pro.df <- cum.pro.df %>% data.frame() %>% 
  group_by(rf.sam) %>% pivot_longer(cols = starts_with("level."), 
                                    names_to = "soil.levels", 
                                    names_prefix = "level.",
                                    values_to = "cum.root.frac") %>%
  mutate(root.frac = cum.root.frac - lag(cum.root.frac, default = 0)) %>%
## after all soil layers with root.frac, at next soil layer, cum.root.frac should be  1, and 0 for the rest
  ungroup(rf.sam) %>%
  mutate(root.frac = if_else(root.frac == -1.0000000000, 0, root.frac))
View(pro.df)
library(plyr)
pro.df <- pro.df %>% mutate(soil.levels = as.numeric(soil.levels),
                            depth = round(plyr::mapvalues(soil.levels, to[-1], from[-1]), 2))

head(pro.df)
detach("package:plyr")
write.csv(pro.df, "data-raw/root.profiles.long.csv", row.names = FALSE)

head(pro.df)
View(pro.df)
root.pro <- pro.df %>% select(-soil.levels, -cum.root.frac) %>%
  pivot_wider(names_from = depth, names_prefix = "depth.", values_from = root.frac)
head(root.pro)
write.csv(root.pro, "data-raw/root.profiles.wide.csv", row.names = FALSE)

####---End--Rooting Profiles for inversion-------------------------------------------

