
# For Conrad catchment observed median Ksat (Lower95%CI, Upper95%CI) 
# for 12.5 cm depth is 38.3 mm/hr (25.4, 51.2, n = 75), while for 60 cm depth 0.7 (0.2, 1.2, n = 40) mm/hr. For reference a storm with 12.5 mm/hr rainfall intensity has a 0.2 probability of occuring in any given rainfall event (Kinner and Stallard 2004).
# In mm/s: for 12.5 cm depth: 0.01063889 (0.007055556, 0.01422222) & for 60 cm depth 0.0001944444 (5.555556e-05, 0.0003333333)
ksat_data = data.frame(depth = c(0.125, 0.60),
                       ksat_median = c(0.01063889, 0.0001944444),
                       ksat_lower = c(0.007055556, 5.555556e-05),
                       ksat_upper = c(0.01422222, 0.0003333333))
# equation for median ksat_obs: y = mx + b; 
# m = (y2-y1)/(x2-x1); b = y- mx
m.median <- c(ksat_data$ksat_median[1] - ksat_data$ksat_median[2])/c(ksat_data$depth[1] - ksat_data$depth[2])
b.median <- ksat_data$ksat_median[1] - m.median*ksat_data$depth[1]
m.median; b.median
# -0.02198831;0.01338743
## see if we get correct y:
m.median*ksat_data$depth[1] + b.median
# 0.01063889
## Equation for lower 95% CI
m.lower <- c(ksat_data$ksat_lower[1] - ksat_data$ksat_lower[2])/c(ksat_data$depth[1] - ksat_data$depth[2])
b.lower <- ksat_data$ksat_lower[1] - m.lower*ksat_data$depth[1]
m.lower; b.lower
# -0.01473684;0.008897661
## see if we get correct y:
m.lower*ksat_data$depth[1] + b.lower
# 0.007055556

# Equation for upper 95% CI
m.upper <- c(ksat_data$ksat_upper[2] - ksat_data$ksat_upper[1])/c(ksat_data$depth[2] - ksat_data$depth[1])
b.upper <- ksat_data$ksat_upper[1] - m.upper*ksat_data$depth[1] 
m.upper; b.upper
# -0.02923976;0.01787719
## see if we get correct y:
m.upper*ksat_data$depth[2] + b.upper
# 0.0003333333

library(tidyverse)
p1 <- ggplot(ksat_data, aes(x = depth)) +
  geom_ribbon(aes(ymin = ksat_lower, ymax = ksat_upper), fill = "grey70") +
  geom_line(aes(y = ksat_median)) + 
  ylab("Ksat Median + 95CI [mm/s]") + xlab("Depth [m]") +
  ggtitle("Observed Ksat data by depth at Conrad Catchment - Godsey et al 2004") +
  geom_text(aes(x = 0.4, y = 0.012, label = "Ksat = -0.022*depth + 0.0134"), size = 6) +
  theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p1
ggsave(file.path("figures/Ksat_profiles/Observed Ksat data by depth at Conrad Catchment - Godsey et al 2004.jpeg"), 
       height = 5, width = 5.5, units='in')

## ksat range required:
# CLM file SoilStateType.F90 line 725
#this%hksat_col(c,lev)  = this%hksat_adj_col(c) * (uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat)

ksat <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86, 
                             4.74, 7.83, 12.93, 21.33, 35.18),
                   #from this%hksat_col(1,1:10) when txr.4 used with hksat_adj = 1:
                   model_ksat = c(0.0013939390371644, 0.0013939390371644, 0.0013939390371644, 
                                  0.00118079650702902, 0.00115943202825438, 0.0011387174232642,
                                  0.00135779423336252, 0.00148280472100791, 0.0015091553933691, 
                                  0.0014569141456897, 0.0014569141456897, 0.0014569141456897,
                                  0.0014569141456897, 0.0014569141456897, 0.0014569141456897)) %>% 
  ## from depth 10 - 15 same ksat values, as same texture and OM
  mutate(ksat_median = m.median*depth + b.median,
         ksat_median = replace(ksat_median, depth > 0.6, ksat_data$ksat_median[2]),
         ksat_lower = m.lower*depth + b.lower,
         ksat_lower = replace(ksat_lower, depth > 0.6, ksat_data$ksat_lower[2]),
         ksat_upper = m.upper*depth + b.upper,
         ksat_upper = replace(ksat_upper, depth > 0.6, ksat_data$ksat_upper[2]))
ksat  
# all depths
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934
p2 <- ggplot(ksat, aes(x = depth)) +
  scale_x_log10() + geom_rug() +
  ylab(expression("Ksat (mm s"^-1*")")) + xlab("Depth (m)") +
  theme(axis.text.y   = element_text(size=12),
        axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  )
p2 +  geom_ribbon(aes(ymin = ksat_lower, ymax = ksat_upper), fill = "grey70") +
  geom_line(aes(y = ksat_median)) + 
  ggtitle("Beyond Observed Ksat data by depth at Conrad Catchment - Godsey et al 2004")
ggsave(file.path("figures/Ksat_profiles/Using Observed Ksat data by depth at Conrad Catchment - Godsey et al 2004.jpeg"), 
         height = 3, width = 3.5, units='in')
# Generating parameter ensembles
## by randomly selecting m and b within the 95% CI range, 
# however m & b cannot be completely independent, so removing those y = mx+b 
# that go beyond the observed intercept for depth 0.6
## changing negative slopes to positive as qunif won't work otherwise
info <- list(par.names = c("m", "b"), min.param = c(-m.lower, b.lower), 
             max.param = c(-m.upper, b.upper)) 
# $par.names
# [1] "m" "b"
# 
# $min.param
# [1] 0.014736843 0.008897661
# 
# $max.param
# [1] 0.02923976 0.01787719

n.sam <- 10000
n.par <- length(info$par.names)
grid <- lhs::randomLHS(n.sam, n.par)
params.1 <- matrix(0, n.sam, n.par)
colnames(params.1) <- info$par.names
## generating ensembles
for (i in 1: n.sam) {
  for (j in 1: n.par) {
    params.1[i, j] <- qunif(grid[i, j], min = info$min.param[j], max = info$max.param[j])
  }
}
summary(params.1)
# m                 b           
# Min.   :0.01474   Min.   :0.008899  
# 1st Qu.:0.01836   1st Qu.:0.011143  
# Median :0.02199   Median :0.013388  
# Mean   :0.02199   Mean   :0.013387  
# 3rd Qu.:0.02561   3rd Qu.:0.015632  
# Max.   :0.02924   Max.   :0.017877
params.1[,"m"] <- -params.1[,"m"]
ksat_boot.1 <- data.frame(depth = ksat$depth)
samp.remove <- vector()
for(i in 1: n.sam){
  col.i <- paste0("ksat.", i)
  ksat_boot.1[, col.i] <- params.1[i, "m"]*ksat_boot.1$depth + params.1[i, "b"]
  # however m & b cannot be completely indeependent, so removing those y = mx+b 
  # that go beyond the observed intercept for depth 0.6
  if (ksat_boot.1[ksat_boot.1$depth == 0.620, col.i] < ksat_data$ksat_lower[2] || ksat_boot.1[ksat_boot.1$depth == 0.620, col.i] > ksat_data$ksat_upper[2] ) {
    samp.remove <- c(samp.remove, i)}
}
ksat_boot.accept <- ksat_boot.1[, -(samp.remove + 1)]
ncol(ksat_boot.accept) -1
## very few profiles would get generated this way
params.accept <- data.frame(params.1[-samp.remove,])
ggplot(params.accept, aes(x = m, y = b)) +
  geom_point() + theme_classic() + xlab("Slope") + ylab("Intercept") +
  geom_text(aes(x = -0.021, y = 0.018, label = "Intercept = -062*Slope + 9.8e-05"), size = 6) +
  ggtitle("Constrain between slope and intercept for Ksat profiles within 95% CI")
ggsave(file.path("figures/Ksat_profiles/Constrain between slope and intercept for Ksat profiles within 95CI.jpeg"), 
       height = 5, width = 5.5, units='in')

## using the constrain between acceptable m and b
lm1 <- lm(b ~ m, data = params.accept)
## generating 10000 ensembles
params.2 <- data.frame(m = params.1[, "m"], b = predict(lm1, newdata = data.frame(m = params.1[, "m"])))

ksat_boot.2 <- data.frame(depth = ksat$depth)
samp.remove.2 <- vector()
for(i in 1: n.sam){
  col.i <- paste0("ksat.", i)
  ksat_boot.2[, col.i] <- params.2[i, "m"]*ksat_boot.2$depth + params.2[i, "b"]
  # however m & b cannot be completely indeependent, so removing those y = mx+b 
  # that go beyond the observed intercept for depth 0.6
  if (ksat_boot.2[ksat_boot.2$depth == 0.620, col.i] < ksat_data$ksat_lower[2] || ksat_boot.2[ksat_boot.2$depth == 0.620, col.i] > ksat_data$ksat_upper[2] ) {
    samp.remove.2 <- c(samp.remove.2, i)}
  ## saturating function beyond the observed data range of depths
  ksat_boot.2[ksat_boot.2$depth < 0.2, col.i] <- params.2[i, "m"]*0.125 + params.2[i, "b"]
  ksat_boot.2[ksat_boot.2$depth > 0.6, col.i] <- params.2[i, "m"]*0.6 + params.2[i, "b"]
}
samp.remove.2 # should be zero
write.csv(ksat_boot.2, file = "data-raw/Ten thousand bootstrapped Ksat profiles within 95CI of observed data.csv", row.names = FALSE)
ksat_boot.long.2 <- gather(ksat_boot.2[,1:201], key = "sample", value = "ksat", -depth)
summary(ksat_boot.long.2)
# depth           sample               ksat          
# Min.   : 0.007   Length:3000        Min.   :0.0004829  
# 1st Qu.: 0.110   Class :character   1st Qu.:0.0006126  
# Median : 1.040   Mode  :character   Median :0.0007297  
# Mean   : 5.936                      Mean   :0.0043582  
# 3rd Qu.: 7.830                      3rd Qu.:0.0090461  
# Max.   :35.180                      Max.   :0.0146341 
p2 %+% ksat_boot.long.2 +
  geom_line(aes(y = ksat, group = sample, color = sample), show.legend = FALSE, size = 0.2) #+
  # ggtitle("Bootstrapped Ksat profiles within observed 95% CI")
ggsave(file.path("figures/Ksat_profiles/Bootstrapped Ksat profiles within observed 95CI_random_two_hundred.jpeg"), 
       height = 3, width = 3.5, units='in')
ggsave(file.path("figures/Ksat_profiles/Bootstrapped Ksat profiles within observed 95CI_random_two_hundred.tiff"), 
       height = 3, width = 3.5, units='in')
### generating only 10 equidistant values of m in the range
few.samp <- 20
m.few <- -seq(from = -max(params.accept$m), to = -min(params.accept$m), length.out = few.samp)
params.3 <- data.frame(m = m.few, b = predict(lm1, newdata = data.frame(m = m.few)))

ksat_boot.3 <- data.frame(depth = ksat$depth)
samp.remove.3 <- vector()
for(i in 1: few.samp){
  col.i <- paste0("ksat.", i)
  ksat_boot.3[, col.i] <- params.3[i, "m"]*ksat_boot.3$depth + params.3[i, "b"]
  # however m & b cannot be completely indeependent, so removing those y = mx+b 
  # that go beyond the observed intercept for depth 0.6
  if (ksat_boot.3[ksat_boot.3$depth == 0.620, col.i] < ksat_data$ksat_lower[2] || ksat_boot.3[ksat_boot.3$depth == 0.620, col.i] > ksat_data$ksat_upper[2] ) {
    samp.remove.3 <- c(samp.remove.3, i)}
  ## saturating funciton beyond the observed data range of depths
  ksat_boot.3[ksat_boot.3$depth < 0.3, col.i] <- params.3[i, "m"]*0.125 + params.3[i, "b"]
  ksat_boot.3[ksat_boot.3$depth > 0.6, col.i] <- params.3[i, "m"]*0.6 + params.3[i, "b"]
}
samp.remove.3 # should be 0
write.csv(ksat_boot.3, file = "data-raw/Few bootstrapped Ksat profiles within 95CI of 
                                        observed data.csv", row.names = FALSE)
ksat_boot.long.3 <- gather(ksat_boot.3, key = "sample", value = "ksat", -depth)
# summary(ksat_boot.long.3)
# depth           sample               ksat          
# Min.   : 0.007   Length:150         Min.   :0.0004962  
# 1st Qu.: 0.110   Class :character   1st Qu.:0.0005888  
# Median : 1.040   Mode  :character   Median :0.0006815  
# Mean   : 5.936                      Mean   :0.0044057  
# 3rd Qu.: 7.830                      3rd Qu.:0.0092268  
# Max.   :35.180                      Max.   :0.0141033
p2 %+% ksat_boot.long.3 +
  geom_line(aes(y = ksat, group = sample, color = sample), show.legend = FALSE) +
  ggtitle("Few Bootstrapped Ksat profiles within observed 95% CI")
ggsave(file.path("figures/Ksat_profiles/Few Bootstrapped Ksat profiles within observed 95CI_random_two_hundred.jpeg"), 
       height = 5, width = 5.5, units='in')
