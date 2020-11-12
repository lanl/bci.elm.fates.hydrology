##---------------------------
# Prepping to run best-fit ensembles on server
# Author: Rutuja Chitra-Tarak
# First written: Oct 27, 2019
##---------------------------

current.folder <- "2019-10-14_5000"
top.few <- 100
params.top.few <- read.csv(file = file.path("results", current.folder, paste0("params.top.few.cond_", top.few, ".csv")), header = TRUE)
params.obj.top.few <- read.csv(file = file.path("results", current.folder, paste0("params.obj.top.few.cond_", top.few, ".csv")), header = TRUE)

paste(params.top.few,  collapse = ",")
# [1] "c(3563, 3535, 2731, 3653, 3444, 1120, 2173, 1611, 3271, 2451, 2495, 2501, 2708, 1881, 3858, 2298, 1017, 1117, 1217, 1317, 1417, 1517, 4185, 2741, 3089, 3157, 2274, 2541, 3124, 4703, 3639, 1867, 2636, 2527, 1754, 1592, 3241, 2031, 2164, 2371, 2029, 2099, 2578, 2693, 3849, 1937, 3005, 3506, 4983, 1929, 2643, 1789, 4142, 2880, 3133, 2929, 3861, 1758, 4047, 2372, 3743, 1793, 2165, 1703, 2688, 1817, 2973, 3145, 4517, 2049, 3762, 3697, 2190, 2714, 1756, 4893, 2196, 1874, 3077, 3303, 2996, 3347, 2833, 1672, \n3588, 3898, 2809, 2766, 2353, 2179, 2666, 3589, 4514, 3361, 2689, 1979, 3450, 2082, 2773, 2805)"
chosen.ranges <- apply(params.obj.top.few[,1:13], 2, range)
save(chosen.ranges, file = "results/chosen.ranges.rda")

