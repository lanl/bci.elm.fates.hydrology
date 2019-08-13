

### Parameters for fsurdat

# double PCT_SAND(nlevsoi, gridcell) ;
# PCT_SAND:long_name = "percent sand" ;
# PCT_SAND:units = "unitless" ;
# double PCT_CLAY(nlevsoi, gridcell) ;
# PCT_CLAY:long_name = "percent clay" ;
# PCT_CLAY:units = "unitless" ;
# double ORGANIC(nlevsoi, gridcell) ;
# ORGANIC:long_name = "organic matter density at soil levels" ;
# ORGANIC:units = "kg/m3 (assumed carbon content 0.58 gC per gOM)" ;
# double FMAX(gridcell) ;
# FMAX:long_name = "maximum fractional saturated area" ;
# FMAX:units = "unitless" ;


###----
### ORGANIC
###----
# Grimm, R., Behrens, T., Marker, M., and Elsenbeer, Helmut. 2008. "Soil organic carbon concentrations and stocks on Barro Colorado Island -- Digital soil mapping using Random Forests analysis." Geoderma. 146 (1-2):102–113. 
# https://doi.org/10.1016/j.geoderma.2008.05.008

# Table 3:
# depth 0–10 10–20 20–30 30–50
# for all obs covering all soil types
# % SOC mean (SD) 5.00 (1.77) 2.13(0.64) 1.53 (0.46) 1.10 (0.35)
# for Ava soil type (Table 1) that represents the BCI 50ha plot; 
# slightly higher concentrations as seen in the maps: Fig 5 & fig 6
# 4.38 (1.88) 2.22 (0.53) 1.65 (0.46) 1.18 (0.35)
# also for Ava soil 
# SOC_Mgperha = 33.27 (15.11), 17.95 (4.91) 13.61 (4.62) 21.80 (7.02)
# Using SOC stock instead because its calculation (eq 1) includes 
# bulk density & stoniness of the soil (data not presented in the paper anywhere)

# Also see Section 3.2.4: 
# In order to compare lateral and vertical SOC distribution patterns 
# we computed the SOC stocks (Fig. 6a–d), since natural pedons include non-soil components 
# such as rocks and pebbles. 
# SOC stocks therefore reflect SOC distributions more realistically.
#  Section 3.1 : In general, high SOC concentrations translated into high carbon stocks. The pale swelling clays (Vertic Luvisol, Acrisol and Vertic Eutric or Alumic Gleysol) constitute an exception in that they have low SOC concentrations but high carbon stocks. 
#  This is due to their comparatively high bulk density and low stoniness. 
# om <- data.frame(depth = c(5, 15, 25, 40), 
#                  stock.soc = c(33.27, 17.95, 13.61, 21.8))
# # Mgperha in 10 cm depth = Kg/m3
# # double ORGANIC(nlevsoi, gridcell) ;
# # ORGANIC:long_name = "organic matter density at soil levels" ;
# # ORGANIC:units = "kg/m3 (assumed carbon content 0.58 gC per gOM)" ;
# v
# om$depth
# # 5 15 25 40 # cm
# om$om
# # 57 31 23 38 # kg/m3

####-----------------
### But according to the BCI soil report (http://biogeodb.stri.si.edu/bioinformatics/bci_soil_map/), see literature/
### Ava red light clays accounting for 72.4% of the plot (page 43)
## AVA soil profile C% is much lower as below
## see BCIAPPBPotsdam0216.docx page 4 in literature/.
# Profile: 	PF01			Soil class: Ava
# Horizon depth	cm  0-6	  6-26	26 – 51/79	51/79– 91/126	91/126-150	150 - 170
# Sample depth      0-6	  6-16	36 - 46	    66 - 76	      136 - 146	  156 - 166
# Organic C	%	      6.99	2.67	1.03	      nd	          nd	        nd

# Granulometry % of fine earth
# Sample depth      0-6	  6-16	36 - 46	    66 - 76	      136 - 146	  156 - 166
# Clay              60.0	64.5	72.0	      69.0	        21.0	      13.0
# Silt              28.0	28.5	22.0	      20.0	        63.0	      73.0
# Sand              11.5	7.0	  6.0	        11.0	        16.0	      14.0

## correspondence:
# for  depths 0.007, 0.03, 0.062, value for 0-6	
# for 0.11, 6-16 
# for 0.21, middle of 6-16	and 36 - 46
# for 0.36 36 - 46
# for 0.62 66 - 76	
# for 1.04, middle of 66 - 76	and 136 - 146	
# for 1.73, 156-166
# for 2.56, 156 - 166 reduced by half of diff between 136 - 146	  156 - 166
om <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                 stock.soc = c(6.99, 6.99, 6.99, 2.67, 1.85, 1.03, 0, 0, 0, 0), # after 0.36, guess work with gradual decrease,
                 PCT_CLAY = c(60.0,	60.0,	60.0,	64.5,	67, 72.0,	69.0, 45, 13.0, 9), # middle of 0
                 PCT_SAND = c(11.5,	11.5,	11.5,	7.0, 6.5, 6.0, 11.0, 13.5,	14.0, 13))
om$om <- round(om$stock.soc/0.58, 1) #in kg/m3
om
# depth stock.soc PCT_CLAY PCT_SAND   om
# 1  0.007      6.99     60.0     11.5 12.1
# 2  0.030      6.99     60.0     11.5 12.1
# 3  0.062      6.99     60.0     11.5 12.1
# 4  0.110      2.67     64.5      7.0  4.6
# 5  0.210      1.85     67.0      6.5  3.2
# 6  0.360      1.03     72.0      6.0  1.8
# 7  0.620      0.00     69.0     11.0  0.0
# 8  1.040      0.00     45.0     13.5  0.0
# 9  1.730      0.00     13.0     14.0  0.0
# 10 2.860      0.00      9.0     13.0  0.0
####-----------------

###-----
### PCT_SAND & PCT_SILT
###-----

# From: data-raw/BCI 50 ha plot soil texture_Matteo.xls
# 
# txr <- readxl::read_xls("data-raw/BCI 50 ha plot soil texture_Matteo.xls")
# head(txr)
# View(txr)
# txr1 <- txr %>% subset(!is.na(txr$Plot))
# head(txr1)
# txr.depth <- txr1 %>% select(Depth, `Bulk density`:Clay) %>% group_by(Depth) %>%
#   summarise_at(vars(`Bulk density`:Clay), mean, na.rm = TRUE)
# head(txr.depth)

# See Sheet 2 : mean_by_depth
# Depth	Bulk density	Sand	Silt	Clay
# cm	fine earth	%	%	%
# g/cm3			
# 
# 0–10	  1.0	12.0	32.6	55.4
# 10–20	  1.0	8.2	  28.8	63.1
# 20–40	  1.0	7.4	  25.7	67.0
# 40–60	  1.2	6.9	  22.6	70.6
# 60–120	1.3	13.1	23.4	63.5
# 120–260	1.3	14.7	26.1	59.2
# 180–260	1.3	15.5	37.2	47.3
# 260–360	 -  22.2	38.7	39.1
# 360–420+				
# And:
om$depth
# 0-10, 10-20, 20-30, 30-50 # cm
om$om
# 57 31 23 38 # kg/m3
# In OM after the first seven values (0-62 cm depth), tsking OM as approximately 10, and then 1, whose value is used for the depths > 3 m
## In the fsurdat file 10 depths are equivalent to the first 10 depths below:
####
nc <- nc_open( "data-raw/DTB4.all.nc", readunlim = FALSE)
# depths
nc$var[['H2OSOI']]$dim[[2]]$vals # 15 depths in m
# [1]  0.007100635  0.027925000  0.062258575  0.118865065  0.212193400  0.366065800  0.619758487
# [8]  1.038027048  1.727635264  2.864607096  4.739156723  7.829766273 12.925320625 21.326469421
# [15] 35.177619934

## the value at 11 cm depth is the middle of 0-10 and 10-20
## from the BCI texture data above, can estimate texture for fsurdat depths:
txr.2 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
           PCT_SAND = c(12, 12, 12, 10, 8.2, 7.4, 6.9, 13.1, 14.7, 15.5),
           PCT_CLAY = c(55.4, 55.4, 55.4, 60, 63.1, 67.0, 70.6, 63.5, 59.2, 47.3),
           # PCT_SILT = c(32.6, 32.6, 32.6, 30, 28.8, 25.7, 22.6, 23.4, 26.1, 37.2),
           #OM = c(57, 57, 57, 44, 31, 23, 38, 30, 10, 1))
           OM = c(12.1, 12.1, 12.1, 4.6, 3.2, 1.8, 0, 0, 0, 0))

## if only texture data from the North Location is used, instead of an average across the four locations:
# Plot	Profile code	Topography	Depth	    Horizon	Bulk density	Sand	Silt	Clay	Textural class
# BCI	North	Well-drained plateau	0–11	    A	      0.92	        10.2	24.8	65.0	Clay
# BCI	North	Well-drained plateau	11–30	    AB	    0.89	        6.0	  22.3	71.7	Clay
# BCI	North	Well-drained plateau	30–56	    Bto1	  1.19	        6.8	  16.1	77.1	Clay
# BCI	North	Well-drained plateau	56–120	  Bto2	  1.4	          8.2	  10.9	80.8	Clay
# BCI	North	Well-drained plateau	120–180	  Bo1	    1.4	          14.6	15.3	70.1	Clay
# BCI	North	Well-drained plateau	180–260+	Bo2	–				

txr.3 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                  PCT_SAND = c(10.2, 10.2, 10.2, 8, 6, 6.8, 7.5, 8.2, 14.6, 15.5),
                  PCT_CLAY = c(65, 65, 65, 60, 66, 71.7, 77.1, 79.0, 70.1, 47.3),
                  OM = c(12.1, 12.1, 12.1, 4.6, 3.2, 1.8, 0, 0, 0, 0))
## Using AVA profile data; because that is where the AVA tower is
# from om

txr.4 <- data.frame(depth = c(0.007, 0.03, 0.062, 0.11, 0.21, 0.36, 0.62, 1.04, 1.73, 2.86),
                  PCT_SAND = c(11.5, 11.5, 11.5, 7.0, 6.5, 6.0, 11.0, 13.5,	14.0, 13),
                  PCT_CLAY = c(60.0, 60.0, 60.0, 64.5, 67, 72.0, 69.0, 45, 13.0, 9),
                  OM = c(12.1, 12.1, 12.1, 4.6, 3.2, 1.8, 0, 0, 0, 0))

### default values:
# txr.1 <-
# PCT_SAND =
#   // PCT_SAND(0,0)
# 42,
# // PCT_SAND(1,0)
# 42,
# // PCT_SAND(2,0)
# 42,
# // PCT_SAND(3,0)
# 41,
# // PCT_SAND(4,0)
# 41,
# // PCT_SAND(5,0)
# 38,
# // PCT_SAND(6,0)
# 37,
# // PCT_SAND(7,0)
# 36,
# // PCT_SAND(8,0)
# 38,
# // PCT_SAND(9,0)
# 36 ;
# 
# PCT_CLAY =
#   // PCT_CLAY(0,0)
# 30,
# // PCT_CLAY(1,0)
# 30,
# // PCT_CLAY(2,0)
# 31,
# // PCT_CLAY(3,0)
# 32,
# // PCT_CLAY(4,0)
# 34,
# // PCT_CLAY(5,0)
# 39,
# // PCT_CLAY(6,0)
# 42,
# // PCT_CLAY(7,0)
# 42,
# // PCT_CLAY(8,0)
# 39,
# // PCT_CLAY(9,0)
# 40 ;
# 
# ORGANIC =
#   // ORGANIC(0,0)
# 34.50443584523,
# // ORGANIC(1,0)
# 35.1579594875661,
# // ORGANIC(2,0)
# 30.220404705696,
# // ORGANIC(3,0)
# 24.7809408727215,
# // ORGANIC(4,0)
# 19.864585793084,
# // ORGANIC(5,0)
# 15.731521316196,
# // ORGANIC(6,0)
# 12.3733677054129,
# // ORGANIC(7,0)
# 9.69339947263121,
# // ORGANIC(8,0)
# 0,
# // ORGANIC(9,0)
# 0 ;
# 
# FMAX = 0.381491014791347 ;