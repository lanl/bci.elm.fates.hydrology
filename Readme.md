Workflow
================
Rutuja Chitra-Tarak
12/2/2019

##### Met data files

Code/0.01\_Checking\_met\_data.R

### Generate parameter files

1.  FATES param files dependencies: Parameter range for roob\_b parameters are generated here: 3.090\_root\_parameters\_b.R Can verify that different root a can give same b: data-raw/estimation of root\_b\_par\_same\_b\_for\_different as.xlsx Other ranges are based on literature or Kunert's data

2.  ELM param files dependencies: fpi\_max ranges are generated from here: \#\# First throughfall from Zimmermann et al and rain data from STRI at Lutz are linked: Code/0.03\_Throughfall.R \#\# As rainfall data is not given in Zimmermann and throughfall data is collected for rainfall events (duration, amount and gap between events), rainfall that composes those events are inferred in Excel. And range of interception found for events &gt; 10 mm: data-raw/throughfall\_rain\_hourly\_for\_inspection\_new.xlsx

3.  Surface data files dependencies: Soil testure data and soil organic content data is gathered and defined here: Code/3.2\_Surfdata\_texture\_options.R But these are not used by the model as it is forced with

<!-- -->

1.  soil characteristic curves (SWC to Matric potential) using Stephan Kupers' data (Yilin's egquation), as well as, bci.hydromet/data-raw/soil\_retention\_curves\_stephan.R
2.  soil hydraulic conductivity, using Godsey & Stallard et al defined here: Code/3.2.5\_Surfdata\_Ksat\_obs\_and\_bootstrapped\_param.R

Surface data for ELM, ELM parameters and FATES paramater files Code/3.1\_Generate\_parameter files.R

These are transfer to the server.

Follow Google Doc 3.0\_Running ELM-FATES with parameter ensembles @ Shared drives/Rutuja\_work/web\_only/ <https://docs.google.com/document/d/1qqbkQGHMG8BMfrUejUlkBU3PvUqPuW0Pc6CBLBQzZZk/edit#heading=h.gjdgxs> Thus run ensembles and extract simulations on the server with Code/2.0\_main.R (that sources Code/1.0\_fun\_extract.R) and transfer back to the desktop storing at a specific dated location such as data-raw/extract/2019-10-14\_5000

Use Code/4.1\_ELM-FATES\_output.R to generate RMSE or Rsq values between simulations and hydrological observations, plot best-fits for individual fluxes and run sensitivity analyses. These outputs are compiled in individual-flux-best-fits.html and Report.html

Outputs from trials such as 4.1\_ELM-FATES\_output\_trial\_leafage\_fdrought.R, 4.1\_ELM-FATES\_output\_trial2.R, 4.1\_ELM-FATES\_output\_trial\_tlai\_reduction.R, are compiled in leafage\_fdrought.html, Effect-of-bbslope---VCmax.html or Max-Rooting-Depth-effect.html

Use Code/4.2\_ELM-FATES\_params\_bestfit.R to generate objective function to choose best-fit simulations that fits all individual fluxes well (QRUNOFF, AET, Soil moisture by depth)

Confirm that individual fluxes are well captured by thus chosen best-fits: Code/4.3\_ELM-FATES\_output\_bestfit.R

Now for the chosen best-fits re-run simulations covering the entire time-period of 1985-2018 Follow Google Doc 4.0\_Running ELM-FATES for best-fit parameter (and Code/5.0\_Prep\_for\_Best-fit\_case\_runs.R to generate par.sam numbers seperated by commas) ensembles<https://docs.google.com/document/d/1q6gC5nPrlcnaxl7WvHHVr-Yzu2xBZlIu4b1OcHdGdOc/edit#heading=h.mn9abqhwo6l5>

Extract data on server and transfer to desktop at location such as data-raw/extract/2019-10-14\_5000/best-fits, Then extract and save full swp, btran time-series to be used with the package. R/6.0\_ELM-FATES\_full best-fit\_swp.R

Plot the entire time-series of SWC, SWP and ELM-FATES generated BTRAN for best-fit simulations. Code/6.1\_Plotting best-fits\_full.R
