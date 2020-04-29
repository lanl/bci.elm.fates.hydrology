Workflow
================
Rutuja Chitra-Tarak
12/2/2019

### 1. Check forcing data

`code/0.01_Checking_met_data.R`

### 2. Generate parameter files

1.  FATES param files dependencies: Parameter range for roob\_b parameters are generated here: 3.090\_root\_parameters\_b.R Can verify that different root a can give same b: data-raw/estimation of root\_b\_par\_same\_b\_for\_different as.xlsx Other ranges are based on literature or Kunert's data

2.  ELM param files dependencies: fpi\_max ranges are generated from here: \#\# First throughfall from Zimmermann et al and rain data from STRI at Lutz are linked: `code/0.03_Throughfall.R` \#\# As rainfall data is not given in Zimmermann and throughfall data is collected for rainfall events (duration, amount and gap between events), rainfall that composes those events are inferred in Excel. And range of interception found for events &gt; 10 mm: data-raw/throughfall\_rain\_hourly\_for\_inspection\_new.xlsx

3.  Surface data files dependencies: Soil testure data and soil organic content data is gathered and defined here: `code/3.2_Surfdata_texture_options.R` But these are not used by the model as it is forced with
    -   soil characteristic curves (SWC to Matric potential) using Stephan Kupers' data (Yilin's equation), as well as, bci.hydromet/data-raw/soil\_retention\_curves\_stephan.R
    -   soil hydraulic conductivity, using Godsey & Stallard et al defined here: `code/3.2.5_Surfdata_Ksat_obs_and_bootstrapped_param.R`
4.  Surface data for ELM, ELM parameters and FATES paramater files `code/3.1_Generate_parameter files.R` These are transferred to the server.

### 3. Run simulations

Follow Google Doc 3.0\_Running ELM-FATES with parameter ensembles @ Shared drives/Rutuja\_work/web\_only/ <https://docs.google.com/document/d/1qqbkQGHMG8BMfrUejUlkBU3PvUqPuW0Pc6CBLBQzZZk/edit#heading=h.gjdgxs>

Thus run ensembles and extract simulations on the server with `code/2.0_main.R` (which in turn sources `code/1.0_fun_extract.R`) and transfer back to the desktop storing at a specific dated location such as data-raw/extract/2019-10-14\_5000

### 4. Calibration and Sensitivity Analyses

1.  Use `code/4.1_ELM-FATES_output.R` to generate RMSE or Rsq values between simulations and hydrological observations, plot best-fits for individual fluxes and run sensitivity analyses. These outputs are compiled in individual-flux-best-fits.html and Report.html

2.  Outputs from trials such as `4.1_ELM-FATES_output_trial_leafage_fdrought.R, 4.1_ELM-FATES_output_trial2.R`, `4.1_ELM-FATES_output_trial_tlai_reduction.R`, are compiled in `leafage_fdrought.html`, `Effect-of-bbslope---VCmax.html` or `Max-Rooting-Depth-effect.html`

3.  Use `code/4.2_ELM-FATES_params_bestfit.R` to generate objective function to choose best-fit simulations that fits all individual fluxes well (QRUNOFF, AET, Soil moisture by depth)

4.  Confirm that individual fluxes are well captured by thus chosen best-fits: `code/4.3_ELM-FATES_output_bestfit.R`

### 4. Re-Run simulations for best-fits for the entire study period

Now for the chosen best-fits re-run simulations covering the entire time-period of 1985-2018 Follow Google Doc 4.0\_Running ELM-FATES for best-fit parameter (and `code/5.0_Prep_for_Best-fit_case_runs.R` to generate par.sam numbers seperated by commas) ensembles<https://docs.google.com/document/d/1q6gC5nPrlcnaxl7WvHHVr-Yzu2xBZlIu4b1OcHdGdOc/edit#heading=h.mn9abqhwo6l5>

### 5. Extract, plot and save Soil Water Potentials as a data-package

1.  Extract data from server and transfer to desktop at the following location `data-raw/extract/2019-10-14_5000/best-fits`,

2.  Then extract and save full swp, btran time-series to be used for generating the data package. `R/6.0_ELM-FATES_full best-fit_swp.R`

3.  Plot the entire time-series of SWC, SWP and ELM-FATES generated BTRAN for best-fit simulations. `code/6.1_Plotting best-fits_full.R`
