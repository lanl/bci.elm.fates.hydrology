Scripts in this repo are a part of the workflow associated with the
following manuscript:

Chitra-Tarak, R, C Xu, S Aguilar, K Anderson-Teixeira, J Chambers, M
Detto, B Faybishenko, RA Fisher, R Knox, C Koven, L Kueppers, N Kunert,
SJ Kupers, NG McDowell, BD Newman, SR Paton, R Pérez, L Ruiz, L Sack, JM
Warren, BT Wolfe, C Wright, SJ Wright, J Zailaa, SM McMahon (2021)
Hydraulically vulnerable trees survive on deep-water access during
droughts in a tropical forest. New Phytologist. https:
//doi.org/10.1111/nph.17464

For the full work-flow and output datasets associated with the above
manuscript see the following archived dataset:

Chitra-Tarak R, Xu C, Aguilar S, Anderson-Teixeira K, Chambers J, Detto
M, Faybishenko B, Fisher R, Knox R, Koven C et al. 2020. Soil water
potentials (1990–2018) from a calibrated ELM-FATES, and rooting depth
analyses scripts, PA-BCI, Panama. 2.0. NGEE Tropics Data Collection.
doi: 10.15486/ngt/1696806

## Workflow Part II: ELM-FATES calibration for Soil Water Potentials

### 1. Check forcing data

`code/1.0_Checking_met_data.R`

### 2. Generate parameter files

1.  FATES param files dependencies:

    -   Parameter range for roob\_b parameters are generated here:
        2.0\_root\_parameters\_b.R
    -   Can verify that different root a can give same b:
        data-raw/estimation of root\_b\_par\_same\_b\_for\_different
        as.xlsx Other ranges are based on literature, or Zailaa et
        al. 2020

2.  ELM param files dependencies: fpi\_max ranges are generated from
    here:

    -   First throughfall from Zimmermann et al and rain data from STRI
        at Lutz are linked: `code/3.0_Throughfall.R`
    -   As rainfall data is not given in Zimmermann and throughfall data
        is collected for rainfall events (duration, amount and gap
        between events), rainfall that composes those events are
        inferred in Excel. Then range of interception is found for
        events &gt; 10 mm:
        data-raw/throughfall\_rain\_hourly\_for\_inspection\_new.xlsx

3.  Surface data files dependencies: Soil texture data and soil organic
    content data is gathered and defined here:
    `code/4.0_Surfdata_texture_options.R` But these are not used by the
    model as it is forced with

    1.  soil characteristic curves (Soil Water Content to Soil Water
        Potential) using Stephan Kupers’ data, as well as,
        bci.hydromet/data-raw/soil\_retention\_curves\_stephan.R

    2.  soil hydraulic conductivity, using Godsey & Stallard et al
        defined here:
        `code/5.0_Surfdata_Ksat_obs_and_bootstrapped_param.R`

4.  Surface data for ELM, ELM parameters and FATES paramater files
    `code/6.0_Generate_parameter files.R` These are transferred to the
    server.

### 3. Run simulations

1.  Follow Google Doc 3.0\_Running ELM-FATES with parameter ensembles @
    Shared drive/Rutuja\_work/web\_only/
    <https://docs.google.com/document/d/1qqbkQGHMG8BMfrUejUlkBU3PvUqPuW0Pc6CBLBQzZZk/edit#heading=h.gjdgxs>

2.  Thus run ensemble members and extract simulations on the server with
    `code/8.0_main.R` (which in turn sources `code/7.0_fun_extract.R`)
    and transfer back to the desktop storing at a specific dated
    location such as data-raw/extract/2019-10-14\_5000

### 4. Calibration and Sensitivity Analyses

1.  Use `code/09.0_ELM-FATES_output.R` to generate RMSE or Rsq values
    between simulations and hydrological observations, plot best-fits
    for individual fluxes and run sensitivity analyses. These outputs
    are compiled in individual-flux-best-fits.html and Report.html

2.  Use `code/10.0_ELM-FATES_params_bestfit.R` to generate objective
    function to choose best-fit simulations that fits all individual
    fluxes well (QRUNOFF, AET, Soil moisture by depth)

3.  Confirm that individual fluxes are well captured by thus chosen
    best-fits: `code/11.0_ELM-FATES_output_bestfit.R`

### 4. Re-Run simulations for best-fits for the entire study period

Now for the chosen best-fits re-run simulations covering the entire
time-period of 1985-2018 Follow Google Doc 4.0\_Running ELM-FATES for
best-fit parameter ensembles @ Shared drive/Rutuja\_work/web\_only/
<https://docs.google.com/document/d/1V2UK_iSdXmkq3jR3TyowW7YllvrC7OVi/edit#heading=h.gjdgxs>
(and `code/12.0_Prep_for_Best-fit_case_runs.R` to generate par.sam
numbers seperated by commas)

### 5. Extract, plot and save Soil Water Potentials as a data-package

1.  Extract data from server and transfer to desktop at the following
    location `data-raw/extract/2019-10-14_5000/best-fits`,

2.  Then extract and save full swp, btran time-series to be used for
    generating the data package.
    `code/13.0_ELM-FATES_full best-fit_swp.R`

3.  Plot the entire time-series of SWC, SWP and ELM-FATES generated
    BTRAN for best-fit simulations.
    `code/14.0_Plotting best-fits_full.R`

## Copyrights

© 2021. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract
89233218CNA000001 for Los Alamos National Laboratory (LANL), which is
operated by Triad National Security, LLC for the U.S. Department of
Energy/National Nuclear Security Administration. All rights in the
program are reserved by Triad National Security, LLC, and the U.S.
Department of Energy/National Nuclear Security Administration. The
Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to
reproduce, prepare derivative works, distribute copies to the public,
perform publicly and display publicly, and to permit others to do so.

## Open Source Redistribution License

This program is open source under the BSD-3 License.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS
IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    groundhog.day = "2021-01-01"
    groundhog.library('rmarkdown', groundhog.day)
    rmarkdown::render("Readme.rmd", output_format = "pdf_document")
    rmarkdown::render("Readme_data.rmd", output_format = "pdf_document")
