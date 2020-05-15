#' @title Soil Water Potential
#' @description Daily Soil Water Potential by depth and parameter realisalition
#' @format A data frame with 10022711 rows and 4 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{depth}}{double COLUMN_DESCRIPTION}
#'   \item{\code{psi}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION} 
#'}
#' @details SWP (MPa)
"psi"
#' @title Mean Soil Water Potential across ensembles by depth
#' @description Mean of Daily Soil Water Potential by depth across parameter ensemble
#' @format A data frame with 137605 rows and 3 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{depth}}{double COLUMN_DESCRIPTION}
#'   \item{\code{psi}}{double COLUMN_DESCRIPTION} 
#'}
#' @details SWP (MPa)
"psi.mean"
#' @title Growth Factor
#' @description Growth factor (0-1, unitless) from Soil Water Potential by depth and par.sam
#' @format A data frame with 10022711 rows and 6 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{depth}}{double COLUMN_DESCRIPTION}
#'   \item{\code{psi}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{psi.2}}{double COLUMN_DESCRIPTION}
#'   \item{\code{gfac}}{double COLUMN_DESCRIPTION} 
#'}
#' @details #' Growth Factor (0-1) from swp
#' From 0 to -0.5MPa f(swp) = 1, 
#' f(swp) = linearly decreases from 0 to 1 f(swp) =  1/(2.42-0.5)*2.42 - 1/(2.42-0.5)*swp
#' f(swp) =  1.260417 - 0.5208333*swp
#' such that at -2.42 MPa f(swp) = 0, the mean TLP in Nobby Kunert's data
"swp.gfac"
#' @title btran
#' @description daily btran (unitless) by depth and par.sam
#' @format A data frame with 899725 rows and 3 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{value}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION} 
#'}
#' @details Unitless
"btran"
#' @title btran mean and 95 CI
#' @description daily mean and CI btran by depth and par.sam
#' @format A data frame with 10585 rows and 4 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{mean}}{double COLUMN_DESCRIPTION}
#'   \item{\code{upper.CI}}{double COLUMN_DESCRIPTION}
#'   \item{\code{lower.CI}}{double COLUMN_DESCRIPTION} 
#'}
#' @details Btran Mean and 95% Confidence Interval (unitless)
"btran.stat"
#' @title Plant Available Water
#' @description Daily Plant Available Water (mm) by depth and par.sam swc - swc_wilt
#' @format A data frame with 10023995 rows and 4 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{depth}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION}
#'   \item{\code{paw}}{double COLUMN_DESCRIPTION} 
#'}
#' @details Plant Available Water (V/V) = swc - swc_wilt (0.17 volumetric) 
#' swc_wilt = water content at wilting point, using Clapp & Hornberger eqn 1978 in FATES
"paw"
#' @title Soil Water Content
#' @description Daily Soil Water Content(V/V) by depth and parameter realisalition
#' @format A data frame with 10023995 rows and 4 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{depth}}{double COLUMN_DESCRIPTION}
#'   \item{\code{swc}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION} 
#'}
#' @details soil water content (V/V)
"swc"
#' @title GPP
#' @description Daily GPP by depth and parameter realisalition
#' @format A data frame with 899725 rows and 3 variables:
#' \describe{
#'   \item{\code{date}}{double COLUMN_DESCRIPTION}
#'   \item{\code{value}}{double COLUMN_DESCRIPTION}
#'   \item{\code{par.sam}}{integer COLUMN_DESCRIPTION} 
#'}
#' @details GPP (gC/m^2/d)
"gpp"