# CRT23
Codes associated with Cattiaux, Ribes and Thompson (2023), Searching for the most extreme temperature events in recent history, submitted to BAMS.

Details:

1. Running the scanning procedure.

scan_mf.R and scan_era5.R are the executable scripts that get the data, perform the analysis and save the results, resp. for Meteo-France (mf) and ERA5 data. Results are stored as "scan_1d" objects, i.e. a list containing all information of a scan at one location with one method (data, trend, fits, p1, etc.).

scan_source.R is sourced by scan_mf.R or scan_era5.R. It contains generic information for the scan: required packages, input / output directories, and a few useful function (e.g. for the computation of p1).

rbase.R is sourced by scan_source.R. It contains basic functions, e.g. the treatment of netcdf files (function myno() imports a netcdf in a list, etc.) or the treatment of time/dates, etc. It also loads several packages that the user will need to install.

wraf_source.R is sourced by scan_source.R. It contains generic information for dealing with WRAF regions (IDs, names..) and functions to get data or families of regions.

forced.response.R is sourced by scan_source.R. It contains functions used for the detrending procedure. See Rigal et al. (2019) and Ribes et al. (2022) for further details.

scan_1d_exe.R is a subroutine of scan_mf.R and scan_era5.R. It performs the scan at one location for one method and returns a scan_1d object.

scan_1d_sub.R is sourced by scan_mf.R or scan_era5.R. It contains functions useful to scan one location, eg. compute.stat.1d() that generates the nday x nduration matrix of p1, p0, etc. values.

2. Visualizing the results.

scan_functions.R contains functions to get scan.1d objects and select the most extreme events (chronology at one location).

scan_figures.R contains functions to plot the results.

CRT23.R is an executable script that reproduces tables and figures of the paper.
