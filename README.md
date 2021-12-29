## Angular Momentum Budget in Spherical Pressure Coordinates from Reanalysis datasets
Script and documentation for deriving hourly Angular Momentum Budgets from reanalysis data on pressure levels.

please cite as:

creator: mhell@ucsd.edu

needed python modules:
glob, xarray, time, shutil, numpy, cdsapi 
# configure cdsapi credentials to download ERA5 data:
https://cds.climate.copernicus.eu/api-how-to 

# File description
## technical_notes.pdf
Derivation of the budget in Pressure coordinates and description of the in and output files:
https://github.com/mochell/AM_budget_pressure_ERA5/blob/main/technical_notes.pdf


## A02_download_file.py
example file of how to download data from the ECMWF server

## A03_cal_budget_worker_optimized.py
calculated the budget as described in technical_notes.pdf
needs:
### module file
tools_AM_budget.py
### contains various constants
config.json        

## example output:
exmpl_repres_dmean_ps_iews_1980-01-01.pdf

exmpl_budget_dmean_1980-01-01.pdf

worker_log_num250_1980-01-01.hist.txt
