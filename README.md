## Angular Momentum Budget in Spherical Pressure Coordinates from Reanalysis datasets
Example script and documentation for deriving hourly Angular Momentum Budgets from reanalysis data on pressure levels

# File description
## technical_notes.pdf
Derivation of the budget in Pressure coordinates and description of the in and output files:
[embed] https://github.com/mochell/AM_budget_pressure_ERA5/blob/main/technical_notes.pdf [/embed]

<object data="https://github.com/mochell/AM_budget_pressure_ERA5/blob/main/technical_notes.pdf" type="application/pdf" width="700px" height="700px">
    <embed src="https://github.com/mochell/AM_budget_pressure_ERA5/blob/main/technical_notes.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="http://yoursite.com/the.pdf">Download PDF</a>.</p>
    </embed>
</object>


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
