# Overview of different scripts used for PWG WRF data analysis in Iversen et al. (2022)

## Scripts to create necessary intermediate files 

Based on monthly wrfout files with hourly resolution
The wrfout files are compressed and temp is added

* create\_mnt\_means\_chunking.py: Creates monthly mean or total files for chosen variables.
* diff\_warm\_past.py: Calculates difference between past and warm period (season mean and monthly), based on files created with the above script. 
* cut\_out\_tile.py: Cuts out a tile from wrfout files based on grid indices, and a time period.
* interp\_to\_P.py: Interpolates wrf data from model levels to presure levels. Avoid large data files by reading in a tile created with the above script. 
* create\_mnt\_means\_plevs.py: Creates monthly mean or total files for chosen variables, for wrfout files which are interpolated to pressure levels. 
* create\_lat\_mean\_plevs.py: Create lattitudinally averaged files for chosen variables based on a tile defined by grid indices

## Scripts for plotting:
* plot\_Xvars\_Yyears.py (FIG 1): Subplot figure where vars on x axis, years on y axis, inkl windrose and NAO index.
* plot\_vertcross3.py (FIG 2): Plot vertical cross sections (longitude - pressure levels) for chosen lattitudinally averaged (over a tile) variables and winter season
* plot\_vertcross\_case\_PR.py (FIG 3): Same as above, but for a case period, including process rates extracted from Thompson&Eidhammer microphysics scheme from WRF.
* plot\_de\_desi\_RH.py (FIG 4): Plot change in e (vapor pressure) and esi (saturation vapor pressure wrt. ice) for given changes in water vapor mixing ratio (Q) and temperature (T)
