# Cuts out a tile from wrfout files based on grid indices, and a time period

import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
import os, fnmatch
from resource import *
import matplotlib.pyplot as plt
import calendar

dask.config.set(**{'array.slicing.split_large_chunks': False})

y = 2015
m = '02'
d = 1; dd1 = 7; dd2 = 9
t_from = '2015-02-17 18'
t_to = '2015-02-19 16'
#X = [20,185]; Y = [200,285] #Stølsheimen
X = [140,220]; Y = [427,460] # Trondelag
periods = ['past','warm']
xychunk = 16
#fix = 'Stlshmn' 
fix = 'Trn'

#Dirs
dirs = {}
#dirs['past'] = '/mnt/elephant/WRFsimulations/KVT_NSF02/rawdata_feb2015'
#dirs['past'] = '/mnt/elephant/WRFsimulations/KVT_NSF02/archive2/postpost'
dirs['past'] = '/mnt/elephant/WRFsimulations/ICEBOX_PGW/archive2/20150216/past'
#dirs['warm'] = '/mnt/elephant/WRFsimulations/ICEBOX_PGW/archive2'
#dirs['warm'] = '/mnt/elephant/WRFsimulations/ICEBOX_PGW/archive2/postpost'
dirs['warm'] = '/mnt/elephant/WRFsimulations/ICEBOX_PGW/archive2/20150216/warm'
outdir = '/work/users/emiive/PGW/tilefiles'


for period in periods:

    wrf = xr.open_mfdataset('%s/wrfout_d02_%i-%s-%i[%i-%i]*.nc4' %(dirs[period],y,m,d,dd1,dd2),\
            chunks={'Time':720, 'west_east':xychunk, 'south_north':xychunk, \
              'soil_layers_stag':4, 'west_east_stag':xychunk, 'south_north_stag':xychunk})

    wrf = wrf.assign_coords(Time=wrf.times)
    wrf = wrf.sel(Time=slice(t_from, t_to))
    #wrf = wrf.sel(bottom_top=slice(0,10), bottom_top_stag=slice(0,11))

    tile = wrf.sel(south_north=slice(Y[0],Y[1]),west_east=slice(X[0],X[1]),\
                   south_north_stag=slice(Y[0],Y[1]),west_east_stag=slice(X[0],X[1]))
    
    if d:
        delayed_obj = tile.to_netcdf('%s/Tile_%s_%i_%s_%i%i-%i%i_%s.nc' %(outdir,period,y,m,d,dd1,d,dd2,fix), compute=False)
    else:
        delayed_obj = tile.to_netcdf('%s/Tile_%s_%i_%s.nc' %(outdir,period,y,m), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

