# Create lattitudinally averaged files for chosen variables based on a tile defined by grid indices

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

y = sys.argv[1]
y = int(y)
period = sys.argv[2]

plot = False
#matplotlib.use('Agg')  # To not open a plot window
domain = 'd02'
#varlist = list(['TEMPERATURE','RR','RH','QCLOUD','QRAIN','QSNOW','QGRAUP','QVAPOR','zerocross'])  
varlist = list(['QCLOUD','QSNOW','QGRAUP','QICE'])
#varlist = list(['TEMPERATURE'])
m1 = '2' # '1[%s]' %mnt1
m2 = '12' 
fix = 'DJF'
X = [20,185]; Y = [200,285] #For transect mean

#Dirs
wrfdir = '/work/users/emiive/PGW/plevfiles'
outdir = ('./files/%s' %period) 
figdir = ('./figures')


########### The main routine ##################

listOfFiles = os.listdir(wrfdir)
listOfFiles.sort() 
for v in varlist:
  i=0  
  for entry in listOfFiles:
      if fnmatch.fnmatch(entry,'%i_1[%s]_%s_%s_plevs.nc' %(y,m1,period,v)) \
      or fnmatch.fnmatch(entry,'%i_0[%s]_%s_%s_plevs.nc' %(y+1,m2,period,v)):
         i=i+1
         wrf = xr.open_mfdataset('%s/%s' %(wrfdir,entry))
                 #chunks={'south_north': 100, 'west_east': 100})

         tile = wrf.sel(south_north=slice(Y[0],Y[1]),west_east=slice(X[0],X[1]))
         dslat1 = tile.mean('south_north')

         if i==1:
            dslat = dslat1
         else:
            dslat = xr.concat([dslat,dslat1], dim="times") 

  if v == 'TEMPERATURE':
      dslat['TEMPERATURE'] = dslat['TEMPERATURE']-273
      V = 'T'
  else:
      V = v

  delayed_obj = dslat.to_netcdf('%s/%i_%s_%s_latmean_plevs_%s.nc' %(outdir,y,period,V,fix), compute=False)
  with ProgressBar():
      results = delayed_obj.compute()

