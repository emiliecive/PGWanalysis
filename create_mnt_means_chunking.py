# Create yearly and monthly mean files for chosen vars

import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr
import dask
from dask.diagnostics import ProgressBar
import os, fnmatch
from resource import *
#import matplotlib.pyplot as plt
import calendar
import plot_maps3

dask.config.set(**{'array.slicing.split_large_chunks': False})

y = sys.argv[1]
y = int(y)
period = sys.argv[2]

plot = False
#matplotlib.use('Agg')  # To not open a plot window
domain = 'd02'
#varlist = list(['T','RR','RH','QCLOUD','QRAIN','QSNOW','QGRAUP','QVAPOR','zerocross'])  
#varlist = list(['QCLOUD','QRAIN','QSNOW','QGRAUP'])
#varlist = list(['SR','rr','RH','QVAPOR'])
varlist = list(['QVAPOR'])
m1 = '012' # '1[%s]' %mnt1
m2 = '1234' 

#Dirs
if period == 'warm':
   dirr = '/mnt/elephant/WRFsimulations/ICEBOX_PGW'
elif period == 'past':
   dirr = '/mnt/elephant/WRFsimulations/KVT_NSF02'
wrfdir = ('%s/archive2/postpost' %dirr)
outdir = ('./files/%s' %period) 
figdir = ('./figures')


########### The main routine ##################

listOfFiles = os.listdir(wrfdir)
listOfFiles.sort() 
i=0
for entry in listOfFiles:
   if fnmatch.fnmatch(entry,'wrfout_%s_%i-1[%s].nc4' %(domain,y,m1)) \
   or fnmatch.fnmatch(entry,'wrfout_%s_%i-0[%s].nc4' %(domain,y+1,m2)):
      i=i+1
      wrf1 = xr.open_mfdataset('%s/%s' %(wrfdir,entry), drop_variables={'T','XTIME'},\
              chunks={'Time':720, 'bottom_top':10, 'south_north': 32, 'west_east': 32,\
              'soil_layers_stag':4, 'bottom_top_stag':10, 'west_east_stag':32, 'south_north_stag':32})
      wrf1 = wrf1.sel(bottom_top=slice(0,10),bottom_top_stag=slice(0,10))  
      if i==1:
          wrf = wrf1
      else:
          wrf = xr.concat([wrf,wrf1], dim="Time") 

wrf = wrf.assign_coords(Time=wrf.times)

for v in varlist:
  if v == 'T':
    T = wrf['TEMPERATURE']-273
    T.attrs['units'] = 'C'
    T = T.groupby('Time.month').mean('Time')
    #T.to_netcdf('%s/%i_%s_T_mean.nc' %(outdir,y,period))
    delayed_obj = T.to_netcdf('%s/%i_%s_T_mean.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'RR':
    RR = wrf['RAINNC'].diff("Time")  # Get hourly precip rate from accumulated time series
    RR = RR.groupby('Time.month').sum('Time')
    delayed_obj = RR.to_netcdf('%s/%i_%s_RR_sums.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'SR':
    SR = wrf['SNOWNC'].diff("Time")  # Get hourly precip rate from accumulated time series
    SR = SR.groupby('Time.month').sum('Time')
    delayed_obj = SR.to_netcdf('%s/%i_%s_SR_sums.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'FF':
    FF = wrf['FF']
    FF = FF.groupby('Time.month').mean('Time')
    delayed_obj = FF.to_netcdf('%s/%i_%s_FF_mean.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'RH':
    T = wrf['TEMPERATURE']
    P = wrf['PRESSURE']
    q = wrf['QVAPOR']
    RH = (0.263*P*q)/(np.exp((17.67*(T-273.15))/(T-29.65)))
    RH.attrs['description'] = 'Relative humidity'
    RH.attrs['units'] = '%'
    RH = RH.groupby('Time.month').mean('Time')
    delayed_obj = RH.to_netcdf('%s/%i_%s_RH_mean.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if 'Q' in v:
    #qvar = cloud_water(wrf,v)  # Convert to kg/m3
    qvar = wrf[v]
    qvar = qvar.groupby('Time.month').sum('Time')
    delayed_obj = qvar.to_netcdf('%s/%i_%s_%s_sums.nc' %(outdir,y,period,v),compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'P':
    P = wrf['PRESSURE']
    P = P.groupby('Time.month').mean('Time')
    delayed_obj = P.to_netcdf('%s/%i_%s_P_mean.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'Z':
    try:
        Z = wrf['Z']
    except:
        g = 9.81
        P = wrf['PRESSURE']
        T = wrf['TEMPERATURE']
        P0 = 10132500
        R = 8.3143
        M = 0.02896
        Z = -((R*T)/(M*g))*np.exp(P/P0)
        Z.attrs['description'] = 'Height'
    Z.attrs['units'] = 'm'
    Z = Z.groupby('Time.month').mean('Time')
    delayed_obj = Z.to_netcdf('%s/%i_%s_H_mean.nc' %(outdir,y,period), compute=False)
    with ProgressBar():
        results = delayed_obj.compute()

  if v == 'zerocross':
    T = wrf['TEMPERATURE'][:,0]-273
    Tdiff = T.diff("Time",label='upper')
    T = T[1:]
    da = xr.full_like(T,fill_value=0)
    mask = ((Tdiff>0) & (T==0))
    da = xr.where(mask, 1, da)
    zerocross = da.groupby('Time.month').sum('Time')
    delayed_obj = zerocross.to_netcdf('%s/%i_%s_%s_sums.nc' %(outdir,y,period,v),compute=False)
    with ProgressBar():
        results = delayed_obj.compute()


##### Plotting #######
if plot:
    plot_maps3.plot_mnt_tots(varlist,y,period)
