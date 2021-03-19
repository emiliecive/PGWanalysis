# Interpolate wrf data from model levels to presure levels
# Will create very large data files, so smart to cut out the relevant area first using
# cut_out_tile.py

import sys
import numpy as np
import pandas as pd
import xarray as xr
from metpy.interpolate import log_interpolate_1d
from metpy.units import units
import os, fnmatch
from resource import *
#import matplotlib.pyplot as plt
#from cartopy import crs
#import cartopy.feature as cfeature
#from wrf import get_cartopy
#from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
#                 cartopy_ylim, latlon_coords, ll_to_xy, CoordPair)

domain = 'd02'
wrfdir = '/work/users/emiive/PGW/tilefiles'
outdir = '/work/users/emiive/PGW/plevfiles'

def interp_p(plevs,y,m,period,v):
   """
   Interpolates one variable and saves to file
   """

   #if period == 'warm':
   #   dirr = '/mnt/elephant/WRFsimulations/ICEBOX_PGW'
   #elif period == 'past':
   #   dirr = '/mnt/elephant/WRFsimulations/KVT_NSF02'
   #wrfdir = ('%s/archive2/postpost' %dirr)

   #filename = ('wrfout_%s_%s-%s.nc4' %(domain,y,m))
   filename = ('Tile_%s_%s_%s.nc' %(period,y,m)) 
   wrf = xr.open_dataset('%s/%s' %(wrfdir,filename))
   wrf = wrf.sel(bottom_top=slice(0,10),bottom_top_stag=slice(0,10))
   P = wrf['PRESSURE']

   var_p = np.empty([P.shape[0],len(plevs),P.shape[-2],P.shape[-1]])

   for i in range(len(wrf.Time)): # DIVIDE AND CONQUER ON TIME DIM
       print(i)
       Pi = units.Quantity(wrf['PRESSURE'][i].values, 'Pa')
       P_int = np.expand_dims(Pi,axis=0)
       var = wrf[v][i].values
       var_int = np.expand_dims(var,axis=0)
       var_p[i] = log_interpolate_1d(plevs, P_int, var_int, axis=1)

   ds = xr.Dataset({v: (['times','plevs','south_north','west_east'], var_p),},
           coords={'times': (['times'], wrf['times'].values),
                    'plevs': (['plevs'], plevs),
                    'XLONG': (['south_north','west_east'], wrf.XLONG.values),
                    'XLAT': (['south_north','west_east'], wrf.XLAT.values),},)

   print('saving ds') 
   ds.to_netcdf('%s/%s_%s_%s_%s_plevs.nc' %(outdir,y,m,period,v))


def interp_p_varlist(plevs,y,m,period,varlist,fix):
   """
   Interpolates a list of variables and saves to one file
   """

   filename = ('Tile_%s_%s_%s*%s.nc' %(period,y,m,fix))
   wrf = xr.open_mfdataset('%s/%s' %(wrfdir,filename))
   P = wrf['PRESSURE']

   ds = xr.Dataset()

   for v in varlist:
       var_p = np.empty([P.shape[0],len(plevs),P.shape[-2],P.shape[-1]])
       for i in range(len(wrf.Time)): # DIVIDE AND CONQUER ON TIME DIM
           print(i)
           Pi = units.Quantity(wrf['PRESSURE'][i].values, 'Pa')
           P_int = np.expand_dims(Pi,axis=0)
           var = wrf[v][i].values
           var_int = np.expand_dims(var,axis=0)
           var_p[i] = log_interpolate_1d(plevs, P_int, var_int, axis=1)

       ds[v] = xr.DataArray(var_p,
                    coords={'times': (['times'], wrf['times'].values),
                    'plevs': (['plevs'], plevs),
                    'XLONG': (['south_north','west_east'], wrf.XLONG.values),
                    'XLAT': (['south_north','west_east'], wrf.XLAT.values),},
                    dims=['times','plevs','south_north','west_east'],)

   print('saving ds')
   ds.to_netcdf('%s/%s_%s_%s_%s_plevs.nc' %(outdir,y,m,period,fix))
