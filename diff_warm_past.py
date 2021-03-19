# Find difference between historic and warm period (season mean and monthly)
# Based on monthly mean files created with create_mnt_means_chunking.py  

import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch


#varlist = list(['T','RR','QSNOW','QCLOUD','QRAIN'])  
varlist = list(['QVAPOR'])
#varlist = list(['RR','SR'])
#years = [2009,2012,2014,2015,2019] 
years = [2012]
months = False #[11,12,1,2,3,4]
season = [12,1,2]
perc = 0
percdiffseason = 0
totaldiff = 0

# Dirs:
outdir = ('./files/')


for v in varlist:

   if (v == 'T') or (v == 'FF') or (v == 'RH'):
      stat = 'mean'
   else:
      stat = 'sums'

   i=0
   for y in years:

      i=i+1
      warm = xr.open_dataset('%s/warm/%s_warm_%s_%s.nc' %(outdir,y,v,stat))
      past = xr.open_dataset('%s/past/%s_past_%s_%s.nc' %(outdir,y,v,stat))

      if months:
          warm = warm.sel(month=months)
          past = past.sel(month=months)

      diff = warm - past
      diff.to_netcdf('%s/%i_%s_diff_warm_past.nc' %(outdir,y,v))
      if perc:
         percdiff = (diff/past)*100
         percdiff.to_netcdf('%s/%i_%s_percdiff_warm_past.nc' %(outdir,y,v))


      # Season mean/total:
      if percdiffseason or totaldiff:

         warm = warm.sel(month=season)
         past = past.sel(month=season)
         if (v == 'T') or (v == 'FF'):
             warm = warm.mean('month')
             past = past.mean('month')
         else:
             warm = warm.sum('month')
             past = past.sum('month')
      
      # Percentage change for season:
      if percdiffseason:
         diff = warm - past
         percdiff = (diff/past)*100
         percdiff.to_netcdf('%s/%i_%s_percdiff_warm_past_DJF.nc' %(outdir,y,v))

      # Difference as averaged over all years
      if totaldiff:
          if i==1:
              warm1 = warm
              past1 = past
          else:
              warm1 = xr.concat([warm1,warm],'year') 
              past1 = xr.concat([past1,past],'year')
          if i == len(years):
              warm = warm1.mean('year')
              past = past1.mean('year')
              diff = warm-past
              diff.to_netcdf('%s/%s_diff_warm_past.nc' %(outdir,v))
              percdiff = (diff/past)*100
              percdiff.to_netcdf('%s/%s_percdiff_warm_past.nc' %(outdir,v))
