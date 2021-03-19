# Plot a vertical cross section (longitude - pressure levels) for a lattitudinally averaged tile and a winter season
# subplots for chosen variables, a column for past and one for diff (warm - past)

import sys
import numpy as np
import xarray as xr
import os, fnmatch
from resource import *
from matplotlib import cm
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from cartopy import crs
import geopandas as gpd
#matplotlib.use('Agg')  # To not open a plot window
#import calendar

#Debug
from IPython import embed

domain = 'd02'
varlist = list(['T','QVAPOR','QICE','QSNOW','QGRAUP','QCLOUD','QRAIN','WC','RR'])
y = sys.argv[1]
y = int(y)
months = [12,1,2]
fix = 'DJF'
savefig = 1
perc = 0.6 # Level for percentile
X = [20,185]; Y = [200,285] #For transect mean
periods = ['past','warm']

#Dirs
dirr = '/mnt/elephant/WRFsimulations/ICEBOX_PGW'
wrfdir = ('%s/archive2/postpost' %dirr)
datadir = '/work/users/emiive/PGW/plevfiles'
outdir = ('../files')
figdir = ('./figures')

data_proj = crs.Stereographic()
map_proj = crs.PlateCarree()
boundary = gpd.read_file('../boundaryshp/ne_50m_admin_0_boundary_lines_land.shp')
HGT = xr.open_dataset('%s/wrfout_d02_2019-10.nc4' %wrfdir)
HGT = HGT['HGT']


def plot_tile(HGT,x,y):

    xlen = np.arange(X[0],X[-1]+1)
    ylen = np.arange(Y[0],Y[-1]+1)
    plt.figure(figsize=(7,8))
    hgt = HGT.plot(cmap='gist_earth_r',vmin=0,vmax=1800)
    plt.plot(xlen, np.full((len(xlen)),Y[0]), color='black')
    plt.plot(xlen, np.full((len(xlen)),Y[1]), color='black')
    plt.plot(np.full((len(ylen)),X[0]), ylen, color='black')
    plt.plot(np.full((len(ylen)),X[1]), ylen, color='black')
    plt.show(block=False)


# Get the variables
Vars = {}
for period in periods:
    Vars[period] = {}
    Qds = xr.Dataset()

    for v in varlist:

        if v == 'RR':
            ds = xr.open_dataset('%s/%s/%i_%s_%s_sums.nc' %(outdir,period,y,period,v))
            da = ds['RAINNC'].sel(month=months,south_north=slice(Y[0],Y[1]),west_east=slice(X[0],X[1]))
            da = da.sum('month')

        elif v != 'WC':
            ds = xr.open_mfdataset('%s/*_%s_%s_plevs.nc' %(datadir,period,v))
            ds['plevs'] = ds['plevs']/100 # hPa

            if 'Q' in v:
                da = ds[v].sum('times')
                da = da*1000 # Convert to g/kg

                # Latmean_plev files for calc. of warm clouds (WC)
                try:
                    ds = xr.open_dataset('%s/%s/%i_%s_%s_latmean_plevs_%s.nc' %(outdir,period,y,period,v,fix))
                    ds['plevs'] = ds['plevs']/100
                    Qds[v] = ds[v]*1000 # g/kg
                except:
                    pass
            elif v == 'T':
                da = ds['TEMPERATURE']-273.15
                da = da.mean('times')
                da.rename('T')
            #elif 'RH' in v:
            #    da = ds[v].where(ds[v] > 100).count('times') / len(ds[v].times)
            #    da = da.where(da > 0, drop=True)
            else:
                da = ds[v].mean('times')

        if v != 'WC':
            dalat = da.mean('south_north',skipna=True)
            Vars[period][v] = dalat

    if 'WC' in varlist:
        # Calculate frequency of warm clouds (WC)
        clouds = Qds['QCLOUD'].where(Qds['QCLOUD'] > 0).count('times')
        Qds = Qds.where(Qds['QSNOW'] < 1e-5)
        Qds = Qds.where(Qds['QGRAUP'] < 1e-5)
        Qds = Qds.where(Qds['QICE'] < 1e-7)
        Qda = Qds['QCLOUD']
        Vars[period]['WC'] = Qda.where(Qda > 0).count('times') / clouds 


# Get lon-labels for x-axis
lontile = ds['XLONG'] #.sel(south_north=slice(Y[0],Y[1]),west_east=slice(X[0],X[1]))
xlon = lontile[np.round(len(lontile.south_north)/2).astype(int)] # Choosing middel tile lat for xlon

# Calc. diff warm - past
diff = {}
for v in varlist:
    diff[v] = Vars['warm'][v] - Vars['past'][v]
    diff[v] = diff[v].where(diff[v] != 0) 

# Set plotting limits
cmap = {}; vmin = {}; vmax = {}; levels = {}; extend = {}; lab={}; diffmin={}; diffmax={}
for v in varlist:
    cmap[v] = 'gist_earth_r'
    vmin[v] = 0
    levels[v] = 11
    lab[v]='g/kg'
vmin['QICE']=0; vmax['QICE']=0.01; levels['QICE']=21;
diffmin['QICE']=-0.003; diffmax['QICE']=0.003
cmap['WC']='YlOrRd'; vmin['WC']=0; vmax['WC']=1; levels['WC']=11; extend['WC']='neither'; lab['WC']='f'
diffmin['WC'] = -0.25; diffmax['WC'] = 0.25
cmap['T'] = 'bwr'; vmin['T'] = -10; vmax['T'] = 10; lab['T'] = 'C'
diffmin['T'] = -2.5; diffmax['T'] = 2.5
cmap['W'] = 'bwr_r'; vmin['W'] = -1.2; vmax['W'] = 1.2
#cmap['RHw'] = 'YlOrRd';vmin['RHw'] = 0; vmax['RHw'] = 0.5; lab['RHw'] = 'f'
#diffmin['RHw'] = -0.15; diffmax['RHw'] = 0.15
#cmap['RHi'] = 'YlOrRd'; vmin['RHi'] = 0; vmax['RHi'] = 0.5; lab['RHi'] = 'f'
#diffmin['RHi'] = -0.15; diffmax['RHi'] = 0.15
vmax['RHw'] = 100; lab['RHw'] = '%'
diffmin['RHw'] = -6; diffmax['RHw'] = 6
vmax['RHi'] = 100; lab['RHi'] = '%'
diffmin['RHi'] = -6; diffmax['RHi'] = 6

##################################################
# FIGURE
##################################################
axlen = len(varlist)+1
hrat = [1]*len(varlist)
hrat.append(0.5)
fig, ax = plt.subplots(ncols=2,nrows=axlen,sharex=True,figsize=(10,13),\
        gridspec_kw={'height_ratios': hrat, 'hspace': 0.3}) #, 'top':0.9, 'bottom':0.1})

for j in range(2):
    for v,i in zip(varlist,range(len(varlist))):
    
      if j==0:
          latmean = Vars['past'][v]
      else:
          latmean = diff[v]
          cmap[v] = 'bwr_r';  cmap['T'] = 'bwr'; cmap['WC'] = 'bwr'

      if v == 'RR':
          latmean.plot(ax=ax[i,j])
          pos = ax[0,j].get_position()
          pos2 = ax[i,j].get_position()
          ax[i,j].set_position([pos.x0,pos2.y0,pos.width,pos2.height])
          ax[i,0].set_ylabel('mm')
          if j==1:
              ax[i,j].set_ylim([-250,50])
      else:
          if (j==0) and (v in ['T','RHw','RHi','QNICE']):
              latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],\
                                    vmin=vmin[v],vmax=vmax[v],levels=levels[v],\
                                    cbar_kwargs={'label':lab[v]})
          elif (j==1) and (v in ['T','RHw','RHi','QICE','WC']):
              latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],levels=levels[v],\
                                    vmin=diffmin[v],vmax=diffmax[v],\
                                    cbar_kwargs={'label':lab[v]})
          else:
              latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],levels=levels[v],\
                                    cbar_kwargs={'label':lab[v]})
          ax[i,0].set_ylabel('hPa')
          ax[i,j].invert_yaxis()
          ax[i,1].set_yticklabels('')
          
      ax[i,j].grid()
      ax[i,j].set_xlabel('')
      ax[i,1].set_ylabel('')
      ax[i,j].set_title(v)

for j in range(2):
    hgt = HGT.sel(south_north=slice(Y[0],Y[1]),west_east=slice(X[0],X[1]))
    hgt = hgt.mean('south_north')
    ax[-1,j].plot(hgt,label='HGT')
    pos = ax[0,j].get_position()
    pos2 = ax[-1,j].get_position()
    ax[-1,j].set_position([pos.x0,pos2.y0,pos.width,pos2.height])
    #ax[-1,j].legend()
    ax[-1,j].set_xticks(np.arange(0,175,25))
    ax[-1,j].set_xticklabels(xlon[np.arange(0,175,25)].values.round(1))
    ax[-1,j].grid()
    ax[-1,j].set_xlabel('lon')
    ax[-1,j].set_title('HGT')
    ax[-1,0].set_ylabel('m')
    ax[-1,1].set_yticklabels('')

#plt.subplots_adjust(hspace=0.11)

if savefig:
    plt.savefig('%s/%i_vertcross3_warm_past_%s.png' %(figdir,y,fix), 
                bbox_inches='tight', dpi=120)
else:
    plt.show(block=False)

