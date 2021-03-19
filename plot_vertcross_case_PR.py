# Plot a vertical cross section (longitude - pressure levels) for a lattitudinally averaged tile and a case period
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

y = 2015
m = '02'
d1 = '17'; d2 = '19'
slice_time = 1
t_from = '2015-02-17 20'
t_to = '2015-02-18'
fix = 'Stlshmn'
#fix = 'Trn'
savefig = 1 
plev = 1
#varlist = list(['W','T','RHw','RHi','QVAPOR','QNICE','QSNOW','QGRAUP','QCLOUD','QRAIN','RR'])
varlist = list(['T','QVAPOR','RHi100','RHi125','PRI_INU','PRI_IDE','QICE','PRS_IAU','PRS_SCW','QSNOW',\
                'QGRAUP','QCLOUD','QRAIN','RR'])
periods = ['past','warm']

#Dirs
plevdir = '/work/users/emiive/PGW/plevfiles/'
tiledir = '/work/users/emiive/PGW/tilefiles/' 
figdir = ('./figures')

# Get the variables
Vars = {}
for period in periods:

    wrf = xr.open_dataset('%s/Tile_%s_%i_%s_%s-%s_%s.nc' %(tiledir,period,y,m,d1,d2,fix))
    if slice_time:
        wrf = wrf.sel(Time=slice(t_from,t_to))
    if plev:
        hgt = wrf['HGT'][0]
        RAINNC = wrf['RAINNC']
        wrf = xr.open_dataset('%s/%i_%s_%s_%s_plevs.nc' %(plevdir,y,m,period,fix))
        wrf = wrf.rename({'times':'Time'})
        if slice_time:
            wrf = wrf.sel(Time=slice(t_from,t_to))

    Vars[period] = {}

    for v in varlist:
        if v == 'RR':
            da = RAINNC.diff("Time")
            da = da.sum('Time')
        elif v == 'T':
            da = wrf['TEMPERATURE']-273.15
            da = da.mean('Time')
            da = da.rename('T')
        elif 'RH' in v:
            eps = 0.622
            Q = wrf['QVAPOR']
            T = wrf['TEMPERATURE']-273.15
            try:
                P = wrf['PRESSURE']
            except:
                P = wrf.plevs
            e = Q/(Q+eps) * P
            if 'RHw' in v:
                es = 6.112*np.exp(17.62*T/(243.12+T))*100
            elif 'RHi' in v:
                es = 6.112*np.exp(22.46*T/(272.62+T))*100
            da = e/es *100
            if v == 'RHi125':
                da = da.where(da > 125).count('Time') / len(da.Time)
            else:
                da = da.where(da > 100).count('Time') / len(da.Time) 
            da = da.rename(v)
            #da = da.mean('Time')
        else:
            da = wrf[v].mean('Time')
        if 'Q' in v:
            da = da*1000 # g/kg
        dalat = da.mean('south_north',skipna=True)
        dalat = dalat.where(dalat != 0)
        try:
            dalat['plevs'] = dalat['plevs']/100 # hPa
        except:
            pass
        Vars[period][v] = dalat


# Calc. diff warm - past
diff = {}
for v in varlist:
    diff[v] = Vars['warm'][v] - Vars['past'][v]
    diff[v] = diff[v].where(diff[v] != 0)

# Set plotting limits
cmap = {}; vmin = {}; vmax = {}; levels = {}; extend = {}; lab={}
for v in varlist:
    cmap[v] = 'gist_earth_r'
    vmin[v] = 0
cmap['T'] = plt.get_cmap('Blues_r'); cmap['T'].set_over('mistyrose'); vmin['T'] = -48; vmax['T'] = 0
cmap['W'] = 'bwr_r'; vmin['W'] = -1.2; vmax['W'] = 1.2
#vmin['RHw'] = 85; vmax['RHw'] = 100
#vmin['RHi'] = 85; vmax['RHi'] = 100
vmax['RHw'] = 0.8
vmax['RHi100'] = 0.7
vmax['RHi125'] = 0.3
vmax['QNICE'] = 1e5
vmax['QICE'] = 4e-3
vmax['PRI_INU'] = 2.5e-7
#vmax['PRI_INU'] = 1e-4
cmap['PRI_IDE'] = 'RdBu'; vmin['PRI_IDE'] = -0.5e-3; vmax['PRI_IDE'] = 0.5e-3
#cmap['PRI_IDE'] = 'RdBu'; vmin['PRI_IDE'] = -0.1; vmax['PRI_IDE'] = 0.1
vmax['PRS_IAU'] = 0.8e-3
#vmax['PRS_IAU'] = 0.2
vmax['PRI_IHM'] = 5e-5
vmax['PRS_SCW'] = 0.02
vmax['QSNOW'] = 0.24

vmin['diff'] = {}; vmax['diff'] = {};
vmin['diff']['T'] = -2; vmax['diff']['T'] = 2
vmin['diff']['QVAPOR']=-0.35; vmax['diff']['QVAPOR']=0.35
vmin['diff']['RHi100'] = -0.2; vmax['diff']['RHi100'] = 0.2
vmin['diff']['RHi125'] = -0.2; vmax['diff']['RHi125'] = 0.2
vmin['diff']['RHw'] = -0.2; vmax['diff']['RHw'] = 0.2
#vmin['diff']['RHi'] = -8; vmax['diff']['RHi'] = 8
#vmin['diff']['RHw'] = -8; vmax['diff']['RHw'] = 8
vmin['diff']['QICE']=-1e-3; vmax['diff']['QICE']=1e-3
vmin['diff']['PRI_INU']=-0.5e-7; vmax['diff']['PRI_INU']=0.5e-7
vmin['diff']['PRI_IDE']=-0.3e-3; vmax['diff']['PRI_IDE']=0.3e-3
vmin['diff']['PRS_IAU']=-0.2e-3; vmax['diff']['PRS_IAU']=0.2e-3
vmin['diff']['PRI_IHM']=-3e-5; vmax['diff']['PRI_IHM']=3e-5
vmin['diff']['PRS_SCW']=-0.01; vmax['diff']['PRS_SCW']=0.01
vmin['diff']['QSNOW']=-0.07; vmax['diff']['QSNOW']=0.07
vmin['diff']['QCLOUD']=-0.1; vmax['diff']['QCLOUD']=0.1
vmin['diff']['QRAIN']=-0.02; vmax['diff']['QRAIN']=0.02

##################################################
# FIGURE
##################################################
cols = 2
axlen = len(varlist)+1
hrat = [1]*len(varlist)
hrat.append(0.5)
fig, ax = plt.subplots(ncols=cols,nrows=axlen,sharex=True,\
                       gridspec_kw={'height_ratios': hrat},figsize=(9,12))

for v,i in zip(varlist,range(len(varlist))):
  for j in range(cols):

      if j==0:
          latmean = Vars['past'][v]
      elif (j==1) and (cols==3):
          latmean = Vars['warm'][v]
      else:
          latmean = diff[v]
          cmap[v] = 'bwr_r'; cmap['T'] = 'bwr'

      if v == 'RR':
          latmean.plot(ax=ax[i,j])
          pos = ax[0,j].get_position()
          pos2 = ax[i,j].get_position()
          ax[i,j].set_position([pos.x0,pos2.y0,pos.width,pos2.height])
          ax[i,0].set_ylabel('mm')
          #if j==cols-1:
              #ax[i,j].set_ylim([-110,5])
              #ax[i,j].set_ylim([-10,60])
      else:
          if j==(cols-1):
              try:
                  latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],extend='both',\
                                        vmin=vmin['diff'][v],vmax=vmax['diff'][v])
              except:
                  latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],extend='both')
          else:
              try:
                  latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],extend='both',\
                                        vmin=vmin[v],vmax=vmax[v])
              except:
                  latmean.plot.contourf(ax=ax[i,j], cmap=cmap[v],extend='both')

          if plev:
              ax[i,0].set_ylabel('hPa')
              ax[i,j].invert_yaxis()
          else:
              ax[i,0].set_ylabel('ML')
          ax[i,1].set_yticklabels('')

      ax[i,j].grid()
      ax[i,j].set_xlabel('')
      ax[i,1].set_ylabel('')
      #ax[i,j].set_title(v)

hgt = hgt.mean('south_north')
for j in range(2):
    ax[-1,j].plot(hgt,label='HGT')
    pos = ax[0,j].get_position()
    pos2 = ax[-1,j].get_position()
    ax[-1,j].set_position([pos.x0,pos2.y0,pos.width,pos2.height])
    ax[-1,j].legend(fontsize='xx-small')
    #ax[-1,j].set_xticks(np.arange(0,175,25))
    #ax[-1,j].set_xticklabels(xlon[np.arange(0,175,25)].values.round(1))
    ax[-1,j].grid()
    ax[-1,j].set_xlabel('lon')
    #ax[-1,j].set_title('HGT')
    ax[-1,0].set_ylabel('m')
    ax[-1,1].set_yticklabels('')

if savefig:
    #plt.savefig('%s/Vertcross_%i_%s_%s-%s_%s_PR.png' %(figdir,y,m,d1,d2,fix),\
    plt.savefig('%s/Vertcross_%i_%s_18_%s_PR.png' %(figdir,y,m,fix),\
                bbox_inches='tight', dpi=120)
else:
    plt.show(block=False)
