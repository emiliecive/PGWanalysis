# FIGURE 2 IN PAPER
# Subplot figure where vars on x axis, years on y axis, inkl windrose and NAO index

import sys
import numpy as np
import pandas as pd
import xarray as xr
import os, fnmatch
from resource import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from cartopy import crs
import geopandas as gpd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import windrose
#matplotlib.use('Agg')  # To not open a plot window
#import cartopy.feature as cfeature
#from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
#                 cartopy_ylim, latlon_coords, ll_to_xy, CoordPair)
import calendar

#Debug
from IPython import embed

varlist = list(['T','RR','RR_diff'])
months = [12,1,2]
fix = 'DJF'
years = np.array([2009,2012,2014,2015,2019])
savefig = 0
domain = 'd02'

# For plotting tile
plot_tile = 1
X = [20,185]; Y = [200,285] #Stolsheimen

# For windrose:
lon=120; lat = 243

# For NAO:
File = './monthly_NAO_index.txt'
index =  pd.read_table(File,names=['year','month','val'],parse_dates=[[0,1]],date_parser=lambda arg: pd.to_datetime(arg, format='%Y %m'))
winter = index[(index.year_month.dt.month<months[-1]+1) | (index.year_month.dt.month>months[0]-1)]
pos = winter[winter.val>0]
neg = winter[winter.val<0]

#Dirs
dirr = '/mnt/elephant/WRFsimulations/ICEBOX_PGW'
wrfdir = ('%s/archive2/postpost' %dirr)
outdir = ('./files')
figdir = ('./figures')

data_proj = crs.Stereographic()
map_proj = crs.PlateCarree()
boundary = gpd.read_file('./boundaryshp/ne_50m_admin_0_boundary_lines_land.shp')
HGT = xr.open_dataset('%s/wrfout_d02_2019-10.nc4' %wrfdir)
HGT = HGT['HGT']

# Figure settings:
stat = {}; lab = {}; cmap = {}; vmin = {}; vmax = {}
stat['T'] = 'mean'; lab['T'] = '$\degree$C'; cmap['T'] = 'RdBu_r'; vmin['T'] = -10; vmax['T'] = 10
stat['RR'] = 'sums'; lab['RR'] = 'mm'; cmap['RR'] = 'gist_earth_r'; vmin['RR'] = 0; vmax['RR'] = 800
lab['RR_diff'] = 'mm'; cmap['RR_diff'] = 'RdBu'; vmin['RR_diff'] = -150; vmax['RR_diff'] = 150


fig = plt.figure(figsize=(10,13))
gs = fig.add_gridspec(len(years),len(varlist)+2,width_ratios=[1/4.5,1/4.5,1/4.5,1/4.5,0.5/4.5],wspace=0.1)
#fig, ax = plt.subplots(ncols=len(varlist)+1,nrows=len(years),figsize=(10,13),subplot_kw={'projection': data_proj})

for v,j in zip(varlist,range(len(varlist))):

    if 'diff' in v:
        ds = xr.open_mfdataset('%s/*_%s_warm_past.nc' %(outdir,v), concat_dim='year')
    else:
        ds = xr.open_mfdataset('%s/past/*_past_%s_%s.nc' %(outdir,v,stat[v]), concat_dim='year')

    ds = ds.assign_coords(year=years)
    if stat == 'sums':
        ds = ds.sel(month=months).sum('month')
    else:
        ds = ds.sel(month=months).mean('month')
    var = list(ds.keys())[0]

    ######## Plotting #########

    for i in range(len(years)):
        ax = fig.add_subplot(gs[i,j], projection = data_proj)

        if len(ds[var].shape) == 3:
           Map = ds[var].sel(year=years[i])\
                        .plot.pcolormesh('XLONG','XLAT',ax=ax,
                                         transform=map_proj,
                                         cmap=cmap[v],add_colorbar=False,
                                         vmin=vmin[v], vmax=vmax[v])
        elif len(ds[var].shape) == 4:
           Map = ds[var].sel(bottom_top=0,year=years[i])\
                        .plot.pcolormesh('XLONG','XLAT',ax=ax,
                                         transform=map_proj,
                                         cmap=cmap[v],add_colorbar=False,
                                         vmin=vmin[v],vmax=vmax[v])

        ax.set_title('')
        #ax.text(-0.2,0.5,years[i]+1,verticalalignment='center',size=11,rotation=90,transform=ax[i,0].transAxes)
        ax.coastlines('50m', linewidth=0.6)
        ax.add_geometries(boundary.geometry, crs = crs.PlateCarree(),\
                               facecolor='none', edgecolor='k', linewidth=0.8)
        ax.set_xlim(1.447*10**5,1.817*10**6)
        ax.set_ylim(7.02*10**6,9.4*10**6)

        if j == 0:
            ax.text(-0.2,0.5,years[i]+1,verticalalignment='center',size=10,rotation=90,transform=ax.transAxes)

        if i == 0:
            ax.set_title(v)
        elif i == len(years)-1:
            # Add colorbar
            axins = inset_axes(ax,width="100%",height="7%",loc='lower center',\
                         bbox_to_anchor=(0, -0.2, 1, 0.7),bbox_transform=ax.transAxes)
            cbar = fig.colorbar(Map,cax=axins,shrink=0.7,orientation='horizontal',extend='both')
            cbar.set_label(lab[v])

            # Add tile
            if (v == 'RR_diff') and plot_tile:
                xlen = np.arange(X[0],X[-1]+1)
                ylen = np.arange(Y[0],Y[-1]+1)
                ax.plot(xlen, np.full((len(xlen)),Y[0]), color='black', transform=map_proj)
                ax.plot(xlen, np.full((len(xlen)),Y[1]), color='black', transform=map_proj)
                ax.plot(np.full((len(ylen)),X[0]), ylen, color='black', transform=map_proj)
                ax.plot(np.full((len(ylen)),X[1]), ylen, color='black', transform=map_proj)

# Add windroses:
for y,i in zip(years,range(len(years))):
    ax1 = fig.add_subplot(gs[i,-2], projection='windrose')
    wrf = xr.open_mfdataset(('%s/wrfout_%s_%s-12.nc4' %(wrfdir,domain,y)) and \
                            ('%s/wrfout_%s_%s-0[1-2].nc4' %(wrfdir,domain,y+1)),\
                             drop_variables={'T','XTIME'})
    dd = wrf['DD']
    ff = wrf['FF']
    dd = dd[:,9,lat,lon].values
    ff = ff[:,9,lat,lon].values

    ax1.bar(dd, ff, normed=True, opening=0.6, bins=np.arange(0,30,5),nsector=8,
            edgecolor='black', cmap=cm.gray_r)
    ax1.set_yticks(np.arange(10, 40, step=10))
    ax1.set_yticklabels(np.arange(10, 40, step=10))
    if i == 0:
        ax1.set_title('WS / WD')
    if i == len(years)-1:
        ax1.set_legend(bbox_to_anchor=(0,-0.5),fontsize=3,markerscale=0.8)

# Add NAO index:
for y,i in zip(years,range(len(years))):
    ax2 = fig.add_subplot(gs[i,-1])

    for m in months:
        if m > 4:
            posplot = pos[(pos.year_month.dt.year == y) & (pos.year_month.dt.month == m)]
            negplot = neg[(neg.year_month.dt.year == y) & (neg.year_month.dt.month == m)]
        else:
            posplot = pos[(pos.year_month.dt.year == y+1) & (pos.year_month.dt.month == m)]
            negplot = neg[(neg.year_month.dt.year == y+1) & (neg.year_month.dt.month == m)]
        ax2.bar(negplot.year_month,negplot.val,width=15,align='edge',color='blue')
        ax2.bar(posplot.year_month,posplot.val,width=15,align='edge',color='red')
    ax2.set_ylim([-2, 2.4])
    ax2.set_xticklabels('')
    ax2.grid()
    if i == 0:
        ax2.set_title('NAO index')
    elif i == len(years)-1:
        ax2.set_xlabel('Dec, Jan, Feb')

    if savefig:
        plt.savefig('%s/Xvars_Yyears_windrose_NAO.png' %(figdir), bbox_inches='tight', dpi=120)
    else:
        plt.show(block=False)
