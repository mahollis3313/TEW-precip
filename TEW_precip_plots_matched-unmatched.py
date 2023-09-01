#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 16:47:35 2021

@author: margaret

Let's make a plot
It ended up making way more sense to just get the precip, and so now we get to add up the annual files and plot things!
"""
# import packages
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import netCDF4
from netCDF4 import Dataset, num2date

level = int(sys.argv[1])

#paths to the relevant files

#the standard files for all of the TEWs
#NHTEWfile = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_precip_accumulated/precip_noTC_{level}_NH_accumulated.nc'
#SHTEWfile = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_precip_accumulated/precip_noTC_{level}_SH_accumulated.nc'
#annualfile = '/data/deluge/scratch/IMERG/TEW_Precip_updated/Annual_precip/annual_precip.nc'
NHTEWfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_annual/precip_noTC_{level}_NH_accumulated.nc'
SHTEWfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_annual/precip_noTC_{level}_SH_accumulated.nc'
annualfile = '/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/global_annual/annual_precip.nc'


#the vertically collocated precip now
#NHmatchfile = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_matched_precip_accumulated/precip_noTC_matched_{level}_NH_accumulated.nc'
#SHmatchfile = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_matched_precip_accumulated/precip_noTC_matched_{level}_SH_accumulated.nc'
NHmatchfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_collocated/precip_noTC_matched_{level}_NH_accumulated.nc'
SHmatchfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_collocated/precip_noTC_matched_{level}_SH_accumulated.nc'

#add up the precip


#for the annual files
NHTEW  = Dataset(NHTEWfile, mode='r')
SHTEW  = Dataset(SHTEWfile, mode='r')
annual = Dataset(annualfile, mode='r')
    
preciplats = annual['lat'][:]
preciplons = annual['lon'][:]

allprecip = np.nansum(annual['precip'], axis=0)

TEWprecip = (np.nansum(NHTEW['precip'], axis=0) + np.nansum(SHTEW['precip'], axis=0))[0,:,:]
avgTEWprecip = np.nanmean(NHTEW['precip'][:] + SHTEW['precip'], axis=0)[0,:,:]
ratio = TEWprecip/allprecip

#and for the vertically collocated TEWs
NHmatch = Dataset(NHmatchfile, mode='r')
SHmatch = Dataset(SHmatchfile, mode='r')

#precip lats and lons are same

matchprecip = (np.nansum(NHmatch['precip'], axis=0) + np.nansum(SHmatch['precip'], axis=0))[0,:,:]

avgMatchPrecip = np.nanmean((NHmatch['precip'][:] + SHmatch['precip']), axis=0)[0,:,:]
matchRatio     = matchprecip/allprecip
matchToTEW     = matchprecip/TEWprecip

#calculate unmatched precip
unmatchedPrecip = (np.nansum(((NHTEW['precip'][:] - NHmatch['precip'][:]) + (SHTEW['precip'][:] - SHmatch['precip'][:])), axis=0))[0,:,:]

avgUnmatchPrecip = (np.nanmean(((NHTEW['precip'][:] - NHmatch['precip'][:]) + (SHTEW['precip'][:] - SHmatch['precip'][:])), axis=0))[0,:,:]

unmatchRatio = unmatchedPrecip/allprecip
unmatchtoTEW = unmatchedPrecip/TEWprecip



#plotting time

#some setup stuff
plotsfolder = '/home/mhollis/Precip_stuff/Figures/PublicationUpdate/MERRA2'

#a/b panel labels
if level == 700:
    panel = 'a'
elif level == 850:
    panel = 'b'

#accumulated precip from all waves
fig = plt.figure(dpi=300, figsize=(8,4))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
levels=(0.00001, 1, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, TEWprecip, levels=levels, cmap='RdYlGn_r', norm=norm) #, colors=('lightgreen','green','gold','darkorange','red','darkred'))
plt.colorbar(shrink=0.62, pad=0.08)

ax.set_ylim(-70, 70)
grid = ax.gridlines(draw_labels=True)
grid.top_labels = False
grid.left_labels = False
grid.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.title(f'Accumulated Precipitation, mm, \n from TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/TEWprecip{level}hPa.png', bbox_inches = "tight", dpi = 300)


#average annual precip from all waves
fig1 = plt.figure(dpi=300, figsize=(8,4))
ax1 = plt.axes(projection=ccrs.PlateCarree())
ax1.coastlines()
levels=(0.00001, 0.1, 1, 5, 10, 50, 100, 250, 500, 1000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, avgTEWprecip, levels=levels, norm=norm, cmap='RdYlGn_r') #, colors=('lightgreen','green','gold','darkorange','red','darkred'))
plt.colorbar(shrink=0.62, pad=0.08)

ax1.set_ylim(-70, 70)
grid1 = ax1.gridlines(draw_labels=True)
grid1.top_labels = False
grid1.left_labels = False
grid1.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid1.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.title(f'Average Annual Precipitation, mm, \n from TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/avgTEWprecip{level}hPa.png', bbox_inches = "tight", dpi = 300)


#fraction of total annual precip from all waves
fig2 = plt.figure(dpi=300, figsize=(8,4))
ax2 = plt.axes(projection=ccrs.PlateCarree())
ax2.coastlines()
#levels = (0.0000000000001, 0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50, 1)
#levels = (0.00000001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.30, 0.50, 1)
#levels = np.linspace(0, 0.3, 16) #I have tried multiple linspaces it's all meh
#levels = np.nanpercentile(ratio, np.linspace(50.45,100,101))
#levels = np.geomspace(0.005,1,10)
levels = (0.0000000000001, 0.02, 0.04, 0.06, 0.08, 0.1 , 0.12, 0.14, 0.16, 0.18, 0.2 ,
       0.22, 0.24, 0.26, 0.28, 0.3)
norm = mcolors.BoundaryNorm(levels, 256)
map2 = ax2.contourf(preciplons, preciplats, ratio, levels=levels, norm=norm, cmap='YlGnBu', extend='max')
fig2.colorbar(map2, shrink=0.62, pad=0.08)

ax2.set_ylim(-70, 70)
grid2 = ax2.gridlines(draw_labels=True)
grid2.top_labels = False
grid2.left_labels = False
grid2.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid2.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#ax2.set_title(f'Fraction of Total Annual Precipitation \n from TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/Fractionprecip{level}hPa.png', bbox_inches = "tight", dpi = 300)



#and the vertically collocated TEWs

#accumulated precip from collocated waves
fig3 = plt.figure(dpi=300, figsize=(8,4))
ax3 = plt.axes(projection=ccrs.PlateCarree())
ax3.coastlines()
levels=(0.00001, 1, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, matchprecip, levels=levels, cmap='RdYlGn_r', norm=norm)
plt.colorbar(shrink=0.62, pad=0.08)

ax3.set_ylim(-70, 70)
grid3 = ax3.gridlines(draw_labels=True)
grid3.top_labels = False
grid3.left_labels = False
grid3.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid3.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.set_title(f'Accumulated Precipitation, mm, \n from Vertically Collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/TEWprecipMatch{level}hPa.png', bbox_inches = "tight", dpi = 300)


#average annual precip from collocated TEWs
fig4 = plt.figure(dpi=300, figsize=(8,4))
ax4 = plt.axes(projection=ccrs.PlateCarree())
ax4.coastlines()
levels=(0.00001, 0.1, 1, 5, 10, 50, 100, 250, 500, 1000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, avgMatchPrecip, levels=levels, norm=norm, cmap='RdYlGn_r')
plt.colorbar(shrink=0.62, pad=0.08)

ax4.set_ylim(-70, 70)
grid4 = ax4.gridlines(draw_labels=True)
grid4.top_labels = False
grid4.left_labels = False
grid4.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid4.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.set_title(f'Average Annual Precipitation, mm, \n from Vertically Collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/avgTEWprecipMatch{level}hPa.png', bbox_inches = "tight", dpi = 300)


#fraction of total annual precip from collocated TEWs
fig5 = plt.figure(dpi=300, figsize=(8,4))
ax5 = plt.axes(projection=ccrs.PlateCarree())
ax5.coastlines()
#levels = (0.00000001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.30, 0.50, 1)
#levels = np.linspace(0, 0.5, 21)
levels = (0.0000000000001, 0.02, 0.04, 0.06, 0.08, 0.1 , 0.12, 0.14, 0.16, 0.18, 0.2 ,
       0.22, 0.24, 0.26, 0.28, 0.3)
norm = mcolors.BoundaryNorm(levels, 256)
map5 = ax5.contourf(preciplons, preciplats, matchRatio, levels=levels, norm=norm, cmap='YlGnBu', extend='max')
fig5.colorbar(map5, shrink=0.62, pad=0.08)

ax5.set_ylim(-70, 70)
grid5 = ax5.gridlines(draw_labels=True)
grid5.top_labels = False
grid5.left_labels = False
grid5.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid5.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#ax5.set_title(f'Fraction of Total Annual Precipitation \n from Vertically Collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/FractionprecipMatch{level}hPa.png', bbox_inches = "tight", dpi = 300)



#and the unmatched TEWs

#accumulated precip from unmatched TEWs
fig6 = plt.figure(dpi=300, figsize=(8,4))
ax6 = plt.axes(projection=ccrs.PlateCarree())
ax6.coastlines()
levels=(0.00001, 1, 10, 50, 100, 500, 1000, 2500, 5000, 10000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, unmatchedPrecip, levels=levels, cmap='RdYlGn_r', norm=norm)
plt.colorbar(shrink=0.62, pad=0.08)

ax6.set_ylim(-70, 70)
grid6 = ax6.gridlines(draw_labels=True)
grid6.top_labels = False
grid6.left_labels = False
grid6.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid6.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.set_title(f'Accumulated Precipitation, mm, \n from Non-collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/TEWprecipUnmatch{level}hPa.png', bbox_inches = "tight", dpi = 300)

#annual average precip from unmatched TEWs
fig7 = plt.figure(dpi=300, figsize=(8,4))
ax7 = plt.axes(projection=ccrs.PlateCarree())
ax7.coastlines()
levels=(0.00001, 0.1, 1, 5, 10, 50, 100, 250, 500, 1000)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, avgUnmatchPrecip, levels=levels, norm=norm, cmap='RdYlGn_r')
plt.colorbar(shrink=0.62, pad=0.08)

ax7.set_ylim(-70, 70)
grid7 = ax7.gridlines(draw_labels=True)
grid7.top_labels = False
grid7.left_labels = False
grid7.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid7.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.set_title(f'Average Annual Precipitation, mm, \n from Non-collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/avgTEWprecipUnmatch{level}hPa.png', bbox_inches = "tight", dpi = 300)


#fraction of total annual precip from unmatched TEWs
fig8 = plt.figure(dpi=300, figsize=(8,4))
ax8 = plt.axes(projection=ccrs.PlateCarree())
ax8.coastlines()
#levels = (0.00000001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.30, 0.50, 1)
#levels = np.linspace(0, 0.5, 21)
levels = (0.0000000000001, 0.02, 0.04, 0.06, 0.08, 0.1 , 0.12, 0.14, 0.16, 0.18, 0.2 ,
       0.22, 0.24, 0.26, 0.28, 0.3)
norm = mcolors.BoundaryNorm(levels, 256)
map8 = ax8.contourf(preciplons, preciplats, unmatchRatio, levels=levels, norm=norm, cmap='YlGnBu', extend='max')
fig8.colorbar(map8, shrink=0.62, pad=0.08)

ax8.set_ylim(-70, 70)
grid8 = ax8.gridlines(draw_labels=True)
grid8.top_labels = False
grid8.left_labels = False
grid8.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid8.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#ax8.set_title(f'Fraction of Total Annual Precipitation \n from Non-collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/FractionprecipUnmatch{level}hPa.png', bbox_inches = "tight", dpi = 300)



#and the other plots for inter-TEW ratios
#fraction of TEW precip from vertically collocated TEWs
fig9 = plt.figure(dpi=300, figsize=(8,4))
ax9 = plt.axes(projection=ccrs.PlateCarree())
ax9.coastlines()
#levels = (0.00000001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.30, 0.50, 1)
#levels = np.linspace(0,1,11)
levels = np.linspace(0,1,21)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, matchToTEW, levels=levels, norm=norm, cmap='YlGnBu')
plt.colorbar(shrink=0.62, pad=0.08)

ax9.set_ylim(-70, 70)
grid9 = ax9.gridlines(draw_labels=True)
grid9.top_labels = False
grid9.left_labels = False
grid9.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid9.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.title(f'Fraction of Total TEW Precipitation \n from Vertically Collocated TEWs at {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/MatchToTEW{level}hPa.png', bbox_inches = "tight", dpi = 300)


#fraction of TEW precip from unmatched TEWs
fig10 = plt.figure(dpi=300, figsize=(8,4))
ax10 = plt.axes(projection=ccrs.PlateCarree())
ax10.coastlines()
#levels = (0.00000001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.30, 0.50, 1)
#levels = np.linspace(0,1,11)
levels = np.linspace(0,1,21)
norm = mcolors.BoundaryNorm(levels, 256)
plt.contourf(preciplons, preciplats, unmatchtoTEW, levels=levels, norm=norm, cmap='YlGnBu')
plt.colorbar(shrink=0.62, pad=0.08)

ax10.set_ylim(-70, 70)
grid10 = ax10.gridlines(draw_labels=True)
grid10.top_labels = False
grid10.left_labels = False
grid10.xlocator = mticker.FixedLocator(np.arange(-180, 180, 60))
grid10.ylocator = mticker.FixedLocator([-60, -40, -20, 0, 20, 40, 60])

#plt.title(f'Fraction of Total TEW Precipitation \n from TEWs not Vertically Collocated, {level} hPa')
plt.text(-200, 65, f'{panel})')
plt.savefig(f'{plotsfolder}/maps/UnmatchToTEW{level}hPa.png', bbox_inches = "tight", dpi = 300)
