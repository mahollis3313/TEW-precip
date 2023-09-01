#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 23:15:21 2021

@author: margaret
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
import netCDF4
from netCDF4 import Dataset, num2date

#NH850 = sorted(glob.glob('/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_*_850_NH.nc'))
#SH850 = sorted(glob.glob('/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_*_850_SH.nc'))
NH850 = sorted(glob.glob('/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_*_850_NH.nc'))
SH850 = sorted(glob.glob('/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_*_850_SH.nc'))
               
#NH700 = sorted(glob.glob('/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_*_700_NH.nc'))
#SH700 = sorted(glob.glob('/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_*_700_SH.nc'))
NH700 = sorted(glob.glob('/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_*_700_NH.nc'))
SH700 = sorted(glob.glob('/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_*_700_SH.nc'))


def hemitext(hemi):
    if hemi == 1:
        hemitxt = 'northern'
    elif hemi == 2:
        hemitxt = 'southern'
    else:
        print('Hemisphere is not valid.')
        sys.exit()
    return hemitxt

def hemishort(hemi):
    if hemi == 1:
        hemishrt = 'NH'
    elif hemi == 2:
        hemishrt = 'SH'
    else:
        print('Hemisphere is not valid.')
        sys.exit()
    return hemishrt

# variance shamelessly stolen from stackoverflow, and adapted to fit into my composite function
# https://stackoverflow.com/questions/15638612/calculating-mean-and-standard-deviation-of-the-data-which-does-not-fit-in-memory
def online_variance(data):
    n = 0
    mean = 0
    M2 = 0

    for x in data:
        n = n + 1
        delta = x - mean
        mean = mean + delta/n
        M2 = M2 + delta*(x - mean)

    variance = M2/(n - 1)
    return variance



def compositePrecip(filelist, latlons=False):
    #important variables that will be same file to file
    important  = filelist[0]
    openimport = Dataset(important, mode='r')
    
    lats = openimport['lat'][:]
    lons = openimport['lon'][:]
    
    openimport.close()
    
    #now make variables for things we need for averaging
    numrecs = 0
    allprecip = np.zeros((lats.size, lons.size))
    
    #and things for the variance
    n = 0
    varmean = np.zeros((lats.size, lons.size))
    M2 = 0
    
    for file in filelist:
        #netcdf things
        waves = Dataset(file, mode='r')
        precip = waves['precip'][:]
        
        #things for the mean
        allprecip += np.nansum(precip, axis=0)
        numrecs += waves.dimensions['record'].size
        
        #things for the variance
        for time in range(precip.shape[0]):
            x = precip[time,:,:]
            n = n + 1
            delta = x - varmean
            varmean = varmean + delta/n
            M2 = M2 + delta*(x- varmean)
        
        waves.close()
        
    mean = allprecip / numrecs
    variance = M2/(n - 1)
    
    if latlons == True:
        return mean, numrecs, variance, lats, lons
    else:
        return mean, numrecs, variance

mean_NH850, numrecs_NH850, var_NH850, lats, lons = compositePrecip(NH850, latlons=True)
mean_SH850, numrecs_SH850, var_SH850 = compositePrecip(SH850)
mean_NH700, numrecs_NH700, var_NH700 = compositePrecip(NH700)
mean_SH700, numrecs_SH700, var_SH700 = compositePrecip(SH700)

print('Update: Completed composites of all precipitation')

def compositeMatched(filelist, level, hemi, latlons=False):
    #important variables that will be same file to file
    important  = filelist[0]
    openimport = Dataset(important, mode='r')
    
    lats = openimport['lat'][:]
    lons = openimport['lon'][:]
    
    openimport.close()
    
    #now make variables for things we need for averaging
    numrecs = 0
    allprecip = np.zeros((lats.size, lons.size))
    
    #and things for the variance
    n = 0
    varmean = np.zeros((lats.size, lons.size))
    M2 = 0
    
    if level == 850:
        lev=0
    elif level == 700:
        lev=1
    else:
        print('invalid level. exiting')
        sys.exit()
    
    for file in filelist:
        waves = Dataset(file, mode='r')
        year = num2date(waves['time'][0], units=waves['time'].units).year
        
        matchindices = np.loadtxt(f'/home/mhollis/track_stats/CollocationUpdated/match_indices_{year}_{hemishort(hemi)}.csv', delimiter=',').astype('int')
        
        precip = waves['precip'][matchindices[:,lev], :, :]
        
        #things for the mean
        allprecip += np.nansum(precip, axis=0)
        numrecs += len(matchindices)
        
        #things for the variance
        for time in range(precip.shape[0]):
            x = precip[time,:,:]
            n = n + 1
            delta = x - varmean
            varmean = varmean + delta/n
            M2 = M2 + delta*(x- varmean)
        
        waves.close()
        
    mean = allprecip / numrecs
    variance = M2/(n - 1)
    
    if latlons == True:
        return mean, numrecs, variance, lats, lons
    else:
        return mean, numrecs, variance


match_NH850, matchrecs_NH850, matchvar_NH850 = compositeMatched(NH850, 850, 1)
match_SH850, matchrecs_SH850, matchvar_SH850 = compositeMatched(SH850, 850, 2)
match_NH700, matchrecs_NH700, matchvar_NH700 = compositeMatched(NH700, 700, 1)
match_SH700, matchrecs_SH700, matchvar_SH700 = compositeMatched(SH700, 700, 2)

print('Update: Completed composites of precipitation associated with matched TEWs')

def compositeUnmatched(filelist, level, hemi, latlons=False):
    #important variables that will be same file to file
    important  = filelist[0]
    openimport = Dataset(important, mode='r')
    
    lats = openimport['lat'][:]
    lons = openimport['lon'][:]
    
    openimport.close()
    
    #now make variables for things we need for averaging
    numrecs = 0
    allprecip = np.zeros((lats.size, lons.size))
    
    #and things for the variance
    n = 0
    varmean = np.zeros((lats.size, lons.size))
    M2 = 0
    
    
    for file in filelist:
        waves = Dataset(file, mode='r')
        year = num2date(waves['time'][0], units=waves['time'].units).year
        
        indices = np.loadtxt(f'/home/mhollis/track_stats/CollocationUpdated/nomatch_indices_{level}hPa_{year}_{hemishort(hemi)}.csv', delimiter=',').astype('int')
        
        precip = waves['precip'][indices, :, :]
        
        #things for the mean
        allprecip += np.nansum(precip, axis=0)
        numrecs += len(indices)
        
        #things for the variance
        for time in range(precip.shape[0]):
            x = precip[time,:,:]
            n = n + 1
            delta = x - varmean
            varmean = varmean + delta/n
            M2 = M2 + delta*(x- varmean)
        
        waves.close()
        
    mean = allprecip / numrecs
    variance = M2/(n - 1)
    
    if latlons == True:
        return mean, numrecs, variance, lats, lons
    else:
        return mean, numrecs, variance

nomatch_NH850, nomatchrecs_NH850, nomatchvar_NH850 = compositeUnmatched(NH850, 850, 1)
nomatch_SH850, nomatchrecs_SH850, nomatchvar_SH850 = compositeUnmatched(SH850, 850, 2)
nomatch_NH700, nomatchrecs_NH700, nomatchvar_NH700 = compositeUnmatched(NH700, 700, 1)
nomatch_SH700, nomatchrecs_SH700, nomatchvar_SH700 = compositeUnmatched(SH700, 700, 2)

print('Update: Completed composites of precipitation associated with unmatched TEWs')

# # # # # # # #
# PLOTTING TIME
# # # # # # # #

def letterFromCompType(level, hemi):
    if level == 850 and hemishort(hemi) == 'NH':
        return 'a)'
    elif level == 700 and hemishort(hemi) == 'NH':
        return 'b)'
    elif level == 850 and hemishort(hemi)== 'SH':
        return 'c)'
    elif level == 700 and hemishort(hemi) == 'SH':
        return 'd)'

def plotComposite(composite, level, hemi, numrecs, comptype):    
    fig = plt.figure(figsize=(7,5),dpi=300)
    ax = plt.axes()
    levels = (0.000000000000001, .05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9)
    #norm = mcolors.BoundaryNorm(levels, 256)
    #print(levels)
    plt.contourf(lons, lats, composite, levels=levels, cmap='RdYlGn_r', extend='max')
    #set extent would be nice to do so that it's -5 to 5 on both axes
    plt.colorbar()
    ax.set_xlim(-6,6)
    ax.set_ylim(-5,5)
    ax.set_xlabel('Distance from TEW Center (degrees)')
    ax.set_ylabel('Distance from TEW Center (degrees)')
    #ax.set_title(f'Composite TEW precipitation rate, mm/hr \n for {comptype} TEW points in the {hemitext(hemi)} hemisphere at {level} hPa')
    ax.text(-7, 4.5, f'{letterFromCompType(level, hemi)}')
    ax.text(-4.9, -4.9, s=f'n={numrecs}')
    
    ax.grid()
    
    #plt.savefig(f'/home/mhollis/Precip_stuff/Figures/PublicationUpdate/IMERG/composites/composite_{comptype}_{hemishort(hemi)}_{level}hPa_maxcolor.png', bbox_inches = "tight", dpi = 300)
    plt.savefig(f'/home/mhollis/Precip_stuff/Figures/PublicationUpdate/MERRA2/composites/composite_{comptype}_{hemishort(hemi)}_{level}hPa_maxcolor.png', bbox_inches = "tight", dpi = 300)

plotComposite(mean_NH850, 850, 1, numrecs_NH850, 'all')
plotComposite(mean_SH850, 850, 2, numrecs_SH850, 'all')
plotComposite(mean_NH700, 700, 1, numrecs_NH700, 'all')
plotComposite(mean_SH700, 700, 2, numrecs_SH700, 'all')

plotComposite(match_NH850, 850, 1, matchrecs_NH850, 'collocated')
plotComposite(match_SH850, 850, 2, matchrecs_SH850, 'collocated')
plotComposite(match_NH700, 700, 1, matchrecs_NH700, 'collocated')
plotComposite(match_SH700, 700, 2, matchrecs_SH700, 'collocated')

plotComposite(nomatch_NH850, 850, 1, nomatchrecs_NH850, 'unmatched')
plotComposite(nomatch_SH850, 850, 2, nomatchrecs_SH850, 'unmatched')
plotComposite(nomatch_NH700, 700, 1, nomatchrecs_NH700, 'unmatched')
plotComposite(nomatch_SH700, 700, 2, nomatchrecs_SH700, 'unmatched')

print('Update: Completed plots of composite means')

#variance gets its own plotting function because colorbars

def plotVariance(composite, level, hemi, numrecs, comptype):    
    fig = plt.figure(figsize=(7,5),dpi=300)
    ax = plt.axes()
    levels = (0.000000000000001, .25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5)
    #norm = mcolors.BoundaryNorm(levels, 256)
    #print(levels)
    plt.contourf(lons, lats, composite, levels=levels, cmap='RdYlGn_r', extend='max')
    #set extent would be nice to do so that it's -5 to 5 on both axes
    plt.colorbar()
    ax.set_xlim(-6,6)
    ax.set_ylim(-5,5)
    ax.set_xlabel('Distance from TEW Center (degrees)')
    ax.set_ylabel('Distance from TEW Center (degrees)')
    #ax.set_title(f'Composite variance of TEW precipitation rate, \n for {comptype} TEW points in the {hemitext(hemi)} hemisphere at {level} hPa')
    ax.text(-7, 4.5, f'{letterFromCompType(level, hemi)}')
    ax.text(-4.9, -4.9, s=f'n={numrecs}')
    
    ax.grid()
    
    plt.savefig(f'/home/mhollis/Precip_stuff/Figures/PublicationUpdate/IMERG/composites/compvar_{comptype}_{hemishort(hemi)}_{level}hPa_maxcolor.png', bbox_inches = "tight", dpi = 300)
    plt.savefig(f'/home/mhollis/Precip_stuff/Figures/PublicationUpdate/MERRA2/composites/compvar_{comptype}_{hemishort(hemi)}_{level}hPa_maxcolor.png', bbox_inches = "tight", dpi = 300)

plotVariance(var_NH850, 850, 1, numrecs_NH850, 'all')
plotVariance(var_SH850, 850, 2, numrecs_SH850, 'all')
plotVariance(var_NH700, 700, 1, numrecs_NH700, 'all')
plotVariance(var_SH700, 700, 2, numrecs_SH700, 'all')

plotVariance(match_NH850, 850, 1, matchrecs_NH850, 'collocated')
plotVariance(match_SH850, 850, 2, matchrecs_SH850, 'collocated')
plotVariance(match_NH700, 700, 1, matchrecs_NH700, 'collocated')
plotVariance(match_SH700, 700, 2, matchrecs_SH700, 'collocated')

plotVariance(nomatch_NH850, 850, 1, nomatchrecs_NH850, 'unmatched')
plotVariance(nomatch_SH850, 850, 2, nomatchrecs_SH850, 'unmatched')
plotVariance(nomatch_NH700, 700, 1, nomatchrecs_NH700, 'unmatched')
plotVariance(nomatch_SH700, 700, 2, nomatchrecs_SH700, 'unmatched')

print('Update: Completed plots of composite variance')
print('All done.')