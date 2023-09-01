#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 00:11:32 2022

@author: margaret

Code for some categorical verification of precipitation

Do I want to have a threshold for separating convective and stratiform?
Or just do other stats for that?
"""

import numpy as np
from netCDF4 import Dataset, num2date
import scipy
import glob
import matplotlib.pyplot as plt
import sys
from math import ceil

#functions I wrote or shamelessly stole and modified
#let's make some circles to go around the track points
def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1 = map(np.radians, [lon1, lat1])
    lon2, lat2 = map(np.radians, [lon2, lat2])
    return 6371 * (np.arccos(np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)))

#and a function for getting indices for TEW point and center of the mask
#shamelessly stolen from https://stackoverflow.com/questions/30180241/numpy-get-the-column-and-row-index-of-the-minimum-value-of-a-2d-array
#but modified to give back indices as ints
def find_min_idx(x):
    k = x.argmin()
    ncol = x.shape[1]
    return int(k/ncol), k%ncol

level = int(sys.argv[1])
hemi = int(sys.argv[2])

# # # # # # # # #
# lists of files
# # # # # # # # #

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

#IMERG files
IMERGfiles = sorted(glob.glob(f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_*_{level}_{hemishort(hemi)}.nc'))


#MERRA2 files
MERRA2files = []

for year in range(2001,2021):
    MERRA2files.append(f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_{year}_{level}_{hemishort(hemi)}.nc')
    


# Now, for each file, we want to do this
def hits_misses_false_alarms(IMERG_file,MERRA_file):
    
    IMERG = Dataset(IMERG_file, mode='r')
    MERRA = Dataset(MERRA_file, mode='r')
    
    IMERGlon = IMERG.variables['lon'][:]
    IMERGlat = IMERG.variables['lat'][:]
    IMERGlons, IMERGlats = np.meshgrid(IMERGlon, IMERGlat)    
    
    MERRAlon = MERRA.variables['lon'][:]
    MERRAlat = MERRA.variables['lat'][:]
    MERRAlons, MERRAlats = np.meshgrid(MERRAlon, MERRAlat)
    
    records = MERRA.variables['time'][:]
    
    box_height = ceil((MERRAlat[1]-MERRAlat[0])/(IMERGlat[1]-IMERGlat[0]))
    box_width = ceil((MERRAlon[1]-MERRAlon[0])/(IMERGlon[1]-IMERGlon[0]))
            
    #figure out hits, misses, false alarms
    
    hitgrid   = np.zeros(MERRAlons.shape)
    missgrid  = np.zeros(MERRAlons.shape)
    falsegrid = np.zeros(MERRAlons.shape)
    
                
    #find the MERRA point closest to each IMERG point
    for i in range(len(MERRAlon)):
        if i % 10 == 0:
            print(f'working on lon {i}')
        for j in range(len(MERRAlat)):
            lon = MERRAlon[i]
            lat = MERRAlat[j]
        
            #create the box around the point of interest
            dist = great_circle(lon, lat, IMERGlons, IMERGlats)
        
            point = np.asarray(find_min_idx(dist))
            
            
            #then at this point for each time
            for t in range(len(records)):
                
                
                IMERGprecip = IMERG.variables['precip'][t,
                                                        int(point[0]-box_height/2):ceil(point[0]+box_height/2)+1,
                                                        int(point[1]-box_width/2):ceil(point[1]+box_width/2+1)]

                IMERGprecip = np.nansum(IMERGprecip)
                
                MERRAprecip = MERRA.variables['precip'][t,j,i]
                
                
                
                #assign it hit, miss, false alarm
                if IMERGprecip == 0. and MERRAprecip == 0.:
                    continue
                elif IMERGprecip == 0.and MERRAprecip != 0.:
                    falsegrid[j,i] += 1
                elif IMERGprecip != 0. and MERRAprecip == 0.:
                    missgrid[j,i] += 1
                elif IMERGprecip != 0. and MERRAprecip != 0.:
                    hitgrid[j,i] += 1
    
            
        
    print(f'Done with {IMERG_file}')
    return hitgrid, missgrid, falsegrid

MERRA = Dataset(MERRA2files[0], mode='r')
MERRAlon = MERRA['lon'][:]
MERRAlat = MERRA['lat'][:]
MERRAlons, MERRAlats = np.meshgrid(MERRAlon, MERRAlat)

hitgrids   = np.zeros(MERRAlons.shape)
missgrids  = np.zeros(MERRAlons.shape)
falsegrids = np.zeros(MERRAlons.shape)


for IMERGfile, MERRA2file in zip(IMERGfiles, MERRA2files):
    hitgrid_time, missgrid_time, falsegrid_time = hits_misses_false_alarms(IMERGfile, MERRA2file)
    
    hitgrids   += hitgrid_time
    missgrids  += missgrid_time
    falsegrids += falsegrid_time     


biasgrid = (hitgrids + falsegrids) / (hitgrids + missgrids)
PODgrid  = hitgrids / (hitgrids + missgrids)
FARgrid  = falsegrids / (hitgrids + falsegrids)

#report total scores across the whole box
totalBias = (np.nansum(hitgrids) + np.nansum(falsegrids)) / (np.nansum(hitgrids) + np.nansum(missgrids))
totalPOD  = (np.nansum(hitgrids)) / (np.nansum(hitgrids) + np.nansum(missgrids))
totalFAR  = (np.nansum(falsegrids)) / (np.nansum(hitgrids) + np.nansum(falsegrids))

print(f'The total bias score is {totalBias}')
print(f'The total probability of detection is {totalPOD}')
print(f'The total false alarm rate is {totalFAR}')


def letterFromCompType(level, hemi):
    if level == 850 and hemishort(hemi) == 'NH':
        return 'a)'
    elif level == 700 and hemishort(hemi) == 'NH':
        return 'b)'
    elif level == 850 and hemishort(hemi)== 'SH':
        return 'c)'
    elif level == 700 and hemishort(hemi) == 'SH':
        return 'd)'

#plots of the probability of detection per pixel?
def plotGrid(grid, level, hemi, gridtype):    
    fig = plt.figure(figsize=(7,5),dpi=300)
    ax = plt.axes()
    levels = (0.000000000000001, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .55,
              .6, .65, .7, .75, .8, .85, .9, .95, 1)
    #norm = mcolors.BoundaryNorm(levels, 256)
    #print(levels)
    plt.contourf(MERRAlons, MERRAlats, grid, levels=levels, cmap='viridis')
    #set extent would be nice to do so that it's -5 to 5 on both axes
    plt.colorbar()
    ax.set_xlim(-6,6)
    ax.set_ylim(-5,5)
    ax.set_xlabel('Distance from TEW Center (degrees)')
    ax.set_ylabel('Distance from TEW Center (degrees)')
    #ax.set_title(f'{gridtype} for MERRA-2 relative to IMERG \n {hemitext(hemi)} hemisphere at {level} hPa')
    ax.text(-7, 4.5, f'{letterFromCompType(level, hemi)}')
    
    
    ax.grid()
    
    plt.savefig(f'/home/mhollis/Precip_stuff/Figures/PublicationUpdate/SpatialVerification/{gridtype}_{hemishort(hemi)}_{level}hPa.png',
                bbox_inches = "tight", dpi = 300)
    
plotGrid(FARgrid, level, hemi, 'False Alarm Ratio')
plotGrid(PODgrid, level, hemi, 'Probability of Detection')
