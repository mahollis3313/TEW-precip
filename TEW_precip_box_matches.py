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
import netCDF4
from netCDF4 import Dataset, num2date

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
    return int(round(k/ncol)), k%ncol



# # # # #
# This is where the level, year, hemisphere information will be taken in
# # # # #
year = int(sys.argv[1])
level = int(sys.argv[2])
hemi = int(sys.argv[3])

if hemi == 1:
    hemitext = 'northern'
    hemishort = 'NH'
elif hemi == 2:
    hemitext = 'southern'
    hemishort = 'SH'
else:
    print('Hemisphere is not valid.')
    sys.exit()


# Things that should be easily accessible to change but are hard coded
    # numbers either about the dataset or picked for the dataset
    # because I want this code to be adaptable to other datasets
    # but am too lazy to figure out a way to do it all with attributes
dx = 20 #130 #20 
dy = 20 #100 #20
xresolution = 0.625 #0.1 #0.625
yresolution = 0.5 #0.1 #0.5

# get the indices for the matched points
indices = np.loadtxt(f'/home/mhollis/track_stats/CollocationUpdated/match_indices_{year}_{hemishort}.csv', delimiter=',').astype('int')
if level == 850:
    levindex = 0
elif level == 700:
    levindex=1
else:
    print('invalid level')
    sys.exit()

# open the track for the year
trackfile = f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.{hemishort}.NoTC.nc'
tracks = Dataset(trackfile, mode='r')

# get the dates, tracklats, tracklons
# dates should be at the matched indices only
times = tracks['time'][:]
tracktimeunits = tracks['time'].units
dates = num2date(times, units=tracktimeunits)

tracklats = tracks['latitude'][:]
tracklons = tracks['longitude'][:]




#then loop through the dates
#there should be a counter for the loop iterations because of setting up files
#and just because it's nice sometimes

count = 0

for index in indices[:,levindex]:
    
    date = dates[index]
    #assemble the precip filepath from the current date
    month = date.month
    day   = date.day
    hour  = date.hour
    
    #adds a zero to the month if single digit, otherwise converts to string
    if month < 10:
        monthtxt = f'0{month}'
    else:
        monthtxt = str(month)
    
    #open the precip file if it was not left open from the previous iteration
    
    #precipfile = f'/data/deluge/scratch/IMERG/IMERG_v6b_precip_3hr/imerg.v6b.precipitationCal.8x.{year}{monthtxt}.nc'
    precipfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_3hr/precip_3hr.MERRA2.{year}{monthtxt}.nc'
    
    #open said precip file
    precip = Dataset(precipfile, mode='r')
    
    
    #get the precip for the date/time of the track
    
    #for IMERG
    #calprecip = precip['precipitationCal'][((day-1)*8+(hour/3))]
    
    #for MERRA-2
    calprecip = precip['PRECTOTCORR'][((day-1)*8+(hour/3))]
    calprecip = calprecip * 3600
    
    #precip lons and lats are same from file to file I sure hope so just need to do this the once
    if count == 0:
        #LATS ARE WEIRD AND GO FROM -90 TO 90 (a problem when I assumed that I just had a nice map ready to go no changes needed)
        preciplat = precip['lat'][:]
        preciplon = precip['lon'][:]
        
        preciplons, preciplats = np.meshgrid(preciplon, preciplat)
    
    
    #get the track point at this time
    tracklon = tracks['longitude'][index]
    tracklat = tracks['latitude'][index]
    
    
    #get precip in a 500km radius, or more accurately the mask for it
    great_circle_dist = great_circle(tracklon, tracklat, preciplons, preciplats)
    
    precipmask = np.greater(great_circle_dist, 500)
    
    waveprecipmasked = np.ma.array(calprecip, mask=precipmask, fill_value=(0))

    
    #add accumulated TEW rain to the appropriate variable
    if count == 0:
        annualprecip  = waveprecipmasked.filled()
    else:
        annualprecip += waveprecipmasked.filled()
    
    
    #check if the next date in the track file is the same month as the current one
    if count+1 == len(indices[:,levindex]):
        break
    else:
        precip.close()

    #progress report on how things are going
    if count % 100 == 0:
        print(f'got precip from record {count} of {len(indices)}')
    
    #increment the counter up by one
    count += 1

#When turning this loose on everything, this total file should be included in the loop
# so that we can cdo and make a nice 38-year climo file when the big loops are done


#second output file is the one for the annual TEW rain
#it will have the TEW total rain
#outfile_accumulated = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_matched_precip_accumulated/precip_noTC_matched_{year}_{level}_{hemishort}_accumulated.nc'
outfile_accumulated = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_collocated/precip_noTC_matched_{year}_{level}_{hemishort}_accumulated.nc'
accumulated_precip_file = netCDF4.Dataset(outfile_accumulated, 'w', format='NETCDF4')

#create dimensions
#follows file.createDimension("dim name", var size or None for unlimited
out_yrs = accumulated_precip_file.createDimension("years", 1)
out_lon = accumulated_precip_file.createDimension("lon", len(preciplon))
out_lat = accumulated_precip_file.createDimension("lat", len(preciplat))
#create variables
#follows file.createVariable("var name", datatype, dimension (leave off if creating scalar)
out_lats    = accumulated_precip_file.createVariable("lat", "f8", ("lat"))
out_lons    = accumulated_precip_file.createVariable("lon", "f8", ("lon"))
out_precip  = accumulated_precip_file.createVariable("precip", "f8", ('years', 'lat', 'lon'))

#variable attributes
#units, long name
out_precip.units = "mm"
#out_precip.long_name = "precipitationCal"
out_precip.long_name = 'PRECTOTCORR'


#the units for the time outputs need to be consistent
#but the numbers will get *really* big if we just use the "minutes since" date units
#so we'll have the output units be the same as the track file

out_lats.units = "degrees_north"
out_lats.long_name = "latitude"

out_lons.units = "degrees_east"
out_lons.long_name = "longitude"

#global attributes
#accumulated_precip_file.description = 'Annually accumulated precipitation within 500km radius of a TEW, from IMERG v6b'
accumulated_precip_file.description = 'Annually accumulated precipitation within 500km radius of vertically collocated TEWs, from MERRA-2 PRECTOTCORR'

#now we can actually save the precipitation and the year
out_yrs           = year
out_lats[:]       = preciplat
out_lons[:]       = preciplon
out_precip[0,:,:] = annualprecip*3

accumulated_precip_file.close() 



