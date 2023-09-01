#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 21:42:41 2021

@author: margaret

So, since my goal is to have TEW precip for a year, here's what I'm going to do:
    This code will *not* have loops for tackling multiple years (those should be enough to implement later)
    
    This code will create a file with all of the rain within a 500km radius of the TEW
    This code will create a file with the annual TEW rain
    
    I need other code to create a file with the annual total rain. Or just a different part of this code. idk
    
    These new files will then be used to make some plots
    
UPDATES 29 APRIL 2022
- so turns out that the x and y resolution variables *did* need to be used further down
    - for those interested, they're used in generating values for out_lons and out_lats
      in the boxes that eventually become the eye plots
      
Updates 1 May 2023
- updated input and output file paths for running updated TEW tracks
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

#functions I wrote
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



# open the track for the year
trackfile = f'/data3/TRACK/TRACKOutput/MERRA2/FilteredNoTCUpdated/{level}hPa/post_filter.{year}.{hemishort}.NoTC.nc'
tracks = Dataset(trackfile, mode='r')

# get the dates, tracklats, tracklons
times = tracks['time'][:]
tracktimeunits = tracks['time'].units
dates = num2date(times, units=tracktimeunits)

tracklats = tracks['latitude'][:]
tracklons = tracks['longitude'][:]



#then loop through the dates
#there should be a counter for the loop iterations because of setting up files
#and just because it's nice sometimes

count = 0
samemonth = False

for date in dates:

    #assemble the precip filepath from the current date
    year  = date.year
    month = date.month
    day   = date.day
    hour  = date.hour
    
    #adds a zero to the month if single digit, otherwise converts to string
    if month < 10:
        monthtxt = f'0{month}'
    else:
        monthtxt = str(month)
    
    #open the precip file if it was not left open from the previous iteration
    if samemonth == False:
        #precip file locations
        #precipfile = f'/data/deluge/scratch/IMERG/IMERG_v6b_precip_3hr/imerg.v6b.precipitationCal.8x.{year}{monthtxt}.nc'
        precipfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_3hr/precip_3hr.MERRA2.{year}{monthtxt}.nc'

        #open said precip file
        precip = Dataset(precipfile, mode='r')
    
    
    #get the precip for the date/time of the track
    
    #for IMERG
    #calprecip = precip['precipitationCal'][((day-1)*8+(hour/3))]
    
    #for MERRA-2
    calprecip = precip['PRECTOTCORR'][((day-1)*8+(hour/3))]
    calprecip = calprecip*3600 #converts from mm/s to mm/hr
    
    #precip lons and lats are same from file to file I sure hope so just need to do this the once
    if count == 0:
        #LATS ARE WEIRD AND GO FROM -90 TO 90 (a problem when I assumed that I just had a nice map ready to go no changes needed)
        preciplat = precip['lat'][:]
        preciplon = precip['lon'][:]
        
        preciplons, preciplats = np.meshgrid(preciplon, preciplat)
    
    
    #get the track point at this time
    tracklon = tracks['longitude'][count]
    tracklat = tracks['latitude'][count]
    
    
    #get precip in a 500km radius, or more accurately the mask for it
    great_circle_dist = great_circle(tracklon, tracklat, preciplons, preciplats)
    
    precipmask = np.greater(great_circle_dist, 500)
    
    waveprecipmasked = np.ma.array(calprecip, mask=precipmask, fill_value=(0))
        
    #make a box around said precip, save to a diff var
    wavecenterindices = find_min_idx(great_circle_dist)
    x1 = int(wavecenterindices[1]-dx/2)
    x2 = int(wavecenterindices[1]+dx/2)
    
    y1 = int(wavecenterindices[0]-dy/2)
    y2 = int(wavecenterindices[0]+dy/2)
    
    #waveprecip = waveprecipmasked.filled()[y1:y2, x1:x2] #so this causes problems when things are at that 0/360 line
    waveprecip = waveprecipmasked.filled().take(range(x1,x2), axis=1, mode='wrap')[y1:y2,:]
    
    #if this is the first time through, set up the output file for the TEW rain
    if count == 0:
        #file is the one indexed like track, but with small boxes for the precip
        #will have the box and the offset (location of the center point)
        #outfile = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TRACK-indexed/precip_noTC_{year}_{level}_{hemishort}.nc'
        outfile = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/precip_noTC_{year}_{level}_{hemishort}.nc'
        track_precip_file = netCDF4.Dataset(outfile, 'w', format='NETCDF4')
        
        #setting up vars I need but didn't already have
        numtracks = tracks.dimensions['tracks'].size
        
        #create dimensions
        #follows file.createDimension("dim name", var size or None for unlimited)
        out_track = track_precip_file.createDimension("tracks", numtracks)
        out_record = track_precip_file.createDimension("record", None)
        out_lon = track_precip_file.createDimension("lon", dx)
        out_lat = track_precip_file.createDimension("lat", dy)
        #create variables
        #follows file.createVariable("var name", datatype, dimension (leave off if creating scalar))
        out_TrackID = track_precip_file.createVariable('TRACK_ID', 'i8', ('tracks'))
        out_firstpt = track_precip_file.createVariable('FIRST_PT', 'i8', ('tracks'))
        out_numpts  = track_precip_file.createVariable('NUM_PTS', 'i8', ('tracks'))
        out_index   = track_precip_file.createVariable('index', 'i8', ('record'))
        out_time    = track_precip_file.createVariable('time', 'f8', ('record'))
        out_TEWlat  = track_precip_file.createVariable('latitude', 'f8', ('record'))
        out_TEWlon  = track_precip_file.createVariable('longitude', 'f8', ('record'))
        out_lats    = track_precip_file.createVariable("lat", "f8", ("lat"))
        out_lons    = track_precip_file.createVariable("lon", "f8", ("lon"))
        out_precip  = track_precip_file.createVariable("precip", "f8", ('record', 'lat', 'lon'))
        
        #variable attributes
        
        #TrackID
        out_TrackID.add_fld_num = 0
        out_TrackID.tot_add_fld_num = 0
        out_TrackID.loc_flags = ""
        out_TrackID.cf_role = "trajectory_id"
        
        #firstpt has no attributes        
        
        #numpts        
        out_numpts.long_name = 'number of obs for this trajectory'
        out_numpts.sample_dimension = 'record'
        
        #index has no attributes
        
        #the units for the time outputs need to be consistent
        #but the numbers will get *really* big if we just use the "minutes since" date units
        #so we'll have the output units be the same as the track file
        out_timeUnits = tracktimeunits
        out_time.starndard_name = "time"
        out_time.long_name = "Time"
        out_time.units = out_timeUnits
        out_time.time_calendar = "gregorian"
        out_time.start = tracks['time'].start
        out_time.step = '3'
        
        #Wave center point lats and lons
        out_TEWlon.standard_name = 'longitude'
        out_TEWlon.long_name = 'TEW Longitude'
        out_TEWlon.units = 'degrees_east'
        
        out_TEWlat.standard_name = 'latitude'
        out_TEWlat.long_name = 'Latitude'
        out_TEWlat.units = 'degrees_north'
        
        #Precip grid lats and lons
        #lats
        out_lats.units = "degrees_north"
        out_lats.long_name = "latitude relative to point closest to TEW center"
        
        #lons
        out_lons.units = "degrees_east"
        out_lons.long_name = "longitude relative to point closest to TEW center"
        
        #precip
        out_precip.units = "mm/hr"
        out_precip.long_name = "precipitationCal"
        
        
        
        #global attributes
        track_precip_file.description = 'MERRA2 PRECTOTCORR precipitation within a 500km radius of the TEW' #IMERG v6b
        
        #save the things that are the same from the track file
        out_TrackID[:] = tracks['TRACK_ID'][:]
        out_firstpt[:] = tracks['FIRST_PT'][:]
        out_numpts[:]  = tracks['NUM_PTS'][:]
        out_index[:]   = tracks['index'][:]
        out_time[:]    = tracks['time'][:]
        out_TEWlon[:]  = tracks['longitude'][:]
        out_TEWlat[:]  = tracks['latitude'][:]
        out_lats[:]    = np.arange(-dy/2,dy/2)*yresolution
        out_lons[:]    = np.arange(-dx/2,dx/2)*xresolution
        
		
    #save the smaller box to the appropriate file
    
    out_precip[count,:,:] = waveprecip
    
    #add accumulated TEW rain to the appropriate variable
    if count == 0:
        annualprecip  = waveprecipmasked.filled()
    else:
        annualprecip += waveprecipmasked.filled()
    
    
    #check if the next date in the track file is the same month as the current one
    if count+1 == len(dates):
        break
    elif month == dates[count+1].month:
        samemonth = True
    else:
        samemonth = False
        precip.close()

    #progress report on how things are going
    if count % 100 == 0:
        print(f'got precip from record {count} of {len(dates)}')
    
    #increment the counter up by one
    count += 1

#When turning this loose on everything, this total file should be included in the loop
# so that there's a nice 38-year climo file when the big loops are done

track_precip_file.close()


#second output file is the one for the annual TEW rain
#it will have the TEW total rain
#outfile_accumulated = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/TEW_precip_accumulated/precip_noTC_{year}_{level}_{hemishort}_accumulated.nc'
outfile_accumulated = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/TEW_annual/precip_noTC_{year}_{level}_{hemishort}_accumulated.nc'
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
out_precip.long_name = "precipitationCal"

#the units for the time outputs need to be consistent
#but the numbers will get *really* big if we just use the "minutes since" date units
#so we'll have the output units be the same as the track file

out_lats.units = "degrees_north"
out_lats.long_name = "latitude"

out_lons.units = "degrees_east"
out_lons.long_name = "longitude"

#global attributes
accumulated_precip_file.description = 'Annually accumulated precipitation within 500km radius of a TEW, from IMERG v6b' #MERRA2 PRECTOTCORR'

#now we can actually save the precipitation and the year
out_yrs           = year
out_lats[:]       = preciplat
out_lons[:]       = preciplon
out_precip[0,:,:] = annualprecip*3

accumulated_precip_file.close()






# =============================================================================
# #I guess if looping through years could then tack on a "oh by the way total annual rainfall"
# if hemi==1 and level == 850:
#     yearfiles = f'/data/deluge/scratch/IMERG/IMERG_v6b_precip_3hr/imerg.v6b.precipitationCal.8x.{year}{monthtxt}.nc'
#     #f'/data3/Data_Processed/MERRA2-Margaret/precip_3hr/precip_3hr.MERRA2.{year}.nc'
#     
#     allprecip = np.zeros((len(preciplat), len(preciplon)))
#     for filepath in yearfiles:
#         file = Dataset(filepath, mode='r')
#         print(f'Now working on file {filepath}')
#             
#         fileprecip = file['precipitationCal'][:,:,:]
#         
#         for time in range(len(fileprecip)):
#             allprecip += fileprecip[time,:,:]
#         file.close()
#     
#     outfile_all = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/Annual_precip/all_precip_{year}.nc'
#     #f'/data3/Data_Processed/MERRA2-Margaret/precip_3hr/global_annual/all_precip_{year}.nc'
#     all_precip_file = netCDF4.Dataset(outfile_all, 'w', format='NETCDF4')
#     
#     
#     #create dimensions
#     #follows file.createDimension("dim name", var size or None for unlimited)
#     out_year_all = all_precip_file.createDimension("years", 1)
#     out_lon_all  = all_precip_file.createDimension("lon", len(preciplon))
#     out_lat_all  = all_precip_file.createDimension("lat", len(preciplat))
#     #create variables
#     #follows file.createVariable("var name", datatype, dimension (leave off if creating scalar))
#     out_lats    = all_precip_file.createVariable("lat", "f8", ("lat"))
#     out_lons    = all_precip_file.createVariable("lon", "f8", ("lon"))
#     out_precip  = all_precip_file.createVariable("precip", "f8", ('lat', 'lon'))
#     
#     #variable attributes
#     #units, long name
#     out_precip.units = "mm"
#     out_precip.long_name = "precipitationCal"
#     
#     #the units for the time outputs need to be consistent
#     #but the numbers will get *really* big if we just use the "minutes since" date units
#     #so we'll have the output units be the same as the track file
#     
#     out_lats.units = "degrees_north"
#     out_lats.long_name = "latitude"
#     
#     out_lons.units = "degrees_east"
#     out_lons.long_name = "longitude"
#     
#     #global attributes
#     all_precip_file.description = 'Global annual precipitation from IMERG v6b'
#     
#     #now we can actually save the precipitation and the year
#     out_yrs           = {year}
#     out_lats[:]       = preciplat
#     out_lons[:]       = preciplon
#     out_precip[:,:] = allprecip*3
#     
#     all_precip_file.close()
# =============================================================================




# =============================================================================
# ratio = annualprecip/allprecip
# 
# #plotting time
# fig = plt.figure(dpi=300, figsize=(8,5))
# ax = plt.axes(projection=ccrs.PlateCarree())
# ax.coastlines()
# plt.contourf(preciplons, preciplats, annualprecip*3, levels=(0.00001, 10, 50, 100, 250, 500,1000), colors=('lightgreen','green','gold','darkorange','red','darkred'))
# plt.colorbar(shrink=0.25)
# plt.title('Accumulated Precipitation from TEWs \n in the Northern Hemisphere in 2014, mm')
# plt.savefig('2014TEWprecip.png', dpi=300)
# 
# fig2 = plt.figure(dpi=300, figsize=(8,5))
# ax2 = plt.axes(projection=ccrs.PlateCarree())
# ax2.coastlines()
# plt.contourf(preciplons, preciplats, ratio, levels=(0.00001, 0.01, 0.10, 0.25, 0.50, 0.75, 1), colors=('lightgreen','green','gold','darkorange','red','darkred'))
# plt.colorbar(shrink=0.25)
# plt.title('Fraction of Total Annual Precipitation from TEWs \n in the Northern Hemisphere in 2014')
# plt.savefig('2014Fractionprecip.png', dpi=300)
# =============================================================================
