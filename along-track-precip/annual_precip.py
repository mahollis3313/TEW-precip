#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 15:14:29 2021

@author: margaret
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
import netCDF4
from netCDF4 import Dataset, num2date

start = 2001
stop = 2020

years = np.arange(start, (stop+1), step=1)

#getlatlonsfile = '/data/deluge/scratch/IMERG/IMERG_v6b_precip_3hr/imerg.v6b.precipitationCal.8x.200101.nc'
getlatlonsfile = '/data3/Data_Processed/MERRA2-Margaret/precip_3hr/precip_3hr.MERRA2.200101.nc'
getlatlons = Dataset(getlatlonsfile, mode='r')

preciplat = getlatlons['lat'][:]
preciplon = getlatlons['lon'][:]

getlatlons.close()

for year in years:
    #get list of files for a given year
    #yearfiles = sorted(glob.glob(f'/data/deluge/scratch/IMERG/IMERG_v6b_precip_3hr/imerg.v6b.precipitationCal.8x.{year}*.nc'))
    yearfiles = sorted(glob.glob(f'/data3/Data_Processed/MERRA2-Margaret/precip_3hr/precip_3hr.MERRA2.{year}*.nc'))
    
    allprecip = np.zeros((len(preciplat), len(preciplon)))
    for filepath in yearfiles:
        file = Dataset(filepath, mode='r')
        print(f'Now working on file {filepath}')
            
        #for IMERG
        #fileprecip = file['precipitationCal'][:,:,:]
        
        #for MERRA-2
        fileprecip = file['PRECTOTCORR'][:,:,:]
        fileprecip = fileprecip*3600
        
        allprecip += np.nansum(fileprecip, axis=0)
        file.close()
    
    #outfile_all = f'/data/deluge/scratch/IMERG/TEW_Precip_updated/Annual_precip/all_precip_{year}.nc'
    outfile_all = f'/data3/Data_Processed/MERRA2-Margaret/precip_TEW_updated/global_annual/all_precip_{year}.nc'
    all_precip_file = netCDF4.Dataset(outfile_all, 'w', format='NETCDF4')
    
    
    #create dimensions
    #follows file.createDimension("dim name", var size or None for unlimited)
    out_year_all = all_precip_file.createDimension("years", 1)
    out_lon_all  = all_precip_file.createDimension("lon", len(preciplon))
    out_lat_all  = all_precip_file.createDimension("lat", len(preciplat))
    #create variables
    #follows file.createVariable("var name", datatype, dimension (leave off if creating scalar))
    out_lats    = all_precip_file.createVariable("lat", "f8", ("lat"))
    out_lons    = all_precip_file.createVariable("lon", "f8", ("lon"))
    out_precip  = all_precip_file.createVariable("precip", "f8", ('lat', 'lon'))
    
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
    #all_precip_file.description = 'Global annual precipitation from IMERG v6b'
    all_precip_file.description = 'Global annual precipitation from MERRA2 PRECTOTCORR'
    
    #now we can actually save the precipitation and the year
    out_yrs           = {year}
    out_lats[:]       = preciplat
    out_lons[:]       = preciplon
    out_precip[:,:] = allprecip*3
    
    all_precip_file.close()