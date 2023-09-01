;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/kf_filter.ncl"
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
; Main
;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
begin

  diri = "/data/reanalysis/MERRA2/Daily/Precip_Daily_Anom/"

; set the wave parameters
; modify for TD region
  tMin =2
  tMax = 7.5
  kMin = -20
  kMax = -6
  hMin = 8
  hMax = 90
  waveName = "TD"  
  obsPerDay = 1

  hMin@_FillValue = default_fillvalue("integer")
  hMax@_FillValue = default_fillvalue("integer")

  hMin = hMin@_FillValue
  hMax = hMax@_FillValue

;Read in precipitation anomalies (same as those calculated for WK spectra)

  setfileoption("nc", "SuppressClose", False)

  filenames = systemfunc("ls "+diri+"MERRA*.nc")
  a = addfiles(filenames,"r")
  ListSetType(a, "join")

  origData = a[:]->precip(:,{-40:40},:) 
  origData!0 = "time"
  origData!1 = "lat"
  origData!2 = "lon"

  time = ispan(0, dimsizes(origData(:,0,0))-1, 1)
  time!0 = "time"
  time&time = time
  time@units = "days since 2000-01-01"

  origData&time = time

  latitude = a[0]->latitude({-40:40})
  longitude = a[0]->longitude

  dims = dimsizes(origData)
  tdims = dims(0)
  latdims = dims(1)
  londims = dims(2)
  
  flter_data = new((/tdims,latdims,londims/), typeof(origData))

  do n=0, latdims-1
    print("latint = "+n)
    
    filtLat = kf_filter(origData(:,n,:), obsPerDay, tMin, tMax, kMin, kMax, hMin, hMax, waveName) 

    flter_data(:,n,:) = filtLat
 
    delete(filtLat)	
  end do

;;Write the NetCDF output

  copy_VarCoords(origData,flter_data)
  flter_data@wavenumber = (/ kMin, kMax /)
  flter_data@period     = (/ tMin, tMax /)
  flter_data@waveName   = "TDfilter"

  flter_data!1 = "latitude"
  flter_data!2 = "longitude"

  printVarSummary(flter_data)
  printVarSummary(origData)
  printMinMax(flter_data,0)

  dims = dimsizes(flter_data)
  ntime = dims(0)
  nlat = dims(1)
  nlon = dims(2)

;  system("rm MERRA_precip_kf_filter.nc")
  fout  = addfile ("MERRA_precip_kf_filter.nc", "c")

  fAtt               = True            ; assign file attributes
  fAtt@title         = "MERRA-2 TD Wave Filtered Daily Anomalies from Smoothed Climatology"
  fAtt@Conventions   = "None"
  fAtt@creation_date = systemfunc ("date")
  fileattdef(fout, fAtt)               ; copy file attributes

  setfileoption(fout,"DefineMode",True)
  dimNames = (/"time", "latitude", "longitude"/)
  dimSizes = (/ntime, nlat , nlon /)
  dimUnlim = (/-1, False, False/)
  filedimdef(fout,dimNames,dimSizes,dimUnlim)

  filevardef(fout, "time" ,typeof(time),getvardims(time))
  filevardef(fout, "latitude"  ,typeof(latitude),getvardims(latitude))
  filevardef(fout, "longitude"  ,typeof(longitude),getvardims(longitude))
  filevardef(fout, "precipitation"   ,typeof(flter_data)  ,getvardims(flter_data))

  filevarattdef(fout,"time", time)
  filevarattdef(fout,"latitude", latitude)
  filevarattdef(fout,"longitude", longitude)
  filevarattdef(fout,"precipitation", flter_data)

  fout->time    = (/time/)
  fout->latitude    = (/latitude/)
  fout->longitude    = (/longitude/)
  fout->precipitation     = (/flter_data/)

end
