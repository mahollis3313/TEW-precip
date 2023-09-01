;load "$NCARG_ROOT/lib/nc:wqarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

load "/data/shared/NCLcodes/wkSpaceBC.ncl"
 
  var      = "PRECT"
  case     = "MERRA Precip 20 update"
  reanalysis_list = asciiread("/data/shared/NCLcodes/MERRA_Precip_Annual_20_update.txt", -1, "string")
  diro     = "./"        ; output dir: location of plots

  latN     = 15          ; latBound
  latS     = -15       ; symmetric about the equator

  nDayWin  =  80       ; Wheeler-Kiladis [WK] temporal window length (days)
  nDaySkip = -60        ; time (days) between temporal windows [segments]

  opt      = True
  opt@debug= True
  opt@pltType     = "png"
  opt@cnLinesOn   = False   ; turn off contour lines

                     ; specify the plot intervals

  spd      = 1       ; samples per day (6 hour sampling rate)
  ;var      = "OLR"  ; user specified

 f = addfiles(reanalysis_list, "r")
 ListSetType(f, "join") 
  
  lat       = f[0]->latitude
  lon       = f[0]->longitude
  ;time      = f[0]->time
  ;recip    = f[:]->precipitation  ;precip2    = dim_avg_n_Wrap(precip,1)
  precip    = f[:]->precip
  ;DTDT    = dim_avg_n_Wrap(LH_Vavg,1)
  precip@_FillValue = 0
  delete(precip@_FillValue)

  wkSpaceTime (precip(:,:,:), diro, case, var             \
              ,latN, spd, nDayWin, nDaySkip, opt  )

