;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
;************************************************

begin

;************************************************
; open file and read in data
;************************************************
; Read in the file extract OLR

;diri = "/data/shared/NCL-WaveHeating/WKVariance/"
diri = "/data/users/justin/HCproject/Variance/FILTERED/7pt5day/"

DATASET = 4 ; 1 TRMM, 2 IMERG, 3 MERRA-2 Precip, 4 MERRA-2 Heating

;;Grab data

if DATASET.eq.1 then
  PLOT_NAME = "TRMM_filter_plot-anom"
  fili_filt = "TRMM_Precip_kf_filter.nc"
  f_filt = addfile(diri+fili_filt, "r")
  OLR_filter = f_filt->pr(730:,{-31:31},:)           ;mm/hr, 2000-2019 (data starts with 1998)
  time = f_filt->time(730:)
  time@units = "days since 1997-01-01"
end if

if DATASET.eq.2 then
  PLOT_NAME = "IMERG_filter_plot"
  fili_filt = "IMERG_Precip_kf_filter.nc"
  f_filt = addfile(diri+fili_filt, "r")
  OLR_filter = f_filt->precipitationCal(214:7518,{-31:31},:) ;mm/day, 2001-2020 (data runs 6/2006-9/2021)
  OLR_filter = OLR_filter/24                         ;mm/hr
  time = f_filt->time(214:7518)
  time@units = "days since 2000-06-01"
end if

if DATASET.eq.3 then
  PLOT_NAME = "MERRA2_Precip_filter_plot"
  fili_filt = "MERRA_precip_kf_filter.nc"
  f_filt = addfile(diri+fili_filt, "r")
  OLR_filter = f_filt->precipitation(:,{-31:31},:)   ;kg/(m^2*s) or mm/s
  OLR_filter = OLR_filter*3600                       ;mm/hr
  time = f_filt->time
  time@units = "days since 2000-01-01"
end if

if DATASET.eq.4 then
  PLOT_NAME = "MERRA2_DTDT_filter_plot"
  fili_filt = "MERRA_DTDT_kf_filter.nc"
  f_filt = addfile(diri+fili_filt, "r")
  OLR_filter = f_filt->DTDT(:,{-31:31},:)            ;K/day
  OLR_filter = OLR_filter*5/3                        ;mm/hr
  time = f_filt->time
  time@units = "days since 2000-01-01"
end if

printMinMax(OLR_filter,0)

;;Time Coordinate Transform
utc = cd_calendar(time,0)
year = utc(:,0)
month = utc(:,1)
day = utc(:,2)
mind = ind(month.ge.6 .and. month.le.9)

print(year+" "+month+" "+day)

;Calculate the variance of TD filtered precipitation anomalies.
OLR_VAR = (OLR_filter)^2
OLR_VAR = abs(OLR_filter)
print("The anomaly amplitudes are")
printMinMax(OLR_VAR,0)

copy_VarCoords(OLR_filter, OLR_VAR)

test_avg = dim_avg_n_Wrap(dim_avg_n_Wrap(OLR_VAR,2),1)

;print(test_avg)
;printMinMax(test_avg,0)
;exit

;;Choice of Annual or JJAS
;PLOT_PERIOD = "-annual"
PLOT_PERIOD = "-JJAS"
;JJAS_VAR = dim_avg_n_Wrap(OLR_VAR(:,:,:),0)
JJAS_VAR = dim_avg_n_Wrap(OLR_VAR(mind,:,:),0)

print("The time-averaged amplitudes are")
printMinMax(JJAS_VAR,0)

wks = gsn_open_wks("png",PLOT_NAME+PLOT_PERIOD+"_abs")
;gsn_define_colormap(wks,"BlAqGrYeOrRe")
gsn_define_colormap(wks,"GMT_wysiwyg")
  
res                     = True
res@gsnDraw             = True           ; don't draw
res@gsnFrame            = True           ; don't advance frame
res@cnInfoLabelOn       = False           ; turn off cn info label
res@cnFillOn            = True            ; turn on color
res@gsnSpreadColors     = True            ; spread out color table
res@lbLabelBarOn        = True           ; turn off individual cb's

;colors
res@cnLevelSelectionMode 	= "ManualLevels"
res@cnMinLevelValF            	= 0.0
res@cnMaxLevelValF		= 0.30
res@cnLevelSpacingF		= 0.01
; res@cnLevelSelectionMode 	= "ExplicitLevels"
;res@cnLevels = (/-24,-20,-16,-12,-8,-4,4,8,12,16,20,24/)
;res@cnLevels = (/-36,-30,-24,-18,-12,-6,6,12,18,24,30,36/)  

res@mpMinLatF            = -30      ; range to zoom in on
res@mpMaxLatF            =  30.
res@mpMinLonF            =  -180.
res@mpMaxLonF            = 180.
 
;turn contour lines off
res@cnLineLabelsOn = False
res@cnLinesOn   = False

; font sizes on x/y labels
res@tmEqualizeXYSizes = True
res@tmXBLabelFontHeightF = .01

; text info for figure labels
res@gsnRightString = ""
res@gsnLeftStringOrthogonalPosF = .1
res@gsnLeftStringFontHeightF = 15
res@tmXBMajorLengthF = 0.008
res@pmLabelBarOrthogonalPosF = 0.3
res@lbTitleOn        = True                  ; turn on title
res@lbTitlePosition = "Bottom"
;res@lbTitleString    = "(mm hr~S~-1~N~)~S~2~N~"
res@lbTitleString    = "mm hr~S~-1~N~"
res@lbTitleFontHeightF= .013
res@lbLabelStride = 2
res@lbAutoManage=False
res@lbJustification = "CenterCenter"
;res@lbLeftMarginF = 0.0 
;res@lbTopMarginF = 0.9
;res@lbBoxMinorExtentF = 0.5
;res@lbLabelOffsetF = 4
res@lbTitleJust = "BottomCenter"
;res@lbTitleOffsetF = 0.8

if DATASET.eq.1 then
  res@gsnLeftString = "(a) TRMM 3B42"
end if
if DATASET.eq.2 then
  res@gsnLeftString = "(b) IMERG"
end if
if DATASET.eq.3 then
  res@gsnLeftString = "(c) MERRA-2 Precipitation"
end if
if DATASET.eq.4 then
  res@gsnLeftString = "(d) MERRA-2 Latent Heating"
end if

plot = gsn_csm_contour_map_ce(wks,JJAS_VAR,res)
 
end
