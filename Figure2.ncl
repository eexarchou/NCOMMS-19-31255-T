; ----------------------------------------------------------------
; Code script presented in the publication:
;
; "Impact of Equatorial Atlantic Variability on ENSO Predictive Skill"
; by Eleftheria Exarchou, Pablo Ortega, Belen Rodrıguez de Fonseca,
; Teresa Losada, Irene Polo, and Chloe Prodhomme 
; in Nature Communications
;
; NMME data have been downloaded from the IRI database
; (http://iridl.ldeo.columbia.edu/SOURCES/.Models/.NMME/). EUROSIP and
; ERA-interim data have been downloaded from the MARS computing facility
; of ECMWF (https://www.ecmwf.int/en/forecasts/). 
; Data from EC-Earth predictions are available on request from 
; Eleftheria Exarchou. 
;
; HadISSTv1.1 and ERSSTv4 are publicly available
; here:https://www.metoffice.gov.uk/hadobs/hadisst/data/ and 
; https://www1.ncdc.noaa.gov/pub/data/.
;
; The data have been formatted to follow a particular data
; structure that allows the function "Load" by the s2dverification
; package to read the data. All data have been loaded and saved as an R object.      
;
; Additional information:
; Period: 1981 - 2011
; Ensemble members: all available (information on Supplementary Table 1)
;
; Date: 15 Dec 2020
; Author: Eleftheria Exarchou (eleftheria.exarchou@gmail.com)
; ----------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

begin
;----------------------------------------------------------------------
; Data is hard-coded here. It was originally calculated by running
; a compute-intensive script.
;----------------------------------------------------------------------

 exps = (/ \ 
             "cfs_v2"        , "ecmwf/system4",       "cancm4"   ,\       
             "cancm3"        ,   "cm2p5-flor-a06",   "cm2p5-flor-b01" ,  "rsmas-ccsm4"  ,\   
             "rsmas-ccsm3"   ,   "cm2p1"         ,   "cm2p1-aer04"    ,  "echam4p5"     ,\    
             "gmao-062012"   , " ecmwf/system5 " ,   "ncep/system2"   ,\
             "EC-Earth"      ,\    
             "obs", "multi-model mean" /)       

; -- Skill and telecoonection strength indices are originally calculated
;             by running a computationally expensive R-script 

skill_atl3_jja   = asciiread("data1.asc",-1,"float") 
skill_nino3_jja = asciiread("data2.asc",-1,"float") 
skill_nino3_son = asciiread("data3.asc",-1,"float") 
skill_nino3_ndj = asciiread("data4.asc",-1,"float") 
cor_atl3jja_nino3jja  = asciiread("data8.asc",-1,"float") 
cor_atl3jja_nino3son  = asciiread("data9.asc",-1,"float") 
cor_atl3jja_nino3ndj  = asciiread("data10.asc",-1,"float") 
                     
                    
skill_atl3_jja@_FillValue = 0. 
skill_nino3_jja@_FillValue = 0. 
skill_nino3_son@_FillValue = 0. 
skill_nino3_ndj@_FillValue = 0. 

cor_atl3jja_nino3jja@_FillValue = 0. 
cor_atl3jja_nino3son@_FillValue = 0. 
cor_atl3jja_nino3ndj@_FillValue = 0. 

 x1 = skill_atl3_jja
 y1 = skill_nino3_jja

 x2 = skill_atl3_jja
 y2 = skill_nino3_son

 x3 = skill_atl3_jja
 y3 = skill_nino3_ndj
 
 x4 = (/cor_atl3jja_nino3son(0:14)/) 
 y4 = (/skill_nino3_son(0:14)/) 
 
 x5 = (/cor_atl3jja_nino3ndj(0:14)/) 
 y5 = (/skill_nino3_ndj(0:14) /)
 
 rc1    = regline_stats(x1,y1 )          ; slope     
 rc2    = regline_stats(x2,y2 )          ; slope     
 rc3    = regline_stats(x3,y3 )          ; slope     
 rc4    = regline_stats(x4,y4 )          ; slope     
 rc5    = regline_stats(x5,y5 )          ; slope     

 data1       = new ( (/2,dimsizes((x1 ))/), typeof((x1 )))
 data2       = new ( (/2,dimsizes((x1 ))/), typeof((x1 )))
 data3       = new ( (/2,dimsizes((x1 ))/), typeof((x1 )))
 data4       = new ( (/2,dimsizes((x1 ))/), typeof((x1 )))
 data5       = new ( (/2,dimsizes((x1 ))/), typeof((x1 )))
 
 data1(0,:)  = y1
 data2(0,:)  = y2
 data3(0,:)  = y3
 data4(0,:)  = y4
 data5(0,:)  = y5
 
 data1(1,:)  = rc1 * ( x1-rc1@xave) + rc1@yave
 data2(1,:)  = rc2 * ( x2-rc2@xave) + rc2@yave
 data3(1,:)  = rc3 * ( x3-rc3@xave) + rc3@yave
 data4(1,:)  = rc4 * ( x4-rc3@xave) + rc4@yave
 data5(1,:)  = rc5 * ( x5-rc3@xave) + rc5@yave

y11= array_append_record( cor_atl3jja_nino3jja, \ 
array_append_record(cor_atl3jja_nino3son,cor_atl3jja_nino3ndj ,0 ),0) ;; Jun ini 

x11= (/1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,1,1, \ 
       2, 2,2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2,2,2, \
       3, 3,3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,3,3,3,3  /)
x22= (/1, 1,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ,1,1,1,1, \ 
       2, 2,2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,2,2,2,2  /)
      
x33=x11 
x43=x22 

y11@_FillValue = 0. 

  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.011                     ; font smaller. default big
  txres@txFont        = "helvetica"                     ; font smaller. default big
  
  txres2               = True                     ; text mods desired
  txres2@txFontHeightF = 0.02                     ; font smaller. default big
  txres2@txFont        = "helvetica-bold"                     ; font smaller. default big



  colors  = (/"white","black",\ 
             "deeppink2",  "forestgreen","firebrick",\
             "dodgerblue2","darkorange2","forestgreen","deeppink2",\
             "darkorange2","purple","red", "brown", "deeppink2", "blue", "red","blue", "black", "grey", \ ; here 
			 "deeppink2", "forestgreen","firebrick",\                           
             "dodgerblue2","darkorange2","forestgreen","deeppink2",\           
             "darkorange2","purple","red", "brown", "deeppink2", "blue", "red","blue", "black","grey",\  ; here ;  /)
             "deeppink2", "forestgreen","firebrick",\                           
             "dodgerblue2","darkorange2","forestgreen","deeppink2",\          
             "darkorange2","purple","red", "brown", "deeppink2", "blue", "red","blue", "black", "grey"/)

;----------------------------------------------------------------------
; Open an eps  file and define the color table
;----------------------------------------------------------------------
 wks = gsn_open_wks ("eps","Figure2")
  gsn_define_colormap(wks,colors)
  plot = new(6,graphic)
;----------------------------------------------------------------------
; Set some plot resources
;----------------------------------------------------------------------
  res                    = True

  res@gsnDraw            = False       ; Don't draw plot
  res@gsnFrame           = False       ; Don't advance frame

  res@vpWidthF          = 0.5
  res@vpHeightF         = 0.5


;---XY marker resources
  res@xyMarkLineMode     = "Markers"
  res@xyMarkerSizeF      = 0.03
  res@xyMarkerThicknessF = 3.
  res@xyMonoMarkerColor  = False          ; Allow different colors for markers
  res@xyMarkerColors     = colors(2:)     ; It's okay to list more than

  res@tmXBMinorOn  = False
  res@tmYLMinorOn  = False
  res@tmYLAutoPrecision = False
  res@tmYLPrecision = 2
  res@tmYLFormat   = "0@*+^sg"
;---Title resources
;  res@tiMainString       = "30 av Annual % Change P vs Delta T~C~" + \
;                           "               for the WetSouth"
  res@tiMainFontHeightF    = 0.02
  res@tmXBLabelFontHeightF = 0.02
  res@tiXAxisFontHeightF   = 0.03
  res@tiYAxisFontHeightF   = 0.03

  ;res@xyMarkers         = (/ 2,3,4,5,6,7,8,9,12,14,15, 16,4,5,16, 16 /)
  res@xyMarkers          = (/2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16,\  
                             2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16,\ 
                             2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16/)



;; ======  RES2 ================================

  res2                    = True
  res2@gsnDraw            = False       ; Don't draw plot
  res2@gsnFrame           = False       ; Don't advance frame
  res2@tmYLMinorOn  = False

  res2@trXMinF            =   0.                  ; min value on x-axis
  res2@trXMaxF            =   4.                  ; max value on x-axis
  res2@trYMinF            =  -0.55                   ; min value on x-axis
  res2@trYMaxF            =   .1                     ; max value on x-axis
  res2@tmXBMode         = "Manual"
  res2@tmXBTickSpacingF = 1.

  res2@tmXBMode        = "Explicit"	
  res2@tmXBValues      = (/0.,0.5, 1., 1.5,  2., 2.5, 3./)

  res2@tmXBLabels      = (/"", "",  "Nino3.4~B4~JJA ~C~ Atl3~B4~JJA", "", "Nino3.4~B4~SON ~C~ Atl3~B4~JJA", "", "Nino3.4~B4~NDJ ~C~ Atl3~B4~JJA"  /)
  res2@tmXBLabelFontHeightF = 0.02
  res2@vpWidthF          = 0.5
  res2@vpHeightF         = 0.5


;---XY marker resources
  res2@xyMarkLineMode     = "Markers"
  res2@xyMarkerSizeF      = 0.025
  res2@xyMarkerThicknessF = 3.
  res2@xyMonoMarkerColor  = False          ; Allow different colors for markers
  res2@xyMarkerColors     = colors(2:)     ; It's okay to list more than
                                          ; you need here
  res2@xyMarkers          = (/2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16,\  
                              2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16,\ 
                              2,3,4,5,6,7,8,9,12,14,15,16,4,5,2,16, 16/)
;                           "               for the WetSouth"
  res2@tiMainFontHeightF  = 0.02
  res2@tiXAxisFontHeightF  = 0.03
  res2@tiYAxisFontHeightF  = 0.03
  res2@tmXBLabelFontHeightF = 0.02
  res2@gsnMaximize       = True                ; maximize plots
  res2@tiYAxisString      = "Correlation Nino3.4 / ATL3"
 
 ;;; === END of RES2 ============================
  
;  res@trYMinF            =  0.7                   ; min value on x-axis
;  res@trYMaxF            =  1.                     ; max value on x-axis

  res@gsnMaximize       = True                ; maximize plots
;----------------------------------------------------------------------
; Create the scatter plot, but don't draw it yet.
;----------------------------------------------------------------------
  ;; AXIS 
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXTBorderOn = False
  res@tmYRBorderOn = False
  res2@tmXTOn = False
  res2@tmYROn = False
  res2@tmXTBorderOn = False
  res2@tmYRBorderOn = False

  res@tiXAxisString      = "Skill in JJA ATL3"
  res@tiYAxisString      = "Skill in JJA Nino3"
   plot(0) = gsn_csm_xy (wks,transpose((/x1,x1/)), transpose((/y1,y1/)),res)
   gsn_text_ndc(wks,"R = "+sprintf("%5.2f",  rc1@r)+", p = "+sprintf("%5.3f",rc1@pval(1)), 0.22,0.66,txres)
   gsn_text_ndc(wks,"a", 0.1 , 0.87,txres2)

  res@tiYAxisString      = "Skill in SON Nino3"
   plot(1) = gsn_csm_xy (wks,transpose((/x2,x2/)), transpose((/y2,y2/)),res)
   gsn_text_ndc(wks,"R = "+sprintf("%5.2f",  rc2@r)+", p = "+sprintf("%5.4f",rc2@pval(1)), 0.6 , 0.66,txres)
   gsn_text_ndc(wks,"b", 0.43 , 0.87,txres2)
   
   res@tiYAxisString      = "Skill in NDJ Nino3"
   plot(2)  = gsn_csm_xy (wks,transpose((/x3,x3/)), transpose((/y3,y3/)),res)
   gsn_text_ndc(wks,"R = "+sprintf("%5.2f",  rc3@r)+", p = "+sprintf("%5.3f",rc3@pval(1)), 0.93 , 0.66,txres)
   gsn_text_ndc(wks,"c", 0.77 , 0.87, txres2)

   res2@tiYAxisString      = "Correlation Nino3 / ATL3"
   plot(3) = gsn_csm_xy (wks,transpose((/x11,x11/)), transpose((/y11,y11/)),res2)
   gsn_text_ndc(wks,"d", 0.1, 0.55,txres2)

   res@trXMinF            =   -0.55                     ; max value on x-axis
   res@trYMinF            =   .68                     ; max value on x-axis

   res@tiXAxisString      = "Cors Nino3(NDJ)/ATL3(JJA)"
   res@tiYAxisString      = "Skill in SON Nino3"
   res@gsnXRefLine =  avg(cor_atl3jja_nino3ndj (15))
   plot(4) = gsn_csm_xy (wks,transpose((/x4,x4/)), transpose((/y4,y4/)),res)
   gsn_text_ndc(wks,"R = -"+sprintf("%5.2f",  rc4@r)+", p = "+sprintf("%6.4f",rc4@pval(1)), 0.55 , 0.34,txres)
   gsn_text_ndc(wks,"e", 0.43 , 0.55,txres2)
   
   res@tiYAxisString      = "Skill in NDJ Nino3"
   res@gsnXRefLine =  avg( cor_atl3jja_nino3ndj (15))
   plot(5)  = gsn_csm_xy (wks,transpose((/x5,x5/)), transpose((/y5,y5/)),res)
   gsn_text_ndc(wks,"R = -"+sprintf("%5.2f",  rc5@r)+", p = "+sprintf("%5.3f",rc5@pval(1)), 0.82 , 0.34,txres)
   gsn_text_ndc(wks,"f", 0.77 , 0.55,txres2)
   
   delete (res@xyMarkers)
   delete(res@xyMarkerColors) 

   
 resp2                = True                   ; modify the panel plot
 resp2@gsnPanelBottom = 0.17 	; shrink panel plot so that there is extra room for the left y-axis label
 resp2@gsnFrame       = False                  ; don't advance panel plot
 resp2@gsnMaximize    = True                ; maximize plots
 gsn_panel(wks,plot,(/2,3/),resp2)               ; now draw as one plot
 
;; LEGEND

 labels2 = (/ \ 
             "cfs_v2"        , "ecmwf/system4",       "cancm4"   ,\       
             "cancm3"        ,   "cm2p5-flor-a06",   "cm2p5-flor-b01" ,  "rsmas-ccsm4"  ,\   
             "rsmas-ccsm3"   ,   "cm2p1"         ,   "cm2p1-aer04"    ,  "echam4p5"     ,\    
             "gmao-062012"   , " ecmwf/system5 " ,   "ncep/system2"   ,\
             "EC-Earth"      ,\    
             "obs", "multi-model mean" /)       
  
  markers1 = (/  2,3,4,5,6,7,8,9,12,14,15, 16,4,5, 2, 16,16,  16 /)
  ypos1    = (/  0.08,   0.095,  0.105,    0.14,  0.16 ,  0.18 ,  0.2,  0.22,    0.24, \ 
  0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38 , 0.4 , 0.4, 0.4/)


  colors2  = (/"white","black",\ 
            "deeppink2",  "forestgreen","firebrick",\
              "dodgerblue2","darkorange2","forestgreen","deeppink2",\
              "darkorange2","purple","red", "brown", "deeppink2", "blue", "red", "blue",\
			  "black", "grey" /) ; here 


  mkres               = True         ; Marker resources
  mkres@gsMarkerThicknessF = 3.
  mkres@gsMarkerSizeF      = 0.015
  txres@txFontHeightF = 0.015
  txres@txJust        = "CenterLeft"

  ypos1(0 ) = 0.12
  do i = 0,3 ; 14
     mkres@gsMarkerColor      = colors2(2+i)
     mkres@gsMarkerIndex      = markers1(i)
     gsn_polymarker_ndc(wks,           0.14, ypos1(i),mkres)
     gsn_text_ndc      (wks,labels2(i),0.16,ypos1(i),txres)
	 ypos1(i+1) = ypos1(i)  + 0.03  
  end do

  ypos1(4 ) = 0.12
  do i = 4,7 ; 14
     mkres@gsMarkerColor      = colors2(2+i)
     mkres@gsMarkerIndex      = markers1(i)
     gsn_polymarker_ndc(wks,           0.38, ypos1(i),mkres)
     gsn_text_ndc      (wks,labels2(i),0.40,ypos1(i),txres)
	 ypos1(i+1) = ypos1(i)  + 0.03  
  end do
  
  ypos1(8 ) = 0.12
  do i = 8,11 ; 14
     mkres@gsMarkerColor      = colors2(2+i)
     mkres@gsMarkerIndex      = markers1(i)
     gsn_polymarker_ndc(wks,           0.58, ypos1(i),mkres)
     gsn_text_ndc      (wks,labels2(i),0.60,ypos1(i),txres)
	 ypos1(i+1) = ypos1(i)  + 0.03  
  end do

  
  ypos1(12 ) = 0.12
  do i = 12,16 ; 14
     mkres@gsMarkerColor      = colors2(2+i)
     mkres@gsMarkerIndex      = markers1(i)
     gsn_polymarker_ndc(wks,           0.78, ypos1(i),mkres)
     gsn_text_ndc      (wks,labels2(i),0.8,ypos1(i),txres)
	 ypos1(i+1) = ypos1(i)  + 0.03  
  end do

          end
