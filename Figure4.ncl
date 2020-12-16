; ----------------------------------------------------------------
; Code script presented in the publication:
;
; "Impact of Equatorial Atlantic Variability on ENSO Predictive Skill"
; by Eleftheria Exarchou, Pablo Ortega, Belen Rodrıguez de Fonseca,
; Teresa Losada, Irene Polo, and Chloe Prodhomme 
; in Nature Communications
;
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
; This script reads the correlations between JJA ATL3 and the spatial
; fields of SST/uas/vas for one experiment, calculated by the
; computiationally expensive script Cors-Figure4.R and saved in netCDF
; format, and plots them  
;
; Date: 15 Dec 2020
; Author: Eleftheria Exarchou (eleftheria.exarchou@gmail.com)
; ----------------------------------------------------------------

begin
;; -------- ctr ----------------------------------------
        indir="Data"
        c1 = addfile(indir+"correlation_tos_ctr_son-lon.nc","r")
        c2 = addfile(indir+"correlation_uas_ctr_son-lon.nc","r")
        c3 = addfile(indir+"correlation_vas_ctr_son-lon.nc","r")
;; -------- nud1 ----------------------------------------
        
        g1 = addfile(indir+"correlation_tos_nud1_son-lon.nc","r")
        g2 = addfile(indir+"correlation_uas_nud1_son-lon.nc","r")
        g3 = addfile(indir+"correlation_vas_nud1_son-lon.nc","r")
                     
;; -------- OBS ----------------------------------------
        k1 = addfile(indir+"correlation_tos_erainterim_son-lon.nc","r")
        k2 = addfile(indir+"correlation_uas_erainterim_son-lon.nc","r")
        k3 = addfile(indir+"correlation_vas_erainterim_son-lon.nc","r")
        
;; -- file containing significance of cors differences ------ 
           file3 = addfile(indir+"correlation_diff_tos_nud1_son-lon.nc","r")
        file3uas = addfile(indir+"correlation_diff_uas_nud1_son-lon.nc","r")
        file3vas = addfile(indir+"correlation_diff_vas_nud1_son-lon.nc","r")

		tos3_dif = file3->tos(0,:,:)
		uas3_dif = file3uas->uas(0,:,:)
		vas3_dif = file3vas->vas(0,:,:)

		;; -- Lon/lat ---
        lon = c1->lon
        lat = c1->lat
		

		tos3_ctr = c1->tos(0,:,:)
        uas3_ctr = c2->uas(0,:,:)
        vas3_ctr = c3->vas(0,:,:)

		tos3_win = g1->tos(0,:,:)
        uas3_win = g2->uas(0,:,:)
        vas3_win = g3->vas(0,:,:)
		;; -- obs ---
		tos3_obs = k1->tos(0,:,:)
        uas3_obs = k2->uas(0,:,:)
        vas3_obs = k3->vas(0,:,:)
;;create plots;;
        wks =\
        gsn_open_wks("eps","Figure4-third-row") ; send graphics to PNG file
	    plot    = new(4,graphic) 
        gsn_define_colormap(wks,"hotcolr_19lev") ; choose colormap
        
		pres                         = True
        pres@gsnFrame                = False    ; Turn off so we can add titles
		res                      = True
		res@gsnDraw              = False
        res@gsnFrame             = False
		
		res@cnFillOn             = True               ; turn on color for contours
        res@cnLinesOn            = False              ; turn off contour lines
        res@cnLineLabelsOn       = False              ; turn off contour line labels
        res@gsnScalarContour     = True               ; contours desired
      
        res@mpLandFillColor      = "grey"            ; set land to be gray

        res@gsnStringFont              = "helvetica-bold"  ; will be used for gsnCenterString later
        res@gsnCenterStringFontHeightF = 0.04              ; ditto 
		res@gsnCenterString            = "CTR: SON"
		res@gsnLeftString              = ""
		res@gsnRightString             = ""

      
        res@gsnAddCyclic   = False
      
        res@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
        res@cnMinLevelValF       = -1.               ; set min contour level
        res@cnMaxLevelValF       =  1.                 ; set max contour level
        res@cnLevelSpacingF      =  0.2               ; set contour spacing
      
;; ----- Change aspect ratio of plot ----------------------
		res@mpShapeMode  = "FreeAspect"
        res@vpHeightF    = 0.4         
        res@vpWidthF     = 0.9                  

;; ------ Label resources ---------------------------------		
    res@lbLabelBarOn               =False

; ---------- Set resources for when we turn on labelbar later
    res@lbTitleFontHeightF         = .02                 ; make title smaller
  ;  res@lbTitleString              = "~S~o~N~C/~S~o~N~C"            ; title string
 ;   res@lbTitlePosition            = "Bottom"            ; title position
    res@lbTitleDirection           = "Across"            ; title direction
    res@lbLabelStride              = 2
    res@lbLabelFontHeightF         = 0.03
    res@pmLabelBarHeightF          = 0.2
    res@lbOrientation              = "Horizontal"        ; orientation
    res@lbBoxLinesOn               = False     ; Turn off labelbar box lines 

; ---- Tickmarks & Labels ----------------------
    res@tmYLMinorOn          = False
    res@tmXBMinorOn          = False
	res@tmXTLabelsOn         = False
	res@tmXBLabelsOn         = False
	res@tmYLLabelsOn         = False


;; ----- Map selection -------------------------------------
        res@mpMinLatF            =  -30. 
        res@mpMaxLatF            =   30.
        res@mpMinLonF            =  -270. 
        res@mpMaxLonF            =   30.
        res@mpCenterLonF          = 210
 ;; ------ Vector resources --------------------------------  
		res@vcRefLengthF              = 0.03           ; define length of vec ref
        res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
        res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
        res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref
        res@vcRefMagnitudeF           = 1.

		res@vcMinDistanceF            = 0.017
        ;res@vcGlyphStyle              = "CurlyVector" ;  "LineArrow"     ; turn on curly vectors
        res@vcGlyphStyle              = "LineArrow"     ; turn on curly vectors
        res@vcLineArrowColor          = "black"           ; change vector color
        res@vcLineArrowThicknessF     =  0.8               ; change vector thickness
        res@vcVectorDrawOrder         = "PostDraw"        ; draw vectors last


;;  ----------- Start plotting ------------------------
		
        res@cnMinLevelValF       = -1.               ; set min contour level
        res@cnMaxLevelValF       =  1.                 ; set max contour level
        res@cnLevelSpacingF      =  .2               ; set contour spacing
		
		res@lbLabelBarOn               = True
        
		res@tmXBLabelsOn         = True
		res@tmXBLabelFontHeightF = 0.02
	    res@tmXBMode        = "Explicit"	
        res@tmXBValues      = (/-270, -240, -210,  -180, -150, -120, -90, -60, -30, 0, 30  /)
        res@tmXBLabels      = (/"" , "240W" ,"210W" ,"180W" ,"150W" ,"120W" ,"90W" ,"60W" ,"30W" ,"0", "" /)
      

       txres = True 
       txres@txFontHeightF = 0.012
         txres@txFont        = "helvetica-bold"
       gsn_text_ndc(wks,"Correlations SST SON/ATL3 SON", 0.5 , 0.6,txres)
		res@gsnCenterString            = "a) Reanalysis"
        
        plot(0) = gsn_csm_vector_scalar_map_ce(wks,uas3_obs,vas3_obs,tos3_obs,res)        
         pres@gsnPanelLeft = .0 
         pres@gsnPanelRight = 0.24 
         gsn_panel(wks,plot(0),(/1,1/),pres) 
	
        res@gsnCenterString            = "b) CTR"
        plot(1) = gsn_csm_vector_scalar_map_ce(wks,uas3_ctr,vas3_ctr,tos3_ctr,res)        
         pres@gsnPanelLeft = .25 
         pres@gsnPanelRight = 0.49 
         gsn_panel(wks,plot(1),(/1,1/),pres) 

		res@gsnCenterString            = "c) NUD-JJAS"
        plot(2) = gsn_csm_vector_scalar_map_ce(wks,uas3_win,vas3_win,tos3_win,res)        
         pres@gsnPanelLeft = .5
         pres@gsnPanelRight = 0.74 
         gsn_panel(wks,plot(2),(/1,1/),pres) 

		res@gsnCenterString            = "d) NUD-JJAS - CTR"
        res@cnMinLevelValF       = -0.5               ; set min contour level
        res@cnMaxLevelValF       =  0.5                 ; set max contour level
        res@cnLevelSpacingF      =  0.1               ; set contour spacing
        ;plot(11) = gsn_csm_contour_map_ce(wks,tos3_dif,res)        
        res@vcRefMagnitudeF           = .4
        gsn_define_colormap(wks,"CBR_drywet") ; choose colormap
		gsn_reverse_colormap(wks)
        plot(3) = gsn_csm_vector_scalar_map_ce(wks,uas3_dif, vas3_dif, tos3_dif, res)        
         pres@gsnPanelLeft = .75 
         pres@gsnPanelRight = 0.99 
         gsn_panel(wks,plot(3),(/1,1/),pres) 

        frame(wks)   ; now advance the frame!
    end
