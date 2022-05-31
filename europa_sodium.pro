FUNCTION match_scattered_sunlight, p, x=x, y=y, err=err, fit=fit
  common sunlight_fit_common, fitindices

  ; Multiply P[0], Add P[1] and shift P[2] a reference spectrum (X) until it best matches Y
  ;  stop
  ;shift = p[2]
  ;shifted_x = interpolate(x, findgen(N_elements(x)) - p[2]); , MISSING=!values.F_nan
  ;stop
  fit = P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]

  return, (y - fit)/err
end

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

pro europa_sodium

  ; FIRST, get sunlight spectrum all on same row
  
  files_drive = "D:\More Europa Data\Perkins RIPS - April 2018\rips_europa\Processed"
  files_desktop = "C:\Users\elovett\EuropaResearch\More Europa Data\Perkins RIPS - April 2018\rips_europa\Processed"
  files = files_desktop
  filter = "*BSCR.fits"
  BSCRfiles = FILE_SEARCH(files, filter, count=n_Files)
  newimg = fltarr(1024,544)
  pixelsofsun = 12

  FOR file = 0, n_Files - 1 DO BEGIN
    array = MRDFITS(BSCRfiles[file], 0, header, /Fscale, /silent)
    array = float(array)
    specchan = array[*,480:*]                                           ; Eliminates imaging channel
    sunlightspec = []
    shifted_sunlight = specchan                                         ; Uses spectral array as a dummy array
    sunlight = fltarr(1024, 2*pixelsofsun + 1)                          ; Creates an array that will hold JUST the sunlight spectrum
    sunmax = []

    FOR i = 0, 1024 - 1 DO BEGIN
      column = specchan[i,*]
      suncol = WHERE(column EQ MAX(column),count)                       ; Finds the sunlight spectrum in each column and puts that row into a 1D array
      suncol = suncol[0]
      sunmax = [sunmax, suncol[0]]
    ENDFOR

    newpos = suncol                                                     ; Shifting the spectrum to match up with the right-most sunlight spectrum location
    polyn = poly_fit(findgen(1024), sunmax,2)                           ; Finds 2nd order polynomial of the sunlight spectrum
    polynomial = poly(findgen(1024), polyn)                             ; Puts polynomial into a form polyn[0] + polyn[1]*X + polyn[2]*X^2
    shifted = polynomial - newpos                                       ; finds the difference between actual position and desired position

    FOR i = 0, 1024 - 1 DO BEGIN
      column = specchan[i,*]
      shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
        + shifted[i], Missing = -1)
      shifted_sunlight[i,*] = shifted_array
      lowlim = newpos - pixelsofsun                                     ; Finds the brightest pixels above Europa's disk are and goes some height pix below that
      upplim = newpos + pixelsofsun                                     ; Finds the brightest pixels above Europa's disk are and goes some height pix above that
      sunlit = shifted_array[lowlim:upplim]                             ; Sunlight portion of each column
      sunlight[i,*] = sunlit
    ENDFOR

    ; SECOND, use the moon's drift slit exposure to get the polynomial for the curvature & straighten:
    moonfile = 'C:\Users\elovett\EuropaResearch\More Europa Data\Perkins RIPS - April 2018\Moon_Surface_DriiftScan.fits'
    ;moonfile = 'D:\More Europa Data\Perkins RIPS - April 2018\Moon_Surface_DriiftScan.fits'         ; this is for the hard drive
    array = MRDFITS(moonfile, 0, header, /Fscale, /silent)
    array = float(array)
    moonspecchan = array[*,480:*]
  
    D1cen = 390.
    D2cen = 631.
    windowwidth = 15.
    D1 = moonspecchan[D1cen - windowwidth:D1cen + windowwidth,*]
    D2 = moonspecchan[D2cen - windowwidth:D2cen + windowwidth,*]
  
    D1_center = fltarr(n_elements(moonspecchan[0,*]))
    D2_center = fltarr(n_elements(moonspecchan[0,*]))
    
    ; D1 LINE FIRST
    for i = 0, n_elements(moonspecchan[0,*]) - 1 do begin
      result = mpfitpeak(findgen(windowwidth*2. + 1), moonspecchan[D1cen - windowwidth:D1cen + windowwidth,i], a, STATUS = STATUS)
      if STATUS GT 1 then D1_center[i] = A[1] else D1_center[i] = !values.F_nan
    endfor
    
    y = findgen(n_elements(moonspecchan[0,*]))
    real = where(finite(D1_center), /NULL)
    COEFF = ROBUST_POLY_FIT(y[real], D1_center[real], 2)
    D1location = poly(findgen(n_elements(moonspecchan[0,*])), coeff) + D1cen-windowwidth
    
    D1polyn = poly_fit(findgen(544), D1location,2)                          ; Finds 2nd order polynomial of the sunlight spectrum
    D1polynomial = poly(findgen(544), D1polyn)                              ; Puts polynomial into a form polyn[0] + polyn[1]*X + polyn[2]*X^2
    
    ;D2 LINE SECOND
    for i = 0, n_elements(moonspecchan[0,*]) - 1 do begin
      result = mpfitpeak(findgen(windowwidth*2. + 1), moonspecchan[D2cen - windowwidth:D2cen + windowwidth,i], a, STATUS = STATUS)
      if STATUS GT 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
    endfor
  
    y = findgen(n_elements(moonspecchan[0,*]))
    real = where(finite(D2_center), /NULL)
    COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
    D2location = poly(findgen(n_elements(moonspecchan[0,*])), coeff) + D2cen-windowwidth
  
    shiftdiff = MEAN(D2location)
    
    rectified_europa = specchan
    
    FOR i = 0, 544 - 1 DO BEGIN
      row = shifted_sunlight[*,i]
      shifted = D2location[i] - shiftdiff
      rectified_array = INTERPOLATE(row, findgen(n_elements(row)) + shifted, Missing = -1)
      rectified_europa[*,i] = rectified_array
    ENDFOR
    
    sunlight_1d = TOTAL(sunlight,2)
    fitindices = [indgen(393), indgen(245)+400, indgen(1024-645)+645]     ; Excludes the sodium emission from Europa
    P_returned = fltarr(4,544)

    FOR i = 0, 544-1 DO BEGIN
      row = rectified_europa[*,i]
      midd = row[400:645, *]
      scaledtot  = TOTAL(row)  / TOTAL(sunlight)
      scaledmid  = TOTAL(midd) / TOTAL(sunlight[indgen(225)+405, *])      ; Gets multiplicative factor
      row_err    = sqrt(abs(row))

      ; Fit a y = A*Gauss_smooth(x,C) + B function to the spectrum, where x is the reference solar spectrum

      p0 = [scaledmid, 0.0, 0.0]                                           ; Guess at initial coefficients
      ;p0 = [-6.97679e-06, 0.429791, 0.0]                                ; Initial guesses for when status = 1.0
      parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p0))
      parinfo.value         = p0
      parinfo[1].fixed      = 0
      parinfo[2].fixed      = 1
      parinfo[2].limited    = [1, 1]
      parinfo[2].limits     = [0.0, 20.]

      fa = {x:sunlight_1d[fitindices], y:row[fitindices], err:row_err[fitindices]}
      p = mpfit('match_scattered_sunlight', [parinfo[0].value, parinfo[1].value, parinfo[2].value], PERROR = err_P, funct=fa, status=status, parinfo=parinfo)
      p_returned[*,i] = [p, status]
      telluric_fit_BFR = scale_fit_sunlight(P, sunlight_1d)             ; Puts it into y = A*shift(Gauss_smooth(x,D),C) + B form

      sub = scaledmid * sunlight_1d                                     ; If you JUST want the multiplicative correction (no scattered sunlight accounted for w mpfit)
      totsubtrd  = row  - telluric_fit_BFR                              ; Change back to row - sub to get just the multiplicative factor
      newimg[*,i] = totsubtrd
      ;      CgDisplay, 1530,1000
      ;      !P.Multi = [0,1,2]
      ;      cgplot, findgen(1024), row;TOTAL(specchan,2)
      ;      cgplot, findgen(1024), telluric_fit_BFR, color='red', /overplot
      ;      if i eq suncol then stop

    ENDFOR
  ; This next part is CHECKING the rectification by plotting the actual line centroids (D1/D2_center). Second plot in window 1 is the difference between the rectified lines (best case = 0)
  ; Window 2 shows 2d rectified rebinned image

  window, 1
  !P.Multi = [0, 1, 2]
  cgplot, findgen(n_elements(newimg[0,*])), D1_center + D1cen-windowwidth, YRANGE=[350.,700.], XRANGE=[0,544], title='Europa D1 and D2 Line Curvature', $
    xtitle='Row (Pixel)', ytitle='Centroid of Absorption Line'
  cgplot, findgen(n_elements(newimg[0,*])), D2_center + D2cen-windowwidth, color='black', /overplot

  cgplot, findgen(n_elements(newimg[0,*])), D2_center - D1_center, /ynozero
  
  !P.Multi = 0
    window, 2
    TV, bytscl(rebin(newimg, 1024, 68), -5.,10.)
    
    
    
    
    CgDisplay, 1500,600
    !P.Multi = [0,1,2]
    cgplot, findgen(1024), TOTAL(specchan,2), title="Raw RIPS Spectrum",  xtitle = "Column (pixels)", ytitle="Total Counts"
    cgplot, findgen(1024), TOTAL(newimg,2)  , title="Scaled Subtracted Sunlit Spectrum", xtitle = "Column (pixels)", ytitle="Total Counts";, yrange=[1000.,4000.],xrange=[0.,1024.]
;    if file EQ 3 then stop
    stop
    
    
  ENDFOR
  
end