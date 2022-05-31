FUNCTION match_scattered_sunlight, p, x=x, y=y, err=err, fit=fit
  common sunlight_fit_common, fitindices

  fit = P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]

  return, (y - fit)/err
end

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

pro europa_master, part=part

  ;Written C. Schmidt BU CSP June 5 2018
  ;Modified by P. Lierle BU CSP Nov 30 2021
  ;Modified by E. Lovett BU CSP May 27 2021
  user = 'Emma'
  case 1 of
   user eq 'carl': Directory = 'C:\Users\schmidtc\Perkins RIPS - April 2018'
   user eq 'Emma': begin
     Directory_mac   = '/DATA/Perkins/Perkins RIPS - April 2018/'
     Directory_PC    = 'C:\Users\elovett\EuropaResearch\More Europa Data\Perkins RIPS - April 2018'
     Directory       = Directory_PC
     stop
   end  
  endcase 


  imaging_statsec = '[183 : 809, 39 : 444]'
  spectra_statsec = '[ 99 : 990,502 : 978]'

  if part eq 0 then begin         ; standard image corrections/calibrations

    ; Make the Bias
    biases = FILE_SEARCH(Directory + 'Bias*', count = n_bias)
    bigarray = fltarr(n_bias,1024,1024)
    for i = 0, n_bias - 1 do begin
      bigarray[i,*,*] = MRDFITS(biases[i], 0, header, /Dscale, /silent )
    endfor
    bias = mean(bigarray, dimension = 1)
    rdnoise = mean(stddev(bigarray, dimension = 1))

    ; Make the Flat
    Na_Flats = FILE_SEARCH(Directory + 'Dome_flat_Na'+strcompress(indgen(3),/rem)+'.fits', count = n_Flats)
    bigarray = fltarr(n_Flats,1024,1024)
    for i = 0, n_Flats - 1 do begin
      bigarray[i,*,*] = MRDFITS(Na_Flats[i], 0, header, /Dscale, /silent )
    endfor
    Na_flatBS = median(bigarray, dimension = 1) - Bias

    ; Still some hot pixels in it...fix them
    MWRFITS, Na_flatBS, directory+'Processed\Na_Flat.BS.fits', header, /create, /silent
    la_cosmic, directory+'Processed\Na_Flat.BS.fits', outsuff = "CR", sigclip = 4.5, $
      statsec = imaging_statsec, readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.
    Na_flatBSCR_imaging = MRDFITS(directory+'Processed\Na_Flat.BSCR.fits', 0, header, /Dscale, /silent )
    la_cosmic, directory+'Processed\Na_Flat.BS.fits', outsuff = "CR", sigclip = 4.5, $
      statsec = spectra_statsec, readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.
    Na_flatBSCR_spectra = MRDFITS(directory+'Processed\Na_Flat.BSCR.fits', 0, header, /Dscale, /silent )
    ;now replace each sections in each channel with the CR subtracted image
    Na_flatBS[*, 480:* ] = Na_flatBSCR_spectra[*, 480:* ]
    Na_flatBS[*, 0:479 ] = Na_flatBSCR_imaging[*, 0:479 ]
    ;Normalize each
    Na_flat_spectra = Na_flatBS / mean(Na_flatBS[99 :990,502 : 978])
    Na_flat_imaging = Na_flatBS / mean(Na_flatBS[183: 809, 39: 444])


    ; Process the spectral frames
    window, 0, xs = 1024, ys = 1024
    ;      filenames = ['Moon_Na_Detection_MinSlit-SubsolarLimb-scan_New.fits', 'Moon_Na_Detection_MinSlit-SubsolarLimb-scan2_New.fits']
    filenames = file_search(directory+'rips_europa\'+'Europa*')
    summed_spectrum = fltarr(1024)
    for i = 0, n_elements(filenames)-1 do begin
      filename = strcompress(filenames[i])
      BS = MRDFITS(filename, 0, header, /Dscale, /silent ) - bias
      Gain = 1. ;float(sxpar(header, 'GAIN'))
      Print, filenames[I], (sxpar(header, 'EXPOSURE'))
      new_filename = STRMID(filename, strlen(directory+'rips_europa\'))
      new_filename = STRMID(new_filename, 0, strpos(new_filename,'.fits'))
      MWRFITS, BS, directory+'rips_europa\Processed\' + new_filename + '.BS.fits', header, /create;, /silent
      la_cosmic, directory+'rips_europa\Processed\' + new_filename + '.BS.fits', outsuff = "CR", sigclip = 4.5, $
        statsec = statsec, readn=RDnoise, gain=Gain ;best parameters are tested, slow n good.
      BSCR = MRDFITS(directory+'rips_europa\Processed\' + new_filename + '.BSCR.fits', 0, header, /Dscale, /silent )
      Dummy = BSCR
      ;BSCR[*, 480:* ] = BSCR[*, 480:* ] / smart_shift(Na_flat_spectra[*, 480:* ], 0, -1.5)
      BSCR[*, 480:* ] = BSCR[*, 480:* ] / smart_shift(Na_flat_spectra[*, 480:* ], 0, -1.5)
      ;BSCR[*, 0:479 ] = BSCR[*, 0:479 ] / Na_flat_imaging[*, 0:479 ]
      MWRFITS, BSCR, directory+'rips_europa\Processed\' + new_filename + '.FF.fits', header, /silent, /create
    endfor

  endif

  if part eq 1 then begin                 ; straightening out the image; includes shifting in the y-direction and rectification using the moon drift slit file

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

      newpos = MEAN(suncol)                                                     ; Shifting the spectrum to match up with the right-most sunlight spectrum location
      polyn = poly_fit(findgen(1024), sunmax,2)                           ; Finds 2nd order polynomial of the sunlight spectrum
      polynomial = poly(findgen(1024), polyn)                             ; Puts polynomial into a form polyn[0] + polyn[1]*X + polyn[2]*X^2
      shifted = polynomial - newpos                                       ; finds the difference between actual position and desired position

      FOR i = 0, 1024 - 1 DO BEGIN
        column = specchan[i,*]
        shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
          + shifted[i], Missing = -1)
        shifted_sunlight[i,*] = shifted_array                             ; This is the straightened image in the y-direction
        lowlim = newpos - pixelsofsun                                     ; Finds the brightest pixels above Europa's disk are and goes some height pix below that
        upplim = newpos + pixelsofsun                                     ; Finds the brightest pixels above Europa's disk are and goes some height pix above that
        sunlit = shifted_array[lowlim:upplim]                             ; Sunlight portion of each column
        sunlight[i,*] = sunlit
      ENDFOR

; SECOND, use the moon's drift slit exposure to get the polynomial for the curvature & straighten:
      moonfile = 'C:\Users\elovett\EuropaResearch\More Europa Data\Perkins RIPS - April 2018\Moon_Surface_DriiftScan.fits'
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

      ; D2 LINE SECOND
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

;   Applies this rectification to all Europa files

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

    ENDFOR
; This next part is CHECKING the rectification by plotting the actual line centroids (D1/D2_center).
; Window 1 shows the rectified lines (should be straight, and second plot should be ~ 0)
; Window 2 shows 2d rectified rebinned image
; Window 3 shows the bias/dark/flat/CR corrected, shifted, rectified, sunlight-subtracted 1D light from Europa 

  window, 1
  !P.Multi = [0, 1, 2]
  cgplot, findgen(n_elements(newimg[0,*])), D1_center + D1cen-windowwidth, YRANGE=[350.,700.], XRANGE=[0,544], title='Europa D1 and D2 Line Curvature', $
    xtitle='Row (Pixel)', ytitle='Centroid of Absorption Line'
  cgplot, findgen(n_elements(newimg[0,*])), D2_center + D2cen-windowwidth, color='black', /overplot

  cgplot, findgen(n_elements(newimg[0,*])), D2_center - D1_center, /ynozero
  
  window, 2
  !P.Multi = 0
  TV, bytscl(rebin(newimg, 1024, 68), -5.,10.)
    
  window, 3
  CgDisplay, 1500,600
  !P.Multi = [0,1,2]
  cgplot, findgen(1024), TOTAL(specchan,2), title="Raw RIPS Spectrum",  xtitle = "Column (pixels)", ytitle="Total Counts"
  cgplot, findgen(1024), TOTAL(newimg,2)  , title="Scaled Subtracted Sunlit Spectrum", xtitle = "Column (pixels)", ytitle="Total Counts";, yrange=[1000.,4000.],xrange=[0.,1024.]
  if file EQ 3 then stop

  ENDFOR
  endif
end