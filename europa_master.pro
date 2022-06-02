FUNCTION match_scattered_sunlight, p, x=x, y=y, err=err, fit=fit
  common sunlight_fit_common, fitindices

  fit = P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]

  return, abs(y - fit)/err
end

FUNCTION scale_fit_sunlight, p, x
  return, P[0]*GAUSS_SMOOTH(x,P[2],/EDGE_TRUNCATE) + P[1]
end

pro europa_master, part=part

  ; Written C. Schmidt BU CSP June 5 2018
  ; Modified by P. Lierle BU CSP Nov 30 2021
  ; Modified by E. Lovett BU CSP May 27 2021
  
    user = 'Carl'                   ; Case sensitive
    case 1 of
     user eq 'Carl': Directory = 'D:\Data\Perkins\Perkins RIPS - April 2018\'
     user eq 'Emma': Directory = 'C:\Users\elovett\EuropaResearch\More Europa Data\Perkins RIPS - April 2018\'
    endcase 
  files = Directory + 'rips_europa\Processed'

  imaging_statsec = '[183 : 809, 39 : 444]'  ; coordinates where light falls the detector (length and format matter here)
  spectra_statsec = '[ 0 : 1023,480 : 959]'  ; coordinates where light falls the detector (length and format matter here)
  i_coords        = long(strsplit(repstr( STRMID(imaging_statsec, 1, 19), ':', ',' ), ',', /extract)) ; make coordinates into an array of 4
  s_coords        = long(strsplit(repstr( STRMID(spectra_statsec, 1, 19), ':', ',' ), ',', /extract)) ; make coordinates into an array of 4

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
        statsec = imaging_statsec, readn=RDnoise, gain=Gain ; best parameters are tested, slow n good.
      Na_flatBSCR_imaging = MRDFITS(directory+'Processed\Na_Flat.BSCR.fits', 0, header, /Dscale, /silent )
      la_cosmic, directory+'Processed\Na_Flat.BS.fits', outsuff = "CR", sigclip = 4.5, $
        statsec = spectra_statsec, readn=RDnoise, gain=Gain ; best parameters are tested, slow n good.
      Na_flatBSCR_spectra = MRDFITS(directory+'Processed\Na_Flat.BSCR.fits', 0, header, /Dscale, /silent )

    ; Now replace each section in each channel with the CR subtracted image
      Na_flatBS[i_coords[0]:i_coords[1], i_coords[2]:i_coords[3]] = Na_flatBSCR_spectra[i_coords[0]:i_coords[1], i_coords[2]:i_coords[3]]
      Na_flatBS[s_coords[0]:s_coords[1], s_coords[2]:s_coords[3]] = Na_flatBSCR_spectra[s_coords[0]:s_coords[1], s_coords[2]:s_coords[3]]
      
    ; Normalize each flat field to unity
      Na_flat_imaging = Na_flatBS / mean(Na_flatBS[i_coords[0]:i_coords[1], i_coords[2]:i_coords[3]])
      Na_flat_spectra = Na_flatBS / mean(Na_flatBS[s_coords[0]:s_coords[1], s_coords[2]:s_coords[3]])

    ; Process the spectral frames
      window, 0, xs = 1024, ys = 1024
            filenames = ['Moon_Na_Detection_MinSlit-SubsolarLimb-scan_New.fits', 'Moon_Na_Detection_MinSlit-SubsolarLimb-scan2_New.fits']
      ;filenames = file_search(directory+'rips_europa\'+'Europa*')
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
      
      ;sky = MRDFITS(directory+'Sky_58.fits', 0, header, /Dscale, /silent ) - bias
     
      ; Hack, Hack, Hack, Hack, Hack, Hack, Hack, Hack!!!!!!!
      ; Data are not flat fielded, unclear yet, if flat fielding helps!!!
                ;yshift       = -1.5
                ;sky          = MRDFITS(directory+'Sky_58.fits', 0, header, /Dscale, /silent ) - bias
                flat         = Na_flat_spectra[*, 480:*]
                s            = size(flat, /dim)    
                shifted_flat = INTERPOLATE(flat, findgen(s[0]), findgen(s[1]) + yshift, Missing = -1, /grid)
        shifted_flat = 1.
      stop
      BSCR[*, 480:* ] = BSCR[*, 480:* ] / shifted_flat

      ; Hack, Hack, Hack, Hack, Hack, Hack, Hack, Hack!!!!!!!
      ; not dealing with the imaging channel flat field yet
      ;BSCR[*, 0:479 ] = BSCR[*, 0:479 ] / Na_flat_imaging[*, 0:479 ]
      MWRFITS, BSCR, directory+'rips_europa\Processed\' + new_filename + '.FF.fits', header, /silent, /create
    endfor

  endif

  if part eq 1 then begin                 ; straightening out the image; includes shifting in the y-direction and rectification using the moon drift slit file

; FIRST, get sunlight spectrum all on same row

    filter = "*BSCR.fits"
    BSCRfiles = FILE_SEARCH(files, filter, count=n_Files)
    newimg = fltarr(1024,544)
    pixelsofsun = 6 ;12

    FOR file = 0, n_Files - 1 DO BEGIN
      array = MRDFITS(BSCRfiles[file], 0, header, /Fscale, /silent)
      array = float(array)
      specchan = array[*,480:*]                                           ; Eliminates imaging channel
      sunlightspec = []
      straightened_spectrum = specchan                                    ; Uses spectral channel as a dummy array, we'll write the straightened spectrum into this array
      sunmax = []

      FOR i = 0, 1024 - 1 DO BEGIN
        column = specchan[i,*]
        suncol = WHERE(column EQ MAX(column),count)                       ; Finds the sunlight spectrum in each column and puts that row into a 1D array
        suncol = suncol[0]
        sunmax = [sunmax, suncol[0]]
      ENDFOR

      newpos = MEAN(sunmax)                                               ; Shifting the spectrum to match up with the mean sunlight spectrum location
      polyn = poly_fit(findgen(1024), sunmax,2)                           ; Finds 2nd order polynomial of the sunlight spectrum
      polynomial = poly(findgen(1024), polyn)                             ; Puts polynomial into a form polyn[0] + polyn[1]*X + polyn[2]*X^2
      shifted = polynomial - newpos                                       ; finds the difference between actual position and desired position

      FOR i = 0, 1024 - 1 DO BEGIN
        column = specchan[i,*]
        shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
          + shifted[i], Missing = -1)
        straightened_spectrum[i,*] = shifted_array                        ; This is the straightened image in the y-direction
      ENDFOR

; SECOND, use the moon's drift slit exposure to get the polynomial for the curvature & straighten:
      moonfile = Directory + 'Moon_Surface_DriiftScan.fits'
      moon     = float(MRDFITS(moonfile, 0, header, /Fscale, /silent))
      moonspecchan = moon[*,480:*]
      
      ; straighten the moonspectrum
        straightened_moonspectrum = specchan
          FOR i = 0, 1024 - 1 DO BEGIN
            column = moonspecchan[i,*]
            shifted_array = INTERPOLATE(column, findgen(n_elements(column)) $ ; Shifts each column so that it is positioned at the last suncol
              + shifted[i], Missing = -1)
            straightened_moonspectrum[i,*] = shifted_array                    ; This is the straightened image in the y-direction
          ENDFOR
      
      D1cen = 390.
      D2cen = 631.
      windowwidth = 15.
      D1 = straightened_moonspectrum[D1cen - windowwidth:D1cen + windowwidth,*]
      D2 = straightened_moonspectrum[D2cen - windowwidth:D2cen + windowwidth,*]
      D1_center = fltarr(n_elements(straightened_moonspectrum[0,*]))
      D2_center = fltarr(n_elements(straightened_moonspectrum[0,*]))

      ; D1 LINE FIRST
        for i = 0, n_elements(straightened_moonspectrum[0,*]) - 1 do begin
          result = mpfitpeak(findgen(windowwidth*2. + 1), straightened_moonspectrum[D1cen - windowwidth:D1cen + windowwidth,i], a, STATUS = STATUS)
          if STATUS GT 1 then D1_center[i] = A[1] else D1_center[i] = !values.F_nan
        endfor
  
        y = findgen(n_elements(straightened_moonspectrum[0,*]))
        real = where(finite(D1_center), /NULL)
        COEFF = ROBUST_POLY_FIT(y[real], D1_center[real], 2)
        D1location = poly(findgen(n_elements(straightened_moonspectrum[0,*])), coeff) + D1cen-windowwidth
        
      ; D2 LINE SECOND
        for i = 0, n_elements(straightened_moonspectrum[0,*]) - 1 do begin
          result = mpfitpeak(findgen(windowwidth*2. + 1), straightened_moonspectrum[D2cen - windowwidth:D2cen + windowwidth,i], a, STATUS = STATUS)
          if STATUS GT 1 then D2_center[i] = A[1] else D2_center[i] = !values.F_nan
        endfor
  
        y = findgen(n_elements(straightened_moonspectrum[0,*]))
        real = where(finite(D2_center), /NULL)
        COEFF = ROBUST_POLY_FIT(y[real], D2_center[real], 2)
        D2location = poly(findgen(n_elements(straightened_moonspectrum[0,*])), coeff) + D2cen-windowwidth
      
      ; Take the average deviation from the mean location, then averaged between both D1 & D2 lines
        mean_shift = MEAN([[D1location - mean(D1location)], $          ; Seems good enough, more inspection would be a good idea though!
                           [D2location - mean(D2location)]], dim=2) 
      
;   Apply this rectification to all Europa files
    rectified_moon   = specchan ; dummy array, about to overwrite
    rectified_europa = specchan ; dummy array, about to overwrite
    FOR i = 0, 544 - 1 DO BEGIN
      test_row = straightened_moonspectrum[*,i]
      test_rectified_array = INTERPOLATE(test_row, findgen(n_elements(test_row)) + mean_shift[i], Missing = -1)
      rectified_moon[*,i] = test_rectified_array
      
      row = straightened_spectrum[*,i]
      rectified_array = INTERPOLATE(row, findgen(n_elements(row)) + mean_shift[i], Missing = -1)
      rectified_europa[*,i] = rectified_array
    ENDFOR   
    ;      Inspect the rectified moon
    ;      window, 0, xs = 1024, ys=1088
    ;      tv, bytscl(rebin(rebin(rectified_moon, 1024, 68), 1024, 1088))

    ;----------------------------------Sunlight Subtraction--------------------------------------
    sunlight    = rectified_europa[*, newpos - pixelsofsun : newpos + pixelsofsun]
    sunlight_1d = TOTAL(sunlight,2)
    fitindices  = [indgen(393), indgen(245)+400, indgen(1024-645)+645]     ; Excludes the sodium emission from Europa
    P_returned  = fltarr(4,544)                                            ; Three coefficients + MPFIT's "Status"
    P_guessed   = Fltarr(3,544)                                            ; Initial Guess that we throw at MPFIT

    FOR i = 0, 544-1 DO BEGIN
      ; generate intial guess for multipliciative scaling
        row         = rectified_europa[*,i]
        guess_scale = median(row[fitindices] / sunlight_1D[fitindices])
        row_err     = sqrt(abs(row))

      ; Fit a y = A*Gauss_smooth(x,C) + B function to the spectrum, where x is the reference solar spectrum
        p0 = [guess_scale, 0.0, 0.0]                                           ; Guess at initial coefficients
        parinfo = replicate({value:0., fixed:0, limited:[0,0], limits:[0.,0.]}, n_elements(p0))
        parinfo.value         = p0
        parinfo[1].fixed      = 1
        parinfo[2].fixed      = 1
        parinfo[2].limited    = [1, 1]
        parinfo[2].limits     = [0.0, 20.]

      fa = {x:sunlight_1d[fitindices], y:row[fitindices], err:row_err[fitindices]}
      p = mpfit('match_scattered_sunlight', p0, PERROR = err_P, functargs=fa, status=status, parinfo=parinfo, /quiet)
      P_guessed[*,i]  = p0
      p_returned[*,i] = [p, status]
      scaled_sunlight = scale_fit_sunlight(P, sunlight_1d)             ; Puts it into y = A*shift(Gauss_smooth(x,D),C) + B form

      sub = guess_scale * sunlight_1d                                  ; If you JUST want the multiplicative correction (no scattered sunlight accounted for w mpfit)
      totsubtrd  = row  - scaled_sunlight                              ; Change back to row - sub to get just the multiplicative factor
      newimg[*,i] = totsubtrd

    ENDFOR
    
; This next part is CHECKING the rectification by plotting the actual line centroids (D1/D2_center).
; Window 1 shows the rectified lines (should be straight, and second plot should be ~ 0)
; Window 2 shows 2d rectified rebinned image
; Window 3 shows the bias/dark/flat/CR corrected, shifted, rectified, sunlight-subtracted 1D light from Europa 

;  window, 1
;    !P.Multi = [0, 1, 2]
;    cgplot, findgen(n_elements(newimg[0,*])), D1_center + D1cen-windowwidth, YRANGE=[350.,700.], XRANGE=[0,544], title='Europa D1 and D2 Line Curvature', $
;      xtitle='Row (Pixel)', ytitle='Centroid of Absorption Line'
;    cgplot, findgen(n_elements(newimg[0,*])), D2_center + D2cen-windowwidth, color='black', /overplot
;    cgplot, findgen(n_elements(newimg[0,*])), D2_center - D1_center, /ynozero
  
  window, 2, xs=1024, ys = 68
    !P.Multi = 0
    TV, bytscl(rebin(newimg, 1024, 68), -5.,10.)
    
  window, 3, xs = 1500, ys = 600
    pos = cglayout([2,2], xgap = 0, ygap = 0, iymargin = 0)
    pos[Where(pos eq 0.52790695)] = .7
    cgplot, findgen(1024), TOTAL(rectified_europa,2), title="Raw RIPS Spectrum", ytitle="Total Counts", /xstyle, pos = pos[*,0], xtickformat = '(A1)'
    cgplot, findgen(1024), TOTAL(newimg,2), title="Scaled Subtracted Sunlit Spectrum", xtitle = "Column (pixels)", ytitle="Total Counts", /xstyle, pos = pos[*,2], /noerase
    cgplot, TOTAL(rectified_europa, 1), findgen(N_elements(TOTAL(rectified_europa, 1))), /xstyle, pos = pos[*,1],/noerase, ytickformat = '(A1)'
    cgplot, [-1.e6, 1.e6], [newpos, newpos], linestyle = 2, color = 'red', /overplot
    cgplot, TOTAL(newimg, 1), findgen(N_elements(TOTAL(newimg, 1))), /xstyle, pos = pos[*,3],/noerase, ytickformat = '(A1)'
    cgplot, [-1.e6, 1.e6], [newpos, newpos], linestyle = 2, color = 'red', /overplot
stop
  ENDFOR
  endif
end