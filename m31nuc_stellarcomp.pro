pro m31nuc_stellarcomp

  loadct, 38
  ; nu pav data
  swspath='/Users/peeters/data/Obs/ISO/SWS/Sources/nuPav/'
  swsfile='com_stand.fits'
  a = rfsws(swspath+swsfile)
  ind=where(a.data.line eq 6)
  a.data[ind].flux= a.data[ind].flux* .97d0
  a=select(a, a.data.line ne 5 or $
     (a.data.line eq 5 and a.data.wave lt 5.3))
  a=select(a, a.data.line ne 6 or $
     (a.data.line eq 6 and a.data.wave gt 5.3 and a.data.wave lt 7.03))
  a=select(a, a.data.line ne 7 or $
     (a.data.line eq 7 and a.data.wave gt 7.03 and a.data.wave lt 12.11))
  ind=where(a.data.line eq 9)
  a.data[ind].flux= a.data[ind].flux* .95d0
 a=select(a, a.data.line ne 9 or $
     (a.data.line eq 9 and a.data.wave gt 12.11 and a.data.wave lt 16.13))
 ind=where(a.data.line eq 10)
  a.data[ind].flux= a.data[ind].flux* .88d0
 a=select(a, a.data.line ne 10 or $
     (a.data.line eq 10 and a.data.wave gt 16.13 and a.data.wave lt 19.5))
 ind=where(a.data.line eq 11)
  a.data[ind].flux= a.data[ind].flux* .85d0
 a=select(a, a.data.line ne 11 or $
     (a.data.line eq 11 and a.data.wave gt 19.5))
nupav_mir = a 
  
  nupav_nir =[2010d0 , 3150d0 , 2680d0] ; Jy
  wave_nir=[1.235, 1.662, 2.159]
    ;;; 2mass
   ;;; j=-.213, h=-1.197, k=-1.526 --> saturated according to Greg.
 ; new from Greg (email July 17, 2015: his 2015 paper on Giants (Catalog of Infrared Observations (CIO, version 5.1; Gezari et al. 2000).)
;  J = -0.250 +/- 0.028
; H = -1.220 +/- 0.071
; K = -1.510 +/- 0.014

  pl, nupav_mir, col=200, yr=[0,3200]
  oplot, wave_nir, nupav_nir, col=200

  ;;; quiescent galaxy template Kaneda 2008
  pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'
  galfile = 'IRS-egal-bgd_Kaneda_June20_2015.txt'
  read_spec, pad+galfile, gal, startline=7, columns=3
  ;; remove overlapping regions between different orders
  ; find where the wavelength decreases
  ;idx_jump = where((gal[0,*]-shift(gal[0,*],-1)) gt 0)
  ; find for each jump, how many overlapping points there are
  ;idx_count = dblarr(n_elements(idx_jump))
  ;for loop=0, n_elements(idx_jump)-2 do begin
  ;   posi=where(gal[0,*] eq gal[0,idx_jump[loop]])
  ;   if n_elements(posi) eq 2 then idx_count[loop] = posi[1]-idx_jump[loop]
  ;end 
  idx_nooverlap = [indgen(73), indgen(13)+76, indgen(158)+102, indgen(6)+268, indgen(103)+278]
  gal = gal[*, idx_nooverlap]
  ;; aperture SL 25"x3.6"; LL 71"x10.2"
  gal_Jy = reform(convert_flux(gal[0,*], gal[1,*] , out='Jy/""', in='MJy/sr')*25d0*3.6d0)
  plot, gal[0,*], gal[1,*]

; NGC0337
pad2 = '/Users/peeters/Dropbox/Dimuthu/nucleus/'
read_spec, pad2+ 'ngc0337_PAHFit_full.tbl', ngc, columns=3, startline=3
wave=reform(ngc[0,*])
flux=reform(ngc[1,*])
unc=reform(ngc[2,*])

; pahfit
restore, filename=pad2+'pahfit_ngc0337.xdr'
cyfit=pahfit_components(wave,fit_ngc0337,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
ccont=dust_tot+stars
ref = flux - ccont

;;; PART I: LARGE APERTURE NUCLEUS
; spectrum nucleus -- large aperture
; check how much % of the flux at 5.4 is due to stellar
; contribution. Use this fraction to normalize the smaller aperture
; spectrum
; ~100% is stellar at 5.4 (upto 8.x um!)

read_spec, pad+ 'nucFLUX',nuc, columns=2, startline=1
nuc_Jy = nuc[1,*] * 1d6 * (30d0*50d0/4.25d10)

; From Pauline's email on June 17, 2015:
;# M31 nucleus photometry - for the LARGE aperture!!!@
;# same 30" x 50" aperture as IRS spectral extraction
;# band  flux density (Jy)
nir = [5.501371568858024474e+00, 8.375078554756264282e+00 ,6.868690818058364123e+00] ; J, H, K
; normalize nupav and nucleus on H (Vega et al. 2010)

tops, 'm31nuc_largeap_stellarcomp.ps'
plot, nuc[0,*], nuc_Jy / nir[1], xr=[0,22], yr=[-.05,1.05], ys=1
oplot, wave_nir, nir/nir[1],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/nupav_nir[1]), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/nupav_nir[1], col=200, ps=6
oplot, nuc[0,*], nuc_Jy / nir[1]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, nuc[0,*])
oplot, nuc[0,*], (nuc_Jy / nir[1]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
ver, 5.4, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus large ap., normalized on H [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed, normalized on H',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
dummy =  (nuc_Jy / nir[1]) - flux_interpol
wave= nuc[0,*]
ind= where(wave gt 5.39 and wave lt 5.4)
print, 'stellar contribution at 5.4um (H) ', (dummy[ind]/ (nuc_Jy / nir[1]))*100, '%'

plot, nuc[0,*], nuc_Jy / nir[0], xr=[0,22], yr=[-.05,1.65], ys=1
oplot, wave_nir, nir/nir[0],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/nupav_nir[0]), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/nupav_nir[0], col=200, ps=6
oplot, nuc[0,*], nuc_Jy / nir[0]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, nuc[0,*])
oplot, nuc[0,*], (nuc_Jy / nir[0]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus large ap., normalized on J [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed, normalized on J',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
dummy =  (nuc_Jy / nir[0]) - flux_interpol
wave= nuc[0,*]
ind= where(wave gt 5.39 and wave lt 5.4)
print, 'stellar contribution at 5.4um (J) ', (dummy[ind]/ (nuc_Jy / nir[0]))*100, '%'

normf_nupav = nupav_mir.data[where(nupav_mir.data.wave gt 5.394 and nupav_mir.data.wave lt 5.399)].flux
normf_nuc = nuc_Jy[where(nuc[0,*] gt 5.39 and nuc[0,*] lt 5.4)]

plot, nuc[0,*], reform(nuc_Jy) / normf_nuc[0], xr=[0,22];, yr=[-.05,1.65], ys=1
oplot, wave_nir, nir/normf_nuc[0],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/normf_nupav), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/normf_nupav[0], col=200, ps=6
oplot, nuc[0,*], nuc_Jy /normf_nuc[0]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, nuc[0,*])
oplot, nuc[0,*], (nuc_Jy / normf_nuc[0]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus large ap., normalized at 5.4 [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
ver, 5.4, linestyle=1, col=30

tops, 'm31nuc_largeap_stellarcomp.ps',/finish



;;; PART II: SMALL APERTURE NUCLEUS
; spectrum nucleus -- small aperture

read_spec, pad+ 'm31nuc_ll2_nucCentre.tbl', ll2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl2_nucCentre.tbl', sl2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl1_nucCentre.tbl', sl1, columns=3, startline=16

ind1 = where(sl2[0,*] lt 7.597)
ind2 = where(sl1[0,*] gt 7.597 and sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537)
wave1 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux1 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc1 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]
; flux in MJy/sr (9"x9" aperture)
; 1sr = 4.25d10 "x"
; to Jy
flux1_Jy = flux1 * 1d6 * (9d0^2/4.25d10)
; test: flux3=flux2
;flux3 = convert_flux(wave1, flux1 , out='Jy/""', in='MJy/sr')*9d0^2
nir1 = [0.617470257083d0, 0.951390649613d0, 0.791547282057d0] ;

;spectrum north position
read_spec, 'm31nuc_ll2_nucUp_correct.tbl', ll2, columns=3, startline=239
read_spec, 'm31nuc_sl2_nucUp_correct.tbl', sl2, columns=3, startline=311
read_spec, 'm31nuc_sl1_nucUp_correct.tbl', sl1, columns=3, startline=311

ind1 = where(sl2[0,*] lt 7.53);7.548)
ind2 = where(sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537 and ll2[0,*] lt 20.978)
wave2 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux2 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc2 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]
; flux in MJy/sr (9"x9" aperture)
; 1sr = 4.25d10 "x"
; to Jy
flux2_Jy = flux2 * 1d6 * (9d0^2/4.25d10)
nir2=[0.276751385896d0, 0.423687757329d0, 0.349367215727d0]

stop
tops, 'm31nuc_smallap_stellarcomp54.ps'
plot, wave1, flux1_Jy / nir1[1], xr=[0,22], yr=[-.05,1.05], ys=1
oplot, wave_nir, nir1/nir1[1],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/nupav_nir[1]), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/nupav_nir[1], col=200, ps=6
oplot, wave1, flux1_Jy / nir1[1]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, wave1)
oplot, wave1, (flux1_Jy / nir1[1]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
ver, 5.4, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus centr ap., normalized on H [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed, normalized on H',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
dummy =  (flux1_Jy / nir1[1]) - flux_interpol
ind= where(wave1 gt 5.39 and wave1 lt 5.4)
xyouts, .5, .75, 'non-stellar contribution at 5.4um (H) '+string((dummy[ind]/ (flux1_Jy / nir1[1]))*100)+ ' %',/normal,  charsize=1.1

plot, wave2, flux2_Jy / nir2[1], xr=[0,22], yr=[-.05,1.05], ys=1
oplot, wave_nir, nir2/nir2[1],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/nupav_nir[1]), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/nupav_nir[1], col=200, ps=6
oplot, wave2, flux2_Jy / nir2[1]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, wave2)
oplot, wave2, (flux2_Jy / nir2[1]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
ver, 5.4, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus north ap., normalized on H [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed, normalized on H',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
dummy =  (flux2_Jy / nir2[1]) - flux_interpol
ind= where(wave2 gt 5.39 and wave2 lt 5.4)
xyouts, .5, .75, 'non-stellar contribution at 5.4um (H) '+string((dummy[ind]/ (flux2_Jy / nir2[1]))*100)+ ' %',/normal,  charsize=1.1

plot, wave2, flux2_Jy / nir2[0], xr=[0,22], yr=[-.05,1.65], ys=1
oplot, wave_nir, nir2/nir2[0],ps=5
norm_nupav =  sh_calcaar(sh_calcaar(nupav_mir, factor=1d0/nupav_nir[0]), smoo=10)
pl, norm_nupav, /o, col=200 
oplot, wave_nir, nupav_nir/nupav_nir[0], col=200, ps=6
oplot, wave2, flux2_Jy / nir2[0]
flux_interpol = interpol(norm_nupav.data.flux, norm_nupav.data.wave, wave2)
oplot, wave2, (flux2_Jy / nir2[0]) - flux_interpol, col=100
hor, 0, linestyle=1, col=30
ver, 5.4, linestyle=1, col=30
xyouts, .5, .9, 'M31 nucleus north ap., normalized on J [Jy]',/normal, charsize=1.1
xyouts, .5, .85, 'nu Pav, SWS, smoothed, normalized on J',/normal, col=200, charsize=1.1
xyouts, .5, .8, 'nucleus- nu pav',/normal, col=100, charsize=1.1
dummy =  (flux2_Jy / nir2[0]) - flux_interpol
ind= where(wave2 gt 5.39 and wave2 lt 5.4)
xyouts, .5, .75, 'non-stellar contribution at 5.4um (J) '+string((dummy[ind]/ (flux2_Jy / nir2[0]))*100)+ ' %',/normal,  charsize=1.1

tops, 'm31nuc_smallap_stellarcomp54.ps',/finish


tops, 'm31nuc_smallaps_stellarcomp.ps'
;; compare nupav & galaxy normalized at 5.39
  norm_nupav = nupav_mir.data[where(nupav_mir.data.wave gt 5.394 and nupav_mir.data.wave lt 5.399)].flux
  pl, sh_calcaar(nupav_mir, factor=1d0/norm_nupav), xr=[5,25],title=' '
  norm_gal = gal_Jy[where(gal[0,*] gt 5.39 and gal[0,*] lt 5.4)]
  oplot, gal[0,*], gal_Jy/norm_gal[0], col=200
  pl, sh_calcaar(sh_calcaar(nupav_mir, fac=1d0/norm_nupav), smoo=10), col=100,/o, thick=2.
  oplot, gal[0,*], smooth(gal_Jy/norm_gal[0], 4), col=250, thick=2.
  xyouts, .5, .9, 'nu Pav, SWS, normalized at 5.4 [Jy]',/normal, charsize=1.1
  xyouts, .5, .85, 'nu Pav, smoothed',/normal, col=100, charsize=1.1
  xyouts, .5, .8, 'quiescent gal. ',/normal, col=200, charsize=1.1
  xyouts, .5, .75, 'quiescent gal., smoothed',/normal, col=250, charsize=1.1
 
 ;; comparenuc & north 
  plot, nuc[0,*], reform(nuc[1,*]) , xr=[5,25], yr=[10, 75], ys=1
  oplot,  wave1, flux1, col=100
  oplot,  wave2, flux2, col=250
  xyouts, .5, .8, 'nucleus (large ap)',/normal, charsize=1.1
  xyouts, .5, .75, 'nucleus center (small ap)',/normal, col=100, charsize=1.1
  xyouts, .5, .7, 'nucleus north (small ap)',/normal, col=250, charsize=1.1

  ;; compare nupav & galaxy & nuc & north normalized at 5.39
  pl, sh_calcaar(nupav_mir, factor=1d0/norm_nupav), xr=[5,25],title=' '
  oplot, gal[0,*], gal_Jy/norm_gal[0], col=200
  oplot, nuc[0,*], reform(nuc_Jy) / normf_nuc[0], col=50
  normf_flux1=flux1_Jy[where(wave1 gt 5.38 and wave1 lt 5.4)]
  normf_flux2=flux2_Jy[where(wave2 gt 5.38 and wave2 lt 5.4)]
  oplot,  wave1, flux1_Jy/ normf_flux1[0], col=100
  oplot,  wave2, flux2_Jy/ normf_flux2[0], col=250
  xyouts, .5, .9, 'nu Pav, SWS, normalized at 5.4 [Jy]',/normal, charsize=1.1
  xyouts, .5, .85, 'quiescent gal. ',/normal, col=200, charsize=1.1
  xyouts, .5, .8, 'nucleus (large ap)',/normal, col=50, charsize=1.1
  xyouts, .5, .75, 'nucleus center (small ap)',/normal, col=100, charsize=1.1
  xyouts, .5, .7, 'nucleus north (small ap)',/normal, col=250, charsize=1.1


  ;; compare nupav & galaxy & nuc & north normalized at 5.39
  pl, sh_calcaar(nupav_mir, factor=1d0/norm_nupav), xr=[5,25],title=' ',/nodata
  oplot,  wave1, flux1_Jy/ normf_flux1[0], col=100
  oplot,  wave2, flux2_Jy/ normf_flux2[0], col=250
  xyouts, .5, .9, 'nucleus center (small ap)',/normal, col=100, charsize=1.1
  xyouts, .5, .85, 'nucleus north (small ap)',/normal, col=250, charsize=1.1


  ; subtract nuc small from others.  
  plot, nuc[0,*],( reform(nuc_Jy) / normf_nuc[0]) - (flux1_Jy/ normf_flux1[0]),  xr=[5,20], yr=[-0.05, .3], ys=1
  oplot,  wave2,( flux2_Jy/ normf_flux2[0])- (flux1_Jy/ normf_flux1[0]), col=250
  xyouts, .5, .9, 'nucleus (large ap) - nucleus center (small ap)',/normal,charsize=1.1
  xyouts, .5, .85, 'nucleus north (small ap)- nucleus center (small ap)',/normal, col=250, charsize=1.1
   xyouts, .5, .8, 'normalized at 5.4 [Jy]',/normal,  charsize=1.1
   hor, 0, linestyle=1, col=30
   dummy =( flux2_Jy/ normf_flux2[0])- (flux1_Jy/ normf_flux1[0])
   ind=where(wave gt 11. and wave lt 11.6)
   max2= max(ref[ind])/ max(dummy[ind])
   oplot, wave, ref/max2, col=kleur('blue')
xyouts, .5, .75, 'NGC0337', col=kleur('blue'),/normal,  charsize=1.1

  
  ;; subtract nupav and galaxy from nucleus spectrum
  ;; nupav = Jy, galaxy = MJy/sr, nucleus: flux1=MJy/sr, flux1_Jy=Jy
  ;; normalize at 5.4um
  ; put on same wavelength grid
  flux_nupav = interpol(nupav_mir.data.flux, nupav_mir.data.wave, wave1)
  flux_gal = interpol(gal[1,*], gal[0,*], wave1)
  ; normalizing factor
  dummy = where(wave1 gt 5.38 and wave1 lt 5.4)
  factor_nupav= flux1_Jy[dummy]/flux_nupav[dummy]
  factor_gal= flux1[dummy]/flux_gal[dummy]
  ; check
  plot, wave1, flux1_Jy, xr=[0,22], yr=[-0.02,.2], ys=1,title='M31 nucleus (small ap) vs. nu Pav (orange) normalized at 5.4um', xtitle='Wavelength', ytitle='Flux density [Jy]'
  oplot, wave1, flux_nupav*factor_nupav[0], col=200
  oplot, wave1, flux1_Jy- flux_nupav*factor_nupav[0], col=100
  hor, 0
stop
  plot, wave1, flux1, xr=[0,22], yr=[-5,90], ys=1, title='M31 nucleus (small ap) vs. quiescent gal (orange) normalized at 5.4um', xtitle='Wavelength', ytitle='Surface Brightness [MJy/sr]'
  oplot, wave1, flux_gal*factor_gal[0], col=200
  oplot, wave1, flux1- flux_gal*factor_gal[0], col=100
  hor, 0


  ;;; PART III: SMALL APERTURE NORTH POSITION NUCLEUS.
;spectrum north position
read_spec, 'm31nuc_ll2_nucUp_correct.tbl', ll2, columns=3, startline=239
read_spec, 'm31nuc_sl2_nucUp_correct.tbl', sl2, columns=3, startline=311
read_spec, 'm31nuc_sl1_nucUp_correct.tbl', sl1, columns=3, startline=311

ind1 = where(sl2[0,*] lt 7.53);7.548)
ind2 = where(sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537 and ll2[0,*] lt 20.978)
wave2 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux2 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc2 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]
; flux in MJy/sr (9"x9" aperture)
; 1sr = 4.25d10 "x"
; to Jy
flux2_Jy = flux2 * 1d6 * (9d0^2/4.25d10)

;; subtract nupav and galaxy from nucleus spectrum
  ;; nupav = Jy, galaxy = MJy/sr, nucleus: flux1=MJy/sr, flux1_Jy=Jy
  ;; normalize at 5.4um
  ; put on same wavelength grid
  flux_nupav = interpol(nupav_mir.data.flux, nupav_mir.data.wave, wave2)
  flux_gal = interpol(gal[1,*], gal[0,*], wave2)
  ; normalizing factor
  dummy = where(wave2 gt 5.38 and wave2 lt 5.4)
  factor1_nupav= flux2_Jy[dummy]/flux_nupav[dummy]
  factor1_gal= flux2[dummy]/flux_gal[dummy]
  dummy = where(wave2 gt 5.28 and wave2 lt 5.31)
  factor2_nupav= flux2_Jy[dummy]/flux_nupav[dummy]
  factor2_gal= flux2[dummy]/flux_gal[dummy]
  ; check
  plot, wave2, flux2_Jy, xr=[0,22], yr=[-0.03,.08], ys=1,title='M31 nucleus-North vs. nu Pav (color) normalized at 5.4um (orange) or 5.3um (green)', xtitle='Wavelength', ytitle='Flux density [Jy]'
  oplot, wave2, flux_nupav*factor1_nupav[0], col=200
  oplot, wave2, flux_nupav*factor2_nupav[0], col=100
  oplot, wave2, flux2_Jy- flux_nupav*factor1_nupav[0]-0.02, col=200
  oplot, wave2, flux2_Jy- flux_nupav*factor2_nupav[0]-0.02, col=100
  hor, -0.02

  plot, wave2, flux2, xr=[0,22], yr=[-5,40], ys=1, title='M31 nucleus-North vs. quiescent gal (color) normalized at 5.4um (orange) or 5.3um (green)', xtitle='Wavelength', ytitle='Surface Brightness [MJy/sr]'
  oplot, wave2, flux_gal*factor1_gal[0], col=200
  oplot, wave2, flux2- flux_gal*factor1_gal[0], col=200
  oplot, wave2, flux_gal*factor2_gal[0], col=100
  oplot, wave2, flux2- flux_gal*factor2_gal[0], col=100
  hor, 0

; compare north & center nucleus

  ; put on same wavelength grid
  flux1_wave2 = interpol(flux1, wave1, wave2)
  ; normalizing factor
  dummy = where(wave2 gt 5.38 and wave2 lt 5.4)
  factor1= flux2[dummy]/flux1_wave2[dummy]
  dummy = where(wave2 gt 5.28 and wave2 lt 5.31)
  factor2= flux2[dummy]/flux1_wave2[dummy]
  ; check
  plot, wave2, flux2, xr=[0,22], yr=[-5,40], ys=1,title='M31 nucleus-North vs. centre (color) normalized at 5.4um (orange) or 5.3um (green)', xtitle='Wavelength', ytitle='Flux density [Jy]'
  oplot, wave2, flux1_wave2*factor1[0], col=200
  oplot, wave2, flux1_wave2*factor2[0], col=100
  oplot, wave2, flux2- flux1_wave2*factor1[0]-0.02, col=200
  oplot, wave2, flux2- flux1_wave2*factor2[0]-0.02, col=100
  hor, -0.02

tops, 'm31nuc_smallaps_stellarcomp.ps',/finish
  
stop

end
