pro fig_m31nuc_norm54

; nu pav
  
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

; quiescent galaxy
 ;;; quiescent galaxy template Kaneda 2008
  pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'
  galfile = 'IRS-egal-bgd_Kaneda_June20_2015.txt'
  read_spec, pad+galfile, gal, startline=7, columns=3
  ;; remove overlapping regions between different orders
   idx_nooverlap = [indgen(73), indgen(13)+76, indgen(158)+102, indgen(6)+268, indgen(103)+278]
  gal = gal[*, idx_nooverlap]
  ;; aperture SL 25"x3.6"; LL 71"x10.2"
  gal_Jy = reform(convert_flux(gal[0,*], gal[1,*] , out='Jy/""', in='MJy/sr')*25d0*3.6d0)
  plot, gal[0,*], gal[1,*]

  ; spectrum nucleus -- large aperture
read_spec, pad+ 'nucFLUX',nuc, columns=2, startline=1
nuc_Jy = nuc[1,*] * 1d6 * (30d0*50d0/4.25d10)

  
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


  norm_nupav = nupav_mir.data[where(nupav_mir.data.wave gt 5.394 and nupav_mir.data.wave lt 5.399)].flux
  norm_gal = gal_Jy[where(gal[0,*] gt 5.39 and gal[0,*] lt 5.4)]
  normf_flux1=flux1_Jy[where(wave1 gt 5.38 and wave1 lt 5.4)]
  normf_flux2=flux2_Jy[where(wave2 gt 5.38 and wave2 lt 5.4)]
  normf_nuc = nuc_Jy[where(nuc[0,*] gt 5.39 and nuc[0,*] lt 5.4)]

  
; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_m31nuc_norm54_log.eps', /portrait,xsize=15,ysize=10,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=1.4
!y.charsize=1.4
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38

pl, sh_calcaar(nupav_mir, factor=1d0/norm_nupav), xr=[5,25],title=' ', ytitle='normalized Flux density [Jy]', position=[.13,.15,.98,.98], /xlog 
  oplot, gal[0,*], gal_Jy/norm_gal[0], col=200
  oplot, nuc[0,*], reform(nuc_Jy) / normf_nuc[0], col=100
   oplot,  wave1, flux1_Jy/ normf_flux1[0], col=50
  oplot,  wave2, flux2_Jy/ normf_flux2[0], col=250
  xyouts, .55, .88, 'NU Pav',/normal, charsize=1.1
  xyouts, .55, .83, 'quiescent galaxy template',/normal, col=200, charsize=1.1
  xyouts, .55, .78, 'M31 nuclear',/normal, col=50, charsize=1.1
  xyouts, .55, .73, 'M31 full',/normal, col=100, charsize=1.1
  xyouts, .55, .68, 'M31 north',/normal, col=250, charsize=1.1
  xyouts, .84, .075, '20',/normal, charsize=1.5

device,/close
set_plot,'x'

end
