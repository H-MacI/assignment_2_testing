pro fig_pahfit_m31_nucleus

pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'

;spectrum nucleus
read_spec, pad+ 'm31nuc_ll2_nucCentre.tbl', ll2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl2_nucCentre.tbl', sl2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl1_nucCentre.tbl', sl1, columns=3, startline=16

ind1 = where(sl2[0,*] lt 7.597)
ind2 = where(sl1[0,*] gt 7.597 and sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537)
wave1 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux1 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc1 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]

;spectrum north position

read_spec, 'm31nuc_ll2_nucUp_correct.tbl', ll2, columns=3, startline=239
read_spec, 'm31nuc_sl2_nucUp_correct.tbl', sl2, columns=3, startline=311
read_spec, 'm31nuc_sl1_nucUp_correct.tbl', sl1, columns=3, startline=311

ind1 = where(sl2[0,*] lt 7.53)
ind2 = where(sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537 and ll2[0,*] lt 20.978)
wave2 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux2 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc2 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]

; get pahfit results
restore, pad+'pahfit_nucUp_correct.xdr' ;fit_nucUp
restore, pad+'pahfit_nucCentre.xdr' ;fit_nucCentre

;; not usefull .... (PAHFIT to nucleus)
plot, wave1, flux1, xr=[5, 22], xs=1, ys=1, yr=[10, 90]
oplot, wave1, fit_nucCentre.final_fit, col=kleur('red')
yfit1=pahfit_components(wave1,fit_nucCentre,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave1,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave1,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave1, cont, COLOR=kleur('magenta'),LINESTYLE=0

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_sp_m31_nucleus.ps', /portrait,xsize=15,ysize=20,/color, bits_per_pixel=8

cleanplot
!x.charsize=1.6
!y.charsize=1.6
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38

; plot
plot, wave2, flux2, xr=[5, 22], xs=1, ys=1, yr=[5, 90], position=[.15,.1,.97,.97],$
   xtitle=textoidl('Wavelength (\mum)'), $
   ytitle=textoidl('Intensity I_{\nu} (MJy/sr)')
oplot, wave2, fit_nucUp.final_fit, col=kleur('red')
yfit2=pahfit_components(wave2,fit_nucUp,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave2,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
  xyouts, 19, 13, 'North', charsize=1.5

oplot, wave1, flux1
ind=where(wave1 gt 7.5d0 and wave1 lt 15.3d0)
;oplot, wave1[ind], fit_nucCentre.final_fit[ind], col=kleur('red')
yfit1=pahfit_components(wave1,fit_nucCentre,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
 cont=dust_tot+stars
 oplot, wave2[ind], cont[ind], COLOR=kleur('magenta'),LINESTYLE=0
 oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
 ;rat=(cont+lines[*,7])
 ;oplot,wave1[ind],rat[ind],COLOR=kleur('blue'),LINESTYLE=0
 xyouts, 19, 25, 'Centre', charsize=1.5

device,/close
set_plot,'x'

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_sp_m31_nucleus_vs2.eps', /portrait,xsize=15,ysize=15,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=1.6
!y.charsize=1.6
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38

; plot
plot, wave2, flux2, xr=[5, 22], xs=1, ys=1, yr=[5, 45], position=[.15,.12,.97,.97],$
   xtitle=textoidl('Wavelength (\mum)'), $
   ytitle=textoidl('Intensity I_{\nu} (MJy/sr)')
oplot, wave2, fit_nucUp.final_fit, col=kleur('red')
yfit2=pahfit_components(wave2,fit_nucUp,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave2,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
  xyouts, 19, 31, 'North', charsize=1.5

device,/close
set_plot,'x'


;; comparison spectrum: NGC0337 HII-type galaxy

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

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_sp_m31_nucleus_vs3.eps', /portrait,xsize=15,ysize=15,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=1.4
!y.charsize=1.4
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38
!p.multi=[0,1,2]

x0=.15
y0=.11
x1=.98
y1=.98
ddy = 0.0
dy = (y1-y0 - 1.*ddy)/2.
xt = 16.2

; plot
plot, wave2, flux2, xr=[5, 22], xs=1, ys=1, yr=[5, 45], position=[x0, dy+y0+ddy, x1, y1], xtickname=[' ', ' ', ' ', ' ',' ', ' ', ' ', ' ']
oplot, wave2, fit_nucUp.final_fit, col=kleur('red')
yfit2=pahfit_components(wave2,fit_nucUp,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave2,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
  xyouts, 17, 40, 'North', charsize=1.5

ind=where(wave2 gt 11. and wave2 lt 11.6)
max1= max(flux2[ind]-cont[ind])
ind=where(wave gt 11. and wave lt 11.6)
max2= max(ref[ind])

plot, wave2, flux2, xr=[5,22], xs=1, ys=1, yr=[-.1, 1.1], position=[x0, y0, x1, y0+dy],/nodata
oplot, wave2, (flux2-cont)/max1
oplot, wave, ref/max2, col=kleur('blue')
xyouts, 17, .95, 'North', charsize=1.5
xyouts, 17, 0.8, 'NGC0337', charsize=1.5, col=kleur('blue')

xyouts, x0+(x1-x0)/2., 0.02, textoidl('Wavelength (\mum)'),/normal, alignment=.5, charsize=1.4
xyouts, 0.04, y0 +dy+ddy +dy/2., textoidl('Intensity I_{\nu} (MJy/sr)') ,/normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.04, y0 + dy/2., textoidl('Normalized Intensity'), /normal, orientation=90, alignment=.5, charsize=1.4


device,/close
set_plot,'x'

cleanplot
; use nuclear center as stellar spectrum and subtract this from norht
; spectrum
stellarcomp =  smooth(flux1, 4)
normf_stellarcomp=stellarcomp[where(wave1 gt 5.38 and wave1 lt 5.4)]
normf_flux2=flux2[where(wave2 gt 5.38 and wave2 lt 5.4)]

flux2_nostellar = flux2 - stellarcomp/normf_stellarcomp[0] * normf_flux2[0]

void=pahfit(wave2,PARINFO=para2,/NO_FIT)

fit_nucUp_nostellar=pahfit(wave2, flux2_nostellar, unc2, parinfo=para2, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+ 'pahfit_nucUp_nostellar_report_correct.txt')

save, fit_nucUp_nostellar, filename=pad+'pahfit_nucUp_nostellar_correct.xdr'

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_sp_m31_nucleus_vs4.eps', /portrait,xsize=15,ysize=15,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=1.4
!y.charsize=1.4
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38
!p.multi=[0,1,2]

x0=.15
y0=.11
x1=.98
y1=.98
ddy = 0.0
dy = (y1-y0 - 1.*ddy)/2.
xt = 16.2

; plot
plot, wave2, flux2_nostellar, xr=[5, 22], xs=1, ys=1, yr=[-1, 12], position=[x0, dy+y0+ddy, x1, y1], xtickname=[' ', ' ', ' ', ' ',' ', ' ', ' ', ' ']
oplot, wave2, fit_nucUp_nostellar.final_fit, col=kleur('red')
yfit2=pahfit_components(wave2,fit_nucUp_nostellar, $
                         DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave2,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
  xyouts, 17, 10, 'North', charsize=1.5
  xyouts, 17, 8.5, 'stellar-subtracted', charsize=1.

ind=where(wave2 gt 11. and wave2 lt 11.6)
max1= max(flux2_nostellar[ind]-cont[ind])
ind=where(wave gt 11. and wave lt 11.6)
max2= max(ref[ind])

plot, wave2, flux2_nostellar, xr=[5,22], xs=1, ys=1, yr=[-.1, 1.1], position=[x0, y0, x1, y0+dy],/nodata
oplot, wave2, (flux2_nostellar-cont)/max1
oplot, wave, ref/max2, col=kleur('blue')
xyouts, 17, .95, 'North', charsize=1.5
xyouts, 17, 0.8, 'NGC0337', charsize=1.5, col=kleur('blue')

xyouts, x0+(x1-x0)/2., 0.02, textoidl('Wavelength (\mum)'),/normal, alignment=.5, charsize=1.4
xyouts, 0.04, y0 +dy+ddy +dy/2., textoidl('Intensity I_{\nu} (MJy/sr)') ,/normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.04, y0 + dy/2., textoidl('Normalized Intensity'), /normal, orientation=90, alignment=.5, charsize=1.4


device,/close
set_plot,'x'

; normalize on H instead of 5.4 ...

nir2=[0.276751385896d0, 0.423687757329d0, 0.349367215727d0]
nir1 = [0.617470257083d0, 0.951390649613d0, 0.791547282057d0] ;

stellarcomp =  smooth(flux1, 4)

flux2_nostellar_nir = flux2 - stellarcomp/nir1[1] * nir2[1]

void=pahfit(wave2,PARINFO=para2,/NO_FIT)

fit_nucUp_nostellar_nir=pahfit(wave2, flux2_nostellar_nir, unc2, parinfo=para2, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+ 'pahfit_nucUp_nostellar_report_correct.txt')

save, fit_nucUp_nostellar_nir, filename=pad+'pahfit_nucUp_nostellar_nir_correct.xdr'

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_sp_m31_nucleus_vs5.eps', /portrait,xsize=15,ysize=15,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=1.4
!y.charsize=1.4
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38
!p.multi=[0,1,2]

x0=.15
y0=.11
x1=.98
y1=.98
ddy = 0.0
dy = (y1-y0 - 1.*ddy)/2.
xt = 16.2

; plot
plot, wave2, flux2_nostellar_nir, xr=[5, 22], xs=1, ys=1, yr=[-1, 12], position=[x0, dy+y0+ddy, x1, y1], xtickname=[' ', ' ', ' ', ' ',' ', ' ', ' ', ' ']
oplot, wave2, fit_nucUp_nostellar_nir.final_fit, col=kleur('red')
yfit2=pahfit_components(wave2,fit_nucUp_nostellar_nir, $
                         DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
   cont=dust_tot+stars
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   
  for i=0,(size(lines,/DIMENSIONS))[1]-1 do begin 
     rat=(cont+lines[*,i])
     oplot,wave2,rat,COLOR=kleur('green'),LINESTYLE=0
  endfor 
   oplot, wave2, cont, COLOR=kleur('magenta'),LINESTYLE=0
  xyouts, 17, 10, 'North', charsize=1.5
  xyouts, 17, 8.5, 'stellar-subtracted', charsize=1.

ind=where(wave2 gt 11. and wave2 lt 11.6)
max1= max(flux2_nostellar_nir[ind]-cont[ind])
ind=where(wave gt 11. and wave lt 11.6)
max2= max(ref[ind])

plot, wave2, flux2_nostellar_nir, xr=[5,22], xs=1, ys=1, yr=[-.1, 1.1], position=[x0, y0, x1, y0+dy],/nodata
oplot, wave2, (flux2_nostellar_nir-cont)/max1
oplot, wave, ref/max2, col=kleur('blue')
xyouts, 17, .95, 'North', charsize=1.5
xyouts, 17, 0.8, 'NGC0337', charsize=1.5, col=kleur('blue')

xyouts, x0+(x1-x0)/2., 0.02, textoidl('Wavelength (\mum)'),/normal, alignment=.5, charsize=1.4
xyouts, 0.04, y0 +dy+ddy +dy/2., textoidl('Intensity I_{\nu} (MJy/sr)') ,/normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.04, y0 + dy/2., textoidl('Normalized Intensity'), /normal, orientation=90, alignment=.5, charsize=1.4


device,/close
set_plot,'x'

stop
end
