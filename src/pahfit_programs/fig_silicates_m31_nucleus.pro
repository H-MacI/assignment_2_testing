function fit_func,x, p

f = p[2]*mygauss(x, peak=p[0], fwhm=p[1],/maxnorm) 
Return, f

end

pro fig_silicates_m31_nucleus

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

; get pahfit results
restore, pad+'pahfit_nucCentre.xdr' ;fit_nucCentre

oplot, wave1, flux1
ind=where(wave1 gt 7.5d0 and wave1 lt 15.3d0)
;oplot, wave1[ind], fit_nucCentre.final_fit[ind], col=kleur('red')
yfit1=pahfit_components(wave1,fit_nucCentre,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
 cont=dust_tot+stars

; get other silicate profiles

; Galactic Center (abs)
read_spec, '/Users/peeters/data/Analyse/Extinction/Silicate_profile/silprofile_chiar.dat', galcen, columns=3, startline=15
; subtract continuum extinction to get silicate profile. (continuum
; extinction is fixed beyond 8um - Chiar & Tielens 2006)
ind=where(galcen[0,*] eq 7.940d0) 
galcenprof = reform(galcen[1,*] - galcen[1,ind[0]])

; ULIRG in Smith et al. 2010 (digitized) - (abs).
read_spec, pad+'iras08572_smith2010.dat', ulirg, columns=2, startline=3

; M81 in Smith et al. 2010 (obtained from Howard Smith) - (em.)
read_spec, pad+'M81best_AARSmZmergedSL1SL2LL.dat', m81, columns=3, startline=1
ind= UNIQ(m81[0,*], SORT(m81[0,*]))
m81w= m81[0,ind]
m81f=m81[1,ind]
m81x=m81[2,ind]
; this is crap -- so take the digitized spectrum in the mean time.
read_spec, pad+'m81_smith2010.dat', m81, columns=2, startline=3
; 2nd try
dummy=sap_rfits(pad+'M81_Howard_test.fits', m81struct)
m81=dblarr(2, n_elements(m81struct.data.wave))
m81[0,*]=m81struct.data.wave
m81[1,*]=m81struct.data.flux

ctptw=[8.36, 8.83, 24.7, 26.45, 27.79]
ctptf=[0.192, 0.195, 0.280, 0.291, 0.311]

plot, m81[0,*], m81[1,*]
cm81= spline(ctptw, ctptf, m81[0,*])
oplot, m81[0,*], cm81 , col=222


; normalization factors
ind=where(wave1 gt 9 and wave1 lt 12)
maxsp = max(smooth(reform(flux1[ind]-cont[ind]),5))
ind=where(galcen[0,*] gt 9 and galcen[0,*] lt 12)
maxgalcen = max(galcenprof[ind])
ind=where(m81[0,*] gt 9 and m81[0,*] lt 11)
maxm81=max(m81[1,ind]-cm81[ind])


plot, wave1, (flux1-cont)/maxsp, xr=[7,21], xs=1
oplot, wave1, smooth(reform(flux1-cont),5)/maxsp, col=kleur('gold')
hor, 0, linestyle=3
oplot, galcen[0,*], galcenprof/maxgalcen, col=kleur('blue')
;oplot, galcen[0,*]+.35, galcenprof/maxgalcen, col=kleur('red')
; exclude ULIRG - similar to gal.cen.
;oplot, ulirg[0,*], ulirg[1,*], col=kleur('green')
oplot, m81[0,*], (m81[1,*]-cm81)/maxm81, col=kleur('magenta')


a =  smooth(reform(flux1-cont),5)
b =  smooth(reform(flux1-cont),10)
ori = flux1-cont
ploterror, wave1, (flux1-cont),unc1, xr=[8,14]
oplot, wave1, a, col=kleur('blue')
oplot, wave1, b, col=kleur('red')
ind=where(wave1 lt 14 and wave1 gt 8)
print, max(a[ind], ia),  max(b[ind], ib), max(ori[ind], io) 
print, wave1[ind[ia]], wave1[ind[ib]], wave1[ind[io]]
ver, wave1[ind[ia]], col=kleur('blue')
ver,  wave1[ind[ib]], col=kleur('red')
ver,  wave1[ind[io]]
; peak pos: 9.9+/- 0.05
ind=where(wave1 gt 9.8 and wave1 lt 10)
maxia=moment(a[ind])
maxib=moment(b[ind])
maxio=moment(ori[ind])

; fit Gauss to left wing
ind=where(wave1 gt 8 and wave1 lt 9.9)
w=wave1[ind] & f=ori[ind] & e=unc1[ind]
start = [9.9d0, 0.1d0, 10d0]
parinfo = replicate({value:0d0, fixed:0, limited:[0,0], $
                limits:[0D0, 0D0]}, 3)
   parinfo(2).limited=[1,0] ; not negative
   parinfo(2).limits=[0d0, 0d0]
   parinfo(0).fixed=1 ; fixed pos
   parinfo(*).value=start
   parl = mpfitfun('fit_func', w, f, e, parinfo=parinfo, perror=perrorl, $
                  yfit=yfit, bestnorm=chisq, /quiet)
   dof=n_elements(w)-2.        ; degree of freedom
   rchisq=chisq/dof
   perrorl=perrorl*sqrt(rchisq)
   ;flux = par[2]
   ;flux_err = perror[2]
   ;peakpos = par[0]
   ;peakpos_err = perror[0]
   fwhml = parl[1]
   fwhml_err = perrorl[1]

; fit Gauss to right wing
ind=where(wave1 gt 9.9 and wave1 lt 12)
w=wave1[ind] & f=ori[ind] & e=unc1[ind]
start = [9.9d0, 0.1d0, 10d0]
parinfo = replicate({value:0d0, fixed:0, limited:[0,0], $
                limits:[0D0, 0D0]}, 3)
   parinfo(2).limited=[1,0] ; not negative
   parinfo(2).limits=[0d0, 0d0]
   parinfo(0).fixed=1 ; fixed pos
   parinfo(*).value=start
   parr = mpfitfun('fit_func', w, f, e, parinfo=parinfo, perror=perrorr, $
                  yfit=yfit, bestnorm=chisq, /quiet)
   dof=n_elements(w)-2.        ; degree of freedom
   rchisq=chisq/dof
   perrorr=perrorr*sqrt(rchisq)
   ;flux = par[2]
   ;flux_err = perror[2]
   ;peakpos = par[0]
   ;peakpos_err = perror[0]
   fwhmr = parr[1]
   fwhmr_err = perrorr[1]


; let's check the gaussians ... looks good.
plot, wave1, ori, xr=[8,15]                                                     
oplot, wave1, parl[2]*mygauss(wave1, peak=parl[0], fwhm=parl[1],/maxnorm) , col=222
oplot, wave1, parr[2]*mygauss(wave1, peak=parr[0], fwhm=parr[1],/maxnorm) , col=333
ver, 9.9

; FWHM = 
print, 'FWHM left:', fwhml/2d0, fwhml_err/2d0
print, 'FWHM right:', fwhmr/2d0, fwhmr_err/2d0
print, 'FWHM tot:',  fwhml/2d0 +  fwhmr/2d0

; lets check it also for M81 ... Smith et al. didn't subtract
; continuum looks like it.
ff=m81[1,*]-cm81
ww=m81[0,*]
ee=f*.1d0 ; no errors in plot 

; fit Gauss to left wing
ind=where(ww gt 8 and ww lt 10.56)
w=ww[ind] & f=ff[ind] & e=ee[ind]
start = [10.56d0, 0.8d0, .1d0]
parinfo = replicate({value:0d0, fixed:0, limited:[0,0], $
                limits:[0D0, 0D0]}, 3)
   parinfo(2).limited=[1,0] ; not negative
   parinfo(2).limits=[0d0, 0d0]
   parinfo(0).fixed=1 ; fixed pos
   parinfo(*).value=start
   parl = mpfitfun('fit_func', w, f, e, parinfo=parinfo, perror=perrorl, $
                  yfit=yfit, bestnorm=chisq, /quiet)
   dof=n_elements(w)-2.        ; degree of freedom
   rchisq=chisq/dof
   perrorl=perrorl*sqrt(rchisq)
   ;flux = par[2]
   ;flux_err = perror[2]
   ;peakpos = par[0]
   ;peakpos_err = perror[0]
   fwhml = parl[1]
   fwhml_err = perrorl[1]

; fit Gauss to right wing
ind=where((ww gt 10.56 and ww lt 11.06) or (ww gt 11.8 and ww lt 12.59))
w=ww[ind] & f=ff[ind] & e=ee[ind]
start = [10.56d0, 0.8d0, .1d0]
parinfo = replicate({value:0d0, fixed:0, limited:[0,0], $
                limits:[0D0, 0D0]}, 3)
   parinfo(2).limited=[1,0] ; not negative
   parinfo(2).limits=[0d0, 0d0]
   parinfo(0).fixed=1 ; fixed pos
   parinfo(*).value=start
   parr = mpfitfun('fit_func', w, f, e, parinfo=parinfo, perror=perrorr, $
                  yfit=yfit, bestnorm=chisq, /quiet)
   dof=n_elements(w)-2.        ; degree of freedom
   rchisq=chisq/dof
   perrorr=perrorr*sqrt(rchisq)
   ;flux = par[2]
   ;flux_err = perror[2]
   ;peakpos = par[0]
   ;peakpos_err = perror[0]
   fwhmr = parr[1]
   fwhmr_err = perrorr[1]


; let's check the gaussians ... looks good for left wing, right
; wing not so well fitted, but fwhm seems fine.
plot, ww, ff, xr=[8,15]                                                     
oplot, ww, parl[2]*mygauss(ww, peak=parl[0], fwhm=parl[1],/maxnorm) , col=222
oplot, ww, parr[2]*mygauss(ww, peak=parr[0], fwhm=parr[1],/maxnorm) , col=333
ver, 10.56

; FWHM = 
print, 'FWHM left:', fwhml/2d0, fwhml_err/2d0
print, 'FWHM right:', fwhmr/2d0, fwhmr_err/2d0
print, 'FWHM tot:',  fwhml/2d0 +  fwhmr/2d0


stop

; make nice fig. for paper

; settings for Fig.
set_plot,'ps'
device,filename=pad+'fig_silicates_m31_nucleus.eps', /portrait,xsize=15,ysize=20,/color, bits_per_pixel=8,/encapsulated

cleanplot
!x.charsize=2.3
!y.charsize=2.3
!x.thick=5
!y.thick=5
!p.thick=5
!p.charthick=3
loadct,38
!p.multi=[0,1,3]

x0=.15
y0=.08
x1=.98
y1=.98
ddy = 0.0
dy = (y1-y0 - 2.*ddy)/3.
xt = 16.2

plot, wave1, flux1, xr=[7, 22], xs=1, ys=1, yr=[5, 65], $
     position=[x0, 2*dy+y0+2*ddy, x1, y1], $
      xtickname=[' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], yminor=2
ind=where(wave1 gt 7.5d0 and wave1 lt 15.3d0)
;oplot, wave1[ind], fit_nucCentre.final_fit[ind], col=kleur('red')
yfit1=pahfit_components(wave1,fit_nucCentre,DUST_CONTINUUM=dusts, $
                         TOTAL_DUST_CONTINUUM=dust_tot,STARLIGHT=stars, $
                         DUST_FEATURES=features, $
                         TOTAL_DUST_FEATURES=features_tot, $
                         LINES=lines,TOTAL_LINES=lines_tot)
 cont=dust_tot+stars
 oplot, wave1[ind], cont[ind], COLOR=kleur('magenta'),LINESTYLE=0
 oplot, wave1, cont, COLOR=kleur('magenta'),LINESTYLE=0
 ;rat=(cont+lines[*,7])
 ;oplot,wave1[ind],rat[ind],COLOR=kleur('blue'),LINESTYLE=0
 xyouts, xt, 55, 'M31-centre', charsize=1.4

plot, m81[0,*], m81[1,*], xr=[7, 22], yr=[0.15, 0.45], ys=1, xs=1,position=[x0, 1*dy+y0+1*ddy, x1, 2*dy+y0+2*ddy], xtickname=[' ', ' ', ' ', ' ',' ', ' ', ' ', ' ']
cm81= spline(ctptw, ctptf, m81[0,*])
oplot, m81[0,*], cm81 , COLOR=kleur('magenta')
xyouts, xt, .4, 'M81-nucleus', charsize=1.4

plot, wave1, (flux1-cont)/maxsp, xr=[7,22], xs=1, position=[x0, y0, x1, y0+dy]
hor, 0, linestyle=3
oplot, galcen[0,*], galcenprof/maxgalcen, col=kleur('blue')
oplot, m81[0,*], (m81[1,*]-cm81)/maxm81, col=kleur('red')
xyouts, xt, 1.3, 'M31-centre', charsize=1.4
xyouts, xt, 1.15, 'M81-nucleus', charsize=1.4, col=kleur('red')
xyouts, xt, 1., 'Galactic Center', charsize=1.4, col=kleur('blue')


xyouts, x0+(x1-x0)/2., 0.01, textoidl('Wavelength (\mum)'),/normal, alignment=.5, charsize=1.4
xyouts, 0.04, y0 + 2.*dy + 2*ddy + dy/2., textoidl('Intensity I_{\nu} (MJy/sr)') ,/normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.03, y0 + dy/2., textoidl('Normalized Intensity,'), /normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.07, y0 + dy/2., textoidl('flux density, \tau'), /normal, orientation=90, alignment=.5, charsize=1.4
xyouts, 0.04, y0 + dy + ddy + dy/2.,textoidl('Flux density F_{\nu} (Jy)'), /normal, orientation=90, alignment=.5, charsize=1.4
device,/close
set_plot,'x'

stop
end
