pro m31nuc_nucUP_redo

pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'

read_spec, 'm31nuc_ll2_nucUp_correct.tbl', ll2, columns=3, startline=239
read_spec, 'm31nuc_sl2_nucUp_correct.tbl', sl2, columns=3, startline=311
read_spec, 'm31nuc_sl1_nucUp_correct.tbl', sl1, columns=3, startline=311
;read_spec, 'm31nuc_sl2_nucUp_ep.tbl', sl2, columns=3, startline=311
;read_spec, 'm31nuc_sl1_nucUp_ep.tbl', sl1e, columns=3, startline=311

ind1 = where(sl2[0,*] lt 7.53);7.548)
ind2 = where(sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537 and ll2[0,*] lt 20.978)
wave2 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux2 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc2 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]

plot, sl2[0,*], sl2[1,*], xr=[5, 25]
oplot, sl1[0,*], sl1[1,*], col=300
oplot, ll2[0,*], ll2[1,*], col=213
oplot, wave2, flux2, col=400

void=pahfit(wave2,PARINFO=para2,/NO_FIT)
; put the T of star a free parameter
;para2[1].fixed=0b
;para2[1].limited=[1d0,1d0]
;para2[1].limits=[4000d0, 6000d0] 
; put the dust features at 6.6 to zero (is used to fill up gap left by
; continuum/noise)
w=where(para2.parname eq 'dust_feature_central_inten[6.7]')
para2[w].fixed=1b
para2[w].value=0d0

fit_nucUp=pahfit(wave2, flux2, unc2, parinfo=para2, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+ 'pahfit_nucUp_report_correct.txt')

save, fit_nucUp, filename=pad+'pahfit_nucUp_correct.xdr'

stop

end
