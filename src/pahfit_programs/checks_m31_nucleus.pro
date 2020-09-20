pro checks_m31_nucleus

pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'

read_spec, pad+ 'm31nuc_ll2_nucCentre.tbl', ll2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl2_nucCentre.tbl', sl2, columns=3, startline=16
read_spec, pad+ 'm31nuc_sl1_nucCentre.tbl', sl1, columns=3, startline=16
;read_spec, pad+ 'm31nuc_sl2_nucCentre_ep.tbl', sl2e, columns=3, startline=311
;read_spec, pad+ 'm31nuc_sl1_nucCentre_ep.tbl', sl1e, columns=3, startline=311

plot, sl2[0,*], sl2[1,*], xr=[5, 25]
oplot, sl1[0,*], sl1[1,*], col=300
oplot, ll2[0,*], ll2[1,*], col=213
;oplot, sl1e[0,*], sl1e[1,*], col=210, linestyle=2
;oplot, sl2e[0,*], sl2e[1,*], col=310, linestyle=2


ind1 = where(sl2[0,*] lt 7.597)
ind2 = where(sl1[0,*] gt 7.597 and sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537)
wave = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]


void=pahfit_sile(wave,PARINFO=para,/NO_FIT) 
; put all PAH betw 10-15 to zero ...
; info on regular expressions: file:///Applications/exelis/idl83/help/online_help/IDL/idl.htm#Creating%20IDL%20Programs/Components%20of%20the%20IDL%20Language/Learning_About_Regular_E.htm#strings_3486979161_298753
w=where(stregex(para.parname,'dust_feature_lambda\[1[01234]\.*',/BOOLEAN))
para[w+1].value=0d0
para[w+1].FIXED=1b
; fix one silicate component to ISM values, peakpos=9.8 and
; range=[1.8, 2.6]
w=where(para.parname eq 'line_lambda[Sil9a]')
para[w].fixed=1b
w = where(para.PARNAME eq 'line_frac_fwhm[Sil9a]')
para[w].limits=[ 1.8D/9.8d, 2.6d/9.8d]
para[w].value=2.2d/9.8d

fit_nucCentre=pahfit_sile(wave, flux, unc, parinfo=para, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+'pahfit_nucCentre_report.txt')


read_spec, 'm31nuc_ll2_nucUp.tbl', ll2, columns=3, startline=16
;read_spec, 'm31nuc_sl2_nucUp.tbl', sl2, columns=3, startline=16
read_spec, 'm31nuc_sl1_nucUp.tbl', sl1, columns=3, startline=16
read_spec, 'm31nuc_sl2_nucUp_ep.tbl', sl2, columns=3, startline=311
;read_spec, 'm31nuc_sl1_nucUp_ep.tbl', sl1e, columns=3, startline=311

plot, sl2[0,*], sl2[1,*], xr=[5, 25]
oplot, sl1[0,*], sl1[1,*], col=300
oplot, ll2[0,*], ll2[1,*], col=213

ind1 = where(sl2[0,*] lt 7.548)
ind2 = where(sl1[0,*] gt 7.548 and sl1[0,*] lt 14.537)
ind3 = where(ll2[0,*] gt 14.537)
wave2 = [reform(sl2[0,ind1]),reform(sl1[0,ind2]),reform(ll2[0,ind3])]
flux2 =  [reform(sl2[1,ind1]),reform(sl1[1,ind2]),reform(ll2[1,ind3])]
unc2 =  [ reform(sl2[2,ind1]),reform(sl1[2,ind2]),reform(ll2[2,ind3])]

void=pahfit(wave,PARINFO=para2,/NO_FIT) 
; put the T of star a free parameter
;para2[1].fixed=0b
;para2[1].limited=[1d0,1d0]
;para2[1].limits=[4000d0, 6000d0] 
; put the dust features at 6.6 to zero (is used to fill up gap left by
; continuum/noise)
w=where(para2.parname eq 'dust_feature_central_inten[6.7]')
para2[w].fixed=1b
para2[w].value=0d0

fit_nucUp=pahfit(wave2, flux2, unc2, parinfo=para2, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+ 'pahfit_nucUp_report.txt')

save, fit_nucUp, filename=pad+'pahfit_nucUp.xdr'
save, fit_nucCentre, filename=pad+'pahfit_nucCentre.xdr'

stop
end
