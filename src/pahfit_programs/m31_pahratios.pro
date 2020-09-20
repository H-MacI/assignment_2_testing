pro m31_pahratios

pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'

; get pahfit results
restore, pad+'pahfit_nucUp_correct.xdr' ;fit_nucUp
restore, pad+'pahfit_nucUp_nostellar_correct.xdr' ;fit_nucUp_nostellar


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

stellarcomp =  smooth(flux1, 4)
normf_stellarcomp=stellarcomp[where(wave1 gt 5.38 and wave1 lt 5.4)]
normf_flux2=flux2[where(wave2 gt 5.38 and wave2 lt 5.4)]
flux2_nostellar = flux2 - stellarcomp/normf_stellarcomp[0] * normf_flux2[0]


mf=pahfit_main_feature_power(fit_nucUp,wave1, flux1, unc1)
mf_nostellar=pahfit_main_feature_power(fit_nucUp_nostellar,wave2, flux2_nostellar, unc2)


index_bands=[0,1,3, 10]  ; 6.2, 7.7complex, 8.6, 17complex
; 11.2complex -- 4
for i=0, 3 do begin
print, mf[index_bands[i]].name+'/11.2um'
print, 'PAHFITmethod  ', mf[index_bands[i]].strength/mf[4].strength
print, 'nuclearmethod  ', mf_nostellar[index_bands[i]].strength/mf_nostellar[4].strength
endfor

yfit2=pahfit_components(wave2,fit_nucUp_nostellar, $
                         DUST_CONTINUUM=dusts2, $
                         TOTAL_DUST_CONTINUUM=dust_tot2,STARLIGHT=stars2, $
                         DUST_FEATURES=features2, $
                         TOTAL_DUST_FEATURES=features_tot2, $
                         LINES=lines2,TOTAL_LINES=lines_tot2)
yfit1=pahfit_components(wave2,fit_nucUp, $
                         DUST_CONTINUUM=dusts1, $
                         TOTAL_DUST_CONTINUUM=dust_tot1,STARLIGHT=stars1, $
                         DUST_FEATURES=features1, $
                         TOTAL_DUST_FEATURES=features_tot1, $
                        LINES=lines1,TOTAL_LINES=lines_tot1)
plot, wave2, fit_nucUp_nostellar.final_fit-(dust_tot2+stars2)
oplot, wave2, fit_nucUP.final_fit-(dust_tot1+stars1), col=kleur('red')

stop
   for i=0,(size(features,/DIMENSIONS))[1]-1 do begin 
     rat=(features[*,i])
     oplot,wave2,rat,COLOR=kleur('blue'),LINESTYLE=0
  endfor   

stop

end
