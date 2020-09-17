pro pahfit_ngc0337

pad = '/Users/peeters/Dropbox/Dimuthu/nucleus/'

read_spec, pad+ 'ngc0337_PAHFit_full.tbl', ngc, columns=3, startline=3
wave=reform(ngc[0,*])
flux=reform(ngc[1,*])
unc=reform(ngc[2,*])

plot, ngc[0,*], ngc[1,*]
void=pahfit(wave,PARINFO=para,/NO_FIT) 

fit_ngc0337=pahfit(wave, flux, unc, parinfo=para, /PLOT_PROGRESS, XSIZE=1200, YSIZE=600, /SCREEN, REPORT=pad+'pahfit_ngc0337_report.txt')

save, fit_ngc0337, filename=pad+'pahfit_ngc0337.xdr'

stop


end

