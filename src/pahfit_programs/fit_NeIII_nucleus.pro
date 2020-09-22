function fit_func,x, p

f = p[2]*mygauss(x, peak=p[0], fwhm=p[1],/maxnorm) + p[3] + p[4] * x 
return, f

end

pro fit_NeIII_nucleus

pad1 = '/Users/peeters/Dropbox/Dimuthu/nucleus/'
pad2='/Users/peeters/data/Obs/SIRTF/M31/'
cube = readfits(pad2+'m31nuc_ll_LL2_cube.fits')
unc = readfits(pad2+'m31nuc_ll_LL2_cube_unc.fits')
read_spec, pad2+'wave_ll2.tbl', dummy, columns=3, startline=239
wave=reform(dummy[0,*])

; convert fluxes to  W/m2/um/sr
s = size(cube)
wave3d = rebin(reform(wave,1,1,s[3]),s[1],s[2],s[3])
cubeW = 1d6 * 3d-12 * cube / wave3d^2d0
uncW = 1d6 * 3d-12 * unc / wave3d^2d0

; to save results
fitNeIII=dblarr(s[1], s[2], 6)

tops, pad1+'fit_NeIII_nucleus.ps'
!p.multi=[0,2,2]

For x_loop = 0, s[1]-1 do begin
 For y_loop = 0, s[2]-1 do begin

   ind = where(wave gt 15. and wave lt 16.2)
   w = wave[ind]
   f = reform(cubeW[x_loop, y_loop, ind])
   e = reform(uncW[x_loop, y_loop, ind])
   area = max(f)
   start = [15.5d0, 0.17d0, area, 0d0, 0d0]
   parinfo = replicate({value:0d0, fixed:0, limited:[0,0], $
                limits:[0D0, 0D0]}, 5)
; no negative fluxes
   parinfo(2).limited=[1,0]
   parinfo(2).limits=[0d0, 0d0]
; fix FWHM
   parinfo(1).fixed=1.
; restrict peak pos
   parinfo(0).limited=[1,1]
   parinfo(0).limits=[15.5, 15.6]
   parinfo(*).value=start
   par = mpfitfun('fit_func', w, f, e, parinfo=parinfo, perror=perror, $
                  yfit=yfit, bestnorm=chisq, /quiet)
   dof=n_elements(wave)-2.        ; degree of freedom
   rchisq=chisq/dof
   perror=perror*sqrt(rchisq)
   plot, w, f, title='pixel ('+string(x_loop)+','+string(y_loop)+')', ys=1
   oplot, w, fit_func(w, par), col=kleur('red')
   fitNeIII[x_loop, y_loop, 0:2]=par[0:2]
   fitNeIII[x_loop, y_loop, 3:5]=perror[0:2]

 endfor
endfor


plot, fitNeIII[*,*,0], ps=4, yr=[15.47, 15.7]
dummy= moment(fitNeIII[*,*,0])
print, dummy
hor, dummy[0], col=kleur('red')
plot, fitNeIII[*,*,1], ps=4
dummy= moment(fitNeIII[*,*,1])
print, dummy
hor, dummy[0], col=kleur('blue')
plot, fitNeIII[*,*,2], ps=4

cleanplot
plotimage, bytscl( fitNeIII[*,*,2])

tops, pad1+'fit_NeIII_nucleus.ps',/finish

loadct, 38
NeIII=readfits(pad2+'m31nuc_LL2_NeIIImap_integrated.fits', hdr)
NeIII = reform(fitNeIII[*,*,2])
writefits, pad2+'m31nuc_LL2_NeIIImap_fittedGaussian.fits',NeIII, hdr

stop
end

