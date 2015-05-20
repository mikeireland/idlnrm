pro rdNextGen,name,w,f,b,teff,logg,mh
;------------------------------------
; this routine reads a gzipped NextGen 
; spectrum into the variables given.
;
;-inputs:
; name: input filename
;-outputs:
; teff,logg,mh: output Teff, log(g) and metallicity
; w: wavelength array, Angstroem
; f: flux array, erg/s/cm^2/cm
; b: Planck function array, erg/s/cm^2/cm
;
; version 1.2
; written by Peter H. Hauschildt, August 1997

spawn,string("gunzip -c "+name+' >rdyng_tmp || rm rdyng_tmp'),ierr

openr,myunit,'rdyng_tmp',/get_lun
readf,myunit,teff,logg,mh
readf,myunit,npoint

w = fltarr(npoint)
f = fltarr(npoint)
b = fltarr(npoint)

readf,myunit,w
readf,myunit,f
readf,myunit,b
close,myunit
free_lun,myunit
spawn,string("rm rdyng_tmp"),ierr
end
