;;Given a wavelength and a magnitude, this function returns the 
;;photon rate in photons/m^2/micron/s.
;;e.g. K-band, 1m diameter mirror mag 8.5-> 27600 photons/s
;; mag2prate(2.2,8.5)*!pi*0.5^2*0.05*0.4/500
;;Black-body for 1 arcsec square pixels:
;; lambda = 20000+findgen(4000)
;; total(planck(lambda,283))*!pi*50.^2/!pi*1e-7/206000.^2/e_p*0.95

function mag2prate,  wavelengths,  mag_vector

;;Flux in jy
jy = mag2jy(wavelengths,  mag_vector)
;;Flux per micron
f_m =  jy*1e-26*2.998e8/(wavelengths*1e-6)^2*1e-6
;;Energy of a photon
e_p = 6.626e-34*2.998e8/(wavelengths*1e-6)
;;Now flux in photons.
return,  f_m/e_p

end
