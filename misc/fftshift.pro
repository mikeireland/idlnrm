;;This function shifts a 1D array some number of pixels that
;;doesn't have to be an integer. For an array that isn't "smooth",
;;i.e. that has power right up to (and beyond) Nyquist, then it
;;won't work properly.

function fftshift, array, fraction

L = n_elements(array)
x = findgen(L)
x = ((x + L/2) mod L) - L/2
ftarray = fft(array,1)
return, fft(ftarray*exp(complex(0,2*!pi*x*fraction/L)),-1)

end
