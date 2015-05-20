;This function will return the refractive index of infrasil
;The dispersion formula is from Wikipedia's Sellmeier equation for
;Fused quartz.

function nge,  lambda

n = 9.28156 + 6.72880*lambda^2/(lambda^2 - 0.44105) + 0.21307*lambda^2/(lambda^2-3870.1)

return,  sqrt(n)

end
