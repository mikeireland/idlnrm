;This function will return the refractive index of infrasil
;The dispersion formula is from Wikipedia's Sellmeier equation for
;Fused quartz.

function ninfrasil,  lambda

;b =  [1.03961212, 2.317923e-1, 1.01046945]
;c =  [6.00069867e-3,  2.00179144e-2,  1.03560653e2]
b =  [0.69616630, 0.40794260, 0.89747940]
c =  [0.068404300, 0.11624140, 9.8961610]

n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i]^2)

return,  sqrt(n)

end
