;This function will return the refractive index of BK7

function nsk4,  lambda

b =  [1.32993741E+00, 2.28542996E-01, 9.88465211E-01]
c =  [7.16874107E-03, 2.46455892E-02, 1.00886364E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
