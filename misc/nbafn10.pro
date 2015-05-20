;This function will return the refractive index of BAFN10
;Also known as N-BAF10, see:
;http://www.acoptics.com/optical_glass.htm

function nbafn10,  lambda

b =  [1.58514950E+00, 1.43559385E-01, 1.08521269E+00]
c =  [9.26681282E-03, 4.24489805E-02, 1.05613573E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
