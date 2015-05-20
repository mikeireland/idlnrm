;This function will return the refractive index of BK7

function nsf10,  lambda
b =  [1.62153902E+00,   2.56287842E-01, 1.64447552E+00]
c =  [1.22241457E-02,   5.95736775E-02, 1.47468793E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
