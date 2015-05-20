;This function will return the refractive index of BK7

function nsf11,  lambda

b =  [1.73759695E+00,   3.13747346E-01, 1.89878101E+00]
c =  [1.31887070E-02,   6.23068142E-02, 1.55236290E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
