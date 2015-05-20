;This function will return the refractive index of SFL6,
;also called N-SF6

function nsfl6,  lambda

b =  [1.77931763E+00,   3.38149866E-01, 2.08734474E+00]
c =  [1.33714182E-02,   6.17533621E-02, 1.74017590E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
