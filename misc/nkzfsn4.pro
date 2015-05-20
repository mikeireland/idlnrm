;This function will return the refractive index of BK7

function nkzfsn4,  lambda

b =  [1.37994218E+00, 1.68496708E-01, 8.74885726E-01]
c =  [8.91159699E-03, 4.05334070E-02, 6.96628237E+01]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
