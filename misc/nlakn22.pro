;This function will return the refractive index of BK7

function nlakn22,  lambda

b =  [1.14229781E+00, 5.35138441E-01, 1.04088385E+00]
c =  [5.85778594E-03, 1.98546147E-02, 1.00834017E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
