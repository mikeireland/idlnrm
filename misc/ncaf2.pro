;This function will return the refractive index of BK7

function ncaf2,  lambda

b =  [0.5675888, 0.4710914, 3.8484723]
c =  [0.050263605,  0.1003909,  34.649040]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i]^2)

return,  sqrt(n)

end
