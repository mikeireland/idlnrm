;This function will return the refractive index of BK7

function nkzfs11,  lambda

b =  [1.3322245,0.28924161, 1.15161734]
c =  [0.0084029848,0.034423972,88.4310532]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
