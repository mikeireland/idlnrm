;This function will return the refractive index of BK7

function nf2,  lambda

b =  [1.34533359e0,2.09073118e-1,9.37357162e-1]
c =  [9.97743871e-3,4.70450767e-2,1.11886764e2]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
