;This function will return the refractive index of BK7

function npk51,  lambda

b =  [1.15610775,0.153229344,0.785618966]
c =  [0.00585597402,0.0194072416,140.537046]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
