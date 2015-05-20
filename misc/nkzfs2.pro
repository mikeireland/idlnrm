;This function will return the refractive index of N-KZFS2

function nkzfs2,  lambda

b =  [1.23697554, 0.153569376,0.903976272]
c =  [0.00747170505, 0.0308053556, 70.1731084]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
