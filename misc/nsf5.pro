;This function will return the refractive index of SF5, also called N-SF5

function nsf5,  lambda

b =  [1.46141885E+00, 2.47713019E-01, 9.49995832E-01]
c =  [1.11826126E-02, 5.08594669E-02, 1.12041888E+02]
b =  [1.52481889E+00, 1.87085527E-01, 1.42729015E+00]
c =  [1.12547560E-02, 5.88995392E-02, 1.29141675E+02]
n = 1.0
for i =  0, n_elements(b)-1 do $
 n +=  b[i]*lambda^2/(lambda^2-c[i])

return,  sqrt(n)

end
