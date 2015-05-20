;This function will return the ordinary ray refractive index
;for YVO4

function nlinbo3,  lambda,  n_e = n_e

n_e =  4.5820 + 0.099169/(lambda^2 - 0.04443) - 0.021950*lambda^2
n_e =  sqrt(n_e)
n = 4.9048+0.11768/(lambda^2-0.04750)-0.027169*lambda^2

return,  sqrt(n)

end
