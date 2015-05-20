;This function will return the ordinary ray refractive index
;for Calcite

function ncalcite,  lambda,  n_e = n_e

n_e =2.18438+0.0087309/(lambda^2-0.01018)-0.0024411*lambda^2
n_e = sqrt(n_e)
n = 2.69705+0.0192064/(lambda^2-0.01820)-0.0151624*lambda^2

return,  sqrt(n)

end
