;This function will return the ordinary ray refractive index
;for a-BBO
;A compensated prism would be:
;n_o = nyvo4(lambda,n_e=n_e)
;devy =  1-n_e/n_o
;n_o = nabbo(lambda,n_e=n_e)
;deva =  1-n_e/n_o
;devcomp = deva + 0.14*devy

function nabbo,  lambda,  n_e = n_e
n_e = 2.3174+0.01224/(lambda^2-0.01667)-0.01516*lambda^2
n_e = sqrt(n_e)
n = 2.7471+0.01878/(lambda^2-0.01822)-0.01354*lambda^2

return,  sqrt(n)

end
