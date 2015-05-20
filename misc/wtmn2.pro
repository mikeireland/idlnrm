;A weighted average routine...
;Now matches the Wikipedia entry, and only calculates the
;error one way - where the formal error is scaled by the square
;root of chi^2.
;
; original                                                MJI  ????
; changed MSE to option;                                  PGT Nov03
; changed to wtmn2					  MJI Jul09

function wtmn2, x, sigma, sdev, chi2=chi2

ix = where(sigma gt 0.)
if (ix[0] eq -1) then begin
 sdev =  0.0
 MSE = 0.0
 wts = 1.
 return,  x[0]
endif
weights=1./sigma[ix]^2
nelt = float(n_elements(x[ix]))
sumwt = total(weights)
mean = total(x[ix]*weights)/sumwt 
e1 = 1./sumwt
chi2 = (nelt gt 1)? total(abs(x[ix]-mean)^2*weights)/(nelt-1):1
sdev=sqrt(e1*chi2)

return, mean
end

