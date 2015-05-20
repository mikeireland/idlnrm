;This function is based on Christian Hummel's 2001 Michelson School
;web document. It's designed to be used directly for fitting...
;Updates from original: vectorized multiple times
;params are... 
; T0: The time of periastron passage
; P:  The period
; a:  The semi-major axis
; e:  The eccentricity
; n:  Capital omega ( an orbital angle )
; w:  Little omega
; i:  The inclination
;NB for differentiation, see mathematica/binary_diff.nb

function binary_noderiv, params, jds, projfact=projfact,  deriv = deriv

t = jds-params[0]
P = params[1]
a = params[2]
e = abs(params[3])
n = params[4]*!pi/180.
w = params[5]*!pi/180.
i = params[6]*!pi/180.

;The mean anomaly 
;(p,t) -> M 
M = 2*!pi*(t mod p)/p
;The eccentric anomaly, M = E - esinE
;(M,e) -> (bE,e) ;Tr2
;1 = dbEdM - e dbEdM cos(bE) 
;0 = dbEde - e*dbEde*cos(bE) - sin(bE)
bE = M+e*sin(M)+e^2/2*sin(2*M)
for k = 0,4 do bE=bE+(M-bE+e*sin(bE))/(1-e*cos(bE))
;The `true anomaly'. With a pi ambiguity,
;nu = 2*atan(sqrt((1+e)/(1-e))*tan(bE/2))
;(bE,e) -> (nu,e) ;Tr3
nu=2*atan(sqrt((1+e)/(1-e))*sin(bE/2), cos(bE/2))
;Derivatives are now for alpha (just an offset from nu)
;From mathematica...

;Offset for nu
;(nu,w) -> alpha ;Tr4
alpha=nu+w

;Final calculations (square brackets mean trivial):
;(alpha,e,i) [+a,n] -> (rho,theta) Tr5
;We have dAd(p,t,e,w), with alpha=A. Only complex for
;e, where we need:
;drhode = drhodA.dAde + drhodnu.dAde + drhode
;dthde = dthdA.dAde + dthde
;Also, drhodp = drhodA.dAdp etc...
;!!! Much of the following can be sped-up by pre-calculation.
rho=a*(1-e^2)/(1+e*cos(nu))*sqrt(cos(alpha)^2+sin(alpha)^2*cos(i)^2)
theta=atan(sin(alpha)*cos(i),cos(alpha))+n


projfact = (cos(alpha))^2 ; + e*cos(w)

return, [[rho], [theta*180/!pi]]
end
