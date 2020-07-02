; JDM 2001Dec01		Returns Vis and Phase of Binary Disk for u,v
; ACR 2020Apr10          
; INPUTS:
; Model of a triple star, each with UD sizes and a ratio.
; Params:
;   params(0) = Separation (mas) 
;   params(1) = Position Angle (degs, 
;               E of N, pointing from primary -> secondary)
;   params(2) = Brightness Ratio of Primary over Secondary
;   params(3) = UD size of primary (mas)
;   params(4) = UD size of secondary (mas)
;   params(5) = separations 2 (mas)
;   params(6) = PA or primary--> tertiary
;   params(7) =  UD size of tertiary
;   params(8) = brightness ratio of primary over tertiary
; u,v : units of Number of wavelengths (i.e., rad-1)
;
; outputs: 
;   visib  -- visibility
;   phases -- phases in degrees (-180,180)
;;This is all computed using the complex notation in the first section of the iss notes

pro triple_disks,u,v,params,visib,phases

delta_dec = mas2rad(params(0)*sin( (params(1)+90)*!dpi/180.))
delta_ra = -1*mas2rad(params(0)*cos( (params(1)+90)*!dpi/180.))
 ;recall +ra points at pa 90
delta_dec2 = mas2rad(params(5)*sin( (params(6)+90)*!dpi/180.))
delta_ra2 = -1*mas2rad(params(5)*cos( (params(6)+90)*!dpi/180.))


x=sqrt(u*u + v*v)
;secondary_flux = 1.0/ (1+params(2))
;primary_flux = 1.0 - secondary_flux

;;set total flux to 1
primary_flux = 1/(1+1/params[2]  + 1/params[8])
secondary_flux = 1./params[2]/(1+1/params[2]  + 1/params[8])
tertiary_flux = 1./params[8]/(1+1/params[2]  + 1/params[8])
;stop
; Primary
diameter=mas2rad(abs(params(3)))
intercept=primary_flux
index=where(x eq 0,count)
if (count gt 0) then x(index)=1e-15

      f=2*beselj(!dpi*diameter*x,1) / (!dpi*diameter*x)
if (params(3) ge 0) then f=abs(intercept*f) else f=abs(intercept/f)
vis_primary = f

; Secondary
diameter=mas2rad(abs(params(4)))
intercept=secondary_flux

      f=2*beselj(!dpi*diameter*x,1) / (!dpi*diameter*x)
if (params(4) ge 0) then f=abs(intercept*f) else f=abs(intercept/f)
vis_secondary=f

;tertiary
diameter = mas2rad(abs(params(7)))
intercept = tertiary_flux
f=2*beselj(!dpi*diameter*x,1) / (!dpi*diameter*x)
if ( params(7) ge 0) then f = abs(intercept*f) else f=abs(intercept/f)
vis_tertiary=f

; Using boden notation (I hope!) from michelson book
phase_factor = dcomplex( cos(-2.0*!dpi*(u*delta_ra+v*delta_dec))  ,$
                         sin(-2.0*!dpi*(u*delta_ra+v*delta_dec)))

phase_factor3 = dcomplex( cos(-2.0*!dpi*(u*delta_ra2+v*delta_dec2)),sin(-2.0*!dpi*(u*delta_ra2+v*delta_dec2)))

complex_vis = vis_primary + vis_secondary*phase_factor + vis_tertiary*phase_factor3

ri2at, double(complex_vis), imaginary(complex_vis), visib, phases
phases= mod360(phases)

;;stop
end
