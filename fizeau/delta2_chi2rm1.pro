;; This function returns the chi-squared
function delta2_chi2rm1, delta2, varx_c=varx_c, x_c=x_c
common loci_parms
varx = lp.varx+delta2[0]
a = get_loci_weights(lp.x, varx,lp.tgt, lp.w)
nx = (size(lp.x))[1]
chi2=0.
varx_c=fltarr(nx)
x_c=fltarr(nx)
for i=0,nx-1 do varx_c[i] = (varx[i,lp.tgt] + total(a^2*varx[i,lp.w]))
for i=0,nx-1 do x_c[i]    =  lp.x[i,lp.tgt] - total(a*lp.x[i,lp.w])
for i=0,nx-1 do chi2 += x_c[i]^2/varx_c[i]
return, chi2/nx - 1
end 
