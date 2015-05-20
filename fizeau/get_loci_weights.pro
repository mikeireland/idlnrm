;;Given a value for the systematic error component, this function
;;finds the optimal LOCI-like weights 
function get_loci_weights, x, varx, tgt, w, cov=cov

N = n_elements(w)
hessian = dblarr(N,N)
c = dblarr(N)
for k=0,N-1 do begin
    c[k] = total(x[*,w[k]]*x[*,tgt]/varx[*,tgt])
    hessian[k,k] += total(varx[*,w[k]]/varx[*,tgt])
    for l=0,N-1 do hessian[l,k] += total(x[*,w[k]]*x[*,w[l]]/varx[*,tgt])
endfor
cov = invert(hessian) ;;NB When used below... there may be a factor of 2 missing!!!
a = cov#c
return, a
end
