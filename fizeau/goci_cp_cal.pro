;;This program completely replaces calibrate_v2_cp.pro, using
;;something similar to the LOCI algorithm (but with global support).
;;Before this is called, the cal4src matrix should be
;;roughly OK, i.e. with very bad data removed. Also, bad_baselines
;;should be pre-flagged.
;;
;;NB it is essential that there are (significantly) less calibrators
;;than the number of independent closure-phases.
;;
;;Inputs...
;; cubeinfo_file: This should contain everything we need!
;; model : this is a 9 element vector for 3 point sources.

pro goci_cp_cal, cubeinfo_file, cp_cal=cp, cp_err_cal=cp_err, bsdir=bsdir, delta2s=delta2s, printit=printit, model=model
common loci_parms, lp

root_dir = !ROOT_DIR
; File and Directory Options
if (not keyword_set(src_for_bad_holes)) then src_for_bad_holes = 0
dir =  ''
pos = strpos(cubeinfo_file,'/',  /reverse_search)
if (pos ne -1) then begin
    dir =  strmid(cubeinfo_file,  0,  pos) + '/'
    cubeinfo_file = strmid(cubeinfo_file,  pos+1)
endif
if n_elements(rho) ne n_elements(theta) then begin
 print, 'rho and theta must have the same number of elements!'
 stop
endif
if not keyword_set(bsdir) then bsdir=dir
print, bsdir
restore, dir +cubeinfo_file
;;Some sanity checking
src =  where(total(clog.cal4src, 1) gt 0,  complement = cal)
ncal = n_elements(cal)
nsrc = n_elements(src)
ncubes = nsrc+ncal
if (cal[0] eq -1) then begin
   print,  '!!! No Calibrators !!! You need a calibrator for this program'
   stop
endif 
if (src[0] eq -1) then begin
    print,  '!!! No Sources !!! You need a target for this program'
    stop
endif 
if (olog.logflag ne 2) then begin
    print,  'Please calibrate the data normally first (e.g. using calibrate_v2_cp).'
    stop
endif
;;Set the bad vectors
bad_holes =  clog.bad_holes
bad_baselines =  clog.bad_baselines
bad_bispect =  clog.bad_bispect

;; Lets get our mask info
mf_filestring = root_dir + '/templates/' + plog.mf_file
restore, mf_filestring

;;Get the important stuff.
aquan = get_used_quantities(bsdir+plog.bs_names[0],  bad_baselines,  bad_bs, root_dir, /nows)
quan = replicate(aquan, ncubes)
cal_quan = replicate(aquan, nsrc)
for i=0,ncubes-1 do quan[i]=get_used_quantities(bsdir+plog.bs_names[i],  bad_baselines,  bad_bs, root_dir, /nows)

;;Find the conversion from phase to closure-phase:
n_ind =  n_baselines-n_holes+1
T = dblarr(n_bispect,n_baselines)
for k=0,n_bispect-1 do begin
    T[k,bs2bl_ix[0,k]] =1/3.
    T[k,bs2bl_ix[1,k]] =1/3.
    T[k,bs2bl_ix[2,k]] =-1/3.
endfor
x     = dblarr(n_ind,ncubes)
varx  = dblarr(n_ind,ncubes)
minvar = dblarr(ncubes)
ix = indgen(n_ind)

;;Finally work on the calibration, by creating matrices etc.
delta2s = dblarr(nsrc)
ncp = n_elements(aquan.cp)
mod_cp_cov = dblarr(ncp,ncp,ncubes)
phi = dblarr(ncp,ncubes)

;;First, find the weighted mean covariance matrix 
used = total(clog.cal4src,2) + total(clog.cal4src,1)
used = where(used gt 1)
nused = n_elements(used)
traces = dblarr(nused)
C_all = quan[used[0]].cp_cov
C_all[*]=0
;;used[4] and used[13] are totally stuffed. Re-analysis???
;;These are actually cubes 4 and 13.
for i=0,n_elements(used)-1 do begin
 traces[i]=trace(quan[used[i]].cp_cov)
 C_all += quan[used[i]].cp_cov;/traces[i]
endfor

;;M is a projection onto the space of linearly independent closure-phases.
M = T#transpose(T);;
D3 = la_eigenql(M,  eigenvectors = U3)
;;NB I think P3 is actually the transpose of what I'd like it
;;to be. So use caution interpreting the following code.
P3=transpose(U3[*,n_bispect-n_ind:*])
;;This next line has to be computed once per pa and sep (but
;;not once per source cube)
D2_notprimed = la_eigenql(P3#C_all#transpose(P3),  eigenvectors = U2)
P4 = transpose(U2)#P3

;;Now we go through each (used?) cube and compute the x and
;;binary contrast variables (takes ~1/3 of the time)
for k=0,ncubes-1 do begin
    x[*,k] = P4#quan[k].cp
    D2 = P4#quan[k].cp_cov#transpose(P4)
    varx[*,k] = D2[ix,ix]
    minvar[k] = 0.9*median(varx[*,k])
    varx[*,k] = varx[*,k] > minvar[k]
endfor

varx_c = varx[*,src]
x_c = x[*,src]

;;OK, Now we've re-arranged the problem into one on the
;;linearly-independent set of closure-phases
for i = 0, nsrc-1 do begin 
    ;;Preliminaries: w is the vector of cals, tgt is the current
    ;;target.
    tgt=src[i]
    w = where(clog.cal4src[*, tgt] eq 1, N)
    if (N lt 1) then begin
        print, 'Need at least 1 calibrator per target!'
        stop
    endif

    ;;Now find the model value of x
    if (keyword_set(model)) then begin
       u1_ = u_save[bs2bl_ix[0,*],i]
       u2_ = u_save[bs2bl_ix[1,*],i]
       u3_ = -u_save[bs2bl_ix[2,*],i]       
       v1_ = v_save[bs2bl_ix[0,*],i]
       v2_ = v_save[bs2bl_ix[1,*],i]
       v3_ = -v_save[bs2bl_ix[2,*],i]
       zeros = replicate(0d0, n_elements(v3_))
       t3data={u1:u1_,u2:u2_,u3:u3_,v1:v1_,v2:v2_,v3:v3_,t3amp:zeros,t3phi:zeros,t3amperr:zeros,t3phierr:zeros}
       m1t3 = binary_t3data(double([model[0:2],0.1,0.1]),t3data=t3data)
       m2t3 = binary_t3data(double([model[3:5],0.1,0.1]),t3data=t3data)
       m3t3 = binary_t3data(double([model[6:8],0.1,0.1]),t3data=t3data)
       modelcp =   (m1t3.t3phi  +m2t3.t3phi  +m3t3.t3phi)*!dpi/180
       model_x = P4#modelcp
    endif else model_x=0

;;Now for the actual calibration, which has to know about target and
;;calibrators (other than just sky rotation). Takes about 1/3 of the
;;time
    xtemp=x ;;xtemp has the model closure-phases removed.
    xtemp[*,tgt]-=model_x
    lp={x:xtemp,varx:varx,tgt:tgt,w:w}
    Delta2 = newton( 0, 'delta2_chi2rm1', tolf=1e-4, /double) 
;; *** Note that the next line is critical, as it is actually where
;;the calibration is done. Kind of opaque coding... ***
    dummy = delta2_chi2rm1(Delta2, varx_c=varx_c1, x_c=x_c1)
    print, dummy
    x_c[*,i]=x_c1+model_x ;;Now add the model back in.
    if (Delta2 lt 0) then begin
      varx_c1 -= Delta2
      Delta2=0
    endif
    varx_c[*,i]=varx_c1
    delta2s[i]=Delta2
 
    ;;At this point, the problem is totally defined by:
    ;;P4, x, varx+delta2

    a = get_loci_weights(xtemp, varx+Delta2, tgt, w, cov=cov)

    cal_quan[i].cp = quan[tgt].cp
    for k=0,n_elements(w)-1 do cal_quan[i].cp -= a[k]*quan[w[k]].cp
;    cal_quan[i].v2[*] = 0
;    for k=0,n_elements(w)-1 do cal_quan[i].v2 += a[k]*quan[w[k]].v2
;    cal_quan[i].v2 = quan[tgt].v2/cal_quan[i].v2
    ;;This next formula doens't respect covariances.
    var = quan[tgt].cp_err^2 + Delta2
    for k=0,n_elements(w)-1 do var += a[k]^2*quan[w[k]].cp_err^2   
    for k=0,n_elements(w)-1 do cal_quan[i].cp_err = sqrt(var)

    if keyword_set(printit) then begin
        print, 'Median Error: ', median(sqrt(varx_c[*,i]))*180/!pi
        print, 'Extra Error: ', sqrt(delta2s[i])*180/!pi
        print, 'Weights: ', a, format='(A,20F6.2)'
    endif
    print, 'Done i=', i+1, ' of ', nsrc, ' target cubes'
endfor

cp=cal_quan.cp
cp_err=cal_quan.cp_err
read_oidata, clog.primary_oifits, oiarray, oitarget,oiwavelength,oivis,oivis2, oit3
for i=0,n_elements(oit3)-1 do *(oit3[i].t3phi) = cp[i]*180/!pi
for i=0,n_elements(oit3)-1 do *(oit3[i].t3phierr) = cp_err[i]*180/!pi
write_oidata, clog.outname+'mrg_goci.oifits',oiarray,oitarget,oiwavelength,oivis,oivis2,oit3
save, P4, x_c, varx_c, file=clog.outname+'_x.idlvar'

end
