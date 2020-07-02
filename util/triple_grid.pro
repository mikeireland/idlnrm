;;This procedure runs a grid search for a three source model (as per
;;binary_grid.pro) but has the utility to window in separation and PA
;;for both companions.
;; ;; Currently does not use visibilities!! This need fixing for
;;future versions.
;;-Aaron C Rizzuto
;;printfile:   print results to a file instead of the screen.
;;NB: -Common blocks can be reset with .reset_session
;; - Need a routine that finds the Hessian matrix inverse...
;;TODO:
;; - This script doesn't work for moderate contrast-ratio
;;   binaries. Needs fixing...?


pro triple_grid,  infile,  printfile = printfile, usevis = usevis_in,  nfr = nfr,  $
 init_crat = init_crat,  maxsep = maxsep, minsep =minsep,  $
 apriori = apriori_in,  proj = proj_in,  cor = cor,  covar = covar,  bestp = bestp,  perror = perror, scale_err=scale_err,$
  force_min_sep=force_min_sep0, t3errscale=t3errscale, significance=significance, $
 redchi2=redchi2, max_snr=max_snr,max2sep=maxsep2,min2sep=minsep2,minpa=minpa,maxpa=maxpa,min2pa=minpa2,max2pa=maxpa2

;;Common block for fitting
common t3block,  t3data,  apriori, vis2data, usevis,  cp_cinv,  proj, proj_err, force_min_sep

if (keyword_set(max_snr) eq 0) then max_snr = 5.0
if (keyword_set(usevis_in) eq 0) then usevis_in = 0
usevis = usevis_in
if (keyword_set(apriori_in) eq 0) then apriori_in = {val:[-1, -1, -1],  err:[0, 0, 0]}
apriori = apriori_in
if (keyword_set(printfile)) then openw,  1,  printfile, width=150 else printfile = ''
;;Convert errors in radians^(-2) to degrees^(-2)
if not keyword_set(cp_cov_in) then cp_cov = [-1.] else cp_cov = cp_cov_in/(!pi/180.)^2
cp_cov_null = cp_cov
cp_cinv = invert(cp_cov)
if not keyword_set(confidence_level) then confidence_level = 0.999
p = printfile
if (p ne '') then print,  '********** Binary fit for file: ',  p,  '**********' $
else print,  '********** Starting Binary fit. No save file **********'
if not keyword_set(nsim) then nsim =long(10000) 
if (keyword_set(init_crat) eq 0) then init_crat =  250
if keyword_set(fix_crat) then init_crat = fix_crat


extract_t3data,  file = infile,  t3data
if (keyword_set(t3errscale)) then t3data.t3phierr *= t3errscale
w = where(t3data.flag)
if (w[0] ne -1) then t3data[w].t3phierr = 1e4
;;If we have the matrix, project the closure-phases onto the
;;linearly dependent set of closure-phases.
if (keyword_set(proj_in)) then begin
 proj = proj_in
 nproj = (size(proj))[2]
 if (cp_cov[0] eq -1) then stop
 cp_cov = transpose(proj)#cp_cov#proj
 ix =  indgen(nproj)
 proj_err = sqrt(cp_cov[ix, ix])
 cp_cov = [-1]
 cp_cinv = [-1]
endif else begin
  proj = [-1]
  proj_err = 0
endelse
extract_vis2data,  file = infile,  vis2data
read_oidata,infile,oiarray,oitarget
target_name=oitarget.target
target_name=target_name[0]

;Enforce a maximum SNR of 5.0 on visibility data.
vis2data.vis2err = sqrt(vis2data.vis2err^2 + (vis2data.vis2data/max_snr)^2)


;;Find the number of rotations on the sky (i.e. independent cubes) by
;;finding how many times the first baseline in the file is repeated,
;;where we define "baseline" as containing uniqu stations indices.
;;WARNING: breaks for a >100 hole mask.
ix = vis2data.sta_index[0] + 100*vis2data.sta_index[1]
nrotations = n_elements(where(ix eq ix[0]))
nv2 = n_elements(vis2data)/nrotations
nt3 = n_elements(t3data)/nrotations
;Number of independent degrees of freedom in the data. The complicated
;expression at the end is (n_holes - 1)
print, "usevis: " + strtrim(usevis,2)
if (p ne '') then printf,1,"usevis: ", usevis
ndf_cp =  nrotations*(nv2 - (3*nt3/nv2+1) )

if (keyword_set(nfr)) then begin
 ;Median single-frame noise in radians
 med_noise =  median(t3data.t3phierr)*sqrt(nfr)/180*!pi
 ;;Weight the number of dependent and independent closure-phases 
 ;;so that the number of degrees of freedom is intermediate between
 ;;these values for and SNR of 1.0. 
 ndf_cp =  (n_elements(t3data)*med_noise + ndf_cp*(1./med_noise))$
       /(                   med_noise +      1./med_noise)
endif

;;currently does nothing!!
if (usevis le 0) then begin
    ndf = ndf_cp
    print, nrotations, (nv2 - (3*nt3/nv2+1) ), ndf, $
      format = '("Degrees of Freedom:",I5," cubes x ",I5," independent clp data = ",I6)'
    if (p ne '') then printf,1,nrotations, (nv2 - (3*nt3/nv2+1) ), ndf, $
      format = '("Degrees of Freedom:",I5," cubes x ",I5," independent clp data = ",I6)'
  endif else begin
   ndf =  ndf_cp + nrotations*(nv2 - 2*(usevis-1))
   print, nrotations,nv2, (nv2 - (3*nt3/nv2+1)), ndf, $
         format = '("Degrees of Freedom:",I5," cubes x (",I5," vis + ",I5," clp data) = ",I6)'
     if (p ne '') then printf,1,nrotations,nv2, (nv2 - (3*nt3/nv2+1)), ndf, $
         format = '("Degrees of Freedom:",I5," cubes x (",I5," vis + ",I5," clp data) = ",I6)'
  endelse


r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./4/maxr)
if (keyword_set(minsep2) eq 0) then minsep2 =  rad2mas(1./4/maxr)
if (n_elements(force_min_sep0) gt 0) then force_min_sep=force_min_sep0 $
else force_min_sep=minsep
;;The next lines are a little hacked-up, but works for the 9h mask,
;;where the baseline separation really defines the field-of-view...
s =  sort(r)
;;stop
if (keyword_set(maxsep) eq 0) then maxsep =  rad2mas(1./2./min([r, sqrt(!pi*r[s[5]]^2/12.)]))
if (keyword_set(maxsep2) eq 0) then maxsep2 = rad2mas(1./2./min([r, sqrt(!pi*r[s[5]]^2/12.)]))
if keyword_set(minpa) eq 0 then minpa=0.0
if keyword_set(maxpa) eq 0 then maxpa=360.
if keyword_set(maxpa2) eq 0 then maxpa2=360.0
if keyword_set(minpa2) eq 0 then minpa2=0.0

if (p ne '') then printf,1,minsep, maxsep,minpa,maxpa, $
                         format = '("Min and max sep/pa 1 in search (mas): ", 4F7.1)'
if (p ne '') then printf,1,minsep2, maxsep2,minpa2,maxpa2, $
                         format = '("Min and max sep/pa 2 in search (mas): ", 4F7.1)'





;;now generate grid parameters and do the grid search, might have to
;;be coarser than binary_grid for speed

;Set the minimum angle to be 1 radian of fringe phase at the max
;baseline...
delsep =  rad2mas(1/2./!pi/maxr) ;1 radian. This is kind-of arbitrarily chosen...
nr   =  round((maxsep-minsep)/delsep) +1;+2
delang = (delsep/maxsep)*180./!pi
nang =  round( (maxpa-minpa)/delang )
delang = (maxpa-minpa)/nang           
nang = nang +1
delsep =  rad2mas(1/2./!pi/maxr) ;1 radian. This is kind-of arbitrarily chosen...
nr2   =  round((maxsep2-minsep2)/delsep) +1;+2
delang2 = (delsep/maxsep2)*180./!pi
nang2 =  round( (maxpa2-minpa2)/delang2 )
delang2 = (maxpa2-minpa2)/nang2
nang2 = nang2+1


params =  fltarr(9, nr, nang,nr2,nang2)
params[3,*,*,*,*] = 0.1
params[4,*,*,*,*] = 0.1
params[7,*,*,*,*] = 0.1
;;stop
cn=-1
print, 'Generating parameter Grid'
for i =  0, nr-1 do for j = 0, nang-1 do for k=0,nr2-1 do for l=0,nang2-1 do begin
; if (i le 3) then  params[0, i, j] = minsep+i*delsep/2. $;separation: finer sampling
; else params[0, i, j] = minsep+(i-2)*delsep ;separation
   if i ne cn then print,i+1,' of ',nr
   cn=i
 params[0,i,j,k,l] = minsep+i*delsep ;separation 1
 params[1,i,j,k,l] = minpa + j*delang ;position angle 1
 params[5,i,j,k,l] = minsep2+k*delsep ;separation 2
 params[6,i,j,k,l] = minpa2 + l*delang2 ;position angle 2
endfor
params[2,*,*,*,*] = init_crat ;;contrast ratio1
params[8,*,*,*,*] = init_crat ;;contrast ratio 2

;;stop
;Now search for the best chi^2 over the grid...We need to adjust the
;contrasts at every grid point because with two companions these
;things are degenerate
chi2_arr =  fltarr(nr, nang,nr2,nang2)
chi2_vis = fltarr(nr, nang,nr2,nang2)
cn=-1

parinfo = replicate({ fixed:0, limited:[0,0],limits:[0.D,0]}, 6)
parinfo[0].fixed=1
parinfo[1].fixed=1
parinfo[3].fixed=1
parinfo[4].fixed=1

;stop
print, 'Running grid search'
for i =  0, nr-1 do for j = 0, nang-1 do for k=0,nr2-1 do for l=0,nang2-1 do begin

   if i ne cn then print,i+1,' of ',nr
  cn=i
;;old school version, we're fitting now on limited set with mpfit (just
;;to make use of parinfo.fixed)
;;  modelt3 = triple_t3data(params[*, i, j,k,l],t3data=t3data)
;;  cpresid = mod360(modelt3.t3phi - t3data.t3phi)
;;  chi2_arr[i,j,k,l] =  total(cpresid^2/t3data.t3phierr^2)
  thisp = params[[0,1,2,5,6,8],i,j,k,l]
  ppp = mpfit('triple_oifits_resid',thisp,parinfo=parinfo,bestnorm=bn,/quiet,maxiter=10)
  chi2_arr[i,j,k,l] = bn
  params[2,i,j,k,l ] =ppp[2]
  params[8,i,j,k,l] = ppp[5]
endfor

w =  where(params[0, *, 0,0,0] lt rad2mas(1./2/maxr))
filter_array =  fltarr(nr, nang,nr2,nang2)
filter_array[*] = 1.0
filter_array[w,*,*,*]  = 100.0

w =  where(params[5,0, 0,*,0] lt rad2mas(1./2/maxr))
if w[0] ne -1 then filter_array[w,*,*,*]  = 100.0



;;if (w[0] ne -1) then for i = 0, nang-1 do filter_array[w, i] =  10.0

dummy =  min((chi2_arr*filter_array), m)
in = array_indices(chi2_arr, m)

;;initialise this structure
apriori = {val:double([ 238.86587248546860  ,     146.51658110575497     ,  9.7279238180774907  ,     88.017738437772877    ,   345.0,   21.990320672245204]),err:double([1e3, 1e3, 1e5,1e3,1e3,1e5])}
;;put in the best grid fit
apriori.val= params[[0,1,2,5,6,8],in[0],in[1],in[2],in[3]]
if apriori.val[2] lt 0 then apriori.val[2] *= -1
if apriori.val[5] lt 0 then apriori.val[5] *= -1


bestp = apriori.val
xixi = 0.1*identity(6)
powell,  bestp, xixi ,  1e-4,  bestchi2,  'triple_oifits_chi2',iter=iter,itmax=200
if (bestp[2] lt 0) then bestp[2] *= -1
;;OK - no idea why this is so, but a negative contrast ratio is
;;     equivalent to a positive one.
print,'Best Grid Solution'
print, bestp
print,'Chi2 ', bestchi2
;;stop

modelt3 = triple_t3data([bestp[0],bestp[1],bestp[2], 0.1, 0.1,bestp[3],bestp[4],0.1,bestp[5]],t3data=t3data)
good = where(t3data.t3phierr lt 60)
;Find the closure-phase error to add in quadrature so that chi^2=1
extra_error =  [0, 0.01 * 10.^(alog10(300.)*findgen(29)/20.) ]
newchi2 =  fltarr(30)
newchi2_null = fltarr(30)
cpresid = double(mod360(t3data.t3phi-modelt3.t3phi))
if (cp_cinv[0] eq -1) then begin
 if (proj[0] eq -1) then begin
  ;;stop
  for i =  0, 29 do newchi2[i] =  total(cpresid^2/(t3data.t3phierr^2+extra_error[i]^2)) 
  for i =  0, 29 do newchi2_null[i] =  total(t3data.t3phi^2/(t3data.t3phierr^2+extra_error[i]^2)) 
  newchi2 = newchi2/n_elements(good)*ndf_cp/(ndf_cp-3) ;The -3 is for the three parameters that are fit.
  newchi2_null = newchi2_null/n_elements(good)
 endif else begin
  for i =  0, 29 do newchi2[i] =  total((cpresid#proj)^2/(proj_err^2+extra_error[i]^2)) 
  for i =  0, 29 do newchi2_null[i] =  total((t3data.t3phi#proj)^2/(proj_err^2+extra_error[i]^2)) 
  newchi2 = newchi2/(nproj-3)
  newchi2_null = newchi2_null/nproj
 endelse
 extra_error_null = interpol(extra_error,  newchi2_null,  1)
 extra_error =  interpol(extra_error,newchi2,1)
endif
if (cp_cinv[0] eq -1) then begin
 if (newchi2[0] lt 1.0) then begin
  extra_error =  0.0 
  print,  'New Chi2: ',  newchi2[0]
  if (p ne '') then printf,  1, 'New Chi2: ',  newchi2[0]
 endif
 if (newchi2_null[0] lt 1.0) then begin
  extra_error_null = 0.
  print,  'New Chi2: ',  newchi2_null[0]
  if (p ne '') then printf,  1, 'New Chi2 (null): ',  newchi2_null[0]
 endif
 print,  extra_error,  extra_error_null, $
   format = '("Extra error for binary and single star soln (degs): ", 2F7.2)' 
 if (p ne '') then printf, 1, extra_error,  extra_error_null, $
   format = '("Extra error for binary and single star soln (degs): ", 2F7.2)'
 ;;Add the extra error in. We used to multiply cperr_null and t3data.t3phierr by sqrt(n_elements(good)/ndf)
 ;;This was a problem... because it didn't distinguish between ndf and ndf_cp
 cperr_null = sqrt(t3data.t3phierr^2+extra_error_null^2)*sqrt(n_elements(good)/ndf_cp)
 if (keyword_set(scale_err)) then begin ;;!!! Not sure this works with usevis !!!
     t3data.t3phierr = t3data.t3phierr*sqrt(newchi2[0]*n_elements(good)/ndf_cp) 
     if (usevis gt 0) then vis2data.vis2err = vis2data.vis2err*sqrt(newchi2[0])
 endif else begin
  t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)*sqrt(n_elements(good)/ndf_cp)
 endelse
 proj_err_null = sqrt(proj_err^2+extra_error_null^2)
 proj_err = sqrt(proj_err^2+extra_error^2)
endif else begin
  ;;OK - really what I want to do is add in extra errors to the
  ;;     covariance matrix.
; newchi2 = transpose(cpresid)#cp_cinv#cpresid/(n_elements(good)-3)
; newchi2_null = transpose(t3data.t3phi)#cp_cinv#t3data.t3phi/(n_elements(good)-3)
  print,  'Inverting and calculating extra error...'
  cor = cov2cor(cp_cov,  sig = sig)
  for i =  0, 29 do begin
   temp_cinv =  invert(double(extra_error[i]^2*cor + cp_cov))
   newchi2[i] =  transpose(cpresid)#temp_cinv#cpresid/(n_elements(good)-3)
   newchi2_null[i] = transpose(t3data.t3phi)#temp_cinv#t3data.t3phi/(n_elements(good)-3)
 endfor
 if (newchi2_null[0] lt 1.0) then begin
  print,  'Chi2 for 1-star solution: ',  newchi2_null[0]
  if (p ne '') then printf,  1, 'Chi2 for 1-star solution: ',  newchi2_null[0]
  cp_cov_null = cp_cov
 endif else begin
  extra_error_null = interpol(extra_error,newchi2_null,1) 
  if (newchi2[0] gt 1) then extra_error =  interpol(extra_error,newchi2,1) else extra_error = 0.
  scale = sqrt(newchi2[0])
  print,  extra_error,  extra_error_null, $
    format = '("Extra error for binary and single star soln (degs): ", 2F7.2)' 
  if (p ne '') then printf, 1, extra_error,  extra_error_null, $
    format = '("Extra error for binary and single star soln (degs): ", 2F7.2)'
  cp_cov += extra_error^2*cor
  cp_cov_null += extra_error_null^2*cor
 endelse
endelse
cp_cinv = invert(double(cp_cov))
cp_cinv_null = invert(double(cp_cov_null))

powell,  bestp,  xixi,  1e-4,  fmin,  'triple_oifits_chi2'

chi2_bin =  triple_oifits_chi2(bestp, sig = perror,  corr = cor,  covar = covar, /print)
redchi2 = chi2_bin/ndf

;;now print some things
printf,1,t3data[0].mjd,format = '("MJD: ", F13.6)'
print,  bestp,                        format = '("Final best solution:   ", 6F9.2)'
if (p ne '') then printf, 1, bestp,   format = '("Final best solution:   ", 6F9.2)'
print,  perror,                       format = '("Errors:                ", 6F9.2)'
if (p ne '') then printf, 1, perror,  format = '("Errors:                ", 6F9.2)'
printf,1,'Covariance Matrix:'
printf,1,covar
close,1

end
