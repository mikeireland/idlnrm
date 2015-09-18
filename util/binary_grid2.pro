;;This program is designed to be a cut-down and simpler version of
;;binary_grid, focussing on diagnostics rather than fancy statistics.
;;Keywords:
;;  infile: The input oifits file. We don't input the
;;    closure-phase covariance matrix to this program - only the data
;;    that is part of the oifits format.
;;  init_crat: The contrast ratio at which the grid search is done. NB
;;    in the high contrast regime, you'll get the same best
;;    solution no matter what init_crat is chosen.
;;  boxsize: The size of the grid box in milli-arcsec
;;  like: An output likelihood map.
;;  x: An output vector giving the axes of the likelihood map in mas.
;;  minsep: The resolution of the grid.
pro binary_grid2,  infile,  init_crat=init_crat, boxsize=boxsize,like=like,x=x, extra_error=extra_error, chi2_arr=chi2_arr, $
  minsep=minsep, plot_chi2=plot_chi2, ps=ps, usevis=usevis, crat_upperlim = crat_upperlim, nsig=nsig, scale_to_null=scale_to_null, $
  plotit=plotit, detection=detection, no_chi2_scaling=no_chi2_scaling, fits_out=fits_out, force_detection=force_detection, fix_crat=fix_crat

if keyword_set(scale_to_null) and keyword_set(no_chi2_scaling) then begin
 print, "ERROR: You can't both scale_to_null and have no_chi2_scaling!"
 stop
endif

if not keyword_set(nsig) then nsig=6.0
if (keyword_set(init_crat) eq 0) then init_crat = 250.
extract_t3data,  file = infile,  t3data
extract_vis2data,  file = infile,  vis2data

;;Deal with "bad" triangles in a simplistic way.
w = where(t3data.flag)
if (w[0] ne -1) then t3data[w].t3phierr = 1e3
 
 if keyword_set(extra_error) then t3data.t3phierr = sqrt(t3data.t3phierr^2 + extra_error^2)
 extract_vis2data,  file = infile,  vis2data
 read_oidata,infile,oiarray,oitarget
 target_name=oitarget.target
 target_name=target_name[0]
 ;Enforce a maximum SNR of 4.0 on visibility data.
 vis2data.vis2err = sqrt(vis2data.vis2err^2 + (0.2*vis2data.vis2data)^2)

;;Find the number of rotations on the sky (i.e. independent cubes) by
;;finding how many times the first baseline in the file is repeated,
;;where we define "baseline" as containing uniqu stations indices.
;;WARNING: breaks for a >100 hole mask.
ix = vis2data.sta_index[0] + 100*vis2data.sta_index[1]
nrotations = n_elements(where(ix eq ix[0]))
nv2 = n_elements(vis2data)/nrotations
nt3 = n_elements(t3data)/nrotations


r =  sqrt(vis2data.u^2 + vis2data.v^2)
maxr =  max(r)
if (keyword_set(minsep) eq 0) then minsep =  rad2mas(1./8/maxr)
nsep=boxsize/minsep

;;Make reduced chi^2 equal to 1 by scaling errors, and make sure that
;;our final chi^2 variable has a minimum of ndf (crudely takes into
;;account linear dependence of closure-phases)
;Number of independent degrees of freedom in the closure-phase data. The complicated
;expression at the end is (n_holes - 1)
ndf_ind =  nrotations*(nv2 - (3*nt3/nv2+1) )

;print, 'Number of independent degrees of freedom: ', ndf_ind

if (not keyword_set(use_vis)) then nv2=0
;;The following just for chi^2 scaling.
ndf = nrotations*(nt3 + nv2)

;print, 'Number of degrees of freedom: ', ndf

chi2_null  = total(t3data.t3phi^2/t3data.t3phierr^2)

;Now search for the best chi^2 over the grid, and make a map of
;5-sigma upper limits
chi2_arr =  fltarr(nsep, nsep)
crat_expectation = fltarr(nsep,nsep)
crat_sigma = fltarr(nsep,nsep) ;;NB **not** taking into account
x=(findgen(nsep)-nsep/2+0.5)*minsep
make_2d, x,x,xx, yy
rho = sqrt(xx^2+yy^2)
theta = atan(-xx,yy)*180/!pi
for i =  0, nsep-1 do for j = 0, nsep-1 do begin
  params = [rho[i,j], theta[i,j], init_crat,0.1,0.1]
  modelt3 = binary_t3data(params,t3data=t3data)
  ;;This has to be a weighted
  crat_expectation[i,j] = total(modelt3.t3phi*t3data.t3phi/t3data.t3phierr^2)/$
                          total(modelt3.t3phi^2/t3data.t3phierr^2)/init_crat
  crat_sigma[i,j] = sqrt( total(modelt3.t3phi^2*t3data.t3phierr^2/t3data.t3phierr^4)/$
  						  total(modelt3.t3phi^2/t3data.t3phierr^2)^2 ) / init_crat
  ;;If we're fixing contrast, then the model *is* the model that came straight from binary_t3data.
  ;;If we're not fixing contrast, then scale the model closure phase by the fitted contrast.
  ;;This is the linear approximation for closure-phase (only valid for moderate contrasts, e.g.
  ;;10:1 or more.
  if not keyword_set(fix_crat) then $
     modelt3.t3phi = modelt3.t3phi*crat_expectation[i,j]*init_crat
  cpresid = mod360(modelt3.t3phi - t3data.t3phi)
  chi2_arr[i, j] =  total(cpresid^2/t3data.t3phierr^2) 
endfor

;;Here we scale the error (sigma) values so that reduced chi-squared is 1.0
;;for either the null model (i.e. single star) or the best fit.
if keyword_set(scale_to_null) then begin
 crat_sigma *= sqrt(chi2_null/ndf_ind)
endif else if not keyword_set(no_chi2_scaling) then begin
 crat_sigma *= sqrt(min(chi2_arr)/ndf_ind)
endif 

crat_upperlim = (crat_expectation + nsig*crat_sigma)<1.0

det_sigma = max(crat_expectation/crat_sigma*(crat_expectation gt 0), ix)
print, infile, " Signficance: ", det_sigma

if (det_sigma gt nsig or keyword_set(force_detection)) then begin
 detection=[rho[ix],  theta[ix], crat_expectation[ix]]
endif else begin
 detection=[-1]
endelse

;ndf =  n_elements(vis2data) -(n_elements(t3data)/float(n_elements(vis2data))*3.+1)
;;The likelihood map. !!! May have errors !!!
if keyword_set(plot_chi2) then begin
  like=chi2_arr/ndf
endif else begin
  chi2_arr /= min(chi2_arr)
  chi2_arr *= ndf
  like = -exp(-(chi2_arr-ndf)/2)
endelse
y=x

if (keyword_set(ps)) then begin
 set_plot, 'ps'
 device, /enc, xsize=11, ysize=11
 !p.thick=2
 !x.thick=2
 !y.thick=2
endif

if (keyword_set(plotit)) then begin
!p.position=[0.18,0.1,0.95,0.8]
 image_cont_deluxe, like, /noc, xv=-x, yv=y, ytit='Dec (mas)', xtit='RA (mas)'
 oplot, [0], [0], psym=2, col=255, symsize=1.5
 
cgColorBar, range=[min(like), max(like)], position=[0.18,0.9,0.95,0.98], charpercent=0.6
;im = image(like, position=[0.1,0.3,0.95,0.95])
;xaxis = AXIS('X', tickvalues=interpol(indgen(nsep), x, [-100,-50,0,50,100]), tickname=['-100','-50','0','50','100']) 
;c = COLORBAR(TARGET=im, $
; POSITION=[0.15,0.1,0.95,0.2], $
; TITLE='Reduced Chi-Square')
endif

if (keyword_set(ps)) then begin
 device, /close
 set_plot, 'x'
endif
 ;ic_jdm, chi2_arr, /noc
 ;stop


if (keyword_set(fits_out)) then begin

if keyword_set(scale_to_null) then begin
	scaled_chi2 = chi2_arr / (chi2_null/ndf_ind)
endif else if not keyword_set(no_chi2_scaling) then begin
	scaled_chi2 = chi2_arr / (min(chi2_arr)/ndf_ind)
endif else scaled_chi2 = chi2_arr

;; Please someone (e.g. Dary) figure out how to put useful header information into
;; this, e.g. comments about which extension is which.
writefits, fits_out, crat_expectation
writefits, fits_out, crat_sigma, /append
writefits, fits_out, scaled_chi2, /append

endif

end