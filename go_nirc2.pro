;; This IDL script shows how to process a block of data, assuming that the flats etc are already made.
fs  =  [868,876,888,896,902] & analysis_string = 'Mar15'
nfs = -[8,4,8,6,8]

;--------
cinfofile = 'cubeinfo'+analysis_string+'.idlvar'

;** Comment out the next 3 lines for re-calibration, after an initial analysis.
qbe_nirc2, '~/data/nirc2/140813/', analysis_string, '~/tel/nirc2/analysis/121202/flat_CH4S_512x512.idlvar',$
 437+indgen(10),fstart = fs, nfs=nfs, /median_sub, ddir_sky='~/data/nirc2/140812/', /trust_object_name
calc_bispect,  cinfofile,    /reset

calibrate_v2_cp,  cinfofile,  cal4src=-1, /reset, /v2div, /see_all
restore,  cinfofile

;** Example uses of binary_grid
binary_grid,  clog.primary_oifits,  printfile = clog.outname+'binary.txt',  nfr = n_elements(olog.frames), $
  maxsep = 240, minsep = 10, init_crat=100, nsim=1e3
binary_grid,  clog.primary_oifits,  printfile = clog.outname+'binary_vis.txt',  nfr = n_elements(olog.frames), $
  maxsep = 160, minsep = 7, init_crat=2, /nosim, /usevis

end
