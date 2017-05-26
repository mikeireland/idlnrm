;; This IDL script shows how to process a block of data, assuming that the flats etc are already made.
;; The output files go in the same directory as this script, and key final outputs have the "analysis_string" in
;; them so you can have multiple analyses in one directory if you like.
;; fs is the 
fs  =  [868,876,888,896,902] & analysis_string = 'Mar15'
nfs = -[8,4,8,6,8]
data_dir = '~/data/nirc2/140813/'
flat     = '~/tel/nirc2/analysis/121202/flat_CH4S_512x512.idlvar'
sky_fnums = 437+indgen(10) ;;File names for dark or sky files.

;-------- Automatic (mostly) from here ----- 
cinfofile = 'cubeinfo'+analysis_string+'.idlvar'

;; ** Comment out the next 3 lines for re-calibration, after an initial analysis.
;; Special options on the second line include:
;; median_sub: good for faint data where the background changes.
;; ddir_sky: data directory to find darks/skies in, in case e.g. you are analysing data
;;    before having taken darks at end of night!
;; trust_object_name: Use the OBJECT rather than TARGNAME to decide on objects. This
;;    only applies to e.g. when dithering and manually changing the object from the command
;;    line of waikoko.
qbe_nirc2, data_dir, analysis_string, flat, dark_fnums, fstart = fs, nfs=nfs, $
    /median_sub, ddir_sky='~/data/nirc2/140812/', /trust_object_name
calc_bispect,  cinfofile,    /reset

;---- Run from here only to re-calibrate, e.g. with a block of data ----

;; Special options include mf_file=mf_g9_kp_pk.idlvar, which seems to calibrate a little 
;;    better than the default matched filter file.
calibrate_v2_cp,  cinfofile,  cal4src=-1, /reset, /v2div, /see_all
restore,  cinfofile

;** Example uses of binary_grid
binary_grid,  clog.primary_oifits,  printfile = clog.outname+'binary.txt',  nfr = n_elements(olog.frames), $
  maxsep = 240, minsep = 10, init_crat=100, nsim=1e3
binary_grid,  clog.primary_oifits,  printfile = clog.outname+'binary_vis.txt',  nfr = n_elements(olog.frames), $
  maxsep = 160, minsep = 7, init_crat=2, /nosim, /usevis

end
