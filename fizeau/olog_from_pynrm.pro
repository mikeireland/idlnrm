;; This procedure creates a cubeinfo file including an olog from a directory
;; or a series of pynrm cubes. 
;; e.g. 
;; olog_from_pynrm, mask='g9'
;; calc_bispect, 'cubeinfo_pynrm.idlvar', mf_file='nirc2/mf_g9_kp_128.idlvar'
;;
;; TODO: Check maskname extraction from pynrm, and enable non-256x256 mf_file inputs in
;; inquire_nirc2
;;
;; Inputs:
;; analysis_dire: The directory for the cubeinfo file. By convention, the cubeinfo file goes in the *same*
;; 			directory as the cubes unless this is specified.
;; tsize: The target angular diameters (positive for calibrators, negative for targets).

pro olog_from_pynrm, mask=mask, cubedir=cubedir, files=files, analysis_dir=analysis_dir, tsize=tsize, analysis_string=analysis_string

;;Sort out the filename
if not keyword_set(analysis_string) then analysis_string='_pynrm'

;;Sort out directories. 
if not keyword_set(cubedir) then cd, current=cubedir
if not keyword_set(analysis_dir) then analysis_dir=cubedir
cd, cubedir, current=initial_dir

;;If no files input, assume that all fits files beginning with "cube" are 
;;cubes to use.
if not keyword_set(files) then spawn, 'ls cube*fits', files


n_cubes = n_elements(files)
f=fltarr(n_cubes,2)             ; this holds frame statistic and err
fl= fltarr(n_cubes,10,2)          ; average and errors.  src#, aperture size, 0=mean/1=error
framestats={xpk:f,ypk:f,totflx:f,pkflx:f,skybgr:f,phot:fl}   ; ,good_frames:ina}

;;Go through the cubes and get first-pass information from the header:
;; frame numbers of all targets.
;; frame numbers of skies.
frames_ix = intarr(2, n_cubes)
frames = []
for i=0,n_cubes-1 do begin
	h = headfits(cubedir + '/' + files[i])
	history = sxpar( h, 'HISTORY' )
	frames_ix[0,i] = n_elements(frames)
	for j=0, n_elements(history)-1 do begin
		if (strpos(history[j], 'Input') lt 0) then continue
		fit_pos = strpos(history[j], '.fit')
		if (fit_pos lt 0) then continue
		;; !!! Warning - Dodgy line !!! (might not always work)
		frames = [frames, fix(strmid(history[j], fit_pos-4, 4))] 
	endfor
	frames_ix[1,i] = n_elements(frames)-1
endfor
;; Initially, all objects are targets.
if not keyword_set(tsize) then tsize = -0.1*findgen(n_cubes)
;; Make the olog and fill in everything we can. As the python code has separate darks/skies
;; for each cube, just set this to -2
olog = make_olog(n_cubes,tsize, frames, -2)
olog.frames_ix = frames_ix
olog.cubedir = cubedir +'/'
;; Fill in the most important parts of olog...
for i=0,n_cubes-1 do begin
	h = headfits(cubedir + '/' + files[i])
	olog.cube_fname[i,0] = files[i]
	fit_pos = strpos(files[i], '.fit')
	if (fit_pos lt 0) then stop
	;; !!! Warning - Dodgy line !!! (might not always work - assumed 3 digit cube extension
	olog.cube_fname[i,1] = strmid(files[i], fit_pos-3, 3)
	olog.cube_sz[i,*] = [sxpar(h,'NAXIS1'), sxpar(h,'NAXIS2'), sxpar(h,'NAXIS3')]
	;;!!! Dodgy - this should use freud. Currently only working for NIRC2
	olog.instrument[i] = 'NIRC2'
	olog.filter[i] = strtrim(sxpar(h, 'FILTER'),2)
	olog.jd[i] = sxpar(h,'MJD-OBS')
	olog.ra[i] = sxpar(h,'RA')
	olog.dec[i] = sxpar(h,'DEC')
	;;A special target name from the pynrm
	olog.source_name[i] = sxpar(h, 'TARGNAME')
	olog.nax1[i] = sxpar(h,'SZX')
	olog.nax2[i] = sxpar(h,'SZY')
	stats = mrdfits(cubedir + '/' + files[i], 1)
	framestats.xpk[i,*]   = [mean(stats[*].xpeak),stdev(stats[*].xpeak)]
	framestats.ypk[i,*]   = [mean(stats[*].ypeak),stdev(stats[*].ypeak)]
	framestats.pkflx[i,*] = [mean(stats[*].max), stdev(stats[*].max)]
	;;Fill in the pa. The strange definition of del_pa assumes the files
	;;are in time order, but is consistent with fizeau/documentation.
	olog.pa[i] = mean(stats[*].pa)
	olog.del_pa[i] = stats[0].pa - stats[n_elements(stats)-1].pa
	;;Fill in the mask unless over-ridden
	if keyword_set(mask) then olog.mask = mask $
	else olog.mask = sxpar(h, 'PMASK')
endfor

save,olog,framestats,file='cubeinfo'+analysis_string+'.idlvar'

cd, initial_dir

end
