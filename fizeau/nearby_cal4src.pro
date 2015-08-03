;;script to make each target only get calibrated by nearby calibrators
;;in time. 

;;Made by ACR 01/2015
;;Updated many times since

;; Inputs:
;; cal4src: the cal4src variable calibrate_v2_cp.pro outputs
;; npos   : the number of nearby calibrator observations on each side of the
;; target observations that you want to use
;; information=information: outputs info on which object is a
;; calibrator (not tested really)
;;
;; Outputs: A new cal4src variable with the new configuration, has
;; same shape as cal4src
;;
;;USE:
;;qbe_nirc2/qbe_conica
;;cal_bispec
;;calibrate_v2_cp with all cals selected top build the basic cal4src array
;;new_cal4src = nearby_cal4src(cal4src,npos)
;;calibrate_v2_cp,.......,cal4src = new_cal4src
;;binary_grid......


function nearby_cal4src,cal4src,npos,information=information
information = strarr(n_elements(cal4src[0,*]))
;;keep for saving purposes
ocal4src = cal4src
for i=0,n_elements(cal4src[0,*])-1 do begin
   ;;Figure out where the targets are, they should be the only part of
   ;;the arrays with something other than zero
   if max(cal4src[*,i]) gt 0 then begin
      information[i] = 'src'
      cals = where(cal4src[*,i] eq 1)
      adist = cals-i
      spotu = where(adist gt 0)
      spotl = where(adist lt 0)
      ucals = cals[spotu]
      lcals = cals[spotl]
      ;;reset cal4src
      cal4src[*,i] = 0
     ; if i eq 9 then stop
      if spotu[0] gt -1 then begin
         sortu = sort(abs(adist[spotu]))
         nearestu = ucals[sortu[0:npos-1]]
         cal4src[nearestu,i] = 1
      endif
      if spotl[0] gt -1 then begin
         sortl = sort(abs(adist[spotl]))
         nearestl = lcals[sortl[0:npos-1]]
         cal4src[nearestl,i] = 1
      endif
     ;if i eq 9 then stop
   endif else information[i] = 'cal'
endfor

;stop
return, cal4src
end
