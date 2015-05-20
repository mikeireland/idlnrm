;######## RELEASE CLIP HERE ######## 
;+
; inquire(paramname,infostr,ix=ix)
;
; inquire will return a default setting (for example, a template) by
;   using the observing log info structure.

; Input:  name         - the name of the thing you want
;         ix           - index for addressing infostr(if needed)
; Returns: infostr     - a structure with everything you wanted to know
; 
; created                                                    PGT 22Oct05
;
; This is really example code to show how to build it at
; present
function inquire_conica,paramname,infostr,ix

camera=infostr.instrument[ix]

case paramname of
  'template': begin
     mask_xtra=''
     if(infostr.mask eq '18Holes') then mask_xtra='-1'
     datasize=fix(infostr.cube_sz[ix,0])
     if(datasize eq 256) then szname='' else szname='_'+strtrim(string(datasize),2)
     ; Tag on a date to allow the rotation to be tuned each run
     date_id=''
     if(strpos(infostr.date[0],'2010-03') ne -1) then date_id='_Mar10' 
     return, 'conica/mf_'+infostr.mask+mask_xtra+'_'+infostr.filter[ix]+'_'+infostr.lyot[ix]+szname+date_id+'.idlvar'
   end 
  'clog': return,  make_clog(infostr)  ; Defaults fine??.
  'n_blocks': return,  10
  'pscale' : if(infostr.lyot eq "S13") then return,13.27 else return,27.15
  'mask':return,  infostr.mask 
 endcase
return,0

end
 

