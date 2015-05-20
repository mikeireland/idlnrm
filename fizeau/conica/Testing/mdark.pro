; 20Feb04
; PGT
;
; USAGE:
; 
;    mdark,start_num,end_num
;

pro mdark,start_num,end_num,flat,bad,cutoff, $
    prefx=prefx,sufx=sufx,datadir=datadir

if (keyword_set(start_num) eq 0) then begin
 print,'USAGE:'
 print,'mdark,start_num,end_num,flat,bad,datadir=datadir'
 return
endif


if (keyword_set(datadir) eq 0) then datadir='./'
if (keyword_set(prefx) eq 0) then prefx='NACO_IMG_CAL049_'
if (keyword_set(sufx) eq 0) then sufx='.fits'

flatnums=[start_num+1,start_num+2]
nflats=end_num-start_num


for i=start_num,end_num do begin
   filename=datadir+prefx+string(i,format="(I4.4)")+sufx
   nflat=float(reform(readfits(filename,head)))
   sz=size(nflat)
   if (i eq start_num) then begin
       flat=fltarr(sz[1],sz[2]) 
       dumbad=flat
   endif



   plothist,nflat,tit='Click on left then right Edge of allowed range'
   cursor,v1,a &wait,.6
   print,' Lower Bound ',v1
   cursor,v2,a & wait,.6
   print,' Upper Bound ',v2
   good=where(nflat gt v1 and nflat lt v2)
   plothist,nflat[good],tit='SECOND CHANCE: Click on left then right Edge of allowed range'
   cursor,v1,a &wait,.6
   print,' Lower Bound 2',v1
   cursor,v2,a & wait,.6
   print,' Upper Bound 2',v2
   good=where(nflat gt v1 and nflat lt v2)
   w=where(nflat lt v1 or nflat gt v2)
   plothist,nflat[good],tit='Final Histogram'

   flat=flat+nflat
   dumbad(w)=dumbad(w)+1
endfor

bad=where(dumbad gt cutoff)

;coadds=sxpar(head,'COADDS')
;flat=flat/float(coadds)
flat=flat/nflats


print,n_elements(bad),' bad pix flagged'

end






