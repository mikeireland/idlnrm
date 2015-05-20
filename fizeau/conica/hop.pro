; For use with CONICA star-hopping to find jump in RA and Dec. 
; Inputs: RA and Dec
; Note this program is a very quick hack-up for Conica_Jun2010.
; It has (at least) three bugs which can cause it to fail completely:
; (1) A dec of -00 23 45 will be treated as +00 23 45. So within minus 
;       one degree of 00 you must input 00 -23 45 
; (2) The code is not clever enough to jump across the zero RA line, and
;       will fail if jumping from RA 23hours to RA 00hours or back. A
;       quick work around is to just add one hour (jump 00 to 01 hours).
; (3) To work properly for a large number of jumps, the code properly should
;       re-calculate its RA and DEC after each jump ... so that after a number
;       of hops, errors build up. Experience Jun2001: don't do more than ~15 hops.
;Optional: A starlist name. Default is root_dir + 'fizeau/nirc2/'
;
; NEW VERSION: updated according to finding from the March 2011 run:
; You acquire on A.
; A->B you use the declinaison of B to calculate the offset
; B->A you also use the declinaison of B to calculate the offset
; B->C you calculate the offset to go from A to C and you calculate the
; offset to go from A to B. Then you subtract the two.


pro hop,  starlist=starlist, $ 
          ra0h,ra0m,ra0s,dec0d,dec0m,dec0s, $
          ra1h,ra1m,ra1s,dec1d,dec1m,dec1s, $
          ra2h,ra2m,ra2s,dec2d,dec2m,dec2s

if(keyword_set(starlist) eq 0) then starlist='star_catalog'

typ=size(ra0h,/type)

if(typ eq 0) then begin
  print,'You can drive this program two ways:'
  print,'hop, 16, 52, 31.89, -42, 59, 30.45,   16, 52, 31.89, -42, 59, 30.45,    16, 52, 20.14,  -38, 01 ,03.13 '
  print,'...OR...'
  print,'hop, starlist="star_catalog","HD 139614","HD 142666","HD 143006"'
  print,""
  print,'where HD 139614 is the target on which acquisition was performed'
  print,"HD 142666 the target on which we are presently"
  print,"and HD 143006 where we want to go"
  goto,enditall
endif

if(typ eq 7) then begin
  a = READSTARLIST(starlist)
  namelist=strtrim(strupcase(a.targ))
  acqname=ra0h
  fromname=ra0m
  toname=ra0s
  acqix=where(namelist eq strupcase(acqname))
  fromix=where(namelist eq strupcase(fromname))
  toix=where(namelist eq strupcase(toname))

  ; in case we have multiply-entered names in catalog
  acqix=acqix[0]
  fromix=fromix[0]
  toix=toix[0]

  if(acqix eq -1) then begin
    print,'Could not find your star ',acqname
    goto,enditall
  endif
  if(toix eq -1) then begin
    print,'Could not find your star ',toname
    goto,enditall
  endif
  if(toix eq -1) then begin
    print,'Could not find your star ',toname
    goto,enditall
  endif

;  ra1h=(a.rah)[fromix] & ra1m=(a.ram)[fromix] & ra1s=(a.ras)[fromix]  
;  dec1d=(a.ded)[fromix] & dec1m=(a.dem)[fromix] & dec1s=(a.des)[fromix]  
;  ra2h=(a.rah)[toix] & ra2m=(a.ram)[toix] & ra2s=(a.ras)[toix]  
;  dec2d=(a.ded)[toix] & dec2m=(a.dem)[toix] & dec2s=(a.des)[toix]  

 ;;The coordinate conversions are already done correctly in READSTARLIST
  ra0ss  = 3600.*(a.raj2000)[acqix]
  dec0ss = 3600.*(a.dej2000)[acqix]
  ra1ss  = 3600.*(a.raj2000)[fromix]
  dec1ss = 3600.*(a.dej2000)[fromix]
  ra2ss  = 3600.*(a.raj2000)[toix]
  dec2ss = 3600.*(a.dej2000)[toix]
endif else begin
 ;;NB This does NOT work between -00 and -01 degrees Dec.
 ;;Surprisingly, this is *not* fixed with "ten.pro" in astrolib, and
 ;;requires a string search for "-"
 ra0ss=ra0s+ra0m*60.+ra0h*3600.
 dec0ss=dec0s+dec0m*60.+abs(dec0d)*3600.
 if (dec0d lt 0) then dec0ss *= -1

 ra1ss=ra1s+ra1m*60.+ra1h*3600.
 dec1ss=dec1s+dec1m*60.+abs(dec1d)*3600.
 if (dec1d lt 0) then dec1ss *= -1

 ra2ss=ra2s+ra2m*60.+ra2h*3600.
 dec2ss=dec2s+dec2m*60.+abs(dec2d)*3600.
 if (dec2d lt 0) then dec2ss *= -1

 ; convert RAs from 24 hours to 360 degrees
 ra0ss = ra0ss*360./24.
 ra1ss = ra1ss*360./24.
 ra2ss = ra2ss*360./24.

endelse

delra=(ra2ss-ra0ss)*cos(dec2ss*!pi/180/3600.) - (ra1ss-ra0ss)*cos(dec1ss*!pi/180/3600.);
        ;this is in seconds of arc
deldec=dec2ss-dec1ss 

nhop=1;
if (sqrt(deldec^2+delra^2) gt 3599) then nhop=long(sqrt(deldec^2+delra^2)/3500)+1;

print,delra,deldec,format='("delta RA = ",F11.2,"    delta Dec = ",F11.2)'
print,nhop,format='("The number of needed hops is : ",I3)'

nrass=indgen(nhop+1)*(ra2ss-ra1ss)/nhop+ra1ss;
ndecss=indgen(nhop+1)*(dec2ss-dec1ss)/nhop+dec1ss;

for n=0,nhop-1 do begin
delra=(nrass(n+1)-ra0ss)*cos(ndecss(n+1)*!pi/180/3600.) - (nrass(n)-ra0ss)*cos(ndecss(n)*!pi/180/3600.);
        ;this is in seconds of arc
deldec=ndecss(n+1)-ndecss(n);
print,n+1,delra,deldec,format='("Hop # ",I3," -> delta RA = ",F11.2,"    delta Dec = ",F11.2)'
endfor

enditall:
end
