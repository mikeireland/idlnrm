pro make_circle,im,x,y,r,sz,pvect=pvect
; Program makes a circle at given location in image.
; Written by John sometime
; PGT - added option to pass out pixel vector PVECT

if (keyword_set(im) eq 0) then begin
  print,'pro make_circle,im,x,y,r,sz'
  return
endif

if (keyword_set(sz) eq 0) then sz=1.0
info=size(im)

xxx=findgen(info(1))
proj=replicate(1,info(1))
xx=xxx#proj
yy=proj#xxx
xx=xx-x
yy=yy-y

IMage=sqrt(xx*xx+yy*yy)
pvect=where(image le r)
image(pvect)=1.0*sz
image(where(image gt r))=0.0
IM=IM+image

;for i=0,info(1)-1 do begin
;for j=0,info(2)-1 do begin
;   r_d=sqrt((float(i)-float(x))^2+(float(j)-float(y))^2)
;   if (r_d le r) then im(i,j)=im(i,j)+size
;endfor
;endfor

end

