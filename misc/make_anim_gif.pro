;Procedure to make an animated gif from a .fits file.

pro make_anim_gif,  filename,  outfile = outfile, nfiles = nfiles

loadct,  13
array =  float(readfits(filename))
if (keyword_set(nfiles) eq 0) then nfiles = (size(array))[3]
if (keyword_set(outfile) eq 0) then outfile = 'anim.gif'
print,  nfiles
for i = 0, nfiles-1 do begin
 frame =  array[*, *, i]
 s =  sort(-frame)
 frame =  frame/mean(frame[s[0:3]])*240
 frame =  byte(frame < 255)
 write_gif, outfile, frame,  /multiple
endfor
write_gif,  /close
;for i = 0, nfiles-1 do begin
; if (i lt 10) then numstr =  '00'+strtrim(i, 2)
; else if (i lt 100) then numstr =  '0'+strtrim(i, 2)
; else numstr =  strtrim(i, 2)
; write_gif,  
;endfor

end
