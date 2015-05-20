;This will makd an mpeg movie using standard settings from a .fits file.

pro make_mpeg,  filename,  scale = scale

array =  float(readfits(filename))^0.5
array =  reverse(array, 2)
if (keyword_set(scale) eq 0) then scale =  2
sz =  size(array)
;newarray =  fltarr(sz[1]*scale,  sz[2]*scale,  sz[3])
;mid =  mpeg_open(dim*scale)
peak =  max(array[*, *, sz[3]/2:*])
array =  array/peak*200
newarray =  rebin(array, sz[1]*scale,  sz[2]*scale,  sz[3])
;for i =  0, (size(array))[3]-1 do newarray[*, *, i]=  rebin(, sz[1]*scale,  sz[2]*scale)
 ;mpeg_put,  mid,  /color, image = im,  frame = i
;mpeg_save,  mid;,  filename = 'output.mpg'
;mpeg_close,  mid
make_MPEG_movie, newarray;, /Color, Table=5

end
