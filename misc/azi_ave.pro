;This is a program to azimuthally average an array

function azi_ave,  a,  av2d = av2d, nelt=nelt, x=x, halfint=halfint

sz =  (size(a))[1]
if keyword_set(halfint) then begin
 ix = findgen(sz)-sz/2.0 + 0.5
 make_2d, ix, ix, xx, yy
 d = sqrt(xx^2 + yy^2)
endif else begin
 d =  shift(dist(sz), sz/2,  sz/2)
else 
nelt =  fltarr(sz)
for i =  long(0),n_elements(d)-1 do nelt[d[i]]++
azi_ave =  fltarr(sz)
for i =  long(0),n_elements(d)-1 do azi_ave[d[i]] += a[i]
azi_ave =  azi_ave/(nelt > 1)
if (arg_present(av2d)) then begin
 av2d =  a
 av2d[*] = azi_ave[d[*]]
endif
x = 0.5 + findgen(sz)

return,  azi_ave

end
