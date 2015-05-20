;;Robust sigma repeatedly gives dumb results, e.g. for a 2-valued distribution.
function robust_sigma_mike, y
ny = n_elements(y)
medy = median(y)
square_devs = (y-medy)^2
s=sort(square_devs)
return, sqrt(square_devs[s[0.68*ny]])
end
