;; This script uses ggv to view images, as image_cont (or tvscl etc.)
;; doesn't seem to work.

pro psim, im

set_plot, 'ps'
device, xsize=10, ysize=10
image_cont, im, /noc
device, /close
spawn, 'display idl.ps &'
set_plot, 'x'

end
