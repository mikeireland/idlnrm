;This function is for fitting to oifits data using mpfit in the grid
;search for triple_grid.pro.
;;This is just the closure phase model though, no visibility as yet

;;p is the paramater vector without disk sizes and is the only input

FUNCTION triple_oifits_resid, p,  hessian = hessian,  covar = covar,  corr = corr,  sig = sig, print=print
  ;;stop
common t3block,  t3data,  apriori, vis2data, usevis,  cp_cinv,  proj,  proj_err, force_min_sep
  delpar =  [0.1, 0.1, 0.01*p[2],0.1,0.1,0.01*p[5]]
  pfull =  [p[0],p[1],p[2], 0.1, 0.1,p[3],p[4],0.1,p[5]]
  modelt3 = triple_t3data(double(pfull),t3data=t3data)
  t3resid = double(mod360(t3data.t3phi-modelt3.t3phi))
  return, t3resid/t3data.t3phierr
END
