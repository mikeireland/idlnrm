;This function will return the ordinary ray refractive index
;for mgf2. From the Javascript code at: http://www.korth.de/calc/mgf2.htm

function nmgf2,  L,  n_e = n_e
 a1o =  0.487551072; 
 a2o =  0.39875034;
 a3o =  2.312035347;
 L1o =  0.04338408;
 L2o =  0.094614424;
 L3o = 23.79360419;
 a1e =  0.413440227; 
 a2e =  0.504974988;
 a3e =  2.490486218;
 L1e =  0.03684262;
 L2e =  0.09076162;
 L3e = 23.77199519;
 no = sqrt(1+a1o*L*L/((L*L)-(L1o*L1o))+a2o*L*L/((L*L)-(L2o*L2o))+a3o*L*L/((L*L)-(L3o*L3o)))
 n_e = sqrt(1+a1e*L*L/((L*L)-(L1e*L1e))+a2e*L*L/((L*L)-(L2e*L2e))+a3e*L*L/((L*L)-(L3e*L3e)))

return, no

end

