;This procedure converts HA and Dec (in degrees) to
;alitutude and azimuth. It also converts the delta_ha and delta_dec
;unit vectors to vectors in a (alt,az) coordinate frame by returning a
;matrix and calculates the angle between vertical and north.
;
;It can be used as a GENERIC coordinate transformation routine, by
;setting poleaz to the azimuth of the hadec coordinate system in the
;hadec coordinage system, and by setting lat to the altitude of the
;pole of the hadec coord system.
;
;ha:  hour angle in degrees (West of North, i.e. ANTICLOCKWISE)
;dec: Declination
;alt: altitude
;az: azimuth (East of North, i.e. also ANTICLOCKWISE)
;par: the angle between vertical (altitude) and North, measured
; in an anticlockwise direction.
;rotm: rotm#[N,W] converts a vector to [alt,az]

pro hadec2altaz_rot,  ha,  dec,  lat,  alt,  az,  par,  $
    rotm = rotm,  poleaz = poleaz

if (keyword_set(poleaz) eq 0) then poleaz =  0

hadec2altaz,  ha,  dec,  lat,  alt,  az
;From http://www.petermeadows.com/html/parallactic.html
sinp =  sin(az*!pi/180.)*cos(lat*!pi/180)/cos(dec*!pi/180)
cosp =  (sin(lat*!pi/180) - sin(dec*!pi/180)*sin(alt*!pi/180))/$
  (cos(dec*!pi/180)*cos(alt*!pi/180))
;print,  sqrt(sinp^2+cosp^2)
par =  atan(sinp, cosp)

rotm = [[cos(par), -sin(par)], $
        [sin(par), cos(par)]]
par = par*180/!pi
az += poleaz

end
