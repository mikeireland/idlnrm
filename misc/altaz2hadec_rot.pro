;This procedure converts Altitude and Azimuth (in degrees) to
;Hour angle and dec (all indegrees). It also converts the delta_alt and delta_az
;unit vectors to vectors in a (HA,dec) coordinate frame by returning a
;matrix and calculates the angle between vertical and north.
;
;It can be used as a GENERIC coordinate transformation routine, by
;setting poleaz to the azimuth of the hadec coordinate system in the
;hadec coordinage system, and by setting lat to the altitude of the
;pole of the hadec coord system.
;
;ha:  hour angle in degrees (West of North, i.e. ANTICLOCKWISE on sky)
;dec: Declination
;alt: altitude
;az: azimuth (East of North, i.e. also ANTICLOCKWISE on sky)
;par: the angle between vertical (altitude) and North, measured
; in an anticlockwise direction. (i.e. parallactic angle)
;rotm: rotm#[alt,az] converts a vector to [N,W]

pro altaz2hadec_rot,  alt,  az,  lat,  ha,  dec,  par,  $
    rotm = rotm,  poleaz = poleaz

if (keyword_set(poleaz) eq 0) then poleaz =  0

az -= poleaz
alt2 = alt
az2  = az
altaz2hadec,  alt2,  az2,  lat,  ha,  dec
;From http://www.petermeadows.com/html/parallactic.html
sinp =  sin(az*!pi/180.)*cos(lat*!pi/180)/cos(dec*!pi/180)
cosp =  (sin(lat*!pi/180) - sin(dec*!pi/180)*sin(alt*!pi/180))/$
  (cos(dec*!pi/180)*cos(alt*!pi/180))
par =  atan(sinp, cosp)

rotm = [[cos(par), sin(par)], $
        [-sin(par), cos(par)]]
par = par*180/!pi

end
