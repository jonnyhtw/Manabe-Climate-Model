      common / param  / ww    , eradsq, grav, tw   , erad , radin, rasq
      common / param  / absvor, dtx2  , dt  , dtsq , rgas , cp   , cappa
      common / param  / o2g   , tr2   , ara , cpp  , vrt3 , rix  , rqcp
      common / param  / gm    , og    , sbc , tfra , tfro , frlxx4
      common / param  / hmf   , krs1  , krs2, rswt1, rswt2, vk0   
      common / param  / z0land, z0sea , pi  , hpi  , twopi, radian
      common / param  / alphao, pialph, pid24
#ifdef sponge
      common / param  / dmptop
#endif
c
c vk0     = von Karman constant
c z0land  = roughness length (land)
c z0sea   = roughness length (sea)
c pi      = 4.0*atan(1)   (3.14159...)
c hpi     = pi/2
c twopi   = pi*2
c radian  = 360/twopi (57...)
c alphao  = empirical constant used in computing O3 absorption
c           used in LWTRAN (0.28)
c pialph  = pi*alphao used in LWTRAN
c pid24   = pi/24 used in ZENITH
c dmptop  = damping coef. for the sponge layer at the top of the atm.
c           it is used in LINEAR and defined in PARMTR.
