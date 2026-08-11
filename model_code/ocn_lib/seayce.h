c====================== include file "seayce.h" ========================
#ifdef seaice
c
c     variables which are associated with the prediction of sea ice
c     (see "ice.F")
c
c     yc     = sea ice thickness in cm for tau-1 & tau time levels.
c     yca    = sea ice thickness in cm for tau+1 time level.
c     ycejam = sea ice jam criterion [cm], above this value advection
c              of sea ice may be set to zero
c     freeze = temperature [deg C] at which water freezes into sea ice
c     altht  = latent heat of freezing for sea water (liquid==>solid)
c              [cal/cm**3]
c     ahyc   = horizontal sea ice diffusion coefficient [cm**2/sec]
c     dtyce  = heat & water contribution to terr & serr from ice changes
c              [same units as terr & serr]
c     stfyc  = heat & water flux into the ocean due to ice changes
c              [cal/cm**2/sec & psu/cm**2/sec]
c
      common /seayce/ yc(imt,jmt,ntau), yca(imt,jmt)
      common /seayce/ ycejam, freeze, altht, ahyc
      common /seayce/ ycknt
      common /seayce/ dtyce(imt,jmt,2)
      common /seayce/ stfyc(imt,jmt,2)
c
#endif
