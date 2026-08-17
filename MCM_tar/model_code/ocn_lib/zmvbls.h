
c====================== include file "zmvbls.h" ========================
c
c     variables used for computing zonal mean diagnostics for tracers 
c     and the u, v and w components of velocity 
#ifdef seaice
c     and sea ice
#endif
c
c     zmt = zonal mean tracers computed row by row and level by level
c           for each tracer
c     zmu = zonal mean of u component of velocity (by rows & levels)
c     zmv = zonal mean of v component of velocity (by rows & levels)
c     zmw = zonal mean of w component of velocity (by rows & levels)
#ifdef seaice
c     zmi = zonal mean of sea ice (by rows)
#endif
c     zwt = zonal weighting (e-w distance) for ocean points at a given 
c           level and j-row for the ts grid.  Length in the east-west
c           direction of the ocean point on the ts grid.
c     zwu = zonal weighting (e-w distance) for ocean points at a given
c           level and j-row for the uv grid.  Length in the east-west
c           direction of the ocean point on the uv grid.
c
      common /zmvbls/ zmt(jmt,km,nt), zmu(jmt,km), zmv(jmt,km)
     $,               zmw(jmt,km)
#ifdef seaice
     $,               zmi(jmt)
#endif
c
      common /zmindx/ zwt(jmt,km), zwu(jmt,km)
c
