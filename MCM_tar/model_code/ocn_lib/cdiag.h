c====================== include file "cdiag.h" =========================
c
c     variables used for computing diagnostics:
c
c     for precision purposes, ektot, dtabs are summed over i for each
c     k & j. tvar, engint, & buoy are summed over i & j for each k.
c     the total integrals are then summed over k and stored into the 
c     zeroth element (k=0)
c
c     energetics are done over the entire ocean domain:
c
c     ektot    = total kinetic energy per unit volume
c     dtabs    = rate of change of tracer per unit volume
c     engint   = internal mode energy integral components per
c                unit volume
c     engext   = external mode energy integral components per
c                unit volume
c     buoy     = buoyancy per unit volume
c
c     term balances are instantaneous breakdowns of all terms in the 
c     momentum & tracer equations. They are averaged over ocean volumes 
c     defined by horizontal and vertical regional masks (see ocean.F):
c     
c
c     termbt   = term balance components for time rate of change of 
c                tracers per unit volume. the total time rate of change
c                is broken down into components as follows:
c                (1) = total time rate of change for the tracer
c                (2) = change due to zonal advection
c                (3) = change due to meridional advection
c                (4) = change due to vertical advection
c                (5) = change due to zonal diffusion
c                (6) = change due to meridional diffusion
c                (7) = change due to vertical diffusion
c                (8) = change due to source term 
#ifdef kdtermb
c                (9) = change due to convective adjustment
c     tfilt    = time rate of change due to filtering and truncation
c                (computed as a residual)
c    stamp1    = timestamp of first sample in term balance time mean
#else
c     tconv    = time rate of change due to convection per unit volume
#endif
c     tvar     = rate of change of tracer variance per unit volume
c     atwv     = average tracer within volume
c     asst     = average sea surface tracer for regional surface areas
c     stflx    = surface tracer flux per unit volume
c
c     termbm   = term balance components for time rate of change of 
c                momentum per unit volume. the total time rate of change
c                is broken down into components as follows:
c                (1) = total time rate of change for the momentum
c                (2) = change due to the pressure gradient
c                (3) = change due to zonal advection
c                (4) = change due to meridional advection
c                (5) = change due to vertical advection
c                (6) = change due to zonal diffusion
c                (7) = change due to meridional diffusion
c                (8) = change due to vertical diffusion
c                (9) = change due to metric terms
c                (10) = change due to coriolis terms
c                (11) = change due to source terms 
c                (12) = change due to surface pressure gradient
c                (13) = change due to metric advection
c     smflx    = surface momentum flux per unit volume
c     avwv     = average velocity within volume
c 
c     zuneng   = vertical average of zonal forcing on row j
c     zuseng   = vertical average of zonal forcing on row j-1
c     zvneng   = vertical average of meridional forcing on row j
c     zvseng   = vertical average of meridional forcing on row j-1
c
c
      common /cdiag/ ektot(0:km,jmt), dtabs(0:km,nt,jmt)
#ifdef kdtermb
     $,              engint(0:km,8), engext(8), tfilt(nt)
#else
     $,              engint(0:km,8), engext(8), tconv(nt)
#endif
     $,              atwv(nt,0:numreg), asst(nt,0:nhreg), avwv(3,numreg)
#ifdef kdtermb
     $,              termbt(0:km,9,nt,0:numreg)
#else
     $,              termbt(0:km,8,nt,0:numreg)
#endif
     $,              termbm(0:km,13,2,numreg)
     $,              tvar(0:km,nt,0:numreg), smflx(2,0:nhreg)
     $,              stflx(nt,0:nhreg), buoy(0:km), zuneng(imt,8) 
     $,              zuseng(imt,8), zvneng(imt,8), zvseng(imt,8)
#ifdef kdtermb
      character * 32 stamp1
      common /cdiagc/ stamp1
#endif
c
c     vbr      = zonal sum of meridional velocity
c     tbrn     = zonal sum of tracer (north row)
c     tbrs     = zonal sum of tracer (south row)
c     ttn      = northward transport of tracer
c         1       2       3       4       5       6      7       8
c       x mean  x eddy  z mean  z eddy   ekman tot adv diffus   total
c
c     tmt      = meridional mass transport
c     ustf     = names & units for surface tracer fluxes
c
#ifdef testnorth
c
c     the ttn, tmt and tbrn arrays are dimensioned such that the 0 index
c     refers to the global totals, 1 represents the Atlantic & Arctic,
c     and 2 refers to the Indo-Pacific.
c
      character * 16 nornam(0:2)
      common /cdiagc/ nornam
cSGI moved to block data
cc    data (nornam(ns),ns=0,2)
cc   $     /'  Global Total  ', ' Atlantic/Arctic', '  Indo/Pacific  '/
c
      common /cdiagm/ normsk(imt,jmt)
      common /cdiaga/ vbr(km), tbrn(km,nt,0:2), tbrs(km,nt)
     $       , ttn(8,jmt,ntmin2,0:2), tmt(jmt,km,0:2)
#else
      common /cdiaga/ vbr(km), tbrn(km,nt),
     $          tbrs(km,nt), ttn(8,jmt,ntmin2), tmt(jmt,km)
#endif
      character*15 ustf(nt,2)
      common /cdiagb/ ustf
c
