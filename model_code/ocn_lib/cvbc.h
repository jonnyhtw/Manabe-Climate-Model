c====================== include file "cvbc.h" ==========================
c
c     vertical boundary condition variables:
c
c     smf = surface momentum flux
c          1 => zonal wind stress (dynes/cm**2)
c          2 => meridional wind stress (dynes/cm**2)
c     bmf = bottom momentum flux
c          1 => zonal bottom drag (dynes/cm**2)
c          2 => meridional bottom drag (dynes/cm**2)
c     stf = surface tracer flux
c          1 => surface heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c                             ( assuming rho*cp = 1 cal/degC/cm**3)
c          2 => surface water flux (ppt/cm**2/sec)
#if defined cpldleg3 && !defined cpld_restore
c     stfadj = surface tracer flux adjustment for heat and water
c          1 => surface heat flux adjustment, terr
c          2 => surface water flux adjustment, serr
#endif
c     btf = bottom tracer flux (for consistency but normally zero!)
c          1 => bottom heat flux (cal/cm**2/sec = cm*degC/sec = ly/sec)
c                             ( assuming rho*cp = 1 cal/degC/cm**3)
c          2 => bottom water flux (ppt/cm**2/sec)
c
#if defined worldocn && defined seasonal
c     smfavg = annual mean wind stresses in x & y directions; sometimes
c              used to force the external mode when DTTS/DTUV is large
c              (dynes/cm**2)
c
#endif
#ifdef ccl4tracer
c     variables use in calculating "stf" for carbon tetrachloride
c
c     CL4air = atmospheric mixing ratio of CCl4 (parts per trillion vol])
c     CL4sol = solubility of CCl4 in seawater
c                              (moles/kg/atmosphere)
c     CL4tvl = transfer (piston) velocity for CCl4 in seawater (cm/sec)
c     CCL4sc = schmidt number for CCl4 in seawater
c 
#endif
#ifdef cfctracers
c     variables use in calculating "stf" values for cfc-11 and cfc-12
c
c     f11air & f12air = atmospheric mixing ratios of cfc-11 and cfc-12
c                              ( ?? )
c     f11sol & f12sol = solubilites of cfc-11 and cfc-12 in seawater
c                              (moles/liter/atmosphere)
c     f11tvl & f12tvl = transfer (piston) velocities for cfc-11 and
c           cfc-12 in seawater (cm/sec)
c
#endif
ccc# if defined worldocn && defined diskless
ccc#  ifdef seasonal
#if defined worldocn 
c
c----------------------------------------------------------------------
c     in this version of the world ocean model, the 12 monthly fields
c     for each surface boundary condition are read in on the first
c     timestep of each run in subroutine "setvbc".
c     they are dimensioned below and end with the suffix "obs"
c     if memory use is a concern, a more complicated scene to read in
c     the surface boundary conditions one a month and only keeping two
c     months' worth memory should be developed.
c----------------------------------------------------------------------
c
c     include monthly surface boundary conditions
c
# if defined cybervbc || defined cpldleg2 || defined cpld_restore
      common /xtvbc/ sstobs (imt,jmt,12), salobs (imt,jmt,12)
# endif
# ifdef cybervbc
      common /xtvbc/ wsxobs (imt,jmt,12), wsyobs (imt,jmt,12)
# endif
# if defined cybervbc || defined cpldleg2 || defined cpldleg3
      common /xtvbc/ wsxavg (imt,jmt),    wsyavg (imt,jmt)
# endif
# if defined cfctracers || defined ccl4tracer || defined cpldleg2 || defined cpld_restore
      common /xtvbc/ yceobs (imt,jmt,12)
# endif
# ifdef seaice
      common /xtvbc/ sno (imt,jmt,12)
# endif
c
# if defined cpldleg2 || defined cpld_restore
      common /xtvbc/ terr (imt,jmt), serr (imt,jmt)
     $,              terrj (imt),    serrj (imt)
     $,              tobs (3)
# endif
# if defined cpldleg3 && !defined cpld_restore
cSGI
      real*4 terr, serr
cSGI
      common /xtvbc/ terr (imt,jmt,2), serr (imt,jmt,2)
     $,              tterr, tserr
# endif
#endif
c
#ifndef multitasking
# define task
#endif
cSGI
cc    task common /cvbc/ smf(imt,2),  bmf(imt,2)
           common /cvbc/ smf(imt,2),  bmf(imt,2)
     $,                  stf(imt,nt), btf(imt,nt)
#if defined cpldleg3 && !defined cpld_restore
     $,                  stfadj(imt,nt)
#endif
#if defined worldocn && defined seasonal
     $,                  smfavg(imt,2)
#endif
#ifdef seaice
     $,                  snowfx(imt)
#if defined cpldleg2 || defined cpld_restore
     $,                  adjyc(imt)
#endif
#endif
#ifdef ccl4tracer
     $,                  CL4air(imt)
     $,                  CCL4sc(imt)
     $,                  CL4sol(imt)
     $,                  CL4tvl(imt)
#endif
#ifdef cfctracers
     $,                  f11air(imt), f12air(imt)
     $,                  f11sc (imt), f12sc (imt)
     $,                  f11sol(imt), f12sol(imt)
     $,                  f11tvl(imt), f12tvl(imt)
#endif
#if defined ccl4tracer || defined  cfctracers
     $,                  wspd(imt)
#endif
c
