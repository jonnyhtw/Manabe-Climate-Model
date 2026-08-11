c====================== include file "param.h" =========================
c
c     main parameter file which sets ocean characteristics:
c
c     imt    = number of grid points in longitudinal direction
c     jmt    = number of grid points in latitudinal direction
c     km     = number of grid points in the vertical direction
c     nt     = number of tracers
c     lseg   = maximum number of longitudinal stream function segments
c     nisle  = number of islands
c     nslabs = number of j rows (slabs) in the memory slab window
c     ntau   = number of time levels in the memory slab window
c     numjpt = number of pointers for internal mode velocities
c              (also used for del**2 fields with biharmonic option)
c     istart = starting longitude for calculations
c              (if "skipland" is set, "istart" is calculated for each
c               longitudinal strip of ocean ... skipping land points)
c     iend   = ending longitude for calculations
c              (if "skipland" is set, "iend" is calculated for each
c               longitudinal strip of ocean ... skipping land points)
c     maxipp = max number of island perimeter points
c     nvar   = number of prognostic variables
c     nvarbh = number of prognostic variables using biharmonic mixing 
c     ntlev  = number of time levels on disk
c     nxres  = number of sub domains of either constant or smoothly
c              varying grid box width (see "coord.h")
c     nyres  = number of sub domains of either constant or smoothly
c              varying grid box height (see "coord.h")
c     nhreg = number of regions in the horizontal used for averaging 
c             tracers.
c     nvreg = number of regions in the vertical used for term balance 
c             calculations. note "nvreg" is not used for tracer averages
c     numreg= total number of regions ( = product of nhreg & nvreg) 
c             used for term balance calculations
c
#ifdef worldocn
c     nptrac= total number of passive tracers added to the model
c     numage, numctc, numcfc = number of idealized age tracers, carbon
c             tetrachloride and cfc tracers carried in the model
c     nttemp= tracer number for potential temperature
c     ntsalt= tracer number for salinity
c     ntage1= tracer number for the 1st idealized age tracer
c     ntccl4= tracer number for carbon tetrachloride
c     ntcfc1= tracer number for cfc-11
c     ntcfc2= tracer number for cfc-12
c
c     parameters for GFDL world ocean model runs
c
# ifdef ocnr15
c     parameters to match ocean with rhomboidal 15 spectral model
c
#ifdef islands
      parameter (imt=98, jmt=42, lseg=7, nisle=1
#else
      parameter (imt=98, jmt=42, lseg=7, nisle=0
#endif
# endif
# ifdef ocnr30
c     parameters to match ocean with rhomboidal 30 spectral model
c
#ifdef islands
      parameter (imt=194, jmt=82, lseg=7, nisle=7
#else
      parameter (imt=194, jmt=82, lseg=7, nisle=0
#endif
# endif
     $,          ntau=2, nxres=1, nyres=1
# ifdef ocnr15
     $,          maxipp=125
# endif
# ifdef ocnr30
     $,          maxipp=600
# endif
     $,          nttemp=1, ntsalt=2
#  ifdef agetracers
     $,          numage=1, ntage1=3
#  else
     $,          numage=0
#  endif
#  ifdef ccl4tracer
     $,          numctc=1, ntccl4=(2+numage+1)
#  else
     $,          numctc=0
#  endif
#  ifdef cfctracers
     $,          numcfc=2, ntcfc1=(2+numage+numctc+1)
     $                   , ntcfc2=(2+numage+numctc+2)
#  else
     $,          numcfc=0
#  endif
     $,          nptrac=(numage+numctc+numcfc), nt = 2+nptrac
# ifdef km12
     $,          km=12
# endif
# ifdef km18
     $,          km=18
# endif
# ifdef km25
     $,          km=25
# endif
#else 
c     parameters for sample MOM run
c
      parameter (imt=92, jmt=60, km=15, nt=2, lseg=5, nisle=2
     $,          ntau=2, nxres=1, nyres=1, maxipp=500
#endif
#ifdef biharmonic
     $,          nslabs=4, numjpt=3
#else
     $,          nslabs=3, numjpt=2
#endif
#if defined worldocn && defined readrmsk 
c     parameters used when reading in GFDL world ocean regional masks
# ifdef ocnr15
     $,          nhreg = 7, nvreg = 2, numreg = nhreg*nvreg
# endif
# ifdef ocnr30
     $,          nhreg = 25, nvreg = 3, numreg = nhreg*nvreg
# endif
#else
     $,          nhreg = 5, nvreg = 2, numreg = nhreg*nvreg
#endif
     $,          nvar=nt+2
#ifdef tcvmix
     $                    +4
# ifdef leq
     $                      +1
# endif
#endif
#ifdef biharmonic
     $,          nvarbh=nt+2
# ifdef tcvmix
     $                    +1
#  ifdef leq
     $                      +1
#  endif
# endif
#endif
#if defined diskless
# ifdef multitasking
     $,          ntlev=3
# else
     $,          ntlev=2
# endif
#else
# if defined twoslabdisks && !defined multitasking
     $,          ntlev=2
# else
     $,          ntlev=3
# endif
#endif
#ifdef cyclic
     $,          imu=imt)
#else
     $,          imu=imt-1)
#endif
      parameter  (imtp1=imt+1, imtm1=imt-1, imtm2=imt-2, imum1=imu-1
     $,          imum2=imu-2, jmtp1=jmt+1, jmtm1=jmt-1, jmtm2=jmt-2
#ifndef skipland
     $,          istart=2, iend=imtm1
#endif
#ifdef symmetry
     $,          jscan=jmtm2+1
#else
     $,          jscan=jmtm2
#endif
     $,          kmp1=km+1, kmp2=km+2, kmm1=km-1
     $,          imtkm=imt*km, nwds=imt*jmt
     $,          nkflds=6
#ifdef seaice
     $                 +ntlev
#endif
     $,          nslab=imt*nvar*km, ntmin2=nt+1/nt)
c
c     add parameter constants
c
#include "pconst.h"
c
