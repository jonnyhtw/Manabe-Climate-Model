      subroutine intcnp
c
c   called from: atmos
c   subroutines called: none
c
c   subroutine intcnp saves off the initial values of the mean
c   of psln and psm in double precision.  It also saves dtx2
c   (time step length * 2) in double precision.
c
c   This is done to improve accuracy of the computation of the
c   robert filter equation in linear.
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
c
      complex vor   , div  , t    , psln, ucos , vcos
      complex vorm  , divm , tm   , psm , ucosm, vcosm
      complex vordt , divdt, tdt  , psdt
      complex psdx  , psdy , zs
c
#include "dpsln.h"
c
#include "nxstps.h"
c
#include "param.h"
c
#include "sfharm.h"
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      apsln = real(psln(1,1,1))
      apsm  = real(psm (1,1,1))
c
      adtx2 = dtx2
c
crjs
c     write(6,400) apsln, apsm, adtx2
c 400 format(' in intcnp, apsln, apsm, adtx2 = ',/,3e32.16)
crjs
      return
      end
