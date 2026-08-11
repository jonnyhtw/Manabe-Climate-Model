       subroutine lismer (windsp, schnum, seamsk, tv)
      implicit real*8 (a-h,o-z)
#if defined ccl4tracer || defined cfctracers
c
c=======================================================================
c
c  this subroutine uses the equations of liss & merlivat to
c  calculate transfer velocities as a function of wind speed
c  and schmidt number.
c
c  reference:
c   P. S. Liss and L. Merlivat in The Role of Air-Sea Exchange in
c   Geochemical Cycling, P. Buat-Menard, ed., pp. 113-127, 
c   (Reidel, Boston, 1986).
c
c  inputs:
c    "windsp" = wind speed [m/sec]
c    "schnum" = Schmidt number [dimensionless]
c    "seamsk" = land/sea mask for t,s j-row [0.0=land, 1.0=sea]
c
c  output:
c    "tv" = transfer velocity [cm/sec]
c  local variables:
c    "hrpsec" = conversion factor for cm/hr to cm/sec [hours per second]
c    "scp"    = Schmidt number divided by 600 (sc' in L & M reference)
c    "tvmin"  = arbitrarily minimum transfer velocity [cm/sec]
c=======================================================================
c
#include "param.h"
c
      dimension windsp(imt), schnum(imt), seamsk(imt), tv(imt)
c
c  set arbitrarily minimum transfer velocity "tvmin" in data
c  statement below.  (note units are [cm/sec])
c
      data tvmin /0.0/
c
c=======================================================================
c
      r600 = c1/600.0
      hrpsec = c1/3600.0
      cm2r3 = -2.0/3.0
c
c  calculate transfer velocities using liss & merlivat relations.
c  allow for arbitrarily set minimum transfer velocity
c  and set transfer velocities to zero over land.
c
      do 100 i=1,imt
        scp = schnum(i)*r600
        if (seamsk(i) .eq. c0) then
          tv(i) = c0
        elseif (windsp(i) .le. 3.60) then
          tv(i) = hrpsec*0.17*windsp(i)*(scp**(cm2r3))
        elseif (windsp(i) .gt. 3.60 .and.  windsp(i) .le. 13.00) then
          tv(i) = hrpsec*(2.85*windsp(i) - 9.65)*(scp**(-p5))
        else
          tv(i) = hrpsec*(5.90*windsp(i) - 49.3)*(scp**(-p5))
        endif
        if (tv(i) .lt. tvmin) tv(i) = tvmin*seamsk(i)
100   continue
c
#endif
      return
      end
