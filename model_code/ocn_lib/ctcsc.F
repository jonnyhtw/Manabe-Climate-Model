      subroutine CTCsc (temp, seamsk, scctc)
      implicit real*8 (a-h,o-z)
#ifdef ccl4tracer
c
c=======================================================================
c
c  this subroutine returns estimates of schmidt numbers.
c  Data for schmidt numbers were obtained from Rik Wanninkof
c  at AOML using the Hayduk and Laudie equation with a 6% correction
c  for salinity.
c
c     schmidt num. = kinematic viscocity / molecular diffusivity
c             k.v. = molecular viscosity / density
c
c  inputs:
c    "temp"   = temperature [deg centigrade]
c    "seamsk" = land/sea mask matching for t,s j-row (1.0=sea, 0.0=land)
c
c  outputs:
c    "scctc" = schmidt number for CCl4
c
c  local variables:
c    "tmin" & "tmax" = minimum and maximum values over which the
c                approximation can be made and is valid [deg c]
c    "sctmin", "sctmax"  = schmidt numbers for CCl4 for the minimum
c                and maximum temperatures
c    "scRikW"  = Schmidt numbers at 35%o salinity (from Rik Wannikof)
c                over the temperature range from "tmin" to "tmax"
c                at one degree C intervals
c    "itemp"   = integer temperature
c    "frtemp"  = frational temperature
c
c=======================================================================
c
#include "param.h"
c
      dimension seamsk(imt), temp(imt), scctc(imt)
      dimension scRikW(0:30)
c
      data (scRikW(N),N=0,30) /
     $    4327.1,
     $    4031.3,  3760.8,  3513.0,  3285.6,  3076.6,
     $    2884.3,  2707.1,  2543.5,  2393.3,  2252.4,
     $    2122.8,  2002.5,  1890.9,  1787.0,  1690.3,
     $    1600.2,  1516.1,  1437.5,  1364.1,  1295.3,
     $    1230.9,  1170.5,  1113.9,  1060.6,  1010.5,
     $     963.38,  918.99,  877.13,  837.64,  800.35/
c
      data tmin,   tmax   /    0.0,  30.0 /
      data sctmin, sctmax / 4327.1, 800.35 /
c
c=======================================================================
c
c  calculate schmidt numbers linearly interpolating from the given data.
c  at land points sc is set to zero.
c
      do 100 i=1,imt
        if (temp(i).le.tmin) then
          scctc(i) = sctmin*seamsk(i)
        elseif (temp(i).ge.tmax) then
          scctc(i) = sctmax*seamsk(i)
        else
          itemp = int(temp(i))
          frtemp = temp(i) - itemp
          scctc(i) = seamsk(i) * (scRikW(itemp)
     $         + frtemp * (scRikW(itemp+1) - scRikW(itemp)))
        endif
100   continue
c
#endif
      return
      end

