      subroutine cfcsol (temp, salt, seamsk, sol11, sol12)
      implicit real*8 (a-h,o-z)
#ifdef cfctracers
c
c  sample call from subroutine setvbc:
c      call cfcsol (t(1,1,jc,nc,1), t(1,1,jc,nc,2), fm(1,1,jc),
c     $             f11sol, f12sol
c
c=======================================================================
c
c  calculate surface seawater solubilities for freons cfc11 & cfc12.
c  computes values for one model j-row at a time.
c  solublities over land are set to be zero.
c  
c  the solubilities are those reported by warner and weiss,
c  in "deep sea research", vol 32, no 12, pp 1485-1497, 1985.
c
c  inputs:
c    "temp"   = temperature (deg centigrade)
c    "salt"   = salinity (ppt)
c    "seamsk" = land/sea mask matching with "temp" and "salt"
c               (1.0=sea, 0.0=land)
c
c  outputs:
c    "sol11" = solubility for freon 11
c    "sol12" = solubility for freon 12
c  if "cfcperkg" is defined when run thru cpp, "sol11" & "sol12" are
c  in units of [mol kg-1 atm-1], otherwise, the default is for
c  volumetric units [mol l-1 atm-1]
c
c  local variables:
c    "sppt" = model salinity units converted to degrees parts per 
c            thousand, for use in the warner & weiss equations
c    "sst"  = model temperature converted to degrees kelvin and
c            divided by 100, for use in the warner & weiss equations
c    "tsq"  = variable "sst" squared
c    "tinv" = reciprocal of variable "sst"
c
c=======================================================================
c
#include "param.h"
c
# ifdef cfcperkg
c  constants for computing cfc solubilities in moist air at a total
c  pressure of 1 atmosphere in gravimetric units
c
      parameter (a1f11=-232.0411, a2f11=322.5546, a3f11=120.4956,
     $                                            a4f11=-1.39165,
     $           b1f11=-0.146531, b2f11=0.093621, b3f11=-0.0160693,
     $           a1f12=-220.2120, a2f12=301.8695, a3f12=114.8533,
     $                                            a4f12=-1.39165,
     $           b1f12=-0.147718, b2f12=0.093175, b3f12=-0.0157340)
# else
c  constants for computing cfc solubilities in moist air at a total
c  pressure of 1 atmosphere in volumetric units
c
      parameter (a1f11=-229.9261, a2f11=319.6552, a3f11=119.4471,
     $                                            a4f11=-1.39165,
     $           b1f11=-0.142382, b2f11=0.091459, b3f11=-0.0157274,
     $           a1f12=-218.0971, a2f12=298.9702, a3f12=113.8049,
     $                                            a4f12=-1.39165,
     $           b1f12=-0.143566, b2f12=0.091015, b3f12=-0.0153924)
# endif
c
      dimension temp(imt), salt(imt), seamsk(imt)
      dimension sol11(imt), sol12(imt)
c
c=======================================================================
c
      do 100 i=1,imt
        sppt = (salt(i) + 0.035)*1000.0
        sst  = (temp(i) + 273.16)*0.01
        tsqr = sst*sst
        tinv = c1/sst
        sol11(i) = exp (a1f11 + a2f11*tinv + a3f11*alog(sst) +
     $       a4f11*tsqr + sppt*(b1f11 + b2f11*sst + b3f11*tsqr) )
     $       * seamsk(i)
        sol12(i) = exp (a1f12 + a2f12*tinv + a3f12*alog(sst) +
     $       a4f12*tsqr + sppt*(b1f12 + b2f12*sst + b3f12*tsqr) )
     $       * seamsk(i)
100   continue
c
#endif
      return
      end

