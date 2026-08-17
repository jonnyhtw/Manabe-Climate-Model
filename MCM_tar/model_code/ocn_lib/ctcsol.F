       subroutine CTCsol (temp, salt, seamsk, solCL4)
      implicit real*8 (a-h,o-z)
#ifdef ccl4tracer
c
c  sample call from subroutine setvbc:
c      call CTCsol (t(1,1,jc,nc,1), t(1,1,jc,nc,2), fm(1,1,jc),
c     $             CL4sol)
c
c=======================================================================
c
c  calculate surface seawater solubilities for CCl4.
c  computes values for one model j-row at a time.
c  solublities over land are set to be zero.
c
c      (author: C.J. Beegle; 1991)
c
c  the solubilities are those reported by Hunter-Smith, Balls and Liss
c  (1983).  Henry's Law constants and the air-sea exchange of various
c  molecular weight halocarbon gases.  Tellus, 35B: 170-176.
c
c  inputs:
c    "temp"   = temperature (deg centigrade)
c    "salt"   = salinity (ppt)
c    "seamsk" = land/sea mask matching with "temp" and "salt"
c               (1.0=sea, 0.0=land)
c
c  outputs:
c    "solCL4" = solubility for carbon tetrachloride
c               incompressible units [mol kg-1 atm-1]
c
c  local variables:
c    "sppt" = model salinity units converted to degrees parts per
c             thousand
c    "sst"  = model temperature converted to degrees kelvin
c    "tinv" = reciprocal of variable "sst"
c    "cw1CL4", "cw2CL4", "dw1CL4", "dw2CL4" = coefficients for the
c             Henry's Law equations for seawater (cw) and distilled
c             water (dw) for solubility for CCl4.
c     "es"  = saturation partial pressure on water in air as a
c             function of temperature (since the air is directly over
c             seawater, assume that the air is saturated). This is
c             Tenet's formula.
c     "R"   = Universal gas constant
c     "rhow"= density of seawater to a mild accuracy.
c
c=======================================================================
c
c
#include "param.h"
c
c      constants for computing carbon tetrachloride solubilities for
c      unionized gas in water.  Formulae are Henry's Law coefficients
c      (H).  The reference gives formulae for distilled water (HdH2O)
c      and 35ppt (HcH2O).  Concentrations for other salinities are 
c      linearly interpolated for the given temperature.
c
      parameter (cw1CL4=-3230., cw2CL4=11.27,
     $           dw1CL4=-2918., dw2CL4=9.77,
     $           rhow=1026.97,
c                             [kg m-3])
     $           R=8.3143   )
c                             [J deg-1 mole-1]
c
      dimension temp(imt), salt(imt), seamsk(imt)
      dimension solCL4(imt)
c
c=======================================================================
c
      do 100 i=1,imt
        sppt = (salt(i) + 0.035)*1000.0
        sst  = (temp(i) + 273.16)
        tinv = c1/sst
        HdH2O = exp(dw1CL4*tinv + dw2CL4)
        HcH2O = exp(cw1CL4*tinv + cw2CL4)
        HCCL4 = HdH2O + (35.0 - sppt)/100.0 * (HcH2O - HdH2O)
c
c      Well, now we have a Henry's law coefficient (which, if you read
c      the paper carefully, is really an inverse Ostwald coefficient).
c      Now come the shinanigans for getting things in the right units
c      so that when the result is multiplied by the atmospheric mole
c      fraction, the results will be in the units [mole kg-1 atmos-1].
c
      es = 6.22*exp( (sst-273.0)/(sst-36.0) * 17.27)
      solCL4(i) = ((1000.0-es)*100.0/(R*sst*HCCL4)/(rhow))*seamsk(i)
c
c              [mbar][100kg m-1 s-2]                1
c                    --------------- * ---------------------------------
c                         [mbar]       [J deg-1 mole-1][kg m2 s-2][ deg]
c                                                      -----------
c                                                          [J]
c
c                    1        [moles]
c               * --------  = -------
c                 [kg m-3]     [kg]
100   continue
#endif
      return
      end


