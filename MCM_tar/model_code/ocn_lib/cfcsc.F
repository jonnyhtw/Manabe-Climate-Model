      subroutine cfcsc (temp, seamsk, sc11, sc12)
      implicit real*8 (a-h,o-z)
#ifdef cfctracers
c
c=======================================================================
c
c  this subroutine returns estimates of schmidt numbers.
c  estimates calculated for cfc11 & cfc12 at temperatures of 0, 5, 10, 
c  15, 20, 25, & 30 degree, have been loaded into array "scest".
c  here approximations for other temperatures are made.  the empirically
c  derived Hayduk-Laudie equation was used in making the "scest"
c  schmidt number estimates.
c     schmidt num. = kinematic viscocity / molecular diffusivity
c             k.v. = molecular viscosity / density
c
c  inputs:
c    "temp"   = temperature [deg centigrade]
c    "seamsk" = land/sea mask matching for t,s j-row (1.0=sea, 0.0=land)
c
c  outputs:
c    "sc11" = schmidt number for cfc-11
c    "sc12" = schmidt number for cfc-12
c
c  local variables:
c    "coef11" & "coef12" = coefficients for computing cfc-11 schmidt
c               numbers from a polynomial of degree 6
c    "tmin & "tmax" = minimum and maximum values over which the 
c               polynomial was calculated and is valid [deg c]
c    "sc11mn", "sc12mn", "sc11mx" & "sc12mx" = schmidt numbers for 
c               cfc-11 and cfc-12 for the minimum and maximum 
c               temperatures used when deriving the polynomial
c
c=======================================================================
c
#include "param.h"
c
      dimension seamsk(imt), temp(imt), sc11(imt), sc12(imt)
      dimension coef11( 7), coef12( 7)
c
      data coef11 /
     $    0.434100E+04,  -.373206E+03,  0.259023E+02,  -.175167E+01,  
     $    0.849122E-01,  -.224272E-02,  0.236452E-04 /
c
      data coef12 /
     $    0.399400E+04,  -.342803E+03,  0.235753E+02,  -.157412E+01,
     $    0.756522E-01,  -.198789E-02,  0.208876E-04 /
c
      data tmin,   tmax   /    0.0,  30.0 /
      data sc11mn, sc11mx / 4341.0, 680.0 /
      data sc12mn, sc12mx / 3394.0, 626.0 /
c
c=======================================================================
c
c  calculate schmidt numbers from a polynomial curve fit.
c  at land points sc is set to zero.
c
      do 100 i=1,imt
        if (temp(i).le.tmin) then
          sc11(i) = sc11mn*seamsk(i)
          sc12(i) = sc12mn*seamsk(i)
        elseif (temp(i).ge.tmax) then
          sc11(i) = sc11mx*seamsk(i)
          sc12(i) = sc12mx*seamsk(i)
        else
          t2 = temp(i)*temp(i)
          t3 = temp(i)*t2
          t4 = t2*t2
          t5 = t2*t3
          t6 = t3*t3
          sc11(i) = (coef11(1) + coef11(2)*temp(i) + coef11(3)*t2 +
     $       coef11(4)*t3 + coef11(5)*t4 +coef11(6)*t5 + coef11(7)*t6)
     $       *seamsk(i)
          sc12(i) = (coef12(1) + coef12(2)*temp(i) + coef12(3)*t2 +
     $       coef12(4)*t3 + coef12(5)*t4 +coef12(6)*t5 + coef12(7)*t6)
     $       *seamsk(i)
        endif
100   continue
c
#endif
      return
      end
