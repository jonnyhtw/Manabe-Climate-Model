      subroutine  CTCair (ccl4yr, aCCL4)
      implicit real*8 (a-h,o-z)
#ifdef ccl4tracer
c
c======================================================================
c
c     This gives the atmospheric mixing ratio in [pptv] for carbon
c     tetrachloride.  Values are from Doug Wallace (personal
c     communication).  The time series runs from 1910.00 to 1990.99
c     Values are considered to be for 1 January of each year, 
c     calculated up to 1974, measured after that through 1988,
c     extrapolated beyond 1987.
c
c  inputs:
c     "ccl4yr" = the real number of years since 1 JAN 1910.
c
c  output:
c     "aCCL4"  = the atmospheric mixing ratio
c
c  internal parameters:
c     "nyrs"   = the length of the time series
c     "CCL4mr" = annual CCL4 mixing ratios for Jan of given year
c                (1910- ??)
c     "CTCyr"  = real year in the CCl4 time series which starts on
c                1 Jan 1910 (or 15 Dec 1909).
c     "iCTCyr" = integer year in the range 1910 to 1991
c     "fCTCyr" = fraction of year
c
c======================================================================
c
#include "param.h"
c
      PARAMETER (numyrs=82)
      DIMENSION CCL4mr(numyrs)
c
c     the 0.02 corresponds to 1 Jan 1910
c
c     [pptv, parts per trillion by volume]
      DATA (CCL4mr(N),N=1,numyrs)/
     $    0.02,  
     $    0.03,   0.06,   0.11,   0.17,   0.24,
     $    0.39,   0.66,   0.97,   1.30,   1.44,
     $    1.58,   1.64,   1.87,   2.16,   2.47,
     $    2.82,   3.24,   3.59,   4.01,   4.76,
     $    5.49,   6.27,   7.08,   7.86,   8.91,
     $   10.13,  11.64,  13.42,  15.01,  16.84,
     $   18.84,  21.24,  24.21,  27.78,  31.62,
     $   34.66,  36.64,  39.68,  42.64,  44.87,
     $   46.83,  46.83,  49.23,  50.56,  52.08,
     $   52.50,  53.42,  53.78,  54.40,  55.85,
     $   57.23,  58.37,  59.15,  60.15,  61.45,
     $   62.64,  63.80,  65.59,  68.26,  71.33,
     $   76.62,  84.86,  92.34,  96.86, 100.66,
     $  105.00, 108.00, 110.00, 112.00, 116.00,
     $  121.00, 123.00, 124.00, 126.00, 127.00,
     $  129.50, 131.00, 132.00, 133.00, 134.50,
     $  136.00/
cc    the following lines extend the extrapolation out to the year 2001
cc    (will need to increase "numyrs" to 92)
cc   $  136.00, 137.50, 139.00, 140.50, 142.00,
cc   $  143.50, 145.00, 146.50, 148.00, 149.50,
cc   $  151.00/
c
c
c======================================================================
c
c     relate the year to the index of the data vector
c
      CTCyr = ccl4yr - 1910.00 + 1.0
      iCTCyr= int(CTCyr)
      fCTCyr= CTCyr-iCTCyr
c
c     estimate the atmospheric mixing ratio of carbon tetrachloride
c     (linearly interpolating with respect to time and checking that
c     the year requested is within the data time range)
c
      if (iCTCyr .le. 0) then
        aCCL4 = c0
      elseif ((iCTCyr+1) .gt. numyrs) then
        write(stdout,9013) ccl4yr
        stop ' CTCair'
      else
        aCCL4 = CCL4mr(iCTCyr) 
     $          + fCTCyr*(CCL4mr(iCTCyr+1) - CCL4mr(iCTCyr))
      endif
c
#endif
      return
9013  format(/' error => beyond end of carbon-tet data: year = ',f9.2)
      end



