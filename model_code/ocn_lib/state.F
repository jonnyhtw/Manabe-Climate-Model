      subroutine state (t, s, rho)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c
c=======================================================================
c
c     state computes one row of normalized densities by using a 3rd
c     order polynomial fit to the knudsen formula, for each level
c     subtract normalizing constants from temperature and salinity
c     and compute polynomial approximation of knudsen density.
c
c     note.. for precision purposes, there is a depth dependent
c     constant subtracted from the density returned by this routine.
c     so... this routine should be used only for horizontal gradients
c     of density.
c
c     inputs:
c
c     t = the input row of temperatures
c     s = the input row of salinities (units: (ppt-35)/1000)
c
c     output:
c
c     rho = normalized densities
c
c=======================================================================
c
c
#include "param.h"
c
      dimension t(imt,km), s(imt,km), rho(imt,km)
cSGI
cc    common /cstate/ to(km), so(km), c(km,9)
#include "cstate.h"
c
c moved to blkdta.F
cinclude "dncoef.h"
c
cSGI
c-----------------------------------------------------------------------
c     statement function
c-----------------------------------------------------------------------
c
      dens (tq, sq, k) = (c(k,1) + (c(k,4) + c(k,7)*sq)*sq +
     $                   (c(k,3) + c(k,8)*sq + c(k,6)*tq)*tq)*tq +
     $                   (c(k,2) + (c(k,5) + c(k,9)*sq)*sq)*sq
c
c-----------------------------------------------------------------------
c
      do 100 k=1,km
        do 90 i=1,imt
          rho(i,k) = dens (t(i,k)-to(k), s(i,k)-so(k), k)
90      continue
100   continue
c
      return
c
c
c
c
      entry statec (t, s, rho, ind)
c
c
c=======================================================================
c
c     statec computes, for one row, the normalized densities by using
c     a 3rd order polynomial fit to the knudsen formula. For
c     purposes of checking vertical stability between adjacent
c     levels, the reference depth for pressure dependence in
c     the knudsen formula must be held constant. that level is
c     determined by "ind".
c
c     inputs:
c
c     t   = the input row of temperatures
c     s   = the input row of salinities (units: (ppt-35)/1000)
c     ind = 1 for comparing levels 1 to 2, 3 to 4, etc.
c           (coefficients for the lower of the 2 levels are used)
c           2 for comparing levels 2 to 3, 4 to 5, etc.
c           (coefficients for the lower of the 2 levels are used)
c
c     output:
c
c     rho = normalized densities
c
c=======================================================================
c
c
      if (ind .lt. 1 .or. ind .gt. 2) then
        write (stderr,99)
        stop '1 statec'
      endif
c
      do 200 l=1,km,2
        if (ind .eq. 1) then
          k = min(l+1,km)
        else
          k = l
        endif
        do 190 i=1,imt
          rho(i,l) = dens (t(i,l)-to(k), s(i,l)-so(k), k)
190     continue
200   continue
c
      do 300 l=2,km,2
        if (ind .eq. 1) then
          k = l
        else
          k = min(l+1,km)
        endif
        do 290 i=1,imt
          rho(i,l) = dens (t(i,l)-to(k), s(i,l)-so(k), k)
290     continue
300   continue
      return
   99 format(/' error => bad "ind" in statec: ind =',i10)
      end
