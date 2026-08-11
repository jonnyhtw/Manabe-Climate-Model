      subroutine filfir (t, f, s, kl, kind, num)
      implicit real*8 (a-h,o-z)
#ifdef firfil
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c     simple finite impulse response filter with [.25, .5, .25] weights
c     (shapiro) modified for assymetric and symmetric boundaries
c
c     inputs:
c             t    = (imt,kl) array of quantity to be filtered along
c                    the first dimension.
c                    note: t(i,k) must be zero where f(i,k) = zero
c                    for this filter to work.
c             f    = (imt,kl) mask of zeroes & ones to indicate land
c                    and ocean. zero indicates a land point
c             s    = (imt,kl) scratch array
c             kl   = number of vertical levels to filter
c             kind = (0,1) = (symmetric, asymmetric) boundary condition
c                    symmetric is appropriate for tracers & vorticity
c                    asymmetric is appropriate for velocities
c             num  = number of times to apply this filter
c
c     outputs:
c             t    = (imt,km) array of filtered quantities
c-----------------------------------------------------------------------
c
#include "param.h"
c
      dimension t(imt,kl), f(imt,kl), s(imt,kl)
c
c=======================================================================
c
      if (kind .eq. 0) then
c
c-----------------------------------------------------------------------
c       apply the filter "num" times using a symmetric (no flux)
c       boundary condition
c-----------------------------------------------------------------------
c
#ifdef cyclic
        do 100 k=1,kl
          t(1,k)   = t(imtm1,k)
          t(imt,k) = t(2,k)
100     continue
#endif
        do 200 n=1,num
c
c         each filtering consists of a double pass
c
          do 190 k=1,kl
            do 170 i=2,imtm1
              s(i,k) = f(i,k)*(p25*(t(i-1,k) + t(i+1,k)) +
     $                 t(i,k)*(c1 - p25*(f(i-1,k) + f(i+1,k))))
170         continue
#ifdef cyclic
            s(1,k)   = s(imtm1,k)
            s(imt,k) = s(2,k)
#endif
            do 180 i=2,imtm1
              t(i,k) = f(i,k)*(p25*(s(i-1,k) + s(i+1,k)) +
     $                 s(i,k)*(c1 - p25*(f(i-1,k) + f(i+1,k))))
180         continue
#ifdef cyclic
            t(1,k)   = t(imtm1,k)
            t(imt,k) = t(2,k)
#endif
190       continue
200     continue
      elseif (kind .eq. 1) then
c
c----------------------------------------------------------------------
c       apply the filter "num" times using an asymmetric (flux)
c       boundary condition
c----------------------------------------------------------------------
c
#ifdef cyclic
        do 300 k=1,kl
          t(1,k)   = t(imtm1,k)
          t(imt,k) = t(2,k)
300     continue
#endif
        do 400 n=1,num
c
c         each filtering consists of a double pass
c
          do 390 k=1,kl
            do 370 i=2,imtm1
              s(i,k) = f(i,k)*(p25*t(i-1,k) + p5*t(i,k) + p25*t(i+1,k))
370         continue
#ifdef cyclic
            s(1,k)   = s(imtm1,k)
            s(imt,k) = s(2,k)
#endif
            do 380 i=2,imtm1
              t(i,k) = f(i,k)*(p25*s(i-1,k) + p5*s(i,k) +p25*s(i+1,k))
380         continue
#ifdef cyclic
            t(1,k)   = t(imtm1,k)
            t(imt,k) = t(2,k)
#endif
390       continue
400     continue
      else
        write (stderr,99) kind
        stop ' filfir '
      endif
c
99    format(/' error => kind =',i10, ' in filfir')
c
#endif
      return
      end
