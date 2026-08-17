      subroutine drycnv(t0)
c
c    called by: convad
c    calls    : nothing
c this subroutine computes a dry adiabatic adjustment.
c
#include "parm.h"
#include "parmkp.h"
#include "parmout.h"
#include "parmleng.h"
c
c
      dimension ta(ix,kx,nhem), t0(ix,kx)
c
#include "auxatm.h"
c
#include "ctlblk.h"
c
#include "grid.h"
c
#include "sgwgt1.h"
c
      dimension test(ix,kx)
c
      equivalence (ta(1,1,1),cgt(1,1,1))
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c here begins code for performing the dry adiabatic adjustment.
c adjustment is performed for those layers where the computed
c temperature lapse rate exceeds the dry adiabatic lapse rate.
c
      ridt = 1.0 / float(idt)
      do 580 iter=1,20
        isum1 = 0
        do 570 k=1,km
         do 76 i=1,ix
          test(i,3) = (dqmqa5(k) - 1.0) * ta(i,k+1,ih)
     *              + (dqmqa5(k) + 1.0) * ta(i,k  ,ih)
     *              + 0.01 + den1(k)
   76    continue
c itest(i,3) = 1 for those layers where the adjustment is to be
c performed and zero for no adjustment.
c
         isum = 0
         do 77 i=1,ix
          if (test(i,3) .lt. 0.0) go to 771
   77    continue
         isum = 0
         go to 772
  771    continue
         isum = 1
  772    continue
         isum1 = isum1 + isum
c
c skip dry adjustment if no points in row at level k require
c adjustment.
c
          if (isum .ne. 0) then
c
c       adjust the points corresponding to itest(i,3)=1
c
            do 78 i=1,ix
             if (test(i,3) .lt. 0.0) then
               ta(i,k  ,ih) = (ta(i,k  ,ih)*dq(k  )
     *            + dq(k+1) * (ta(i,k+1,ih) - den2(k)))
     *                      *  den(k)
               ta(i,k+1,ih) =  ta(i,k  ,ih)*pq(k) + den2(k)
             endif
   78       continue
          endif
  570   continue
c
c leave iterative loop 580 only after all levels and all points are
c dry adiabatically adjusted.
c
        if (isum1 .eq. 0) go to 120
  580 continue
  120 continue
c
c compute rate of change of temperature due to dry adiabatic
c adjustment.
c
      do 573 k=1,kx
        do 573 i=1,ix
          dadadj(i,k) = ridt * (ta(i,k,ih) - t0(i,k))
  573 continue
c
      return
      end
