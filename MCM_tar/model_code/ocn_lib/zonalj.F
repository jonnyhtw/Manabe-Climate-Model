      subroutine zonalj (j)
      implicit real*8 (a-h,o-z)
#ifdef zonalavg
cfpp$ noconcur r
c
c=======================================================================
c
c     compute zonal means for the tracer and u, v, w
#ifdef seaice
c     and sea ice.
#endif
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "cprnts.h"
#include "cregin.h"
#include "ctmngr.h"
#include "docnam.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#include "slabs.h"
#include "switch.h"
#include "zmvbls.h"
c
      if (j .eq. 1) then
c
c
c-----------------------------------------------------------------------
c       initialize various quantities when j = 1
c-----------------------------------------------------------------------
c
        do 100 m=1,nt
          do 90 k=1,km
            do 80 jj=1,jmt
              zmt(jj,k,m) = c0
80          continue
90        continue
100     continue
c
        do 200 k=1,km
          do 190 jj=1,jmt
            zmu(jj,k) = c0
            zmv(jj,k) = c0
            zmw(jj,k) = c0
            zwt(jj,k) = c0
            zwu(jj,k) = c0
190       continue
200     continue
c
#ifdef seaice
        do 300 jj=1,jmt
          zmi(jj) = c0
300     continue
c
#endif
      elseif (j .lt. jmt) then
c 
c-----------------------------------------------------------------------
c        compute zonal means for each tracer and for u, v, w
# ifdef seaice
c        and sea ice
# endif
c-----------------------------------------------------------------------
c
        do 400 m=1,nt
          do 390 k=1,km
            do 380 i=2,imtm1
              if (k .le. kmt(i,j)) then
                zmt(j,k,m) = zmt(j,k,m) + t(i,k,jc,nc,m)*dxt(i)
              endif
380         continue
390       continue
400     continue
c
        do 500 k=1,km
          do 490 i=2,imtm1
            if (k .le. kmt(i,j)) then
              zmw(j,k) = zmw(j,k) + dxt(i)*w(i,k)
              zwt(j,k) = zwt(j,k) + dxt(i)
            endif
490       continue
500     continue
c
        do 600 k=1,km
          do 590 i=2,imtm1
            if (k .le. kmu(i,j)) then
              zmu(j,k) = zmu(j,k) + u(i,k,jc,nc)*dxu(i)
              zmv(j,k) = zmv(j,k) + v(i,k,jc,nc)*dxu(i)
              zwu(j,k) = zwu(j,k) + dxu(i)
            endif
590       continue
600     continue
c
#ifdef seaice
        do 700 i=2,imtm1
          if (kmt(i,j) .ne. 0) then
            zmi(j) = zmi(j) + yc(i,j,nc)*dxt(i)
          endif
700     continue
#endif
c
c-----------------------------------------------------------------------
c
      endif
c
#endif
      return
      end
