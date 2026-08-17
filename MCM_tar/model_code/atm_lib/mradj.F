      subroutine mradj ( rp, ps, zmadj )
c
c         called from: passb
c  subroutines called: none
c        entry points: none
c
c this routine attempts to borrow water from the four closest points
c  in the slab (i-1,i+1,k-1,k+1) inorder to fill the negative mixing
c  ratio points.  if not enough water is avaible no borrowing is
c  performed.
c
#include "parm.h"
#include "parmkp.h"
c
cSGI
      integer*8 itwo
cSGI
c
      logical even, mfirst
c
      real ratio, sum
      dimension rp(ix,kx), ps(ix), zmadj(il)
c     dimension negbit(ixp2,kxp2)
c
#include "auxatm.h"
c
#include "ctlblk.h"
c
      common / mrjcom / dqt(kxp2)
c
#include "sgwgt1.h"
c
      dimension rsave (ix,kx), rptemp(ixp2,kxp2)
      dimension temprp(ix,kx)
      dimension rps(ix), cmadj(ix)
c
      data mfirst / .true. /
c
c  scan and adjust mixing ratios backwards on odd time steps
c
c
cSGI
      itwo = 2
cc    even = mod(ktau,2) .eq. 0
      even = mod(ktau,itwo) .eq. 0
cSGI
c
c  rptemp is used as a work area for rp.
c
      do 100 k=1,kx
        do 100 ip=2,ixp
          rptemp(ip,k+1) = rp(ip-1,k) * ps(ip-1)
  100 continue
c
c  load rps for later use
c
        do 110 i=1,ix
          rps(i) = 1.0 / ps(i)
  110 continue
c
c  load dqt only on first time into mradj
c
      if (mfirst) then
        do 200 k=1,kx
          dqt(k+1) = dq(k)
 200    continue
        dqt(   1) = 0.0
        dqt(kxp2) = 0.0
        mfirst    = .false.
      endif
c
c  rptemp has zeroes in levels 1 & kxp2 with rp placed in between these.
c  thus any point in rp has four neighbor points in rptemp.
c
      do 300 ip=1,ixp2
        rptemp(ip,  1)  = 0.0
        rptemp(ip,kxp2) = 0.0
  300 continue
c
c  set boundaries.
c
      do 500 k=1,kx
        rptemp(  1 ,k+1) = rptemp(ixp,k+1)
        rptemp(ixp2,k+1) = rptemp( 2 ,k+1)
  500 continue
c
c
c  test negative rp points and do adjustment where necessary.
c
      do 650 k=1,kx
       do 601 i=1,ixp
        if (rptemp(i,k+1) .lt. 0.0) go to 771
  601  continue
       ixxx = 0
       go to 772
  771  continue
       ixxx = 1
  772  continue
       if (ixxx .eq. 0) go to 650
        do 600 i=1,ix
          if (even) then
            ip = i + 1
          else
            ip = ixp2 - i
          endif
          if (rptemp(ip,k+1) .lt. 0.0) then
            sum = rptemp(ip,k) * dqt(k) +
     *           (rptemp(ip-1,k+1) + rptemp(ip+1,k+1)) * dqt(k+1) +
     *            rptemp(ip,k+2) * dqt(k+2)
            if (sum .ge. abs(rptemp(ip,k+1) * dqt(k+1))) then
              ratio  = 1.0 + rptemp(ip,k+1) * dqt(k+1) / sum
              rptemp(ip  , k ) = rptemp(ip  , k ) * ratio
              rptemp(ip-1,k+1) = rptemp(ip-1,k+1) * ratio
              rptemp(ip+1,k+1) = rptemp(ip+1,k+1) * ratio
              rptemp(ip  ,k+2) = rptemp(ip  ,k+2) * ratio
              rptemp(ip  ,k+1) = 0.0
c
c reset boundaries of rptemp for cyclic consistency
c
              if (ip .le. 3) then
                rptemp(ixp ,k  ) = rptemp(1,k  )
                rptemp(ixp2,k  ) = rptemp(2,k  )
                rptemp(ixp ,k+1) = rptemp(1,k+1)
                rptemp(ixp2,k+1) = rptemp(2,k+1)
                rptemp(ixp ,k+2) = rptemp(1,k+2)
                rptemp(ixp2,k+2) = rptemp(2,k+2)
              elseif (ip .ge. ix) then
                rptemp(   1,k  ) = rptemp(ixp ,k  )
                rptemp(   2,k  ) = rptemp(ixp2,k  )
                rptemp(   1,k+1) = rptemp(ixp ,k+1)
                rptemp(   2,k+1) = rptemp(ixp2,k+1)
                rptemp(   1,k+2) = rptemp(ixp ,k+2)
                rptemp(   2,k+2) = rptemp(ixp2,k+2)
              endif
            endif
          endif
  600   continue
  650 continue
c
c  replace adjusted rp values in the rp array and compute net filling
c
      do 700 k=1,kx
       do 700 i=1,ix
        temprp(i,k) = rptemp(i+1,k+1) * rps(i)
 700  continue
      do 701 k=1,kx
       do 701 i=1,ix
        fillmr(i,k) = temprp(i,k) - rp(i,k)
        rp    (i,k) = temprp(i,k)
 701  continue
c
c  compute integral of pressure weighted negative mixing ratios only
c   when needed, otherwise store negative infinitity in zmadj
c
 1000 continue
c
      if (.not.balnce) then
        zmadj(jrow) = -999.0e23
      else
        do 810 k=1,kx
         do 810 i=1,ix
          if (rp(i,k) .ge. 0.0) then
            rsave(i,k) = 0.0
          else
            rsave(i,k) = rp(i,k) * ps(i)
          endif
  810   continue
        do 820 i=1,ix
          cmadj(i) = dq(1) * rsave(i,1)
  820   continue
        do 830 k=2,kx
          do 830 i=1,ix
            cmadj(i) = cmadj(i) + dq(k) * rsave(i,k)
  830   continue
        zmadj (jrow) = 0.0
        do 840 i=1,ix
         zmadj(jrow) = zmadj(jrow) + cmadj(i)
  840   continue
      endif
c
      return
      end
