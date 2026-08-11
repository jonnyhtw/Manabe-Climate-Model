      subroutine zonal
      implicit real*8 (a-h,o-z)
#ifdef zonalavg
cfpp$ noconcur r
c
c=======================================================================
c
c     compute zonal means (layer by layer) for the tracers and u, v, w
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
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#include "zmvbls.h"
c
      dimension number(jmt)
c
c 
c-----------------------------------------------------------------------
c        compute zonal means for each tracer and for u, v, w
# ifdef seaice
c        and sea ice
# endif
c-----------------------------------------------------------------------
c
      do 30 m=1,nt
        do 20 k=1,km
          do 10 j=1,jmt
            if (zwt(j,k) .ne. c0) then
              zmt(j,k,m) = zmt(j,k,m)/zwt(j,k)
            endif
10        continue
20      continue
30    continue
c
      do 90 k=1,km
        do 80 j=1,jmt
          if (zwt(j,k) .ne. c0) then
            zmw(j,k) = zmw(j,k)/zwt(j,k)
          endif
80      continue
90    continue
c
      do 200 k=1,km
        do 190 j=1,jmt
          if (zwu(j,k) .ne. c0) then
            zmu(j,k) = zmu(j,k)/zwu(j,k)
            zmv(j,k) = zmv(j,k)/zwu(j,k)
          endif
190     continue
200   continue
c
#ifdef seaice
      do 300 j=1,jmt
        if (zwt(j,1) .ne. c0) then
          zmi(j) = zmi(j)/zwt(j,1)
        endif
300   continue
c
#endif
c-----------------------------------------------------------------------
c    print out zonal mean table
c-----------------------------------------------------------------------
c
#ifdef ocnr15
      iskp = 2
#endif
#ifdef ocnr30
      iskp = 4
#endif
      scale = c1
c
      do 400 j=1,jmt
        number(j) = j
400   continue
c
      do 500 m=1,nt
        if (m .eq. 2) then
          scale = 1000.0
        else
          scale = c1
        endif
        write(6,9001) trname(m), scale, (yt(j),j=2,jmtm1,iskp)
     $,               (number(j),j=2,jmtm1,iskp)
        if (m.eq.3) then
          do 480 k=1,km
            write(6,9021) k, (scale*zmt(j,k,m),j=2,jmtm1,iskp)
480       continue
        else
          do 490 k=1,km
            write(6,9011) k, (scale*zmt(j,k,m),j=2,jmtm1,iskp)
490       continue
        endif
500   continue
c
      scale = 1000.0
      write(6,9001) 'w velocity  ', scale, (yu(j),j=2,jmtm1,iskp)
     $,             (number(j),j=2,jmtm1,iskp)
      do 600 k=1,km
        write(6,9011) k, (scale*zmw(j,k),j=2,jmtm1,iskp)
600   continue
c
c
      scale = c1
      write(6,9001) 'u velocity  ', scale, (yu(j),j=2,jmtm1,iskp)
     $,             (number(j),j=2,jmtm1,iskp)
      do 700 k=1,km
        write(6,9011) k, (scale*zmu(j,k),j=2,jmtm1,iskp)
700   continue
c
      scale = c1
      write(6,9001) 'v velocity  ', scale, (yu(j),j=2,jmtm1,iskp)
     $,             (number(j),j=2,jmtm1,iskp)
      do 800 k=1,km
        write(6,9011) k, (scale*zmv(j,k),j=2,jmtm1,iskp)
800   continue
c
#ifdef seaice
      scale = 0.01
      write(6,9001) 'sea ice     ', scale, (yt(j),j=2,jmtm1,iskp)
     $,             (number(j),j=2,jmtm1,iskp)
      k = 1
      write(6,9011) k, (scale*zmi(j),j=2,jmtm1,iskp)
c
#endif
c-----------------------------------------------------------------------
c
#endif
      return
c
9001  format(/,' zonal means: ',a12,5x,' scale =',e12.6, /,4x,20f6.1,
     $       /,3x,20i6)
9011  format(1x,i2,1x,20f6.2)
9021  format(1x,i2,1x,20f6.0)
c
      end
