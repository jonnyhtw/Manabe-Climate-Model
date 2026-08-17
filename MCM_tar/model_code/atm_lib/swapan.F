      subroutine swapan
#ifdef dirnal
c
c  swapan is called by atmos in diurnal models just before end of run.
c     entry points:       gmtan, zeroac
c     subroutines called: uvphrm (by gmtan)
c  swapan moves diurnal model analysis data from units 32, 33, 53 and 54
c   to units 30, 31, 51 and 52 respectively at the end of the monthly
c   run if the latest data are on the former units.  this occurs on
c   months with an even number of days.
c
c
#include "parm.h"
#include "parmout.h"
#include "parmkp.h"
#include "parmjx.h"
#include "parmdsk.h"
      parameter ( ixil = ix * il )
      parameter ( njm = 2 * jx * mx, njmk = kx * njm )
      parameter ( nhvp = 17 )
#ifdef fixcld
      parameter ( nhv = 13 )
#endif
#ifdef prdcld
      parameter ( nhv = kx + 14 )
#endif
c
      complex vor  , div  , t   , psln
      complex vorm , divm , tm  , psm
      complex vordt, divdt, tdt , psdt
      complex ucos , vcos , psdx, psdy, zs, ek
      complex qdot , somg
      complex qmix , qmixm, qmxdt
#ifdef rad
      complex radten
#endif
c
      dimension ucosg (njmk) , vcosg (njmk), tempg(njmk), qmixg(njmk)
      dimension qdotg (njmk) , omgmt (njmk), radgt(njmk)
      dimension anala (awrds), analb(bwrds)
      dimension accmg (ixil) , field(ixil,nhvp)
      dimension lflag (ixil) , flag (ixil)
c
#include "auxatm.h"
c
#include "ctlblk.h"
c
#include "cloudp.h"
c
#include "grdfld.h"
c
#include "nxstps.h"
c
#include "radten.h"
c
#include "sfharm.h"
c
#include "sqdot.h"
c
      dimension accm(njmk)
c
#include "water.h"
c
      equivalence ( ucosg   (1), ucos  (1,1,1))
      equivalence ( vcosg   (1), vcos  (1,1,1))
      equivalence ( tempg   (1),  t    (1,1,1))
      equivalence ( radgt   (1), radten(1,1,1))
      equivalence ( qmixg   (1), qmix  (1,1,1))
#ifdef dirnal
      equivalence ( omgmt   (1), somg  (1,1,1))
      equivalence ( qdotg   (1), qdot  (1,1,1))
#endif
      equivalence ( accmg   (1), accm      (1))
      equivalence ( field (1,1), ts      (1,1))
      equivalence ( anala   (1), vfu     (1,1))
#ifdef moist
      equivalence ( analb   (1), convr   (1,1))
#endif
#ifdef dry
      equivalence ( analb   (1), convt   (1,1))
#endif
#ifdef prdcld
      equivalence ( lflag   (1), kflags  (1,1))
#endif
c
      data lg / 51 /
      data lh / 52 /
c
      if (nin .eq. 0) return
      nou   = 0
      lgi   = lg + nin
      lgo   = lg + nou
      lhi   = lh + nin
      lho   = lh + nou
      adski = adsk + nin
      adsko = adsk + nou
      bdski = bdsk + nin
      bdsko = bdsk + nou
      rewind lgi
      rewind lgo
      rewind lhi
      rewind lho
      rewind adski
      rewind adsko
      rewind bdski
      rewind bdsko
c
      do 2 lt=1,ltm
        call  rdctl(lgi)
        call wrtctl(lgo, 0)
        call  rdctl(lhi)
        call wrtctl(lho, 0)
        do 3 lv=1,7
          read (lgi) qdot
          write(lgo) qdot
    3   continue
        do 4 lv=1,nhv
          read (lhi) shfx
          write(lho) shfx
    4   continue
    2 continue
c
      do 5 lt=1,ltm
        call  rdctl(adski)
        call wrtctl(adsko, 0)
        call  rdctl(bdski)
        call wrtctl(bdsko, 0)
        do 5 j=1,il
          read (adski) anala
          write(adsko) anala
          read (bdski) analb
          write(bdsko) analb
    5 continue
c
      rewind lgi
      rewind lgo
      rewind lhi
      rewind lho
      rewind adski
      rewind adsko
      rewind bdski
      rewind bdsko
      write(6,20) ktau
      return
c
      entry gmtan
c     ===== =====
c
c  gmtan reads accumulated diurnal model analysis data from disk, adds
c   current data to these gmt fields and writes these sums on another
c   set of disk units.  it is called by atmos near the end of each
c   radiation (gmt output) time step of diurnal model runs.
c
      call uvphrm
      nou = iabs(nin-2)
      lgi = lg + nin
      lhi = lh + nin
      lgo = lg + nou
      lho = lh + nou
      read (lgi) accm
      accm = accm + ucosg
      write(lgo) accm
      read (lgi) accm
      accm = accm + vcosg
      write(lgo) accm
      read (lgi) accm
      accm = accm + tempg
      write(lgo) accm
      read (lgi) accm
      accm = accm + qmixg
      write(lgo) accm
      read (lgi) accm
      accm = accm + qdotg
      write(lgo) accm
      read (lgi) accm
      accm = accm + omgmt
      write(lgo) accm
      read (lgi) accm
      accm = accm + radgt
      write (lgo) accm
c
      do 7 lv=1,nhvp
        if ((lv .ne. 2) .and. (lv .ne. 4) .and. (lv .ne. 12)
     *      .and. (lv .ne. 16)) then
          read (lhi) accmg
          do 577 i=1,ixil
           accmg(i) = accmg(i) + field(i,lv)
  577     continue
          write(lho) accmg
        endif
    7 continue
c
#ifdef prdcld
      do 8 k=1,kx
        read (lhi) accmg
        do 578 i=1,ixil
         accmg(i) = accmg(i) + float(and(1,(lflag(i)/ipower(k))))
  578   continue
        write(lho) accmg
    8 continue
c
      read(lhi) accmg
      do 579 i=1,ixil
       if (lflag(i) .gt. 0) then
         flag(i) = 1.0
       else
         flag(i) = 0.0
       endif
  579 continue
        do 571 i=1,ixil
         accmg(i) = accmg(i) + flag(i)
  571   continue
      write(lho) accmg
#endif
      return
c
      entry zeroac
c     ===== ======
c  zeroac zeroes diurnal model analysis data on disk at the beginning
c   of every run before data accumulation is started.  it is called by
c   atmos near the beginning of this program in diurnal model runs.
c
      rewind lg
      rewind lh
      rewind adsk
      rewind bdsk
      narray(1)    = 1
      qdot         = (0.0,0.0)
      shfx         = 0.0
c
      do 10 lt=1,ltm
        call wrtctl(lg, 0)
        call wrtctl(lh, 0)
        do 17 lv=1,7
          write(lg) qdot
   17   continue
        do 11 lv=1,nhv
          write(lh) shfx
   11   continue
   10 continue
c
      anala = 0.0
      analb = 0.0
c
      do 14 lt=1,ltm
        call wrtctl(adsk, 0)
        call wrtctl(bdsk, 0)
        do 13 j=1,il
          write(adsk) anala
          write(bdsk) analb
   13   continue
   14 continue
c
      rewind lg
      rewind lh
      rewind adsk
      rewind bdsk
      return
c
   20 format(/' analysis data transferred to primary disks at tau=',i7/)
c
#endif
      end
