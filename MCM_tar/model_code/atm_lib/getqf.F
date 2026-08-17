      subroutine getqf ( dayrad )
c
#include "parm.h"
c
#include "grdfld.h"
c
#include "param.h"
c
#include "qa.h"
c
#include "surfo.h"
c
      dimension acoef(ix, 6), bcoef(ix,6)
      dimension qfbar(ix)
      dimension qflux(ix,12)
      dimension sumq (ix)
      dimension sine ( 6,12) , cosine(6,12)
c
      data iunitr / 41 /
      data iunitw / 42 /
c
      rewind iunitw
      factor = dayrad * twopi / 365.0
      do 400 j=1,il
        read ( iunitw )  qfbar, acoef, bcoef
        do 200 i=1,ix
          sumq(i) = qfbar(i)
 200    continue
        do 300 k=1,maxwve
          factk = factor * float(k)
          sink  = sin ( factk )
          cosk  = cos ( factk )
          do 300 i=1,ix
            sumq(i) = sumq(i) + acoef(i,k) * sink
     1                        + bcoef(i,k) * cosk
 300    continue
        do 400 i=1,ix
         if (ts(i,j) .le. 0.0) then
           qbot(i,j) = sumq(i)
         else
           qbot(i,j) = 0.0
         endif
 400  continue
c
      return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      entry qcoef
c     ===== ====
c
      recip = twopi / 365.0
      delta = 365.0 / 12.0 * recip
      do 600 m=1,12
        tm = ( float(m) - 0.50 ) * delta
        do 600 k=1,6
          arg         = float(k) * tm
          sine  (k,m) = sin ( arg )
          cosine(k,m) = cos ( arg )
 600  continue
c
      rewind iunitr
      rewind iunitw
      do 900 j=1,il
        read ( iunitr ) qflux
c
        do 700 i=1,ix
          qfbar(i) = 0.0
          do 700 m=1,12
            qfbar(i) = qfbar(i) + qflux(i,m)
 700    continue
c
        recip = 1.0 / 12.0
        do 710 i=1,ix
         qfbar(i) = qfbar(i) * recip
710   continue
c
        length = ix * 6
        do 720 k=1,6
         do 720 i=1,ix
          acoef(i,k) = 0.0
          bcoef(i,k) = 0.0
  720   continue
        do 800 m=1,12
          do 800 k=1,6
            do 800 i=1,ix
              acoef(i,k) = acoef(i,k) + qflux(i,m) * sine  (k,m)
              bcoef(i,k) = bcoef(i,k) + qflux(i,m) * cosine(k,m)
 800    continue
        recip  = 2.0 / 12.0
        do 810 k=1,6
         do 810 i=1,ix
          acoef(i,k) = acoef(i,k) * recip
          bcoef(i,k) = bcoef(i,k) * recip
  810   continue
        if ( maxwve .eq. 6 )  then
          do 820 i=1,ix
           acoef(i,6) = acoef(i,6) * 0.5
           bcoef(i,6) = 0.0
  820     continue
        endif
        write ( iunitw )  qfbar, acoef, bcoef
 900  continue
c
      return
      end
