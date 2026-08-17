      subroutine fixsst ( dayrad )
c
c          called from : ocean
c         entry points : coef
c   subroutines called : none
c
c  this routine reads monthly mean values for to (mixed layer temp)
c   and yce (sea ice thickness), from which fourier coefficients are
c   computed (entry point coef).  these are used to determine values
c   of to and yce for each day (dayrad).
c
#include "parm.h"
c
#include "grdfld.h"
c
#include "param.h"
c
#include "iohist.h"
c
#include "surfo.h"
c
      common / stemp  / acoef (ix,6,2), bcoef (ix,6,2)
      common / stemp  / sst   (ix,12) , sice  (ix,12)
      common / stemp  / sstbar(ix)    , barice(ix)
      common / stemp  / sums  (ix)    , sumi  (ix)
      common / stemp  / htcont(ix)
      common / stemp  / sine  (6,12)  , cosine(6,12)
c
c set mixed layer depth(heat capicity) and latent heat of freezing
      data co / 5000.0 /, alhfrz / 80.0 /
c
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
      rewind iunitw
      factor = dayrad * twopi / 365.0
c
      do 400 j=1,il
c
c read fourier coefficients for each row
c and synthesize sst and sea ice values
        read ( iunitw )  sstbar, barice, acoef, bcoef
        do 200 i=1,ix
          sums(i) = sstbar(i)
          sumi(i) = barice(i)
 200    continue
        do 300 k=1,maxwve
          factk = factor * float(k)
          sink  = sin ( factk )
          cosk  = cos ( factk )
          do 300 i=1,ix
            sums(i) = sums(i) + acoef(i,k,1) * sink
     1                        + bcoef(i,k,1) * cosk
            sumi(i) = sumi(i) + acoef(i,k,2) * sink
     1                        + bcoef(i,k,2) * cosk
 300    continue
c
c test for sst < tfro and sea ice thickness < 0.
        do 310 i=1,ix
         sums(i) = max(sums(i), tfro)
         sumi(i) = max(sumi(i), 0.0 )
  310   continue
c
c make sst and sea ice consistent by conserving heat content
c  of the mixed layer.
        do 320 i=1,ix
         htcont(i) = (sums(i)-tfro)*co -sumi(i)*alhfrz
         if (htcont(i) .le. 0.0) then
           sumi(i) = -htcont(i)/alhfrz
           sums(i) = tfro
         else
           sumi(i) = 0.0
           sums(i) = tfro + htcont(i)/co
         endif
  320   continue
c
c load to and yce for sea points; otherwise, set to 290 & 0
        do 350 i=1,ix
         if (ts(i,j) .le. 0.0) then
           to (i,j) = sums(i)
           yce(i,j) = sumi(i)
         else
           to (i,j) = 290.0
           yce(i,j) =   0.0
         endif
  350   continue
 400  continue
c
      return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      entry coef
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
        read ( iunitr ) sst , sice
c
        do 700 i=1,ix
          sstbar(i) = 0.0
          barice(i) = 0.0
          do 700 m=1,12
            sstbar(i) = sstbar(i) + sst (i,m)
            barice(i) = barice(i) + sice(i,m)
 700    continue
        recip = 1.0 / 12.0
        do 710 i=1,ix
         sstbar(i) = sstbar(i) * recip
         barice(i) = barice(i) * recip
  710   continue
c
c       length = ix * 6 * 2
        do 720 n=1,2
         do 720 k=1,6
          do 720 i=1,ix
           acoef(i,k,n) = 0.0
           bcoef(i,k,n) = 0.0
  720   continue
        do 800 m=1,12
          do 800 k=1,6
            do 800 i=1,ix
              acoef(i,k,1) = acoef(i,k,1) + sst (i,m) * sine  (k,m)
              bcoef(i,k,1) = bcoef(i,k,1) + sst (i,m) * cosine(k,m)
              acoef(i,k,2) = acoef(i,k,2) + sice(i,m) * sine  (k,m)
              bcoef(i,k,2) = bcoef(i,k,2) + sice(i,m) * cosine(k,m)
 800    continue
        recip  = 2.0 / 12.0
        do 820 n=1,2
         do 820 k=1,6
          do 820 i=1,ix
           acoef(i,k,n) = acoef(i,k,n) * recip
           bcoef(i,k,n) = bcoef(i,k,n) * recip
  820   continue
        if ( maxwve .eq. 6 )  then
           do 830 i=1,ix
            acoef(i,6,1) = acoef(i,6,1) * 0.5
            acoef(i,6,2) = acoef(i,6,2) * 0.5
            bcoef(i,6,1) = 0.0
            bcoef(i,6,2) = 0.0
  830      continue
        end if
        write ( iunitw )  sstbar, barice, acoef, bcoef
 900  continue
c
      return
      end
