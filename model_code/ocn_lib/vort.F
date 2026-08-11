      subroutine vort
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     compute vorticity for relaxation
c
c=======================================================================
c
#include "param.h"
#include "cdiag.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#ifdef multitasking
      dimension tempik(imt,km,1)
#else
#include "temp.h"
#endif
c
c
      do 400 j=2,jmtm1
c
c-----------------------------------------------------------------------
c       form curl of time change in vertically averaged equations
c-----------------------------------------------------------------------
c
#ifdef islands
c       all vorticity pts. are computed so that those needed for the
c       line integral of hole relaxation (immediately adjacent to
c       islands will be defined).
c
        is=2
        ie=imtm1
#else
c       non-active vorticity pts. are kept at 0 for convenience in relax
c       to accomplish this, ztd is computed only between given starting
c       and ending indices.
c
        do 100 l=1,lseg
          is=isz(j,l)
          if(is.eq.0) go to 101
          ie=iez(j,l)
#endif
          do 90 i=is,ie
            ztd(i,j) = ((zu(i,j)*dxu(i) + zu(i-1,j)*dxu(i-1))*csu(j  )
     $             -(zu(i,j-1)*dxu(i) + zu(i-1,j-1)*dxu(i-1))*csu(j-1))
            ztd(i,j)=(((zv(i,j) - zv(i-1,j))*dyu(j  )
     $              +(zv(i,j-1) - zv(i-1,j-1))*dyu(j-1)
     $              -ztd(i,j))*dxt2r(i)*dytr(j))*cstr(j)
90        continue
100     continue
101     continue
c
c-----------------------------------------------------------------------
c       do analysis of external mode forcing on diagnostic timesteps
c-----------------------------------------------------------------------
c
        if (diagts .and. eots) then
          engtmp = c0
          fx = cst(j)*dyt(j)/c2dtsf
          do 200 i=2,imtm1
            tempik(i,1,1) = - p(i,j,1)*ztd(i,j)*fx*dxt(i)
200       continue
          do 300 i=2,imtm1
            engtmp = engtmp +tempik(i,1,1)
300       continue
          engext(1) = engext(1) + engtmp
        endif
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c       filter vorticity at high latitudes
c-----------------------------------------------------------------------
c
        call filz (j)
#endif
400   continue
      return
      end


