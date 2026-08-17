      subroutine nlmix(j)
      implicit real*8 (a-h,o-z)
#ifdef nlhmix
# ifdef microtasking
cfpp$ noconcur r
# endif
c=====================================================================
c
c     compute horizontal mixing coefficients based on Smagorinsky
c     non-linear viscosity as implemented and described in
c     Rosati & Miyakoda(jpo,vol 18,#11,1988).
c
#include "param.h"
#include "cnlmix.h"
#include "coord.h"
#include "cprnts.h"
#include "ctmngr.h"
#include "grdvar.h"
#include "slabs.h"
#include "switch.h"
c=======================================================================
c
c---------------------------------------------------------------------
c     compute tension rate of strain & shearing rate of strain
c     define dten,dshear at t row and u,v column
c---------------------------------------------------------------------
c
      fact1 = p5*cstr(j+1)
      fact2 = dytr(j+1)*cst(j+1)
      do 100 k=1,km
        do 90 i=2,imtm1
            dten(i,k) = ((u(i+1,k,jc,nm) + u(i+1,k,jp,nm))
     $                 - (u(i-1,k,jc,nm) + u(i-1,k,jp,nm)))
     $                 *fact1*dxu2r(i)
     $                 - (csur(j+1)*v(i,k,jp,nm) - csur(j)*v(i,k,jc,nm))
     $                 *fact2
          dshear(i,k) = ((v(i+1,k,jc,nm) + v(i+1,k,jp,nm))
     $                 - (v(i-1,k,jc,nm) + v(i-1,k,jp,nm)))
     $                 *fact1*dxu2r(i)
     $                 + (csur(j+1)*u(i,k,jp,nm) - csur(j)*u(i,k,jc,nm))
     $                 *fact2
90      continue
100   continue
#ifdef cyclic
      do 200 k=1,km
        dten(imt,k)   = dten(2,k)
        dshear(imt,k) = dshear(2,k)
        dten(1,k)     = dten(imtm1,k)
        dshear(1,k)   = dshear(imtm1,k)
200   continue
#endif
c
c---------------------------------------------------------------------
c     compute pure deformation
c---------------------------------------------------------------------
c
      do 300 k=1,km
        do 290 i=1,imt
          dpure(i,k) = sqrt(c2*(dten(i,k)**2 + dshear(i,k)**2))
290     continue
300   continue
c
c     compute effective anisentropic wavenumber of diffusing turbulence
c     for effective 2d isentropic wavenumber set effk & effky
c     to: (c28*sqrt(cs(j)*dxu(i)*dyu(j)))**2*sqrt2r
c
      sqrt2r   = c1/sqrt(c2)
      c28      = 0.28
      effky(j) = (c28*dyu(j))**2*sqrt2r
      do 400 i=1,imt
        effk(i) = (c28*csu(j)*dxu(i))**2*sqrt2r
400   continue
c
      if (diagts) then
c
c---------------------------------------------------------------------
c       print the solution at specified latitudes on diagnostic tsteps
c---------------------------------------------------------------------
c
        do 600 jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. j .and. prlat(jlat) .le. yt(jmt)) then
            do 590 k=1,km
              do 580 i=1,imt
                dkmc(i,k) = effk(i)*dpure(i,k)
580           continue
590         continue
c
c           use pure deformation for diff coeff
c
            hfact = c1
            bah   = c0
            do 570 k=1,km
              do 560 i=2,imt
                dpt1 = p25*(dpure(i,k)  + dpure(i-1,k)
     $              + dpures(i,k) + dpures(i-1,k))
                dkhcx(i,k) = (effk(i)*hfact*dpures(i,k) + bah)
                dkhcy(i,k) = (effky(j)*hfact*dpt1 + bah)
560           continue
570         continue
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
            scl = 1.e7
            write (stdout,9100) 'am ', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (dkmc(1,1), imt, is, ie, ks, ke, scl)
c
            scl = 1.e7
            write (stdout,9100) 'ahx', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (dkhcx(1,1), imt, is, ie, ks, ke, scl)
c
            scl = 1.e7
            write (stdout,9100) 'ahy', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (dkhcy(1,1), imt, is, ie, ks, ke, scl)
          endif
600     continue
      endif
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     $,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     $,', scaling=',1pg10.3)
c
#endif
      return
      end
