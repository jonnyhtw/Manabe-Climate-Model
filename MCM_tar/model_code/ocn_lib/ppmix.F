      subroutine ppmix (j)
      implicit real*8 (a-h,o-z)
#ifdef ppvmix
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c    compute vertical mixing coefficients based on the pacanowski &
c    philander richardson mixing scheme (jpo vol 11, #11, 1981).
c    note: this parameterization was designed for equatorial models
c    and may not do a good job in mid or high latitudes.
c
c     inputs:
c
c     j    = current central slab row
c
c     outputs:
c
c     vertical mixing coefficient arrays "vdc,vvc" in "cvmix.h"
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "cppmix.h"
#include "cprnts.h"
#include "ctask.h"
#include "ctmngr.h"
#include "cvmix.h"
#include "grdvar.h"
#ifdef skipland
#include "index.h"
#endif
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#include "temp.h"
c
      if (j .eq. jsrow) then
c
c-----------------------------------------------------------------------
c       compute density difference across bottoms of "t" boxes for row
c       "jsrow" at time level tau-1.
c-----------------------------------------------------------------------
c
c       compute vertical difference of density.
c
        do 100 ks=1,2
          call statec (t(1,1,jc,nm,1), t(1,1,jc,nm,2), tempik(1,1,1),ks)
          do 90 k=ks,kmm1,2
            do 80 i=1,imt
              drhoc(i,k) = tempik(i,k,1) - tempik(i,k+1,1)
80          continue
90        continue
100     continue
c
c       set no-flux density across bottom
c
        do 200 i=1,imt
          if (kmt(i,j) .ne. 0) then
            drhoc(i,kmt(i,j)) = c0
          endif
200     continue
      endif
c
c-----------------------------------------------------------------------
c     compute richardson number at bottoms of "t" boxes for row j+1
c     at time level tau-1.
c-----------------------------------------------------------------------
c
c     compute vertical difference of density.
c
      do 300 ks=1,2
        call statec (t(1,1,jp,nm,1), t(1,1,jp,nm,2), tempik(1,1,1), ks)
        do 290 k=ks,kmm1,2
          do 280 i=1,imt
            drhop(i,k) = tempik(i,k,1) - tempik(i,k+1,1)
280       continue
290     continue
300   continue
c
c     set no-flux density across bottom
c
      do 400 i=1,imt
        if (kmt(i,j+1) .ne. 0) then
          drhop(i,kmt(i,j+1)) = c0
        endif
400   continue
c
c     compute richardson number on "u" grid box bottoms
c
      fx = -p25*grav
      do 500 k=1,kmm1
        t1 = fx*dzw(k)
        do 470 i=2,imtm1
          riuc(i,k) = t1*gm(i,k+1,jc)*(drhop(i,k) + drhop(i+1,k)
     $               + drhoc(i,k) + drhoc(i+1,k))
     $                /((u(i,k,jc,nm) - u(i,k+1,jc,nm))**2
     $               +  (v(i,k,jc,nm) - v(i,k+1,jc,nm))**2 + epsln)
470     continue
#ifdef cyclic
        riuc(1,k)   = riuc(imtm1,k)
        riuc(imt,k) = riuc(2,k)
#else
        riuc(1,k)   = c0
        riuc(imt,k) = c0
#endif
500   continue
c
      if (j .gt. jsrow) then
c
#ifdef skipland
        do 1400 lstrip=1,lseg
          istart = ist(j,lstrip)
          iend   = iet(j,lstrip)
          if (istart .ne. 0) then
#endif
c
c-----------------------------------------------------------------------
c           compute vertical viscosity coeff on "u" grid box bottoms
c-----------------------------------------------------------------------
c
            do 700 k=1,kmm1
              do 690 i=istart,iend
                t1       = c1/(c1 + c5*riuc(i,k))
                vvc(i,k) = fricmx*t1**2 + bvvc
690           continue
700         continue
c
c           in regions of gravitational instability, reset the vertical
c           mixing coefficient
c
            do 800 k=1,kmm1
              do 790 i=istart,iend
                if (riuc(i,k) .lt. c0) vvc(i,k) = vvclim
790           continue
800         continue
c
c-----------------------------------------------------------------------
c           compute richardson number at bottoms of boxes on "t" grid
c           by averaging
c-----------------------------------------------------------------------
c
            do 900 k=1,kmm1
              do 890 i=istart,iend
                tempik(i,k,1) = p25*(riuc(i,k) + riuc(i-1,k)
     $                             + rium(i,k) + rium(i-1,k))
890           continue
900         continue
c
c-----------------------------------------------------------------------
c           compute vertical diffusivity coeff on "t" grid box bottoms
c-----------------------------------------------------------------------
c
            do 1000 k=1,kmm1
              do 990 i=istart,iend
                t1       = c1/(c1 + c5*tempik(i,k,1))
                vdc(i,k) = fricmx*t1**3 + bvdc
990           continue
1000        continue
c
c           in regions of gravitational instability, reset the vertical
c           mixing coefficient to the limit
c
            do 1100 k=1,kmm1
              do 1090 i=istart,iend
                if (tempik(i,k,1) .lt. c0) vdc(i,k) = vdclim
1090          continue
1100        continue
c
c-----------------------------------------------------------------------
c           approximation for high freq wind mixing near the surface
c-----------------------------------------------------------------------
c
            do 1200 i=istart,iend
              if (vvc(i,1) .lt. wndmix) vvc(i,1) = wndmix
              if (vdc(i,1) .lt. wndmix) vdc(i,1) = wndmix
1200        continue
c
c-----------------------------------------------------------------------
c           set viscosity & diffusivity to zero in land
c-----------------------------------------------------------------------
c
            do 1390 k=1,kmm1
              do 1380 i=istart,iend
                vvc(i,k) = gm(i,k+1,jc)*vvc(i,k)
                vdc(i,k) = fm(i,k+1,jc)*vdc(i,k)
1380          continue
1390        continue
c
c           set viscosity & diffusivity at level km
c
            do 1300 i=istart,iend
              vvc(i,km) = c0
              vdc(i,km) = c0
1300        continue
c
#ifdef skipland
          endif
1400    continue
#endif
      endif
c
c-----------------------------------------------------------------------
c     transfer richardson numbers for use with next row
c-----------------------------------------------------------------------
c
      do 1501 k=1,km
        do 1490 i=1,imt
          drhoc(i,k) = drhop(i,k)
          rium(i,k) = riuc(i,k)
1490    continue
1501  continue
c
      if (diagts .and. eots) then
c
c---------------------------------------------------------------------
c       print the solution at specified latitudes on diagnostic tsteps
c---------------------------------------------------------------------
c
        do 6200 jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. j .and. prlat(jlat) .le. yt(jmt)) then
            is  = indp (prslon(jlat), xt, imt)
            ie  = indp (prelon(jlat), xt, imt)
            ks  = indp (prsdpt(jlat), zt, km)
            ke  = indp (predpt(jlat), zt, km)
            fx  = 1.0e-2
            scl = c1
            write (stdout,9100) 'vvc ', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (vvc(1,1), imt, is, ie, ks, ke, scl)
c
            scl = c1
            write (stdout,9100) 'vdc ', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (vdc(1,1), imt, is, ie, ks, ke, scl)
          endif
6200    continue
      endif
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     $,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     $,', scaling=',1pg10.3)
#endif
      return
      end
