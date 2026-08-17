c
c subroutine to calculate internal ice forces for cavitating fluid
c
      subroutine fint
#if defined ice_dynamics && defined icecpld
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
      dimension fx(ix,iy), fy(ix,iy)
c
c calculate internal ice force from ice pressure (cavitating fluid)
c
      do 10 j = 1, np
      do 10 i = 1, nl
         i0 = i
         if (i0.eq.1) i0=nl+1
c
c calculate latitude difference between centers of grid and the cosine
c
         dphi = abs(arlatc(j+1)-arlatc(j))
         cosp = cos(arlatc(j+1))
c
c note presur is ice pressure and gradient of presure is negative force
c
         fx(i,j)=-(presur(i+1,j+1)-presur(i0,j+1))/(R*dlam*cosp)
         fy(i,j)=-(presur(i+1,j+1)-presur(i+1,j))/(R*dphi)
  10  continue
c
      do 20 j = 1, np-1
      do 20 i = 1, nl
c
c linear interpolation from C-grid to grid center for y-direction force
c
         d1 = abs(arlat(j+1)-arlatc(j+1))*uvmy(i,j)
         d2 = abs(arlatc(j+1)-arlat(j))*uvmy(i,j+1)
         fnum = fy(i,j)*d1 + fy(i,j+1)*d2
         fden = d1 + d2
         if (fden .gt. 0.) then
            fyi(i,j) = fnum/fden
         else
            fyi(i,j) = 0.
         end if
c
c linear interpolation from C-grid to grid center for x-direction force
c note the longitudinal grid spacing is not varying
c
         i1 = i + 1
         if (i1.gt.nl) then i1=1
         d1 = uvmx(i1,j)
         d2 = uvmx(i,j)
         fnum = fx(i,j)*d2 + fx(i1,j)*d1
         fden = d1 + d2
         if (fden .gt. 0.) then
            fxi(i,j) = fnum/fden
         else
            fxi(i,j) = 0.
         end if
c
  20  continue
c
c zero out ice forces when no ice is present
c
      do 30 j = 1, np
      do 30 i = 1, nl
c        fxi(i,j) = fx(i,j)
c        fyi(i,j) = fy(i,j)
         if (hef(i,j,1) .le. 0.) then
            fxi(i,j) = 0.
            fyi(i,j) = 0.
         end if
c        print *, i, j, fxi(i,j), fyi(i,j), hef(i,j,1)
  30  continue
c     stop
c
#endif
      return
      end
