c
c subroutine to input geostrophic wind velocities (b-grid) in ice dynamic
c
      subroutine winds
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c     ihem = 2
c     do j = 1, np+1
c        print *, j, weightin(j), weightin(il-j+1)
c     end do
c     stop
      do 18 j = 1, np
      do 18 i = 1, nl
c
c gairx is meridonal wind speed, gairy is zonal wind speed (in m/sec)
c
         gairx(i,j) = 9.
         gairy(i,j) = 0.
  18  continue
c
c interpolate wind stress vectors from center points to corners
c
      do 19 j = 1, np
      do 19 i = 1, nl
c     do i = 1, nl
c
         i1 = i + 1
c
c cyclic boundary conditions
c
	 if (i1 .gt. nl) then
            i1 = 1
         end if
         j1 = j + 1
c
c calculate area weights
c
         arw = abs(arlatc(j+1)-arlat(j))
         arw1 = abs(arlat(j)-arlatc(j))
c
c re-order gaussian weights to match dynamic variables arrangement
c
         gnw = weightin(np+1-j)
c
         if (j .ne. np) then
            gnw1 = weightin(np+1-j-1)
c
c account for non-existing gaussian weight
c
	 else
            gnw1 = 0.
         end if
c
c re-order southern hemisphere weight arrangement
c
	 if (ihem .eq. 2) then
            if (j .ne. 1) then
               gnw = weightin(il+1-j+1)
            else
               gnw = 0.
            end if
            gnw1 = weightin(il+1-j)
         end if
c
         w11 = arw*gnw
       	 w12 = arw1*gnw1
         w21 = arw*gnw
         w22 = arw1*gnw1
c
         x11 = w11*tauxo(i,j)*heffm(i,j)
         x12 = w12*tauxo(i,j+1)*heffm(i,j+1)
         x21 = w21*tauxo(i+1,j)*heffm(i+1,j)
         x22 = w22*tauxo(i+1,j+1)*heffm(i+1,j+1)
c
         y11 = w11*tauyo(i,j)*heffm(i,j)
         y12 = w12*tauyo(i,j+1)*heffm(i,j+1)
         y21 = w21*tauyo(i+1,j)*heffm(i+1,j)
         y22 = w22*tauyo(i+1,j+1)*heffm(i+1,j+1)
c
         z11 = w11*heffm(i,j)
       	 z12 = w12*heffm(i,j+1)
         z21 = w21*heffm(i+1,j)
         z22 = w22*heffm(i+1,j+1)
c
         xall = x11 + x12 + x21 + x22
         yall = y11 + y12 + y21 + y22
         zall = z11 + z12 + z21 + z22
         pz = z11*z12*z21*z22
c
         if (pz .ne. 0.) then
            tauxin(i,j) = xall/zall
            tauyin(i,j) = yall/zall
         else
            tauxin(i,j) = 0.
            tauyin(i,j) = 0.
         end if
c
c     end do
c     print *, j, gnw, gnw1
  19  continue
c     stop
c
#endif
      return
      end
