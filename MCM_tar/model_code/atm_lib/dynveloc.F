c
c subroutine to separate ice thickness and compactness used in growth
c to two-hemisphere form for ice dynamic
c
      subroutine dynveloc
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c zero out ice variables
c
      do 9 j = 1, np
      do 9 i = 1, nl
      do 9 l = 1, 3
         uice(i,j,l) = 0.
         vice(i,j,l) = 0.
  9   continue
c
c load thickness and compactness
c
      do 19 j = 1, np-1
      do 19 i = 1, nl
      do 19 l = 1, 3
         inew = i
c           if (i .eq. nl-1) then
c              inew = 1
c           end if
c        if (i .eq. nl) then
c      	    inew = 1
c        end if
c
c northern hemisphere
c
         k = j
c
c southern hemisphere
c
         if (ihem .eq. 2) then
            k = il - np + j + 1
         end if
c
c ice velocities
c
         uice(i,np-j,l) = uic(inew,k,l)*uvmx(i,np-j)
         vice(i,np-j,l) = vic(inew,k,l)*uvmy(i,np-j)
c
         if (l .eq. 3) then
c
c wind stresses
c
            tauxin(i,np-j) = taxic(inew,k)*uvmx(i,np-j)
            tauyin(i,np-j) = tayic(inew,k)*uvmy(i,np-j)
c
c ocean currents
c
c           gwatx(i,np-j) = gwxic(inew,k)*uvmx(i,np-j)
c           gwaty(i,np-j) = gwyic(inew,k)*uvmy(i,np-j)
c
         end if
c
  19  continue
c
#endif
      return
      end
