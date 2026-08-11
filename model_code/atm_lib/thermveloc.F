c
c subroutine to store ice velocities for later time step
c
      subroutine thermveloc
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c limit on ice velocities
c
      specap = 0.00001
c
      do 19 j = 1, np-1
      do 19 i = 1, nl
      do 19 l = 1, 3
         inew = i 
c        if (i .eq. 1) then
c           inew = nl
c        end if
c        if (i .eq. 2) then
c           inew = nl
c        end if
c
c northern hemisphere
c
         k = j
c
c southern hemisphere
c
         if (ihem .eq. 2) then
c           k = il + 2 + 1 - j
            k = il - np + j + 1
         end if
c
         uic(i,k,l) = uice(inew,np-j,l)*uvmx(inew,np-j)
         vic(i,k,l) = vice(inew,np-j,l)*uvmy(inew,np-j)
c
c zero out small ice velocities
c
         uspe = abs(uic(i,k,l))
         vspe = abs(vic(i,k,l))
         if (uspe .lt. specap) then
            uic(i,k,l) = 0.
         end if
         if (vspe .lt. specap) then
            vic(i,k,l) = 0.
         end if
  19  continue
c
#endif
      return
      end
