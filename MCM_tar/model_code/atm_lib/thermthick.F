c
c subroutine to put dynamic ice thickness and compactness forms
c to thermodynamic variables used in growth
c
      subroutine thermthick
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c northern hemisphere
c
      if (ihem .eq. 1) then
         jint = 1
         jend = np - 2
      end if
c
c southern hemisphere
c
      if (ihem .eq. 2) then
         jint = 2
         jend = np - 1
      end if
c
c load thickness and compactness (no outflow points) 
c
      do 29 j = jint, jend
      do 29 i = 1, nl
         inew = i  
c        if (i .eq. 1) then
c            inew = nl
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
c           k = il + 1 - j
            k = il - np + j + 1
         end if
c
         heff(i,k,1) = hef(inew,np+1-j,1)*heffm(inew,np+1-j)
         area(i,k,1) = are(inew,np+1-j,1)*heffm(inew,np+1-j)
         presu(i,k) = presur(inew,np+1-j)*heffm(inew,np+1-j)
c
c store ice velocities at grid center
c
         uicec(i,k) = uice(inew,np+1-j,2)*heffm(inew,np+1-j)
         vicec(i,k) = vice(inew,np+1-j,2)*heffm(inew,np+1-j)
c
c store internal ice force
c
         fxicei(i,k) = fxi(inew,np+1-j)*heffm(inew,np+1-j)
         fyicei(i,k) = fyi(inew,np+1-j)*heffm(inew,np+1-j)
c
  29  continue
c
#endif
      return
      end
