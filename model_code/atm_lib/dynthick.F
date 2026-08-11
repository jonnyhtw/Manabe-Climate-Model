c
c subroutine to separate wind stress, ice thickness and compactness 
c to two-hemisphere form for ice dynamics
c
      subroutine dynthick
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
         tauxo(i,j) = 0.
         tauyo(i,j) = 0.
         hef(i,j,1) = 0.
         are(i,j,1) = 0.
         presur(i,j) = 0.
   9  continue
c
      jp = np
c
c avoid reading out of bound for southern hemisphere
c
      if (ihem .eq. 2) then
         jp = np - 1
      end if
c
c load thickness and compactness
c
      do 19 j = 1, jp
      do 19 i = 1, nl
         inew = i
c        if (i .eq. nl-1) then
c           inew = 1
c        end if
c        if (i .eq. nl) then
c           inew = 1
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
c wind stress vectors on ice
c
         tauxo(i,np+1-j) = tauxic(inew,k)*heffm(i,np+1-j)
         tauyo(i,np+1-j) = tauyic(inew,k)*heffm(i,np+1-j)
c
c thickness and compactness and ice internal pressure
c
         hef(i,np+1-j,1) = heff(inew,k,1)*heffm(i,np+1-j)
         are(i,np+1-j,1) = area(inew,k,1)*heffm(i,np+1-j)
         presur(i,np+1-j) = presu(inew,k)*heffm(i,np+1-j)
  19  continue
c
c cyclic boundary conditions
c
      do 29 j = 1, np
         tauxo(nl+1,j) = tauxo(1,j)
         tauyo(nl+1,j) = tauyo(1,j)
         hef(nl+1,j,1) = hef(1,j,1)
         are(nl+1,j,1) = are(1,j,1)
         presur(nl+1,j) = presur(1,j)
  29  continue
c
c set outside grid points
c
      do 39 i = 1, nl+1
         hef(i,1,1) = 0.
         hef(i,np+1,1) = 0.
         are(i,1,1) = 0.
         are(i,np+1,1) = 0.
         presur(i,1) = 0.
         presur(i,np+1) = 0.
  39  continue
c
#endif
      return
      end
