      subroutine albcld(nd, crz)
c
c  This subroutine computes the cloud albedos.
c  This calculation is based on sigma and cloud thickness
c  in the old scheme
c  and sigma, cloud thickness  and latitude in the new scheme.
c     
c       Called from: impt
c       Calls      : nothing
c       Entry      : cldhtb
c
c
#include "parm.h"
c
#include "parmkp.h"
#include "parmjx.h"
#include "parmrad.h"
c
#include "cabd.h"
c
#include "cld.h"
c
#include "cldhts.h"
c
#include "cldals.h"
c
#include "param.h"
c
#include "sgwgt1.h"
c
c
#ifdef fixcld
c load high, middle and low cloud albedos.
      cca(2,1) = zza(2,1)
      cca(2,2) = zza(2,2)
      cca(3,1) = zza(3,1)
      cca(3,2) = zza(3,2)
      cca(4,1) = zza(4,1)
      cca(4,2) = zza(4,2)
#endif
#ifdef prdcld
c
c  Set the cloud albedos, include zenith angle dependence.
c  In the loops below, the cldhm and cldml are sigma levels in the atm.
c   they represent hieghts above the ground and are used to separate
c   low, middle and high clouds.  Cldhl and cldml are compute in 
c   the entry point below (cldhtb).
c  
# ifdef cldht60
      do 3 k=2,nd
       if (kth(k) .ne. kbh(k)) then
c        thick cloud
         cca(k,1) = zza(4,1)
         cca(k,2) = zza(4,2)
         ccb(k,2) = 0.30
       elseif (q(kth(k)) .le. cldhm) then
c        high cloud
         cca(k,1) = zza(2,1)
         cca(k,2) = zza(2,2)
         ccb(k,2) = 0.04
#  ifdef covcld
         cldsw(k) = crz
#  endif
       elseif (q(kth(k)).gt.cldhm .and. q(kth(k)).lt.cldml) then
c        middle cloud
         cca(k,1) = zza(3,1)
         cca(k,2) = zza(3,2)
         ccb(k,2) = 0.30
#  ifdef covcld
         cldsw(k) =  crz
#  endif
       elseif (q(kth(k)) .ge. cldml) then
c        low cloud
c        note: thick cloud and low cloud have the same radiative
c              properties.
         cca(k,1) = zza(4,1)
         cca(k,2) = zza(4,2)
         ccb(k,2) = 0.30
#  ifdef covcld
         cldsw(k) = crz
#  endif
       else
c        problems
         write(6,2005) k,kth,kbh
         stop 'albcld'
       endif
    3 continue
# endif
# ifdef cldht93
c
c  The main difference between the code below and that above is
c  in the assignment of thick cloud albedoes.
c  In the old code, if any two adjacent levels had cloud, it
c  was defined a thick cloud and it was assigned low cloud
c  properties.
c  In the new code, the thick cloud properties on are determined
c  by the height of the base of the thick cloud. So that if there
c  are two adjacent cirrus level clouds, they are still assigned
c  cirrus cloud properties.
c
        do 3 k=2,nd
         if (kth(k) .ne. kbh(k)) then
c          thick cloud
c          determine thick cloud base
c          assign cloud properties as a function of cloud base height
           if (q(kbh(k)) .le. cldhm) then
c            high cloud
             cca(k,1) = zza(2,1)
             cca(k,2) = zza(2,2)
             ccb(k,2) = 0.04
           elseif (q(kbh(k)).gt.cldhm .and. q(kbh(k)).lt.cldml) then
c            middle cloud
             cca(k,1) = zza(3,1)
             cca(k,2) = zza(3,2)
             ccb(k,2) = 0.30
           elseif (q(kbh(k)) .ge. cldml) then
c            low cloud
             cca(k,1) = zza(4,1)
             cca(k,2) = zza(4,2)
#  ifdef new_cld_absor
             ccb(k,2) = 0.40
#  else
             ccb(k,2) = 0.30
#  endif
           else
c            problems
             write(6,2005) k,kth,kbh
             stop 'impt'
           endif
         elseif (q(kth(k)) .le. cldhm) then
c          thin high cloud
           cca(k,1) = zza(2,1)
           cca(k,2) = zza(2,2)
           ccb(k,2) = 0.04
#  ifdef covcld
           cldsw(k) = crz
#  endif
         elseif (q(kth(k)).gt.cldhm .and. q(kth(k)).lt.cldml) then
c          thin middle cloud
           cca(k,1) = zza(3,1)
           cca(k,2) = zza(3,2)
           ccb(k,2) = 0.30
#  ifdef covcld
           cldsw(k) =  crz
#  endif
         elseif (q(kth(k)) .ge. cldml) then
c          thin low cloud
           cca(k,1) = zza(4,1)
           cca(k,2) = zza(4,2)
#  ifdef new_cld_absor
           ccb(k,2) = 0.40
#  else
           ccb(k,2) = 0.30
#  endif
#  ifdef covcld
           cldsw(k) = crz
#  endif
         else
c          problems
           write(6,2005) k,kth,kbh
           stop 'impt'
         endif
    3   continue
# endif
c
      return
     
c
 2005 format(1x,' problems in albcld at level ',i3,/,1x,
     *       'the cloud tops are: ',/,1x,30i4,/,1x,
     *       'the cloud bottoms are:',/,1x,30i4)
c
c
       entry cldhtb(alat)
c      ============

c cldhtb computes the cloud height boundaries for high,
c        middle and low clouds as a function of sigma
c        level.

c  Compute cloud height boundaries to determine low, middle, high
c  cloud properties.
# ifdef cldht60
c  These definitions were determined by Dick and Suki many years
c  ago. (approx. 1960)
      cldhm = 0.4
      cldml = 0.7
# endif
# ifdef cldht93
c  These defintions were determined by Suki, Dick, Mike and Ron
c  in 1993. The difference from those above is that they are
c  a function of latitude.
      cldhp = .7
      cldhe = .4
      cldmp = .85
      cldme = .7
      alatn = abs(alat)
      cldhm = cldhp + (90.0-alatn)*(cldhe-cldhp)/90.0
      cldml = cldmp + (90.0-alatn)*(cldme-cldmp)/90.0
ccc   write(6,4321) jrow, alatn, cldhm, cldml
c4321 format('albcld jrow, alatn, cldhm, cldml =',i5,4e16.8)
c
# endif
#endif
c
      return
      end
