      subroutine cldint(cldobs, kindex, nl)
#ifdef fixcld
c
c  This subroutine computes the indexes for the
c  heights of the cloud tops
c  and bottoms for the fixed cloud model. 
c
c  This routine is called be cldht
c
#include "parm.h"
#include "parmncc.h"
#include "parmkp.h"
      parameter ( lx=ncc+1 )
      parameter (latobs=19)
c
#include "sgwgt1.h"
c
      dimension cldobs(10), kindex(latobs,nc)
      dimension cldlat(latobs)
c
c  Set the global mean atm pressure in millibars
      data psbar / 1013.25 /
c
c
c  Fill in Southern hemisphere cloud heights
      do 10 j=1,10
       cldlat(j) = cldobs(j)
   10 continue
      do 11 j=1,9
       cldlat(j+10) = cldobs(10-j)
   11 continue
c 
c  Start latitude loop to compute index at each latitude.
      do 20 j=1,latobs
c
c  Find first place where the pressure on the model level
c  is greater than the pressure of the cloud height.  Starting
c  form the top of the atm and going down the column.
       do 30 k=1,kx
        alevel = q(k)*psbar
        if (alevel .gt. cldlat(j)) go to 31
   30  continue
       stop 'cldint'
   31  continue
c  k is the index of the first model level below the cloud height.
c
c  compute the pressure half way between the model levels
       prsmid = (q(k)+q(k-1))*psbar*0.5   
c
c  If prsmid is greater than cldlat (cloud height) then the
c  level above is closer to cloud height, otherwise it is the
c  level below.
       if (prsmid .gt. cldlat(j)) then
         kindex(j,nl) = k-1
       else
         kindex(j,nl) = k
       endif
   20 continue
c
#endif
      return
      end
