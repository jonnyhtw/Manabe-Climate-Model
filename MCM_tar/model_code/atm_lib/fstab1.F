      subroutine fstab1 (k, ivar)
#ifdef garratt
c
c introduces Richardson number dependence of  mixing coefficients 
c for vertical diffusion
c
c (based on recommendations of Garratt "The Atmospheric Boundary Layer",
c Cambridge University Press, 1992,  pp 245-246, as supplemented by
c Garratt, personal communication, 1993.
c
c input quantities:
c       ivar   = index for variable to be diffused
c                (1 = u-momentum
c                 2 = v-momentum
c                 3 = heat
c                 4 = water vapor
c          k   = vertical level
c
c output quantities:
c    stbadj(n) = array of adjustments
c
c
#include "parm.h"
c
#include "param.h"
c
#include "parmkp.h"
c
#include "mixlen.h"
c
#include "surfc.h"
c
#include "richrd.h"
c
      dimension x1(ix), x2(ix)
c
c  these are Garratt's recommendations for use in adjusting
c  above-surface-layer K values
      a = 10.0
      c = 0.5*a
      bm = 80.0*fgarr(k-1)*elmix(k-1)**2
      bt = 35.0*fgarr(k-1)*elmix(k-1)**2
c      print *, a, c, bm, bt
c      print *, (ri(i,k-1),i=1,ix)
c
      do 1 i=1,ix
c
        if (ri(i,k-1).lt.0.0) then
            x1(i)=a*ri(i,k-1)
            x2(i)=sqrt(-ri(i,k-1))
            if (ivar.le.2) then
c
c adjustment for momentum
                stbadj(i)=1-x1(i)/(1.+bm*x2(i))
              elseif (ivar.ge.3) then
c
c adjustment for heat and water vapor
                stbadj(i)=1-x1(i)/(1.+bt*x2(i))
            endif
          else
            stbadj(i)=1./((1.+c*ri(i,k-1))**2)
        endif
c
c set floor for adjustment
        stbadj(i) = amax1(stbadj(i),0.1)
c
    1 continue
c
c      print 1414, k, (stbadj(i),i=1,ix)
c 1414 format (' k= ',i3/4(12e10.4/))
c
#endif
      return
      end

