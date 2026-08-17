      subroutine fstab2 (drag, ivar)
#ifdef garratt
c
c introduces Richardson number dependence of surface drag coefficient
c for vertical diffusion
c
c (based on recommendations of Garratt "The Atmospheric Boundary Layer",
c Cambridge University Press, 1992,  pp 243-244
c
c input quantities:
c       ivar   = index for variable to be diffused
c                (1 = u-momentum
c                 2 = v-momentum
c                 3 = heat
c                 4 = water vapor)
c    drag(n)   = array of neutral drag coefficient values
c
c
c   note:  with this formulation, drag is proportional to (Ri)**1/2 as
c     speed approaches zero, which is proportional to 1/u, 
c     so the flux ( C*u*(q-q0), say) remains finite;  this differs
c     from the "gustiness factor" often employed -- there is no need
c     to make umin as large as 1 m/s -- use something like 1.e-06.
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
      dimension drag(ix)
c
c  these are Garratt's recommendations when roughness lengths for 
c     heat and momentum are assumed equal
      a = 10.0
      c = 0.5*a
      bbm = 4.9
      bbt = 2.6

      bbb=bbt/bbm

      do 1 i=1,ix
        bm = drag(i)*a*bbm*exp(vk0/(2.*sqrt(drag(i))))
        bt = bbb*bm
c
        if (rib(i).lt.0.0) then
            x1=a*rib(i)
            x2=sqrt(-rib(i))
            if (ivar.le.2) then
c
c adjustment for momentum
                stbadj(i)=1-x1/(1.+bm*x2)
              elseif (ivar.ge.3) then
c
c adjustment for heat and water vapor
                stbadj(i)=1-x1/(1.+bt*x2)
            endif
          else
            stbadj(i)=1./((1.+c*rib(i))**2)
        endif
c
c set floor for adjustment
        stbadj(i) = amax1(stbadj(i),0.1)
c
    1 continue
c
#endif
      return
      end

