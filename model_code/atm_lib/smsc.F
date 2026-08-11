      subroutine smsc
c
c  called from: pass2
c  subroutines called: none
c
c  smsc computes the components of the ground hydrology, namely
c  soil moisture (soilm), snow depth (snwdpt), runoff (runoff),
c  snowmelt (snwmlt), and glacier melt (ycemlt).
c
#include "parm.h"
c
      dimension tauran(ix),tausno(ix)
c
#include "ctlblk.h"
c
#if defined cpldleg3 || defined cpldleg1
#include "gauss.h"
c
#endif
#include "grdfld.h"
c
#include "grid.h"
c
#include "intgrl.h"
c
#include "param.h"
c
#if defined cpldleg3 || defined cpldleg1
#include "qa.h"
c
#endif
#include "surfc.h"
c
      dimension tsnw(ix),tmps(ix),tmpw(ix)
#if defined cpldleg3 || defined cpldleg1
      dimension snopk(ix),runf(ix)
#endif
c
      equivalence (tauran(1),gpsdx(1,1)), (tausno(1),gpsdy(1,1))
c
c  set water field capacity of the soil. units are cm of liquid water.
      data wfc /15.0/
#if defined cpldleg3 || defined cpldleg1
c  set upper limit for snow depth, excess is glacier runoff.
c   units are cm.
      data snomax /20.0/
#endif
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
c
c  definitions of arrays used below:
c        tmps   = water equivalent net snow from current timestep
c               = (snowfall - snowmelt).
c        tmpw   = net liquid water from current timestep
c               = (rainfall + snowmelt water).
c        ycemlt = water equivalent accumulated glacier melt over land
c                 with no snow cover and glaciers.
c        smelts = water equivalent snowmelt or glacier melt for current
c                 timestep computed in subroutine tsc.
c        soilm  = accumulated soil moisture
c               = (soilm + rain + snowmelt water - evaporation);
c                 soilm .le. wfc for land points; soilm = wfc for sea.
c        snwdpt = water equivalent depth of snow
c               = (snowfall - snowmelt - sublimation).
c        runoff = accumulated excess soil moisture .gt. wfc over land.
c        snwmlt = water equivalent total accumulated snowmelt.
c        altnt  = evaporation(sublimation) in cm from subroutine vertdf.
#if defined cpldleg3 || defined cpldleg1
c        snwrnf = excess snow .gt. snomax over land
#endif
c
c  tsnw=1.0 if snow is present; tsnw=0.0 if no snow.
      do 1 i=1,ix
       if (snwdpt(i,jrow) .gt. 0.0) then
         tsnw(i) = 1.0
       else
         tsnw(i) = 0.0
       endif
    1 continue
c
      do 2 i=1,ix
      tmps(i)=tausno(i)-smelts(i,jrow)*tsnw(i)
      tmpw(i)=tauran(i)+smelts(i,jrow)*tsnw(i)
      soilm (i,jrow)=  soilm(i,jrow)+tmpw(i)-(1.0-tsnw(i))*altnt(i,jrow)
      soilm (i,jrow)=  soilm(i,jrow)*zin(i)+(1.0-zin(i))*wfc
      snwdpt(i,jrow)=(snwdpt(i,jrow)+tmps(i)-tsnw(i)*altnt(i,jrow))
     *                               *zin(i)
      snwmlt(i,jrow)= snwmlt(i,jrow)+smelts(i,jrow)*     tsnw(i) *zin(i)
      ycemlt(i,jrow)= ycemlt(i,jrow)+smelts(i,jrow)*(1.0-tsnw(i))*zin(i)
  2   continue
c  note: altnt in the above loop is really fluxw from vertdf.
c
#if defined cpldleg3 || defined cpldleg1
      do 10 i=1,ix
       runf (i) = 0.0
       snopk(i) = 0.0
   10 continue
#endif
c
      do 3 i=1,ix
       if (soilm(i,jrow) .gt. wfc) then
#if !defined cpldleg3 && !defined cpldleg1
         runoff(i,jrow) = runoff(i,jrow) + (soilm(i,jrow) - wfc)
#endif
#if defined cpldleg3 || defined cpldleg1
         runf  (i)      = soilm(i,jrow) - wfc
#endif
         soilm (i,jrow) = wfc
       endif
    3 continue
c
#if defined cpldleg3 || defined cpldleg1
      do 4 i=1,ix
       if (snwdpt(i,jrow) .gt. snomax) then
         snopk (i)      = snwdpt(i,jrow) - snomax
         snwdpt(i,jrow) = snomax
       endif
    4 continue
c
      do 44 i=1,ix
       snwrnf(i,jrow) = snwrnf(i,jrow) + snopk (i)
   44 continue
c
      do 5 i=1,ix
       runoff(i,jrow) = runoff(i,jrow) + runf(i)
    5 continue
c
      do 6 i=1,ix
       if (ijrnf(i,jrow) .ne. 0) then
         ir     = ijrnf(i,jrow)/100
         jr     = ijrnf(i,jrow)-100*ir
         lsr    = min0(jr,ilp-jr)
         rainwt = 0.5*weight(ls )*cosqr(ls)
         rainwr = 0.5*weight(lsr)*cosqr(lsr)
         wjrt   = rainwt/rainwr
         per(ir,jr) = per(ir,jr)+wjrt*runf (i)
         sno(ir,jr) = sno(ir,jr)+wjrt*snopk(i)
       endif
    6 continue
#endif
c
      return
      end
