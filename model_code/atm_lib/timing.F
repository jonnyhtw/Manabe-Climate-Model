      subroutine timing
c
#include "parm.h"
c
c
c
      data alast/0.0/
c
      call getime(accume,elapsd)
      adelta = accume - alast
      write(6,1) accume,adelta
      alast = accume
c
    1 format(1x,'accumulated c.p.u. time since beginning of run',e12.5,
     *       25x,'elapsed c.p.u. time since last call',e12.5)
c
      return
      end
