c
c subroutine to input geostrophic ocean velocities (b-grid) in ice dynamic
c
      subroutine currs
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "variceco.h"
c
      do 18 j = 1, np
      do 18 i = 1, nl
c
c gwatx is meridonal ocean speed, gwaty is zonal ocean speed (in m/sec)
c
         gwatx(i,j) = 0.*heffm(i,j)
         gwaty(i,j) = 0.*heffm(i,j)
  18  continue
c
#endif
      return
      end
