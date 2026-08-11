c
c subroutine to calculate the number of latitude grid points needed
c This is the first of three subroutines for setting up the boundary
c conditions and the necessary constants in the ice code
c
      subroutine bnds
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "surfo.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c northern hemisphere
c
      jinit = 1
      je = iy
      js = 1
c
c southern hemisphere
c
      if (ihem .eq. 2) then
         jinit = il
         je = iy + 1
         js = -1
      end if
c
c cycle through to find the sea ice boundary
c
      jlimit = 0
      do j = jinit, je, js
      do i = 1, ix
c
c jl is counter of number of grid points 
c
         jl = j
         if (ihem .eq. 2) then
            jl = il - j + 1
         end if
         if (heff(i,j,1) .gt. 0.) then
            jlimit = jl
         end if
      end do
      end do
c
c set lower limit on ice boundary when no ice is present
c
      if (jlimit .lt. 1) then
         jlimit = 1
         print *, 'no ice in the hemisphere, ihem is ', ihem
      end if
c
#endif
      return
      end
