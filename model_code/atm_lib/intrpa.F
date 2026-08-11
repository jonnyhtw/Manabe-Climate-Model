      subroutine intrpa(ain,dout,lgrid,zin,seamku)
#if defined cpldleg3 || defined cpldleg1 || defined ice_dynamics
c
c   interpolate surface boundary data from atmos to ocean grids
c
#include "parm.h"
c
      parameter ( imt=ix*2+2, jmt=il+2 )
c
#include "intbc.h"
c
      dimension ain(ix,il), dout(imt,jmt), zin(ix,il)
      dimension seamku(imt,jmt)
      dimension aini(imt-1,il), zinmk(imt,il)
c
c  description of variables
c     ain......array interpolated from (atmos)
c     dout.....array interpolated to (ocean)
c     lgrid....=0 for t,s grid; =1 for u,v grid
c     zin......land sea mask for the ain grid (1.0=land)
c     seamku...land sea mask for the uv  grid (1.0=sea )
c     zinmk....land sea mask for the aini(ts) grid (1.0=land)
c
c.......................................................................
c
c  east-west interpolation from atmos grid to ocean grid.
c   assumes atmos gridbox is centered at Greenwich meridian,
c       and ocean gridbox interface   at Greenwich meridian.
c
c                    i=     1       2       3  .....
c           atmos grid: |   x   |   x   |   x   |......
c           ocean grid: | o | o | o | o | o | o |.......
c                    i=   2 | 3   4   5   6   7 .....
c                           |
c                           |<-- 0 degrees longitude
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      do 10 j=1,jmt
        do 10 i=1,imt
          dout(i,j) = 0.0
   10 continue
c
      do 20 j=1,il
        do 20 i=1,ix
          aini(2*i-1,j) = ain(i,j)
          aini(2*i  ,j) = ain(i,j)
   20 continue
c
      if (lgrid .ne. 0) then
        do 30 j=1,il
          do 30 i=1,ix
            zinmk(2*i-1,j) = zin(i,j)
            zinmk(2*i  ,j) = zin(i,j)
   30   continue
        do 40 j=1,il
          aini (imt-1,j) = ain(1,j)
          zinmk(imt-1,j) = zin(1,j)
   40   continue
      endif
c
c   interpolation for t,s points
c
      if (lgrid .eq. 0) then
c   reverse j for ocean grid
        do 50 j=1,il
          jr = il - j + 2
          do 50 i=1,ix*2
            dout(i+1,jr) = aini(i,j)
   50   continue
c
        return
      endif
c
c   interpolation for u,v points
c
      do 70 j=2,jmt-2
        do 80 l=1,il
          if (ay(l) .ge. oyv(j)) then
            js = l
            if (js .eq. 1) js = 2
            jn = js - 1
            go to 90
          endif
   80   continue
   90   continue
        c  = oyv(j) - ay(jn)
        d  = ay(js) - oyv(j)
        jn = il - jn + 1
        js = il - js + 1
        c  = c * wta(js)
        d  = d * wta(jn)
c   e-w interpolation for u,v assumes ocean grid is 2x atmos grid
        do 100 i=2,imt-1
          seamku(i,j) = 0.0
          cval1 = c * aini(i  ,js) * (1.0 - zinmk(i  ,js))
          cval2 = c * aini(i-1,js) * (1.0 - zinmk(i-1,js))
          dval1 = d * aini(i  ,jn) * (1.0 - zinmk(i  ,jn))
          dval2 = d * aini(i-1,jn) * (1.0 - zinmk(i-1,jn))
          cc1   = c *                (1.0 - zinmk(i  ,js))
          cc2   = c *                (1.0 - zinmk(i-1,js))
          dd1   = d *                (1.0 - zinmk(i  ,jn))
          dd2   = d *                (1.0 - zinmk(i-1,jn))
          tval1 = cc1 + cc2 + dd1 + dd2
          tval2 = cc1 * cc2 * dd1 * dd2
          if (tval2 .ne. 0.0) then
            dout  (i,j) = (cval1 + cval2 + dval1 + dval2) / tval1
            seamku(i,j) = 1.0
          endif
  100   continue
   70 continue
c
      do 110 i=2,imt-1
        seamku(i,jmt-1) = 0.0
  110 continue
#endif
c
      return
      end
