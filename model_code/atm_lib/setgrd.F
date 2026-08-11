      subroutine setgrd
#if defined cpldleg3 || defined cpldleg1 || defined ice_dynamics
c
#include "parm.h"
c
#include "param.h"
c
      parameter ( jh=(3-nhem)*il/2 )
      parameter ( imt=2*ix+2, jmt=il+2 )
c
      dimension oyh(jmt+1), wth(iy)
c
#include "gauss.h"
c
#include "intbc.h"
c
#include "astro.h"
c
      dimension oy(jmt)
c
c----------------------------------------------------------------------------
c
      write(6,11) weight
c
c   compute ocean weights
c
      do 10 j = 1, iy
        wta(   j  ) = weight(j)
        wta(il+1-j) = weight(j)
   10 continue
      do 20 j = 1, il
        wto(j+1) = wta(j)
   20 continue
      wto( 1 ) = wto(  2  )
      wto(jmt) = wto(jmt-1)
c
c   compute half latitudes for gaussian grid, equator to pole
c
      swt = 0.0
c
c   set equator
c
      oyh(jh+2) = 0.0
      do 200 j = 1, jh
        l   = jh - j + 1
        swt = swt + wta(l)
        if (swt .gt. 1.0) swt = 1.0
        wth(j) = asin(swt)
  200 continue
      do 210 j = 1, jh
        oyh(  j  +1) = -wth(jh-j+1)
        oyh(jmt+1-j) =  wth(jh-j+1)
  210 continue
      oyh(jmt+1) =  oyh( jmt ) + (oyh(jmt) - oyh(jmt-1))
      oyh(  1  ) = -oyh(jmt+1)
c
c   specify t,s gridpoints at midpoint of half latitudes
c
      do 50 j = 1, jmt
        oy (j) = 0.5 * (oyh(j) + oyh(j+1))
        oyv(j) = oyh(j+1)
   50 continue
c
      do 60 j = 1, il
        ay(j) = oy(j+1)
   60 continue
c
c   write latitude data
c
      write(6,61)
      do 70 j = jmt,1,-1
        oyd  = oy ( j ) * radian
        oyvd = oyh(j+1) * radian
        write(6,71) j, oy(j), oyv(j), oyd, oyvd
   70 continue
c
   11 format(/,1x,'latitude weights from setgrd',/,(1x,10f12.9))
   61 format(/,8x,'latitude',7x,'latitude',7x,'latitude',7x,'latitude',
     *       /,8x,'of t,s (rad.)',2x,'of u,v (rad.)',2x,'of t,s (deg.)',
     *         2x,'of u,v (deg.)')
   71 format(1x,i5,4e15.6)
#endif
c
      return
      end
