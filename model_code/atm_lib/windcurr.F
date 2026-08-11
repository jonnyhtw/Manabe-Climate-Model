c
c subroutine to interpolate and/or convert wind stress and ocean currents
c to usable form for ice code
c
      subroutine windcurr
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
      parameter (imt=ix*2+2, jmt=il+2)
c
c location of ocean currents
c
#include "ascomi.h"
c
c input surface temperature to be used in land mask
c
#include "grdfld.h"
c
#include "varice.h"
c
#include "variceco.h"
c
      dimension seln(ix,il), seamku(imt,jmt)
      dimension txic(imt,jmt),tyic(imt,jmt)
      dimension uao(ix,il), vao(ix,il)
c
c cgsms is conversion from mks units to cgs unit
c
      data ccgsms /100./
c
c land-sea mask for intrpa subroutine 
c
      do 10 j = 1, il
      do 10 i = 1, ix
         if (ts(i,j) .gt. 0.0) then
            seln(i,j) = 1.
         else
            seln(i,j) = 0.
         end if
  10  continue
c
      do 20 j = 1, jmt
      do 20 i = 1, imt
         txic(i,j) = 0.
         tyic(i,j) = 0.
c
c input ocean currents and convert units from cm/s to m/s
c
         uao(i,j) = biceo(i,j,1)/ccgsms
         vao(i,j) = biceo(i,j,2)/ccgsms
c
  20  continue
c
c interpolate wind stresses to ocean (u,v) grid
c
      lgri = 1
      call intrpa(tauxic,txic,lgri,seln,seamku)
      call intrpa(tauyic,tyic,lgri,seln,seamku)
c
c shift, invert, and arrange the variables to match the input for ice code
c
      do 30 j = 1, il
      do 30 i = 1, ix
         inew = 2*i + 1
         jnew = il - j + 1
c
c wind stress components
c
         taxic(i,j) = txic(inew,jnew)
         tayic(i,j) = tyic(inew,jnew)
c
c ocean velocity components
c
         gwxic(i,j) = uao(inew,jnew)
         gwyic(i,j) = vao(inew,jnew)
c
  30  continue
c
c     print *, 2*ix+1,taxic(2*ix+1,2),taxic(2*ix+1,3),taxic(2*ix+1,3)
c     do i = 1, ix
c     print *,i, taxic(i,2),taxic(i,3),taxic(i,4)
c     print *,i, tauxic(i,2),tauxic(i,3),tauxic(i,4)
c     do j = 1, il
c     print *, j, ts(1,j),seln(1,j)
c     print *, j, tayic(1,j),tayic(2,j),tayic(3,j)
c     print *, j, tauyic(1,j),tauyic(2,j),tauyic(3,j)
c     do i = 1, ix
c     print *, i,j,taxic(i,j),tauxic(i,j),tayic(i,j),tauyic(i,j)
c     end do
c     end do
c     stop
c
#endif
      return
      end
