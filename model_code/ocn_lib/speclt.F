      subroutine speclt (dyt)
      implicit real*8 (a-h,o-z)
#if defined ocnr15 || defined ocnr30
c
c=======================================================================
c     this subroutine is taken from the GFDL Climate Group's atmospheric
c     model.  It combines elements of the coupled model's "gaussg" and
c     "setgrd" subroutines to calculate the gaussian half-latitudes,
c     which should correspond to u,v gridpoint locations.
c
c    varables:
c      nzero = number of points between pole and equator
c      f     = cosine of the colatitudes [radians]
c      wts   = weights for the gaussian quadrature
c      alat  = the colatitudes [degrees]
c=======================================================================
c
#include "param.h"
c
c------------------------------------------------------
csun  implicit double precision (a-h,o-z)
csun  real dyt
c------------------------------------------------------ 
c
      parameter ( nzero  = jmtm2/2 )
      parameter ( nzermx = 100 )
c
      dimension f(nzero), wt(nzero)
c
      dimension dyt(jmt)              
c
      dimension wth(jmt), wta(jmtm2), wto(jmt)              
      dimension oyh(jmtp1), oyv(jmt)
c
c=======================================================================
c
      if (nzero .gt. nzermx) then
        write (stdout,99)
        write (stderr,99)
        stop ' gaussg '
      endif
c
c      set local constants
c
csun(precision is machine dependent
      xlim  = 1.0e-12
      xlim  = 1.0e-10
csun  xlim  = 1.0e-6
      ir    = nzero + nzero
      fi    = ir
      fi1   = fi + c1 
      pi    = c4 * atan(c1)
      pir2  = p5 * pi
      fn    = pir2 / nzero
c
      do 100 i=1,nzero
        wt(i) = i - p5 
100   continue
c
      do 200 i=1,nzero
        f(i) = sin(wt(i) * fn + pir2 )
200   continue
c
      dn  = fi  / sqrt(c4*fi*fi  - c1)
      dn1 = fi1 / sqrt(c4*fi1*fi1 - c1)
      a   = dn1 * fi
      b   = dn  * fi1
      irp = ir + 1
      irm = ir - 1
c
      do 300 i=1,nzero
        kount = 0
c
201     continue
        kount = kount + 1
        call ordleg (g , f(i), ir )
        call ordleg (gm, f(i), irm)
        call ordleg (gp, f(i), irp)
        gt    = (f(i) * f(i) - c1 ) / (a * gp - b * gm)
        ftemp = f(i) - g * gt
        gtemp = f(i) - ftemp
        f(i)  = ftemp
        if (abs(gtemp) .gt. xlim) go to 201
300   continue
c
      do 400 i=1,nzero
        a = c2 * (c1 - f(i)*f(i))
        call ordleg (b, f(i), irm)
        b     = b * b * fi * fi
        wt(i) = a * (fi - p5 ) / b
400   continue
c
c   (the following code is adapted from the climate model's
c   "setgrd" subroutine)
c
c   compute ocean weights
c  
      do 500 j=1,nzero
        wta(j) = wt(j)
        wta(jmtm2+1-j) = wt(j)
500   continue           
c
      DO 600 J=1,jmtm2   
        wto(j+1) = wta(j)
600   continue           
c
      wto( 1 ) = wto(  2  )
      wto(jmt) = wto(jmtm1)
c
c   compute half latitudes for gaussian grid, equator to pole
c   (these should be the same as the location of the u,v points)
c
      swt = c0 
c
c   set equator
c
      oyh(nzero+2) = c0  
c
      do 700 j = 1, nzero
        l   = nzero - j + 1   
        swt = swt + wta(l)
        if (swt .gt. c1 ) swt = c1   
        wth(j) = asin(swt)
700   continue           
c
      DO 800 j = 1, nzero  
        oyh(  j  +1) = -wth(nzero-j+1) 
        oyh(jmtp1-j) =  wth(nzero-j+1)
800   continue           
c
      oyh(jmtp1) =  oyh( jmt ) + (oyh(jmt) - oyh(jmtm1))  
      oyh(  1  ) = -oyh(jmtp1)           
      deg = 180.0 / pi         
C                  
      do 900 j = 1, jmt
        oyv(j) = oyh(j+1) * deg
900   continue           
c      
      do 1000 j=2,jmt
        dyt(j) = oyv(j) - oyv(j-1)
1000  continue
c
      dyt(1) = dyt(2)
c
   99 format(/ ' error => execution halted in gaussg'/
     $      ' increase size of array maximum'/)
c
#endif
      return
      end
