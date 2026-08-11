      subroutine uvphrm
c
c  uvphrm is called by atmos.  this subroutine converts spectral vortic-
c    ity and divergence to spectral velocity components times cosine
c    of the latitude (ucos and vcos) for the current time step and for
c    of the previous time step for use by vertdf.  spectral
c    log pressure gradients are also computed.
c
#include "parm.h"
#include "parmjx.h"
c
      complex zc(mx,jx,kx), zp(mx,jx,kx)
      complex vor, div, ucos, vcos, t, psdx, psdy, psln, gp, zs
      complex vorm, divm, tm, psm, vordt, divdt, tdt, psdt, ek
      complex ucosm, vcosm, udum, vdum
      complex iz
c
      data iz / (0.0, 1.0) /
c
      dimension gp(mx,jx,kx)
c
#include "dxordy.h"
c
#include "nxstps.h"
c
#include "sfharm.h"
c
#include "sftuvm.h"
c
      equivalence (gp(1,1,1), vcos (1,1,1))
      equivalence (zc(1,1,1), vordt(1,1,1)), (zp(1,1,1),divdt(1,1,1))
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
#ifdef triang
c
c  Zero ucos, etc. fields for triangular truncation computations
c
      do 772 k=1,kx
       do 772 n=1,jx
        do 772 m=1,mx
         ucos (m,n,k) = (0.0,0.0)
         vcos (m,n,k) = (0.0,0.0)
         ucosm(m,n,k) = (0.0,0.0)
         vcosm(m,n,k) = (0.0,0.0)
  772 continue
      do 773 n=1,jx
       do 773 m=1,mx
        psdx(m,n,1) = (0.0,0.0)
        psdy(m,n,1) = (0.0,0.0)
  773 continue
c
#endif
      if (nhem .le. 1) then
c
c  in a hemispheric model, zero all spectral velocity components and
c    log pressure gradients so that uncomputed spherical harmonic
c    components will not be left indefinite or incorrect.
c
        do 19 k=1,kx
          do 19 n=1,jx
            do 19 m=1,mx
              ucos(m,n,k) = (0.0,0.0)
              vcos(m,n,k) = (0.0,0.0)
   19   continue
        do 20 k=1,kx
          do 20 n=1,jx
            do 20 m=1,mx
              ucosm(m,n,k) = (0.0,0.0)
              vcosm(m,n,k) = (0.0,0.0)
   20   continue
        do 21 n=1,jx
          do 21 m=1,mx
            psdx(m,n,1) = (0.0,0.0)
            psdy(m,n,1) = (0.0,0.0)
   21   continue
c
c  if truncation limit jtrun is odd, zero even index jtrun1 values of
c    divergence which may not be initialized in hemispheric models to
c    prevent multiplication by indefinites in loops 3 and 9 below.
c
        if (mod(jx,2) .ne. 0) then
          do 17 k=1,kx
            do 17 m=1,mx
              div (m,jtrun1,k) = (0.0,0.0)
              divm(m,jtrun1,k) = (0.0,0.0)
   17     continue
        else
c
c  if truncation limit jtrun is even, zero odd index jtrun1 values of
c    vorticity which may not be initialized in hemispheric models to
c    prevent multiplication by indefinites in loops 8 and 13 below.
c
          do 23 k=1,kx
           do 23 m=1,mx
             vor (m,jtrun1,k) = (0.0,0.0)
             vorm(m,jtrun1,k) = (0.0,0.0)
   23     continue
        endif
      endif
c
c  compute zonal derivative of del square vorticity (zp) and compute
c    zonal derivative of del square divergence (zc) in spectral space.
c    (see november 4, 1981 notes.)
c
      do 1 k=1,kx
        do 1 n=jin,jx,jin
          do 1 m=1,mx
            zp(m,n,k) = uvdxm(m,n) * (iz * vor(m,n,k))
    1 continue
      do 2 k=1,kx
        do 2 n=1,jx,jin
          do 2 m=1,mx
            zc(m,n,k) = uvdxm(m,n) * (iz * div(m,n,k))
    2 continue
c
c  compute meridional derivatives of del squares of vorticity and diver-
c    gence (stream function and velocity potential respectively) at up-
c    per and lower limits of wave numbers and add to zonal derivatives
c    zc and zp to obtain ucos and vcos.
c
      do 8 k=1,kx
        do 8 m=1,mx
          ucos( m,1,k) = zc(m,1,k) - uvdyp(m, 1) * vor(m,   2  ,k)
          ucos(m,jx,k) =             uvdym(m,jx) * vor(m,jtrun1,k)
    8 continue
      if (jin .ne. 2) then
        do 3 k=1,kx
          do 3 m=1,mx
            vcos(m,1,k) = zp(m,1,k) + uvdyp(m,1) * div(m,2,k)
    3   continue
      endif
      do 33 k=1,kx
        do 33 m=1,mx
          vcos(m,jx,k) = -uvdym(m,jx) * div(m,jtrun1,k)
   33 continue
c
c  compute meridional derivatives of del squares of vorticity and diver-
c    gence by equation 11 of november 9, 1981 notes and add to zonal
c    derivatives to obtain ucos and vcos.
c
      do 4 k=1,kx
        do 4 n=2,jtrun1,jin
          do 4 m=1,mx
            vcos(m,n,k) = -uvdym(m,n) * div(m,n-1,k)
     *                  +  uvdyp(m,n) * div(m,n+1,k)
     *                  +  zp(m,n,k)
    4 continue
      do 11 k=1,kx
        do 11 n=ji,jtrun1,jin
          do 11 m=1,mx
            ucos(m,n,k) = uvdym(m,n) * vor(m,n-1,k)
     *                  - uvdyp(m,n) * vor(m,n+1,k)
     *                  + zc(m,n,k)
   11 continue
c
c  compute zonal derivative of del square of divergence for the previous
c    time step.
c
      do 12 k=1,kx
        do 12 n=1,jx,jin
          do 12 m=1,mx
            zc(m,n,k) = (iz * divm(m,n,k)) * uvdxm(m,n)
   12 continue
c
c  differentiate in the meridional and zonal directions log p* in spec-
c    tral space at meridional wave numbers at upper and lower limits.
c
      do 5 m=1,mx
        psdx(m,1,1)  = psdxm(m) * (iz * psln(m,1,1))
        psdx(m,jx,1) = 0.0
        psdy( m,1,1) = psdyp(m, 1) * psln(m,   2  ,1)
        psdy(m,jx,1) =-psdym(m,jx) * psln(m,jtrun1,1)
    5 continue
c
c  compute meridional and zonal derivatives of log p* in spectral space
c    for the other wave numbers.
c
      do 6 n=2,jtrun1,jin
        do 6 m=1,mx
          psdy(m,n,1) = -psdym(m,n) * psln(m,n-1,1)
     *                + psdyp(m,n) * psln(m,n+1,1)
    6 continue
      do 15 n=ji,jtrun1,jin
        do 15 m=1,mx
          psdx(m,n,1) = psdxm(m) * (iz * psln(m,n,1))
   15 continue
c
c  compute zonal derivative of del square of vorticity for the previous
c    time step.
c
      do 7 k=1,kx
        do 7 n=jin,jx,jin
          do 7 m=1,mx
            zp(m,n,k) = (iz * vorm(m,n,k)) * uvdxm(m,n)
    7 continue
c
c  compute ucosm and vcosm for the previous time
c    step by the same procedure used for computing ucos and vcos above.
c
c  first set ucosm and vcosm at upper and lower meridional limits.
c
      do 13 k=1,kx
        do 13 m=1,mx
          ucosm(m, 1,k) = zc(m,1,k) - uvdyp(m, 1) * vorm(m,   2  ,k)
          ucosm(m,jx,k) =             uvdym(m,jx) * vorm(m,jtrun1,k)
   13 continue
      if (jin .ne. 2) then
        do 9 k=1,kx
          do 9 m=1,mx
            vcosm(m,1,k) = zp(m,1,k) + uvdyp(m,1) * divm(m,2,k)
    9   continue
      endif
      do 99 k=1,kx
        do 99 m=1,mx
          vcosm(m,jx,k) = -uvdym(m,jx) * divm(m,jtrun1,k)
   99 continue
c
c  compute ucosm and vcosm at other meridional waves numbers by
c    equations 3 and 5 of the november 4, 1981 notes.
c
      do 14 k=1,kx
        do 14 n=ji,jtrun1,jin
          do 14 m=1,mx
            ucosm(m,n,k) = uvdym(m,n) * vorm(m,n-1,k)
     *                   - uvdyp(m,n) * vorm(m,n+1,k)
     *                   +    zc(m,n,k)
   14 continue
      do 10 k=1,kx
        do 10 n=2,jtrun1,jin
          do 10 m=1,mx
            vcosm(m,n,k) = -uvdym(m,n) * divm(m,n-1,k)
     *                   +  uvdyp(m,n) * divm(m,n+1,k)
     *                   +     zp(m,n,k)
   10 continue
c
      return
      end
