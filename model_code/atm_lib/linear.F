      subroutine linear
c
c   called from: atmos
c   subroutines called: kernal and impvar, in implicit models only.
c   entry points: tderiv, robtq.
c
c   subroutine linear computes the horizontal diffusion terms and
c   adds them to the tendency terms in the equations of motion.
c   also, radiation tendencies are added to the temperature tendency.
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
c
      complex vor   , div  , t  , psln
      complex vorm  , divm , tm , psm
      complex vorp  , divp , tp , psp
      complex vordt , divdt, tdt, psdt
      complex ucos  , vcos ,psdx, psdy, geopot, zs, ek
      complex qmix  , qmixm, qmxdt, qp
#ifdef rad
      complex radten
#endif
c
c
      dimension geopot (mx,jx,kx)
      dimension vorp   (mx,jx,kx)
      dimension psp    (mx,jx,kx)
c
#include "ctlblk.h"
c
#include "dispat.h"
c
#ifdef consrv_ps
c
#include "dpsln.h"
c
#endif
c
#include "lincoef.h"
c
#include "nxstps.h"
c
#include "param.h"
c
#include "radten.h"
c
#include "sgwgt1.h"
c
#include "sfharm.h"
c
#include "tsmoth.h"
c
      dimension divp(mx,jx,kx), tp(mx,jx,kx)
      dimension qp(mx,jx,kx)
c
#include "water.h"
c
      equivalence (geopot(1,1,1), divp(1,1,1))
      equivalence (vorp  (1,1,1), divp(1,1,1))
      equivalence (psp   (1,1,1), divp(1,1,1))
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      if (implct) then
c
c  in implicit case perform geopotential build-up in grmult and add
c   zeroes to divdt
c
        do 1 k=1,kx
          do 1 n=1,jx
            do 1 m=1,mx
               geopot(m,n,k) = (0.0,0.0)
    1   continue
      else
c
c  in explicit case perform geopotential build-up in spectral domain
c
        do 2 n=1,jx
          do 2 m=1,mx
            geopot(m,n,kx) = grav * zs(m,n,1) + extkx * t(m,n,kx)
     *                     + extkxm * t(m,n,km)
    2   continue
        geopot(1,1,kx) = geopot(1,1,kx)
     *                 + cmplx(sqrt(2.0) * tfra * (extkx + extkxm),0.0)
        do 4 kk=1,km
          k = kp - kk
          do 3 n=1,jx
            do 3 m=1,mx
              geopot(m,n,k-1) = geopot(m,n,k)
     *                        + factq( k ) * t(m,n, k )
     *                        + factr(k-1) * t(m,n,k-1)
    3     continue
          geopot(1,1,k-1) = geopot(1,1,k-1)
     *                    + cmplx(sqrt(2.0) * tfra
     *                            * (factq(k) + factr(k-1)),0.0)
    4   continue
      endif
c
c  compute horizontal diffusion terms and add them to the tendencies.
c  see december 4 1981 notes: linear horizontal diffusion in
c  spectral model.
c
      do 13 k=1,kx
        do 13 n=1,jx,jin
          do 13 m=1,mx
            tdt(m,n,k) = tdt(m,n,k) + txy(m,n) * tm(m,n,k)
c
c   add del square of kinetic energy, ek, (see gekcos computed in
c   grmult) to the spectral divergence tendency, divdt.
c   del square of geopotential, geopot, is also added
c   to divdt in explicit models.
c
            divdt(m,n,k) = divdt(m,n,k)
     *                   +  elcf(m,n) * (0.5 * ek(m,n,k) +geopot(m,n,k))
     *                   +   vxy(m,n) * divm(m,n,k)
   13 continue
      do 14 k=1,kx
        do 14 n=jin,jx,jin
          do 14 m=1,mx
            vordt(m,n,k) = vordt(m,n,k) + vxy(m,n) * vorm(m,n,k)
   14   continue
        do 15 k=1,kx
          do 15 n=1,jx,jin
            do 15 m=1,mx
              qmxdt(m,n,k) = qmxdt(m,n,k) + qxy(m,n) * qmixm(m,n,k)
   15   continue
c
#ifdef sponge
c
c Radiation condition at top level of model.
c Intended to slow winds at top level.
c Time constant, dmptop, is stored in param.h
c and defined in PARMTR.
c
      do 40 n=1,jx,jin
        do 40 m=1,mx
           divdt(m,n,1) = divdt(m,n,1) - div(m,n,1)/dmptop
   40 continue
c
      do 41 n=jin,jx,jin
        do 41 m=1,mx
           vordt(m,n,1) = vordt(m,n,1) - vor(m,n,1)/dmptop
   41 continue
c
#endif
#ifdef rad
c  add radiation tendencies to tdt field
c
      do 16 k=1,kx
        do 16 n=1,jx,jin
          do 16 m=1,mx
             tdt(m,n,k) = tdt(m,n,k) + radten(m,n,k)
   16 continue
#endif
      if (implct) then
        if (resetm)  call kernal
        call impvar
      endif
      return
c
      entry tderiv
c     ===== ======
c
c  perform time extrapolation of spectral variables and do
c  robert time smoothing. see notes december 1 1981: leapfrog
c  time differencing and robert time smoothing performed in tderiv.
c  tderiv is called by atmos immediately after linear is called.
c
c  set single precision vble to get-a-round CRI slowdown for
c  double precision computations
      salpha = alpha
      som2al = om2al
c
      do 24 k=1,kx
        do 24 n=1,jx,jin
          do 24 m=1,mx
            divp(m,n,k) = divm(m,n,k) + dtx2 * divdt(m,n,k)
              tp(m,n,k) =   tm(m,n,k) + dtx2 *  tdt (m,n,k)
            divm(m,n,k) = salpha * (divm(m,n,k) + divp(m,n,k))
     *                  + som2al *   div(m,n,k)
              tm(m,n,k) = salpha *    tm(m,n,k)
     *                  + som2al *    t (m,n,k)
             div(m,n,k) = divp(m,n,k)
               t(m,n,k) =   tp(m,n,k)
              qp(m,n,k) =         qmixm(m,n,k) + dtx2  * qmxdt(m,n,k)
           qmixm(m,n,k) = salpha *qmixm(m,n,k) + som2al *qmix (m,n,k)
            qmix(m,n,k) = qp(m,n,k)
   24 continue
c
c
      do 26 k=1,kx
        do 26 n=jin,jx,jin
          do 26 m=1,mx
            vorp(m,n,k) = vorm(m,n,k) + dtx2 * vordt(m,n,k)
            vorm(m,n,k) = salpha * (vorm(m,n,k) + vorp(m,n,k))
     *                  + som2al * vor(m,n,k)
            vor (m,n,k) = vorp(m,n,k)
   26 continue
      do 27 n=1,jx,jin
        do 27 m=1,mx
          psp (m,n,1) = psm(m,n,1) + dtx2 * psdt(m,n,1)
          psm (m,n,1) = salpha * (psm(m,n,1) + psp(m,n,1))
     *                + som2al * psln(m,n,1)
          psln(m,n,1) = psp(m,n,1)
   27 continue
#ifdef consrv_ps
c
c  The following code conserves ln(ps) over 10**9 timesteps
c  or longer on 64 bit machines.
c
      apsdt = real(psdt(1,1,1))
c
      apsp = apsm + adtx2*apsdt
      apsm = alpha*(apsm+apsp) + om2al*apsln
c
      psm(1,1,1) = cmplx(sngl(apsm),0.0)
#endif
c
c  The following statement will force ps(1,1,1) to be conserved.
c  This is true for decadal time scales.
c
      psln(1,1,1) = psm(1,1,1)
c
      return
c
      entry robtq
c     ===== =====
c
c  robtq is called from atmos after the convective adjustment is done.
c  it completes the robert time smoothing of temperature and
c  mixing ratio. see december 1 1981 notes.
c
c  set single precision vble to get-a-round CRI slowdown for
c  double precision computations
      salpha = alpha

      do 30 k=1,kx
        do 30 n=1,jx,jin
          do 30 m=1,mx
            tm   (m,n,k) =    tm(m,n,k) + salpha *    t(m,n,k)
            qmixm(m,n,k) = qmixm(m,n,k) + salpha * qmix(m,n,k)
   30 continue
c
c
      return
      end
