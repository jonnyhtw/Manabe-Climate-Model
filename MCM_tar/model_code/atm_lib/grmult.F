      subroutine grmult
c
c   called from: passa
c   subroutines called: vertdf (or vrtdfb for boundary layer models)
c
c   this subroutine computes in grid domain the nonlinear terms (such as
c   ut,vt,etc.) in the equations of motion, the thermodynamic equation
c   and the water vapor equation in moist models. vertical velocity and
c   log pressure tendency are computed as well as temperature and mixing
c   ratio tendencies due to adiabatic vertical motions. horizontal ad-
c   vection terms for momentum, heat and water vapor are also computed.
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
c
      complex ps
c
c
      dimension gt    (ix,kx)
      dimension ushear(ix,kp), vshear(ix,kp), tshear(ix,kp)
      dimension qshear(ix,kp), gdumk (ix,kx), gstore(ix,kx)
      dimension scnvr (ix,kx), cmvft (ix   ), cmvfrd(ix   )
      dimension gubar (ix   ), gvbar (ix   ), gtat0 (ix   )
      dimension gdivbr(ix   ), gtat1 (ix   )
c
c
#include "auxatm.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grid.h"
c
#include "gridw.h"
c
#include "guvm.h"
c
#include "intgrl.h"
c
#include "mrstr.h"
c
#include "param.h"
c
#include "pstar.h"
c
#include "sgwgt1.h"
c
#include "vrtcom.h"
c
#include "wimplt.h"
c
#include "zsrvar.h"
c
c   c4a=7/12, c4b=1/12, weights for mixing ratio 4th order differencing
c
      data c4a / 0.58333333 /, c4b / 0.083333333 /
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
#ifdef dirnal
      do 101 i=1,ix
       gwdot(i,kp) = 0.0
  101 continue
#endif
      secsqr      = 1.0 / cosqr(ls)
      costh       =  sqrt(cosqr(ls))
c
c
c   convert model celsius temperatures to kelvin
c
      do 5 k=1,kx
       do 5 i=1,ix
        gt(i,k) = cgt(i,k,ih) + tfra
    5 continue
c
c   compute grid values of vertical mean velocity (u,v) and divergence
c
      do 8 i=1,ix
       gtat0 (i) = -58.0
       gtat1 (i) = +15.0
       gubar (i) = 0.0
       gvbar (i) = 0.0
       gdivbr(i) = 0.0
  8   continue
      do 9 k=1,kx
       do 9 i=1,ix
        gubar (i) = gubar (i) + dq(k) * gucos(i,k,ih)
        gvbar (i) = gvbar (i) + dq(k) * gvcos(i,k,ih)
        gdivbr(i) = gdivbr(i) + dq(k) * gdiv (i,k,ih)
    9 continue
c
c   compute vertical shear of u,v,t for vertical advection
c
      do 10 k=2,kx
       do 10 i=1,ix
        ushear(i,k) = gucos(i,k,ih) - gucos(i,k-1,ih)
        vshear(i,k) = gvcos(i,k,ih) - gvcos(i,k-1,ih)
        tshear(i,k) = cgt  (i,k,ih) - cgt  (i,k-1,ih)
   10 continue
      do 11 i=1,ix
       ushear(i,1) = 0.0
       vshear(i,1) = 0.0
       qshear(i,1) = 0.0
c   tshear(i,1) is multiplied by gwdot(i,1)=0. below
       tshear(i,1)  = 0.5  * (cgt  (i, 1,ih) - gtat0(i))
       ushear(i,kp) = 0.25 * (gucos(i,kx,ih) - gucos(i,km,ih))
       vshear(i,kp) = 0.25 * (gvcos(i,kx,ih) - gvcos(i,km,ih))
       qshear(i,kp) = 0.0
c   tshear(i,kp) is multiplied by gwdot(i,kp)=0. below
       tshear(i,kp) = 0.5 * (gtat1(i) - cgt(i,kx,ih))
   11 continue
c
c
c   compute mean mixing ratio for vertical advection
      do 6 i=1,ix
        qshear(i, 2) = 0.5 * (gqmix(i, 1,ih) + gqmix(i, 2,ih))
        qshear(i,kx) = 0.5 * (gqmix(i,km,ih) + gqmix(i,kx,ih))
    6 continue
      do 7 k=3,km
        do 7 i=1,ix
          qshear(i,k) = c4a * (gqmix(i,k  ,ih) + gqmix(i,k-1,ih))
     *                - c4b * (gqmix(i,k+1,ih) + gqmix(i,k-2,ih))
  7   continue
c
c   compute explicit tendency of log pstar, convert log pstar to pstar
c
      do 12 i=1,ix
       gpseq (i,ih) = -gubar(i) * gpsdx(i,ih) - gvbar(i) * gpsdy(i,ih)
       gvbrps(i,ih) = exp(gps(i,ih))
 12   continue
c
c   compute log pressure gradient times wind component deviations from
c   vertical mean for later use below
c
      do 14 k=1,kx
        do 14 i=1,ix
          gdumk(i,k) = gpsdx(i,ih) * (gucos(i,k,ih) - gubar(i))
     *               + gpsdy(i,ih) * (gvcos(i,k,ih) - gvbar(i))
 14   continue
c
c   compute sigma dot vertical velocity
c
      do 15 k=1,km
        do 15 i=1,ix
          gwdot(i,k+1) = gwdot(i,k) - dq(k) * (gdumk(i,k) * secsqr
     *                 - gdivbr(i) + gdiv(i,k,ih))
 15   continue
c
#ifndef impcit
c
c   compute total tendency of log pressure in non-implicit models
      do 16 i=1,ix
        gpseq(i,ih) = gpseq(i,ih) - gdivbr(i) * cosqr(ls)
   16 continue
c
c   use gstore for temporary storage array.  gstore is redefined below
c
      do 17 k=1,kx
        do 17 i=1,ix
          gstore(i,k) = gdumk(i,k) - gdivbr(i) * cosqr(ls)
   17 continue
c
#endif
#ifdef dirnal
c
c   use gstore for temporary storage array.  gstore is redefined below
c
      do 17 k=1,kx
        do 17 i=1,ix
          gstore(i,k) = gdumk(i,k) - gdivbr(i) * cosqr(ls)
   17 continue
c
c   compute omega (dp/dt) vert velocity for diurnal model analysis data
c
      do 23 k=1,kx
       do 23 i=1,ix
        omg(i,k,ih) = gvbrps(i,ih)*(0.5*(gwdot(i,k) + gwdot(i,k+1))
     *              + q(k)*secsqr*gstore(i,k))
   23 continue
#endif
c
c
c   compute longitude-height slabs of the following quantities:
c    agrid,bgrid:  components of horizontal momentum flux
c    gekcos:       kinetic energy
c    gtu,gtv:      components of horizontal heat flux
c    gqmxeq:       tendency of mixing ratio from vertical transport only
c    gqmixu,gqmixv:components of horizontal moisture flux
 20   continue
      do 25 k=1,kx
        do 25 i=1,ix
          agrid(i,k,ih) = gvor(i,k,ih) * gucos(i,k,ih)
     *                  + rgas * gt(i,k) * gpsdy(i,ih)
     *                  + odq(k) * (gwdot(i,k+1) * vshear(i,k+1)
     *                  + gwdot(i,k) * vshear(i,k))
          bgrid(i,k,ih) = gvor(i,k,ih) * gvcos(i,k,ih)
     *                  - rgas * gt(i,k) * gpsdx(i,ih)
     *                  - odq(k) * (gwdot(i,k+1) * ushear(i,k+1)
     *                  + gwdot(i,k) * ushear(i,k))
 25   continue
      do 28 k=1,kx
        do 28 i=1,ix
          gekcos(i,k,ih) = gucos(i,k,ih)**2 + gvcos(i,k,ih)**2
          gtu   (i,k,ih) = cgt(i,k,ih)      * gucos(i,k,ih)
          gtv   (i,k,ih) = cgt(i,k,ih)      * gvcos(i,k,ih)
 28   continue
c
      do 29 k=1,kx
        do 29 i=1,ix
        gqmxeq(i,k,ih)=-2.*odq(k)*(gwdot(i,k+1)*qshear(i,k+1)
     *             -gwdot(i,k)*qshear(i,k)-gqmix(i,k,ih)*(gwdot(i,k+1)
     *             -gwdot(i,k)))+gdiv(i,k,ih)*gqmix(i,k,ih)
   29 continue
      do 26 k=1,kx
        do 26 i=1,ix
          gqmixu(i,k,ih) = gqmix(i,k,ih) * gucos(i,k,ih)
          gqmixv(i,k,ih) = gqmix(i,k,ih) * gvcos(i,k,ih)
   26 continue
c
#ifndef impcit
c   compute conversion and vertical contribution to temperature tendency
c   for explicit case
c
      do 30 k=1,kx
        do 30 i=1,ix
          scnvr(i,k)=dq(k)*(cappa*gstore(i,k)*secsqr+akapsg(k)*
     *          (gwdot(i,k)+gwdot(i,k+1)))*gt(i,k)*gvbrps(i,ih)
          gteq(i,k,ih)=cappa*gt(i,k)*gstore(i,k)
     *            +cosqr(ls)*(cgt(i,k,ih)*gdiv(i,k,ih)
     *            +akapsg(k)*gt(i,k)*(gwdot(i,k)+gwdot(i,k+1))-odq(k)*
     *            (gwdot(i,k+1)*tshear(i,k+1)+gwdot(i,k)*tshear(i,k)))
 30   continue
#endif
#ifdef impcit
c   compute conversion and vertical contribution to temperature tendency
c   and correct agrid and bgrid for implicit case
c
      do 37 k=1,kx
        do 37 i=1,ix
          gstore(i,k) = gt(i,k) - tmean(k)
   37 continue
      do 32 k=1,km
        do 32 i=1,ix
          tshear(i,k+1) = gstore(i,k+1) - gstore(i,k)
          gwexpt(i,k+1) = gwexpt(i,k) - secsqr * dq(k) * gdumk(i,k)
   32 continue
      do 33 k=1,kx
        do 33 i=1,ix
          agrid(i,k,ih)=agrid(i,k,ih)-rgas*tmean(k)*gpsdy(i,ih)
          bgrid(i,k,ih)=bgrid(i,k,ih)+rgas*tmean(k)*gpsdx(i,ih)
          gteq(i,k,ih)=cappa*gt(i,k)*gdumk(i,k)+cosqr(ls)*(cgt(i,k,ih)
     *           *gdiv(i,k,ih)+akapsg(k)*(gstore(i,k)*(gwdot(i,k+1)
     *           +gwdot(i,k))+tmean(k)*(gwexpt(i,k+1)+gwexpt(i,k)))
     *           -cappa*gstore(i,k)*gdivbr(i)-odq(k)*(gwdot(i,k+1)
     *           *tshear(i,k+1)+gwdot(i,k)*tshear(i,k)
     *           +gwexpt(i,k+1)*tmnshr(k+1)+gwexpt(i,k)*tmnshr(k)))
          scnvr(i,k)=dq(k)*(cappa*(gdumk(i,k)*secsqr-gdivbr(i))
     *       +akapsg(k)*(gwdot(i,k)+gwdot(i,k+1)))*gt(i,k)*gvbrps(i,ih)
 33   continue
#endif
c
c   compute zonal mean of conversion before subr. vertdf uses /vtemp/
c
      zmcnvr(jrow) = 0.0
      do 70 k=1,kx
       do 70 i=1,ix
        zmcnvr(jrow) = zmcnvr(jrow) + scnvr(i,k)
   70 continue
c
c   add contributions for vertical diffusion
c
c
      call vertdf
c
#ifdef dirnal
      do 102 i=1,ix
       gwdot(i,kp) = gvbrps(i,ih)
  102 continue
#endif
c
      do 41 k=1,kx
        do 41 i=1,ix
          agrid(i,k,ih) = agrid(i,k,ih) - vfv(i,k) * costh
          bgrid(i,k,ih) = bgrid(i,k,ih) + vfu(i,k) * costh
          gteq (i,k,ih) = gteq (i,k,ih) + vft(i,k) * cosqr(ls)
 41   continue
c
      do 46 k=kv,kx
        do 46 i=1,ix
          gqmxeq(i,k,ih) = gqmxeq(i,k,ih) + vfr(i,k)
 46   continue
c
c   compute zonal means of t(9), u(3), vft, vfr for printing by zmint
c
      zmt9(jrow) = 0.0
      zmu3(jrow) = 0.0
      do 75 i=1,ix
       zmt9(jrow) = zmt9(jrow) + rix*gt(i,kx)
       zmu3(jrow) = zmu3(jrow) + rix*gucos(i,3,ih)/costh
   75 continue
      do 55 i=1,ix
       cmvft (i) = dq(1) * vft(i,1) * gvbrps(i,ih)
#ifdef moist
       cmvfrd(i) = dq(1) * vfr(i,1) * gvbrps(i,ih)
#endif
   55 continue
      do 56 k=2,kx
        do 57 i=1,ix
          cmvft (i) = dq(k) * vft(i,k) * gvbrps(i,ih) + cmvft(i)
#ifdef moist
          cmvfrd(i) = dq(k) * vfr(i,k) * gvbrps(i,ih) + cmvfrd(i)
#endif
   57   continue
   56 continue
      zmvft(jrow) = 0.0
#ifdef moist
      zmvfr(jrow) = 0.0
#endif
      do 76 i=1,ix
       zmvft(jrow) = zmvft(jrow) + rix*cmvft(i)
#ifdef moist
       zmvfr(jrow) = zmvfr(jrow) + rix*cmvfrd(i)*dt
#endif
   76 continue
c
      return
      end
