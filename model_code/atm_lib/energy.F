      subroutine energy
c
c  this subroutine is called by atmos.  it calls subroutine chksum.
c    energy computes in the spectral domain global integrals of model
c    variables and diagnostic quantities.  (see sept. 12, 1978 notes)
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
c
      complex vor  , div  ,   t  , psln
      complex vorm , divm ,   tm , psm
      complex vordt, divdt,   tdt, psdt
      complex ucos , vcos , psdx , psdy , zs , ek  , ps
      complex qmix , qmixm, qmxdt
c
#include "ctlblk.h"
c
#include "dispat.h"
c
#include "gauss.h"
c
#include "intgrl.h"
c
#include "nxstps.h"
c
#include "onlint.h"
c
#include "param.h"
c
#include "pstar.h"
c
#include "sfharm.h"
c
#include "sgwgt1.h"
c
#include "water.h"
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c  compute mean global integrals of:
c     zeke = variance of wind based on divergence and vorticity
c             (check on integral eke)
c     znke = zonal kinetic energy
c     entr = variance of vorticity
c     divv = variance of divergence
c     diss = del fourth linear horizontal diffusion
c
c  zero integrals initially for accumulation of contributions to them
c
      eke  = 0.0
      zeke = 0.0
      znke = 0.0
      entr = 0.0
      divv = 0.0
      diss = 0.0
c
c  derive global mean surface pressure from spherical harmonic p*(0,0)
c
      psbare = real(ps(1,1,1)) / sqrt(2.0)
c
c  vertically integrate global mean variance of wind derived from
c    spherical harmonic ek(0,0) at all kx levels.
c
      do 23 k=1,kx
        eke = eke + real(ek(1,1,k)) * dq(k)
   23 continue
      eke = 0.5 * eke / sqrt(2.0)
c
c  compute zonal wave number zero odd function contributions
c
      do 5 k=1,kx
        do 5 n=2,jtrun1,jin
c
c  variable coef (1/(l(l+1))**2) in statements below converts squares of
c    vorticity and divergence to squares of velocity potential (chi**2)
c    and stream function (psi**2) respectively (see nov. 4, 1981 notes).
c
c  complex conjugate multiplications below compute variances of chi and
c    psi which when added together yield variance of wind.
c    see subroutine spint comments and notes of nov. 24, 1981 and sept.
c    13, 1983 for information on global integration in spectral space.
c
c  arcane function of coef below results from del fourth operator
c    in equation for horizontal diffusion.
c
          coef = 1.0 / float(n * (n-1))
          rlv  = real(vor(1,n,k) * conjg(vor(1,n,k))) * dq(k)
          entr = entr + rlv
          diss = diss + coef * rlv * ((1.0 / coef) - 2.0) ** 4
          znke = znke + coef * rlv
          zeke = zeke + coef * rlv
    5 continue
c
c  compute zonal wave number zero even function contributions
c
      do 2 k=1,kx
        do 2 n=ji,jtrun1,jin
          coef = 1.0 / float(n * (n-1))
          rld  = real(div(1,n,k) * conjg(div(1,n,k))) * dq(k)
          divv = divv + rld
          diss = diss + coef * rld * ((1.0 / coef) - 2.0) ** 4
          znke = znke + coef * rld
          zeke = zeke + coef * rld
    2 continue
c
c  compute non zero zonal wave number odd function contributions
c
      do 6 k=1,kx
        do 6 n=jin,jtrun1,jin
          do 6 m=2,mx
            coef = 1.0 / float((n + isc*(m-1) - 1) * (n + isc*(m-1)))
            rlv  = 2.0 * real(vor(m,n,k) * conjg(vor(m,n,k))) * dq(k)
            entr = entr + rlv
            diss = diss + coef * rlv * ((1.0 / coef) - 2.0) ** 4
            zeke = zeke + coef * rlv
  6   continue
c
c  compute non zero zonal wave number even function contributions
c
      do 3 k=1,kx
        do 3 n=1,jtrun1,jin
          do 3 m=2,mx
            coef = 1.0 / float((n + isc*(m-1) - 1) * (n + isc*(m-1)))
            rld  = 2.0 * real(div(m,n,k) * conjg(div(m,n,k))) * dq(k)
            divv = divv + rld
            diss = diss + coef * rld * ((1.0 / coef) - 2.0) ** 4
            zeke = zeke + coef * rld
    3 continue
c
c  multiply by functions of physical constants such as earth radius
c    and g, spectral normalization factors and by mean global surface
c    pressure to convert integrals to cgs units.
c
      zeke  = 0.25 * eradsq * zeke
      znke  = 0.25 * eradsq * znke * psbare * og
      eddy  = psbare * og * zeke - znke
      entr  = sqrt(0.5 * entr)
      divv  = sqrt(0.5 * divv)
      diss  = 0.25 * (viscxy * radin * rasq) ** 2 * diss
      diss  = 2.0 * sqrt(zeke * diss) * psbare / grav
      tbare = 0.0
      am    = 0.0
      tke   = 0.0
      qbare = 0.0
c
c  compute pressure-weighted global integrals in spectral space for
c    temperature, kinetic energy, zonal wind and mixing ratio.
c    computation is similar to that used in subroutine spint.
c  for m = 0 waves
c
      do 7 k=1,kx
        do 7 n=1,jtrun1,jin
          rtb   = real(ps (1,n,1) * conjg( t  (1,n,k))) * dq(k)
          rle   = real(ps (1,n,1) * conjg( ek (1,n,k))) * dq(k)
          rtm   = real(ps (1,n,1) * conjg(ucos(1,n,k))) * dq(k)
          tbare = tbare + rtb
          tke   = tke   + rle
          am    = am    + rtm
          rqb   = real(ps (1,n,1) * conjg(qmix(1,n,k))) * dq(k)
          qbare = qbare + rqb
    7 continue
c
c  for m.gt.0 waves
c
      do 8 k=1,kx
        do 8 n=1,jtrun1,jin
          do 8 m=2,mx
            rtb   = 2.0 * real(ps (m,n,1) * conjg( t  (m,n,k))) * dq(k)
            rle   = 2.0 * real(ps (m,n,1) * conjg( ek (m,n,k))) * dq(k)
            rtm   = 2.0 * real(ps (m,n,1) * conjg(ucos(m,n,k))) * dq(k)
            tbare = tbare + rtb
            tke   = tke   + rle
            am    = am    +rtm
            rqb  = 2.0 * real(ps(m,n,1) * conjg(qmix(m,n,k))) * dq(k)
            qbare = qbare + rqb
  8   continue
c
c  process integrals to convert to proper units.
c
      tbare =  0.5 * tbare / psbare + tfra
      qbare =  0.5 * qbare / grav
      tke   =  0.5 * o2g * tke
      hmf   =  float(jin)
      cnvr  = -hmf * o2g * cp * cnvr
      diss  =  diss - hmf * o2g * tvdis
      am    =  erad * o2g * am
c
c  compute angular momentum of resting atmosphere from p*(0,0) and
c    p*(0,2) spherical harmonics (see notes of december 14, 1981).
c
      amrat = (eradsq * ww * sqrt(2.)   / (3.0 * grav))
     *      * (ps(1,1,1) - ps(1,3,1) / sqrt(5.0))
      drag  = 0.5 * hmf * erad * drvor
      trq2  = 0.5 * hmf * erad * zmts
      trqt  = drag + trq2
      enthl = cp * psbare * tbare / grav
      ekpot = tke + enthl
      angm  = am + amrat
      dissp = diss
      cnver = cnvr
      write(6,111) psbare, tbare,  tke,  am  ,  angm,
     *              enthl, ekpot,  eke,  zeke,  divv,
     *               entr, qbare, drag, cnver, dissp,
     *               trq2,  trqt, znke,  eddy
c
      call chksum(psbare,19)
c
      return
c
 111  format(/' psbar =',e14.7,3x,'tbar  =',e14.7,3x,'kener =',e14.7,3x,
     *       'rangm =',e14.7,3x,'aangm =',e14.7/' enthlp=',e14.7,3x,
     *       'ekpot =',e14.7,3x,'ekvar =',e14.7,3x,'divort=',e14.7,3x,
     *       'divar =',e14.7/' vrtvar=',e14.7,3x,'qbar  =',e14.7,3x,
     *       'sfdrag=',e14.7,3x,'cnver =',e14.7,3x,'dissp =',e14.7/,
     *       ' mntrq =',e14.7,3x,'torque=',e14.7,3x,'zoneke=',e14.7,3x,
     *       'eddyke=',e14.7/)
c
      end
