      subroutine lwcool(tj, rj, psj, icf)
c
c  compute net longwave cooling rates
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  routine :   lwcool
c  called by : impt
c  calls :     lwtran
c  purpose :   compute net longwave cooling rates
c
c  inputs:
c     **passed** tj             celsius temperature at each grid point
c     **passed** rj             mixing ratio of h2o at each grid point
c     **passed** psj            surface pressure at each ground point
c     / grdfld / ts             ground temperature at each ground point
c     / trans  / transm         longwave transmission coefficient
c     / trans  / transz         correction term for top of atmosphere
c     / trans  / absrbx         quarter-layer absorption coefficient
c     / cld    / cldlw          fractional cloud cover
c                                  (fixed cloud model)
c     / cldpm  / ccover         fractional cloud cover
c                                  (cloud prediction model)
c
c  outputs:
c     / auxatm / ztufxt         upward flux at top of atmosphere
c     / auxatm / zufxnt         net flux at ground
c     / grdfld / slwdp          downward flux at ground
c     / grid   / grad           cooling rate at each grid point
c
c
c  procedure:
c     the temperature at each grid point is converted from celsius to
c     absolute (kelvin); from the temperatures planck blackbody
c     radiances are computed for each grid point.  the longwave
c     transmission coefficients (recomputed every ldisk timesteps) are
c     then combined with the blackbody radiances to determine clear
c     fluxes (i.e., the flux which would occur in the absence of
c     clouds) between all pairs of levels.  dqx(i,k,l), for example, is
c     the clear-path downward flux from vertical level l to level k.
c     then, from the known distribution of clouds, the actual upward and
c     downward fluxes at each level are computed.  these are used to
c     compute cooling rates as well as some diagnostic quantities.
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
#include "parm.h"
#include "parmkp.h"
#include "parmrad.h"
c
      parameter ( sigma = 5.673e-5, c1 = 3.740e-5, c2 = 1.4385 )
      parameter ( cp1 = 0.24, cp2 = 4.1867e7 * cp1 )
c
#include "auxatm.h"
c
#include "cld.h"
c
#include "ctlblk.h"
c
#include "grid.h"
c
#include "grdfld.h"
c
#include "param.h"
c
#include "planck.h"
c
#include "sgwgt1.h"
c
#include "trans.h"
c
      dimension dqx (ix,0:kp,0:kp), uqx (ix,0:kp,0:kp)
      dimension dfxc(ix,0:kp),      ufxc(ix,0:kp)
      dimension psigl(ix,2:kx)
c
      dimension tj(ix,kx,nhem), rj(ix,kx,nhem), psj(ix,nhem)
c
#ifdef prdcld
      integer kloud (ix,0:kp)
      integer kldtop(ix,2:kx), kldmid(ix,2:kx), kldbot(ix,2:kx)
c
      dimension numtop(2:kx), nummid(2:kx), numbot (2:kx)
      dimension ufxkb (ix),   uqxkt (ix),   ufxb(ix,2:kx)
      dimension dfxkt (ix),   dqxkb (ix),   dqxb(ix,2:kx)
      dimension ppfkb (ix),   ppfkt (ix),   ppfb(ix,2:kx)
#endif
      dimension cdxa  (kx),   cuxa(0:kx)
      dimension dfx(ix,0:kp), ufx(ix,0:kp)
c
      equivalence (dfx, dqx(1,0,0)), (ufx, uqx(1,0,kp))
c
      data cdxa / km * 0.5, 1.0 /, cuxa / -1.0, kx * -0.5 /
c
c ----------------------------------------------------------------------
c  compute clear (cloud-free) fluxes between all pairs of levels
c ----------------------------------------------------------------------
c
c ****** compute pbb (planck blackbody radiance), and d(pbb)/dk
c
      do 5 k=1,kx
        do 5 i=1,ix
          tdegk(i,k) = tj(i,k,ih) + tfra
    5 continue
c
      do 6 i=1,ix
        tdegk(i, 0) = tdegk(i,1)
        tdegk(i,kp) = abs(ts(i,jrow))
    6 continue
c
      do 7 k=0,kp
        do 7 i=1,ix
          tcube(i,k) = tdegk(i,k) * tdegk(i,k) ** 2
          pbb  (i,k) = sigma * tcube(i,k) * tdegk(i,k)
    7 continue
c
      do 8 k=0,kx
        do 8 i=1,ix
          dpbb(i,k) = pbb(i,k+1) - pbb(i,k)
    8 continue
c
c ****** call lwtran to recompute trans matrix on full lw timestep
c
      if (.not. ldisk) go to 10
c
      call lwtran (tj, rj, psj)
c
c ****** compute upward and downward fluxes at "flux" levels
c
   10 continue
      do 11 i=1,ix
        dqx(i,0,0) = 0.0
   11 continue
c
      do 12 k=1,kx
        do 12 i=1,ix
          dqx(i,k,k) = cdxa(k) * dpbb(i,k) * absrbx(i,k,2) + pbb(i,k)
   12 continue
c
      do 13 l=kx-1,0,-1
        do 13 k=l+1,kx
          do 13 i=1,ix
            dqx(i,k,l) = dqx(i,k,l+1) - transm(i,k,l) * dpbb(i,l)
   13 continue
c
      do 14 k=1,kx
        do 14 i=1,ix
          dqx(i,k,0) = dqx(i,k,0) - transz(i,k) * pbb(i,0)
   14 continue
c
      do 15 k=0,kx
        do 15 i=1,ix
          uqx(i,k,k+1) = cuxa(k) * dpbb(i,k) * absrbx(i,k,1)
     $                   + pbb(i,k+1)
   15 continue
c
      do 16 l=2,kp
        do 16 k=0,l-2
          do 16 i=1,ix
            uqx(i,k,l) = uqx(i,k,l-1) + transm(i,k,l-1) * dpbb(i,l-1)
   16 continue
c
c ----------------------------------------------------------------------
c  compute actual fluxes from cloud distribution
c ----------------------------------------------------------------------
c
#ifdef fixcld
      do 20 m=2,ncc
        kb = kbh(m)
        do 20 k=kb,kx
          do 20 i=1,ix
            dfx(i,k) = dfx(i,k)
     $                 * (1.0 - cldlw(m)) + dqx(i,k,kb) * cldlw(m)
   20 continue
c
      do 30 m=ncc,2,-1
        kt = kth(m)
        do 30 k=0,kt-1
          do 30 i=1,ix
            ufx(i,k) = ufx(i,k)
     $                 * (1.0 - cldlw(m)) + uqx(i,k,kt) * cldlw(m)
   30 continue
c
      if (nc .lt. 3) go to 80
c
      do 50 m=2,ncc
        kt = kth(m)
        kb = kbh(m)
        if (kb .eq. kt) go to 50
        do 40 k=kt,kb-1
          pfactr = (qmh(k) - qmh(kt-1)) / (qmh(kb) - qmh(kt-1))
          do 40 i=1,ix
            dfxc(i,k) = dfx(i,kt-1)
     $                  + pfactr * (dqx(i,kb,kb) - dfx(i,kt-1))
            ufxc(i,k) = uqx(i,kt-1,kt)
     $                  + pfactr * (ufx(i,kb) - uqx(i,kt-1,kt))
   40   continue
        do 45 k=kt,kb-1
          do 45 i=1,ix
            dfx(i,k) = cldlw(m) * dfxc(i,k)
     $                 + (1.0 - cldlw(m)) * dfx(i,k)
            ufx(i,k) = cldlw(m) * ufxc(i,k)
     $                 + (1.0 - cldlw(m)) * ufx(i,k)
   45   continue
   50 continue
#endif
c
#ifdef prdcld
c ****** cloud prediction version of computation
c
      do 105 k=0,kp
       do 105 i=1,ix
        if (ccover(i,k) .ne. 0) then
          kloud(i,k) = 1
        else
          kloud(i,k) = 0
        endif
  105 continue
      do 110 k=2,kx
       do 110 i=1,ix
        kldtop(i,k) = kloud(i,k) * (1 - kloud(i,k-1))
        kldmid(i,k) = kloud(i,k) *      kloud(i,k+1)
        kldbot(i,k) = kloud(i,k) * (1 - kloud(i,k+1))
  110 continue
c
      do 20 kb=2,kx
        do 20 k=kb,kx
         do 20 i=1,ix
          if (kldbot(i,kb) .eq. 1)
     *      dfx(i,k) = dfx(i,k) * (1.0-ccover(i,kb))
     *               + dqx(i,k,kb) * ccover(i,kb)
   20 continue
c
      do 30 kt=kx,2,-1
        do 30 k=0,kt-1
         do 30 i=1,ix
          if (kldtop(i,kt) .eq. 1)
     *      ufx(i,k) = ufx(i,k) * (1.0-ccover(i,kt))
     *               + uqx(i,k,kt) * ccover(i,kt)
   30 continue
c
      do 35 k=2,kx
        numtop(k) = 0
        numbot(k) = 0
        nummid(k) = 0
   35 continue
      do 36 k=2,kx
       do 36 i=1,ix
        numtop(k) = numtop(k) + kldtop(i,k)
        numbot(k) = numbot(k) + kldbot(i,k)
        nummid(k) = nummid(k) + kldmid(i,k)
   36 continue
c
      do 40 k=kx,2,-1
       do 115 i=1,ix
        if (kldbot(i,k) .eq. 1) then
          ufxkb(i) = ufx(i,k)
          dqxkb(i) = dqx(i,k,k)
          ppfkb(i) = qmh(k)
        endif
  115  continue
        if (nummid(k) .eq. 0) goto 40
       do 120 i=1,ix
        if (kldmid(i,k) .eq. 1) then
          ufxb(i,k) = ufxkb(i)
          dqxb(i,k) = dqxkb(i)
          ppfb(i,k) = ppfkb(i)
        endif
  120  continue
  40  continue
c
      do 50 k=2,kx
       if (nummid(k) .eq. 0) goto 50
         do 125 i=1,ix
          if (kldtop(i,k) .eq. 1) then
            dfxkt(i) = dfx(i,k-1)
            uqxkt(i) = uqx(i,k-1,k)
            ppfkt(i) = qmh(k-1)
          endif
  125    continue
         do 130 i=1,ix
          if (kldmid(i,k) .eq. 1) then
            psigl(i,k) = (qmh(k) - ppfkt(i)) / (ppfb(i,k) - ppfkt(i))
            dfxc (i,k) = dfxkt(i) + psigl(i,k)*(dqxb(i,k) - dfxkt(i))
            ufxc (i,k) = uqxkt(i) + psigl(i,k)*(ufxb(i,k) - uqxkt(i))
          endif
  130    continue
   50 continue
c
      do 135 k=2,kx
       do 135 i=1,ix
        if (kldmid(i,k) .eq. 1) then
          dfx(i,k) = dfx(i,k)*(1.0-ccover(i,k)) + dfxc(i,k)*ccover(i,k)
          ufx(i,k) = ufx(i,k)*(1.0-ccover(i,k)) + ufxc(i,k)*ccover(i,k)
        endif
  135 continue
#endif
c
c ----------------------------------------------------------------------
c  compute net cooling rates from fluxes
c ----------------------------------------------------------------------
c
c 1.43306e-6 converts ergs / (cm**2*sec) to ly/min
c
   80 continue
c
c  on the first pass (icf=1), save everything, because clouds
c  are the real clouds. On the second pass (icf=2), save just
c  selected variables, because the clouds=0.
      if (icf .eq. 1) then
        do 85 i=1,ix
         ztufxt(i)     = ufx(i,0) * 1.43306e-06
         zufxnt(i)     = (ufx(i,kx) - dfx(i,kx)) * 1.43306e-06
         slwdp(i,jrow) = -dfx(i,kx) * 1.43306e-6
   85   continue
c
c 1.1574e-5 = 1.0 / 86400.0
c
        cof = grav * 86400.0 / cp2
        cof = cof * 1.1574e-5
c
        do 90 k=1,kx
         do 90 i=1,ix
          grad(i,k,ih) = cof / (dq(k) * psj(i,ih))
   90   continue
c
        do 95 k=1,kx
         do 95 i=1,ix
          grad(i,k,ih) = grad(i,k,ih)
     $       * (ufx(i,k) - ufx(i,k-1) + dfx(i,k-1) - dfx(i,k))
   95   continue
#ifdef cloud_forcing
      else
        do 850 i=1,ix
         ztufxtc(i)     = ufx(i,0) * 1.43306e-06
         zufxntc(i)     = (ufx(i,kx) - dfx(i,kx)) * 1.43306e-06
  850   continue
#endif
      endif
c
      return
      end
