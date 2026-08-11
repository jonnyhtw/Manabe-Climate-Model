      subroutine tsc
c
c  tsc is called by vertdf and performs surface heat balance
c  computations to obtain surface temperature (over land and
c  sea ice) and surface heat fluxes.
c
c  subroutine called: esndif
c
#include "parm.h"
#include "parmjx.h"
#include "parmkp.h"
c
      dimension tj(ix,kx,nhem), psj(ix,nhem)
      dimension test  (ix), ft1  (ix), ft2  (ix)
      dimension smelt (ix), rhok (ix), rsts (ix), itest(ix)
      dimension bkp   (ix), temp (ix),t9rqc9(ix), abvlw(ix)
      dimension aa    (ix), bb   (ix), cc   (ix), dd   (ix)
      dimension pst   (ix), drsts(ix), aat3 (ix), temp3(ix)
      dimension tsp   (ix), esdif(ix), htbal(ix), qf   (ix)
      dimension it1   (ix), it2  (ix), cpprh(ix), tgnd (ix)
      dimension rkx   (ix)
      dimension shf   (ix)
#ifdef wet_ep_exp
      dimension test3(ix)
      dimension aat03(ix)
      dimension f1(ix), f2(ix), df1dts(ix), df1dt0(ix), df2dts(ix)
      dimension df2dt0(ix), deti(ix)
      dimension bb1(ix), bb2(ix), ccc(ix), ddd(ix)
#endif
c
#include "albedo.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "grid.h"
c
#include "intgrl.h"
c
#include "mrstr.h"
c
#include "qa.h"
c
#include "param.h"
c
#include "sgwgt1.h"
c
#include "surfc.h"
c
#ifdef ice_dynamics
c
c use ocean temperature for heat flux on open water for ice point
c
#include "surfo.h"
c
c variables for ice model used in climate model
c
#include "varice.h"
#endif
c
#ifdef wet_ep_exp
      common / wet_ep / ts0(ix)
#endif
c
#include "zmflux.h"
c
      equivalence ( it1(    1), ft1(    1)), ( it2(  1), ft2   (  1))
      equivalence ( tj (1,1,1), cgt(1,1,1)), ( psj(1,1), gvbrps(1,1))
#ifdef ice_dynamics
c
c specific humidity for water at sea ice freezing point (cfi)
c
      data es271 / 5265.1 /
c
#endif
c
#ifdef hydrol
      data es273 / 6107.8 /
#endif
      data epspk / 1.7065 /
c  1.7065 (0.005cal/sec/cm*60sec/min)/(0.00293cal/sec/cmsq*60sec/min)
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      rcosls = 1.0 / sqrt(cosqr(ls))
      q60rg  = 60.0 * q(kx) / rgas
#ifdef ice_dynamics
c
c use normalized (by compactness) ice thickness, store old ice (cfi)
c
      do i = 1, ix
         pkiceo(i) = pkice(i)
         pkice(i) = hice(i,jrow)
      end do
c
c arlaic is the concentration of land or sea ice (cfi)
c
      do i = 1, ix
         arlaic(i) = area(i,jrow,2)
c
c set land/ice compactness to one for land points (cfi)
c
         if (zin(i) .eq. 1.) then
            arlaic(i) = 1.
         end if
      end do
#endif
c
c  set it1 = 1 where sea ice is present; in swamp models, it1=0 all pts
c
c  the following if test was changed to get-a-round a compiler error
c   on the sun.
      if (swamp) then
        do 72 i=1,ix
         it1(i) = 0
   72   continue
      else
        do 71 i=1,ix
         if (pkice(i) .gt. 0.0) then
           it1(i) = 1
         else
           it1(i) = 0
         endif
   71   continue
      endif
c
c  compute variables in surface heat
c  balance (see october 12, 1981 notes)
c
        do 1 i=1,ix
#ifdef ice_dynamics
c
c rstsow is open water mixing ratio when ice is present (cfi)
c
          rstsow(i) = (0.622*es271)/(psj(i,ih)-es271)
c
#endif
#ifdef dry
          it1  (i)      = 0
          qf   (i)      = 0.0
#endif
#ifdef moist
          bkp  (i)      = 0.3 * float(it1(i)) / (pkice(i) + epspk)
# if defined cpldleg3 || defined cpldleg1
c turn off heat conduction thru the ice if ice thicker than 2 meters
	  bkp  (i)      = min(bkp(i),2000.0)
#endif
#endif
    1   continue
c
c  ft3 = 1.0 where there is snow cover or ocean.
c
      do 2 i=1,ix
        abvlw(i) = rcosls * sqrt(gucos(i,kx,ih)**2 + gvcos(i,kx,ih)**2)
    2 continue
c
      do 3 i=1,ix
        abvlw(i) = psj(i,ih) * amax1(100.0, abvlw(i))
    3 continue
c
      do 5 i=1,ix
        rhok (i) = q60rg * cd(i,3) * abvlw(i) / (tj(i,kx,ih) + tfra)
#ifdef ice_dynamics
c
c account for open water at ice points (cfi)
c
        rhokow(i) = q60rg*cdow(i,3)*abvlw(i)/(tj(i,kx,ih)+tfra)
c
#endif
        aa   (i) = sbc
        cpprh(i) = rhok (i) * cpp
#ifdef ice_dynamics
c
c add in ice compactness dependency in upper level temperature (cfi)
c
        bb(i) = cpprh(i) * (1.0 - ekd(i)*cdfi(i,3)*arlaic(i))
#else
        bb   (i) = cpprh(i) * (1.0 - ekd(i))
#endif
#ifdef wet_ep_exp
        bb1  (i) =  cpprh(i) + bkp(i)
        bb2  (i) = -cpprh(i) * ekd(i)
#endif
#ifdef moist
        bb   (i) = bb   (i) + bkp(i)
#endif
#ifdef moist
        dd   (i) = hlt  (i) * dw (i) * rhok(i)
#endif
#ifdef dry
        dd   (i) = 0.0
#endif
#ifdef moist
#ifdef ice_dynamics
c
c different surface drag over ice and open water at ice points (cfi)
c
        cc   (i) = slwdp(i,jrow) - cpprh(i) * fkd (i)/cdio(i,3)
     *                           - dd   (i) * fkdr(i)/cdio(i,3)
c
c additional terms for upper level because of ice compactness (cfi)
c
     >       - cpprh(i)*ekd(i)*cdfo(i,3)*tfro*(1.-arlaic(i))
     >       - dd(i)*ekdr(i)*cdfo(i,4)*rstsow(i)
     >         *(1.-arlaic(i))
c
#else
        cc   (i) = slwdp(i,jrow) - cpprh(i) * fkd (i)
     *                           - dd   (i) * fkdr(i)
#endif
        cc   (i) = cc   (i) - bkp(i) * tfro
#ifdef ice_dynamics
c
c add in ice compactness dependency in upper level mixing ratio (cfi)
c
        dd   (i) = dd(i) * (1.0 - ekdr(i)*cdfi(i,4)*arlaic(i))
c
#else
        dd   (i) = dd   (i) * (1.0 - ekdr(i))
#endif
# ifdef wet_ep_exp
        ccc  (i) = slwdp(i,jrow) - cpprh(i) * fkd (i) - bkp(i) * tfro
        ddd  (i) = hlt  (i) * rhok(i)
        ccc  (i) = ccc  (i) - ddd(i) * fkdr(i)
        ddd  (i) = ddd  (i) * (1.0 - ekdr (i))
# endif
#endif
    5   continue
c
c  ft2 = 1. for land or sea ice; ft2 = 0. for open ocean; for swamp
c   models, ft2 = 1.0 everywhere. (meaning load t* from Newton iter)
c
      if (swamp) then
        do 6 i=1,ix
          ft2(i) = 1.0
    6   continue
      else
        do 7 i=1,ix
          ft2(i)  = amax0(ifix(zin(i)), it1(i))
    7   continue
      endif
#ifdef ice_dynamics
c
c let temperature over ice to be initial guess in iteration (cfi)
c
      do i = 1, ix
      if (pkice(i).gt.0. .and. zin(i) .eq. 0. .and. .not.swamp) then
         ts(i,jrow) = -tsice(i,jrow)
      end if
      end do
#endif
c
c  solve tstar eqn using Newton iteration method (do loop 10)
c    oct. 12, 1981 notes.
c
c
c
#ifdef wet_ep_exp
      do  9 i=1,ix
        ts0(i) = ts(i,jrow)
    9 continue
#endif
c
      do 10 iter=1,20
        isum = 0
c
c  get t* values for the row in degrees k and degrees c
c
        do 11 i=1,ix
          tsp (i) = abs(ts(i,jrow))
#ifdef moist
# ifdef wet_ep_exp
          ts0 (i) = abs(ts0(i))
          tgnd(i) = ts0(i) - tfra
# else
          tgnd(i) = tsp(i) - tfra
# endif
#endif
   11   continue
c
c  compute saturation vapor pressure and its derivative with
c   respect to temperature for each t* value
c
#ifdef moist
        call esndif(tgnd, ess, esdif, ix)
#endif
c
      do 12 i=1,ix
c
c  compute saturation mixing ratio
c
#ifdef moist
        pst (i) = 1.0 / (psj(i,ih) - ess(i))
        rsts(i) = 0.622 * ess(i) * pst(i)
#endif
#ifdef dry
        rsts(i) = 0.0
#endif
c
c  compute d(saturation mixing ratio)/dt
c
#ifdef moist
        drsts(i) = 0.622 * psj(i,ih) * esdif(i) * pst(i) * pst(i)
#endif
#ifdef dry
        drsts(i) = 0.0
#endif
c
c  use previous iterative value of tstar as initial guess a plug into
c   newton equation.
c
        aat3 (i) = aa (i) * tsp(i) * tsp(i) * tsp(i)
#ifdef wet_ep_exp
        aat03 (i) = aa (i) * ts0(i) * ts0(i) * ts0(i)
        f1    (i) = (aat3(i) + bb(i)) * tsp(i) + cc(i) + dd(i) * rsts(i)
        f2    (i) = (aat03(i) + bb1(i)) * ts0(i) + bb2(i) * tsp (i)
     *            + ccc(i) + ddd(i) * rsts(i)
        df1dts(i) = 4. * aat3(i) + bb(i)
        df1dt0(i) = dd(i) * drsts(i)
        df2dts(i) = bb2(i)
        df2dt0(i) = 4. * aat03(i) + bb1(i) + ddd(i) * drsts(i)
        deti  (i) = 1./(df1dts(i)*df2dt0(i) - df1dt0(i) * df2dts(i))
        temp3 (i) = tsp(i) - deti(i) * (df2dt0(i)*f1(i) -
     *              df1dt0(i)*f2(i))
        ts0   (i) = ts0(i) - deti(i) * (df1dts(i)*f2(i) -
     *              df2dts(i)*f1(i))
        test3 (i) = -deti(i) * (df1dts(i)*f2(i) - df2dts(i)*f1(i))
#else
        temp3(i) = tsp(i)
     *          - ((aat3(i) + bb(i)) * tsp(i) + cc(i) + dd(i) * rsts(i))
     *          / (4.0 * aat3(i) + bb(i) + dd(i) * drsts(i))
#endif
c
c  compute difference between previous and current values of t*
c
        test(i) = temp3(i) - tsp(i)
  12  continue
c
c  set itest = 1 if the previous and current values of t* differ
c   by more than .0008 degrees. For wet_ep models, the test
c   is also 0.0008 degrees
c
       do 74 i=1,ix
        if (abs(test(i))*ft2(i) .gt. 0.0008) go to 771
#ifdef wet_ep_exp
        if (abs(test3(i))*ft2(i) .gt. 0.0008) go to 771
#endif
   74  continue
       isum = 0
       go to 772
  771  continue
       isum = 1
  772  continue
c
c  load computed t* for land and sea ice
c
        do 75 i=1,ix
         ts(i,jrow) = ft2(i)*temp3(i) + (1.0-ft2(i))*ts(i,jrow)
#ifdef wet_ep_exp
         ts0(i)     = ft2(i)*ts0  (i) + (1.0-ft2(i))*ts(i,jrow)
#endif
   75   continue
c
c  set t* < 0 over sea ice
c
        do 76 i=1,ix
         ts(i,jrow) = (1.0-float(it1(i))*2.0)*ts(i,jrow)
#ifdef wet_ep_exp
         ts0(i)     = (1.0-float(it1(i))*2.0)*ts0(i)
#endif
   76   continue
c
c  stop iterating when all itest = 0  (ie, isum = 0)
c
        if (isum .eq. 0) go to 14
   10 continue
c
c  print diagnostic and terminate run if more than 20
c   iterations are required
c
      write(6,100) jrow,(i,test(i),ts(i,jrow),tsp(i),i=1,ix)
      stop 1234
   14 continue
c
c  make sure sea values in swamp models are set correctly
      if (swamp) then
        do 77 i=1,ix
         if (zin(i).eq.0.0 .and. ts(i,jrow).gt.0.0) then
           ts(i,jrow) = -ts(i,jrow)
         endif
   77   continue
      endif
c
#ifdef hydrol
      units = dtx2 / 9600.0
#endif
c
c  compute snow/ice melt and ocean heat flux (do loops 15, 16, 17)
c
      do 15 i=1,ix
#ifdef moist
        temp  (i)      = abs(ts(i,jrow))
#ifdef ice_dynamics
c
c calculate upper level temperature and mixing ratio (cfi)
c
        t9rqc9(i) = ekd(i)*cdfi(i,3)*temp(i) + fkd(i)/cdio(i,3)
        rkx(i) = rsts(i)*cdfi(i,4)*ekdr(i) + fkdr(i)/cdio(i,4)
c
c account for ice compactness dependency (cfi)
c
	if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
c
c teff is effective surface temperature when ice is present (cfi)
c
        teff = area(i,jrow,2)*cdfi(i,3)*temp(i)
     >         + (1.-area(i,jrow,2))*cdfo(i,3)*tfro
        t9rqc9(i) = ekd(i)*teff + fkd(i)/cdio(i,3)
c
c reff is effective surface mixing ratio when ice is present (cfi)
c
        reff = area(i,jrow,2)*cdfi(i,4)*rsts(i)
     >         + (1.-area(i,jrow,2))*cdfo(i,4)*rstsow(i)
        rkx(i) = reff*ekdr(i) + fkdr(i)/cdio(i,4)
c
	end if
c
#else
        t9rqc9(i)      = ekd  (i) * temp(i) + fkd (i)
        rkx   (i)      = rsts (i) * ekdr(i) + fkdr(i)
#endif
        dd    (i)      = hlt  (i) * dw  (i) * rhok(i)
#endif
        bb    (i)      = cpprh(i)
#ifdef moist
        bb    (i)      = bb   (i) + bkp (i)
#endif
#ifdef hydrol
        smelts(i,jrow) = 0.0
#endif
   15 continue
c
c  compute melt for snow-covered land points and glacier ice
c  note that negative albedo indicates glacier ice
c
#ifdef hydrol
      do 16 i=1,ix
        if (temp(i) .gt. tfra .and. zin(i) .eq. 1.0) then
          if (snwdpt(i,jrow).gt.0.0 .or.
     *        (snwdpt(i,jrow).le.0.0 .and. alb(i,jrow).lt.0.0)) then
            ts    (i,jrow) = tfra
#ifdef wet_ep_exp
            ts0   (i)      = tfra
#endif
#ifdef ice_dynamics
c
c correct for cd in diffusion (cfi)
c
            t9rqc9(i) = ekd(i)*cdfi(i,3)*tfra + fkd(i)/cdio(i,3)
#else
            t9rqc9(i)      = ekd(i) * tfra + fkd(i)
#endif
            rsts  (i)      = (0.622 * es273) / (psj(i,ih) - es273)
#ifdef ice_dynamics
c
c correct for cd in diffusion (cfi)
c
            rkx(i) = rsts(i)*cdfi(i,4)*ekdr(i) + fkdr(i)/cdio(i,4)
#else
            rkx   (i)      = rsts(i) * ekdr(i) + fkdr(i)
#endif
            smelts(i,jrow) = dd(i) * (rkx(i) - rsts(i))
     *                     + bb(i) * (t9rqc9(i) - tfra)
     *                     - slwdp(i,jrow) - aa(i) * frlxx4
          endif
        endif
   16 continue
#endif
c
c  compute  ocean heat flux everywhere; store only open water points
c   for swamp models; qf=0.0, because ft2=1.0 eventho zin=0.0
c
#ifdef moist
      do 17 i=1,ix
        qf(i) = dd(i) * (rkx(i) - rsts(i))
     *        + cpprh(i) * (t9rqc9(i) - temp(i))
     *        - slwdp(i,jrow) - aa(i) * (temp(i) * temp(i))**2
        qf(i) = (1.0 - ft2(i)) * qf(i)
   17 continue
#endif
c
c  compute ocean heat flux by conduction for points with sea ice
c
#ifdef hydrol
      do 18 i=1,ix
#ifdef ice_dynamics
      qficeatm(i) = 0.
#endif
      if (pkice(i) .gt. 0.0 .and. zin(i) .eq. 0.0 .and. .not.swamp) then
        if (temp(i) .le. tfra) then
          qf(i) = bkp(i) * (temp(i) - tfro)
#ifdef ice_dynamics
c
c calculate heat flux to ice for thick ice (cfi)
c
          qficeatm(i) = qf(i)
#endif
        else
c
c  compute ocean heat flux for points with melting sea ice
c
          ts(i,jrow) = -tfra
#ifdef wet_ep_exp
          ts0   (i)  = -tfra
#endif
#ifdef ice_dynamics
c
c account for ice compactness dependency for upper level variables (cfi)
c
          teff = area(i,jrow,2)*cdfi(i,3)*tfra 
     >           + (1.-area(i,jrow,2))*cdfo(i,3)*tfro
          t9rqc9(i) = teff*ekd(i) + fkd(i)/cdio(i,3)
#else
          t9rqc9(i)  =  tfra *ekd(i) + fkd(i)
#endif
          rsts  (i)  = (0.622 * es273) / (psj(i,ih) - es273)
#ifdef ice_dynamics
          reff = area(i,jrow,2)*cdfi(i,4)*rsts(i) 
     >           + (1.-area(i,jrow,2))*cdfo(i,4)*rstsow(i)
          rkx(i) = reff*ekdr(i) + fkdr(i)/cdio(i,4)
#else
          rkx   (i)  = rsts(i) * ekdr(i) + fkdr(i)
#endif
          qf    (i)  = dd  (i)* (rkx (i) - rsts(i))
     *               + cpprh(i) * (t9rqc9(i) - tfra)
     *               - slwdp(i,jrow) - aa(i) * frlxx4
#ifdef ice_dynamics
c
c calculate heat flux to ice for thick ice (cfi)
c
          qficeatm(i) = qf(i)
#endif
          endif
        endif
#ifdef ice_dynamics
c
        qficeatm(i) = qficeatm(i)*(1.-zin(i))
#endif
   18 continue
#ifdef ice_dynamics
c
c calculate open water heat flux to ocean for ice points (cfi)
c
      do i = 1, ix
         qfowatm(i) = qf(i)
         if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
         ddow = hltow(i)*dw(i)*rhokow(i)
         qfowatm(i) = ddow*(rkx(i)-rstsow(i))
     >          + rhokow(i)*cpp*(t9rqc9(i)-tfro)
     >          - slwdp(i,jrow) - aa(i)*tfro*tfro*tfro*tfro
c
         end if
         qfowatm(i) = qfowatm(i)*(1.-zin(i))
      end do
c
c average temperature weighted by heat fluxes or by compactness (cfi)
c
      do i = 1, ix
      if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
         tempi = abs(ts(i,jrow))
c
c weighted by heat fluxes (cfi)
c
         qicew = abs(qficeatm(i))*area(i,jrow,2)
         qoww = abs(qfowatm(i))*(1.-area(i,jrow,2))
         qtotal = qicew + qoww
         if (qtotal .ne. 0.) then
         ths(i,jrow) = tempi*qicew/qtotal + tfro*qoww/qtotal
         else
         ths(i,jrow) = -ts(i,jrow)
         end if
c
c set surface temp to negative for sea ice point (cfi)
c
c        ts(i,jrow) = -1.*ths(i,jrow)
c
c weighted by compactness only (cfi)
c
         tts(i,jrow)=tempi*area(i,jrow,2)+tfro*(1.-area(i,jrow,2))
c
c        ts(i,jrow) = -1.*tts(i,jrow)
c
      end if
      end do
c
c compare heat flux over each grid (cfi)
c
c     write (6, 44)
c  44 format ('checking')
c     if (iday .eq. 2 .or. ktau .eq. 1989040) then
c     do i = 1, ix
c     if (pkice(i) .gt. 0.) then
c     print *, i,jrow,qfowatm(i),qficeatm(i),tts(i,jrow),pkice(i)
c     print *, '-----------------'
c     write (6,49) i,jrow,qfowatm(i),qficeatm(i),-ts(i,jrow)
c    >            ,-tts(i,jrow),arlaic(i),area(i,jrow,2)
   49 format (2i5,6f11.5)
c     end if
c     end do
c     end if
c
c  calculate sensible and latent heat flux over open water at ice 
c  point and do a check (cfi)
c
      do i = 1, ix
      if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
c
c sensible heat flux (cfi)
c
         shfow(i,jrow) = rhokow(i)*cpp*(tfro-t9rqc9(i))
c
c latent heat flux (cfi)
c
         ddow = hltow(i)*dw(i)*rhokow(i)
         altntow(i,jrow) = ddow*(rstsow(i)-rkx(i))
c
c upward longwave radiation (cfi)
c
         uplwow = aa(i)*tfro*tfro*tfro*tfro
c
         htbalow = - slwdp(i,jrow) - uplwow - shfow(i,jrow)
     >             - altntow(i,jrow)
         htbalow = htbalow - qfowatm(i)
         if (abs(htbalow) .gt. 0.005) then
         write (6,59)
   59    format ('heat flux not balanced for open water at ice!')
         stop
         end if
      end if
      end do
c
#endif
#endif
c
c  compute and accumulate sensible and latent heat flux values over the
c   restart interval and perform a heat balance.
c
      do 78 i=1,ix
       test(i) = abs(ts(i,jrow))
   78 continue
      do 19 i=1,ix
       shf  (i)      = cpprh(i) * (test(i) - t9rqc9(i))
       shfx (i,jrow) = shfx(i,jrow) + shf(i)
#ifdef dirnal
       shfx (i,jrow) = shf(i)
#endif
#ifdef moist
       altnt(i,jrow) = dd (i) * (rsts(i) - rkx(i))
#endif
       temp (i)      = aa (i) * ts(i,jrow) * ts(i,jrow)
     *                        * ts(i,jrow) * ts(i,jrow)
       htbal(i)      = - slwdp (i,jrow) - temp(i) - shf (i)
#ifdef moist
     *                 - altnt (i,jrow)
#endif
#ifdef hydrol
     *                 - smelts(i,jrow)
#endif
       htbal(i)      = htbal(i) - qf(i)
   19 continue
      do 79 i=1,ix
       test(i) = abs(htbal(i))
   79 continue
c
c  set itest = 1 where abs(htbal) > .005
c
      isum = 0
       do 80 i=1,ix
        if (test(i) .gt. 0.005) go to 781
   80  continue
       isum = 0
       go to 782
  781  continue
       isum = 1
  782  continue
c
#ifdef ice_dynamics
c
c store heat flux summed over the atmospheric loop (cfi)
c
      do i = 1, ix
      if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
         qfice(i,jrow) = qfice(i,jrow) + qficeatm(i)
      end if
      qfow(i,jrow) = qfow(i,jrow) + qfowatm(i)*(1.0-zin(i))
      end do
#endif
c
c  store ocn heat flux for use by ocn model
c
      do 400 i=1,ix
#ifdef ice_dynamics
c
c include ice and open water heat fluxes in ocean heat flux (cfi)
c
        qft = qf(i)
        if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
        qft=qficeatm(i)*area(i,jrow,2)+qfowatm(i)*(1.0-area(i,jrow,2))
        end if
        qf(i) = qft
        qao(i,jrow) = qao(i,jrow) + qft*(1.0-zin(i))
c
#else
        qao(i,jrow) = qao(i,jrow)+qf(i)*(1.0-zin(i))
#endif
  400 continue
#ifdef ice_dynamics
c
c store temperature with drag used in vertical diffusion (cfi)
c
      do i = 1, ix
      tempi = abs(ts(i,jrow))
      tscd(i,jrow) = tempi*cd(i,3)
      if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
	 tscd(i,jrow) = area(i,jrow,2)*cd(i,3)*tempi
     >                  + (1.-area(i,jrow,2))*cdow(i,3)*tfro
      end if
      end do
c
c restore ice thickness (cfi)
c
      do i = 1, ix
         pkice(i) = pkiceo(i)
      end do
c
c change surface temperature by weighting open water and sea ice (cfi)
c
      do i = 1, ix
c
c initialize ice temperature to melting point of fresh water (cfi)
c
      tsice(i,jrow) = tfra
c
      if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
c
c store ice temperature (cfi)
c
         tsice(i,jrow) = -ts(i,jrow)
c
c weighting by heat fluxes (cfi)
c
c        ts(i,jrow) = -ths(i,jrow)
c
c weighting by temperature (cfi)
c
         ts(i,jrow) = -tts(i,jrow)
      end if
      end do
#endif
c
c   here insert print code
      if (mod(ktau+1,ioprt) .eq. 0) then
        if (jrow .eq. 1) write(6,315)
        iprt = 48
        write(6,320) qao(iprt,jrow), qf(iprt), slwdp(iprt,jrow),
     *    shfx(iprt,jrow) ,ts(iprt,jrow), zin(iprt), altnt(iprt,jrow),
     *    pkice(iprt) ,smelts(iprt,jrow), cd(iprt,3), abvlw(iprt),
     *    tj(iprt,kx,ih), psj(iprt,ih), rsts(iprt), rkx(iprt),
     *    dd(iprt), ekdr(iprt), fkdr(iprt), cpprh(iprt), t9rqc9(iprt),
     *    ekd(iprt), fkd(iprt)
  315   format(/1x,'      qao and components from tsc'//
     * '    qao         qf       slwdp       shfx        ts        zin',
     * '       altnt      pkice      smelts       cd       abvlw   '/
     * '     tj        psj        rsts       rkx         dd    ',
     * '    ekdr       fkdr      cpprh      t9rqc9      ekd    ',
     * '    fkd    '/)
  320   format (2(1x,11e11.5/))
      endif
c
c
c  convert melt from ly/min to cm/timestep
c
#ifdef hydrol
      do 20 i=1,ix
        smelts(i,jrow) = units * smelts(i,jrow)
   20 continue
#endif
c
c  compute hemispheric integrals of sensible and latent heat
c   fluxes on radiation time steps for print out
c
      if (lrad) then
        do 81 i=1,ix
         zmsnht(jrow) = zmsnht(jrow) + rix *   shf(i)
         zmltht(jrow) = zmltht(jrow) + rix * altnt(i,jrow)
         zmqao (jrow) = zmqao (jrow) + rix *    qf(i)
   81   continue
      endif
c
      if (isum .eq. 0) then
c
        return
c       ======
c
      else
c
c  print diagnostic and terminate run
c    if heat balance error exists
c
      do 21 i=1,ix
        write(6,101) i,jrow,htbal(i),slwdp(i,jrow),
     *                                    shfx(i,jrow),ts(i,jrow),zin(i)
#ifdef moist
     *                                   ,altnt(i,jrow),qf(i),pkice(i)
#endif
#ifdef hydrol
     *                                   ,smelts(i,jrow)
#endif
     *              ,cd(i,3),abvlw(i),tj(i,kx,ih),psj(i,ih)
#ifdef moist
     *              ,rsts(i),rkx(i),dd(i),ekdr(i),fkdr(i)
#endif
     *              ,cpprh(i),t9rqc9(i),ekd(i),fkd(i)
   21 continue
      call aquit (101)
c     ==== ====
c
      stop 101
      endif
c
  100 format(1x,'*** maximum iterations in tsc ***  row=',i3/
     *       (1x,i5,3e16.7))
  101 format('0heat balance error in tsc'//1x,2i5,9e13.6/11x,9e13.6/11x,
     *       4e13.6)
c
      end
