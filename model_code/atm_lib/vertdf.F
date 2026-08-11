      subroutine vertdf
c
c  this subroutine is called by grmult.  vertdf computes tendencies of
c    u, v, t (and r if a moist model) due to vertical diffusion alone
c    and stores them in the first 3 (4) fields of common block /auxatm/.
c    it calls subroutine tsc (tscs in annual mean model) to make the t*
c    computation consistent with the vertical heat flux computed here.
c     for boundary layer models, vertical diffusion is computed by a
c     static stability-dependent mixing length formula based on the
c     mellor-yamada theory.
c
#include "parm.h"
#include "parms.h"
#include "parmkp.h"
      parameter ( iter = 6 )
      parameter ( ndtjv = ndt * jv, ndtjvl = ndtjv * lnd )
c
c
#ifdef moist
      dimension vf(ix,kx,4)
#endif
#ifdef dry
      dimension vf(ix,kx,3)
#endif
c
#include "albedo.h"
c
#include "auxatm.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "grid.h"
c
#include "guvm.h"
c
#include "intgrl.h"
c
#include "mixlen.h"
c
#include "mrstr.h"
c
#include "param.h"
c
#if defined cpldleg3 || defined cpldleg1
#include "qa.h"
#else
c
#ifdef ice_dynamics
c
c requires sno variable in ice dynamics calculation (cfi)
c
#include "qa.h"
#endif
c
#endif
#include "richrd.h"
c
#include "sgwgt1.h"
c
#include "surfc.h"
#ifdef ice_dynamics
c
c collect wind stress for ice dynamics calculatin (cfi)
c
#include "varice.h"
#endif
c
#include "vrtcom.h"
c
#ifdef wet_ep_exp
      common/ wet_ep / ts0(ix)
#endif
c
      dimension  umj(ix,kx), vmj(ix,kx),  tmj(ix,kx)
      dimension  zph(ix,kx),  zh(ix,kx)
#ifdef garratt
      dimension  dcdzs(ix,km), dtdz(ix,km)
      dimension  dzh(ix,km),   ac(ix,km)
#endif
      dimension sd(ix,kx), se(ix,kx)
      dimension s1(ix), s3(ix), s5(ix), s6(ix), s7(ix)
      dimension f(ix,kx,lmx), prox(ix,kx,lmx), e(ix,kx,lmx)
      dimension w(ix,kx),   x (ix,kx),  y(ix,kx)
      dimension z(ix,kx), flux(ix,lmx)
      dimension vrt5(ix), vrt6(ix), hm(ix,kx,lmx), vrt7(ix)
      dimension dkmh(ix,kp,lmx), s4(kx), sc(ix,kx), s2(ix)
      dimension rmj(ix,kx), dkmin(ix)
      dimension ruvmag(ix), theta(ix,km), test(ix,kx)
      dimension ulow(ix), vlow(ix), tlow(ix), thlow(ix)
      dimension bvfreq(ix), kcrit(ix)
      dimension tdegc(ix), z0(ix)
#if defined cpldleg3 || defined cpldleg1
# ifdef windbc
      dimension swind(ix)
# endif
#endif
c
#include "zsrvar.h"
c
c   rk0 is the characteristic scale of the roughness.
c
      equivalence (  vfu(1,1),   vf(1,1,1))
c
#if defined cpldleg3 || defined cpldleg1
#ifdef windbc
c   eps is used to avoid taking the sqrt of zero
      data eps / 1.0e-7 /
c
# endif
#endif
      data spmin / 100.0 /
c
c     rwk=1./(.75*15.) = 1./11.25 = .088888889
c
      data rwk / 0.088888889 /
#ifdef ice_dynamics
c
c specific humidity for water at sea ice freezing point (cfi)
c
      data es271 / 5265.1 /
c
#endif
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      secsqr = 1.0 / cosqr(ls)
      sekan  = sqrt(secsqr)
      rsekan = 1.0 / sekan
      rdtx2  = 1.0 / dtx2
c
c  zin=1 for land and zero for sea
c
c   compute drag coefficient cd
c   estimate soil moisture factor in moist models without hydrology
c   note:  in boundary layer models rms and cd are recomputed later
c
      do 2 i=1,ix
        dw(i)  = 1.0 - 0.5 * zin(i)
        s1(i)  = 0.0
    2 continue
c
c  set flags s1, s2, s5 and s6 according to the following definitions -
c     s1  1.0=snow                           0.0=no snow
c     s2  1.0=snow or sea ice or glacier     0.0=rest of points
c     s5  1.0=snow or sea point              0.0=no snow, land points
c     s6  1.0=glacier                        0.0=no glacier
c
      do 3 i=1,ix
       if (snwdpt(i,jrow) .gt. 0.0) then
         s1(i) = 1.0
       else
         s1(i) = 0.0
       endif
    3 continue
      do 4 i=1,ix
       if (pkice(i) .gt. 0.0) then
         s2(i) = 1.0
       else
         s2(i) = 0.0
       endif
    4 continue
      do 5 i=1,ix
        s2(i) = (1.0 - zin(i)) * s2(i)
        s2(i) = max(s1(i),s2(i))
    5 continue
      do 6 i=1,ix
       if (alb(i,jrow) .lt. 0.0) then
         s6(i) = 1.0
       else
         s6(i) = 0.0
       endif
    6 continue
      do 7 i=1,ix
        s2 (i) = max(s6(i),s2(i))
        hlt(i) = 600.0 + 80.0 * s2(i)
#ifdef ice_dynamics
c
c get latent heat coefficient for open water in ice points (cfi)
c
	if (pkice(i) .gt. 0.0) then
           hltow(i) = 600.
        else
           hltow(i) = 0.
        end if
c
#endif
        s5 (i) = 1.0 - zin(i) * (1.0 - s1(i))
    7 continue
c
c   recompute soil moisture factor in models with hydrology
c
c
      do 8 i=1,ix
        dw(i) = rwk * soilm(i,jrow)
    8 continue
      do 9 i=1,ix
        dw(i) = min(dw(i), 1.0)
        dw(i) = max(dw(i), 0.0)
    9 continue
      do 10 i=1,ix
        dw(i) = s5(i) + (1.0 - s5(i)) * dw(i)
   10 continue
c
c   compute height above ground of q levels in boundary layer
c
      do 14 k=1,kx
        do 14 i=1,ix
          tmj(i,k) = gtm(i,k,ih) + tfra
   14 continue
      do 13 i=1,ix
        zph(i,kx) = 0.0
   13 continue
c
c zh is temporarily the layer thickness from k+1/2 to k-1/2
      do 15 k=2,kx
        do 15 i=1,ix
          zh(i,k)  = qlgc3(k) * tmj(i,k) * og
   15 continue
c
c zph is the height of level k+1/2 (grid box bottom)
      do 25 k=1,km
        do 25 i=1,ix
          zph(i,kx-k) = zph(i,kp-k) + zh(i,kp-k)
   25 continue
c
c now zh is height of level k (grid box center)
      do 30 k=1,kx
        do 30 i=1,ix
          zh(i,k) = zph(i,k) + factr(k) * tmj(i,k) * og
   30 continue
c
c compute neutral drag coefficient at each point
      do 11 i=1,ix
        z0(i) = z0land * zin(i) + z0sea * (1.0 - zin(i)) 
   11 continue
c
c same neutral drag coefficient is used for momentum, heat,
c and water vapor
      do 122 l=1,lmx
        do 122 i=1,ix
          cd(i,l) = ( vk0 / alog(zh(i,kx)/z0(i)) ) ** 2
  122 continue
c
c  remove cosine from wind components before using
c
      do 20 k=1,kx
        do 20 i=1,ix
          umj(i,k) = sekan * gucosm(i,k,ih)
          vmj(i,k) = sekan * gvcosm(i,k,ih)
          rmj(i,k) =         gqmixm(i,k,ih)
   20 continue
c
#ifdef no_neg_in_vertdf
c  for diffusion calculations, do not allow mixing ratios
c  to be negative; set them to zero instead
      do 4321 k=1,kx
        do 4321 i=1,ix
          if (rmj(i,k).lt.0.0) then
            rmj(i,k) = 0.0
          endif
 4321 continue
#endif
c
#ifdef garratt
c
c  compute richardson numbers (ri and rf) and stability functions
c   (sm and sh).  see notes of june 1976.
c
      do 35 k=kvp,kx
        do 31 i=1,ix
          ac(i,k-1)    = vrt4(k-1) / (tmj(i,k) + tmj(i,k-1))
          dcdzs(i,k-1) = ((umj(i,k-1) - umj(i,k)) ** 2 +
     *                   (vmj(i,k-1) - vmj(i,k)) ** 2) * ac(i,k-1) ** 2
   31   continue
        do 578 i=1,ix
         dcdzs(i,k-1) = max(dcdzs(i,k-1), 1.0e-10)
  578   continue
        do 36 i=1,ix
          dzh(i,k-1)  = zph(i,k-1) - zph(i,k)
          dtdz(i,k-1) = (tmj(i,k-1) - tmj(i,k)) * ac(i,k-1)
          ri(i,k-1)   = (grav * (dtdz(i,k-1) + gm))
     *                / (0.5 * (tmj(i,k) + tmj(i,k-1)) * dcdzs(i,k-1))
   36   continue
   35 continue
#endif
c
c  initialize d (at k-1/2) and e and f (at k=kv-1) to zero
c
      do 100 l=1,lmx
        do 100 k=1,kp
          do 100 i=1,ix
            dkmh (i,k,l) = 0.0
  100 continue
c
      do 160 i=1,ix
       gmzkx(i) = gm * zh(i,kx)
  160 continue
c
      do 161 l=1,lmx
       do 161 i=1,ix
        f(i,kvm,l) = 0.0
        e(i,kvm,l) = 0.0
  161 continue
c
c   store variables in hm fields and convert temperature to potential
c
      do 200 k=1,kx
        do 200 i=1,ix
          hm(i,k,1) = umj(i,k)
          hm(i,k,2) = vmj(i,k)
          hm(i,k,3) = gtm(i,k,ih) + gm * zh(i,k)
          hm(i,k,4) = rmj(i,k)
  200 continue
c
#ifdef ice_dynamics
c
c cdow is drag coefficient for open water in ice points (cfi)
      do l = 1, lmx
        do i = 1, ix
         cdow(i,l) = cd(i,l)
        end do
      end do
c
c cdio is multiplier for diffusion quantities (cfi)
c
      do l = 1, lmx
        do i = 1, ix
c
c set cdio to original drag value for non-ice points (cfi)
c
           cdio(i,l) = cd(i,l)
           if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
             cdio(i,l) = cd(i,l)*area(i,jrow,2) +
     $                   cdow(i,l)*(1.-area(i,jrow,2))
           end if
        end do
      end do
c
c calculate the fractions of drag coefficients for ice points (cfi)
c
      do l = 1, lmx
        do i = 1, ix
c
c set the ratios to one for non-ice points (cfi)
c
          cdfi(i,l) = 1.
          cdfo(i,l) = 1.
c
	  if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
            cdfi(i,l) = cd(i,l)/cdio(i,l)
            cdfo(i,l) = cdow(i,l)/cdio(i,l)
          end if
        end do
      end do
c
#endif
      if (.not.bndry) then
c
c preliminaries for computation of diffusion coefficients
        do 400 k=kvp,kx
          do 400 i=1,ix
            sc(i,k-kv) = tmj(i,k-1) + tmj(i,k)
            sc(i,k-kv) = 1.0 / (sc(i,k-kv) * sc(i,k-kv) ** 2)
            sd(i,k-kv) = (umj(i,k) - umj(i,k-1)) ** 2
            se(i,k-kv) = (vmj(i,k) - vmj(i,k-1)) ** 2
  400   continue
c
c compute diffusion coefficients at half levels in free atmosphere
        do 410 l=1,lmx
          do 410 k=kvp,kx
#ifdef garratt
c
c determine stability adjustment
            call fstab1(k,l)
#endif
            do 410 i=1,ix
              dkmh(i,k,l)  = vrt2(k-1) * sc(i,k-kv)
     *                     * sqrt(sd(i,k-kv) + se(i,k-kv))
#ifdef garratt
              dkmh(i,k,l)  = dkmh(i,k,l) * stbadj(i)
#endif
  410   continue
#ifdef garratt
c
c compute bulk Richardson number for adjustment of surface drag
        do 415 i=1,ix
          tsfc = abs(ts(i,jrow)) - tfra
          tgrad = hm(i,kx,3) - tsfc
          spdsqr = umj(i,kx)**2 + vmj(i,kx)**2
          spdsqr = max(spdsqr,1.e-10)
          rib(i) = grav*zh(i,kx)*tgrad/((hm(i,kx,3)+tfra)*spdsqr)
  415   continue
c
c determine stability adjustment and adjust drag coefficients
        do 425 l=1,lmx
          call fstab2(cd(1,l),l)
          do 420 i=1,ix
            cd(i,l) = cd(i,l) * stbadj(i)
  420     continue
  425   continue
#endif
c
c compute ground flux coefficient
        do 500 l=1,lmx
          do 500 i=1,ix
            s3(i)        = vrt3 * cd(i,l) / tmj(i,kx)
            dkmh(i,kp,l) = s3(i) * sqrt(umj(i,kx)**2 + vmj(i,kx)**2)
c
c set lower limit for wind speed when sensible heat and evaporation are
c calculated; function min is used because dkmh < 0 (from negative vrt3)
            dkmin(i)     = spmin * s3(i)
            if (l.ge.3) then
              dkmh(i,kp,l) = min(dkmh(i,kp,l),dkmin(i))
            endif
  500   continue
c
#ifdef ice_dynamics
c
c for lead parameterization in sea ice so diffused quantities are cd x (cfi)
c
	do l = 1, lmx
          do i = 1, ix
            dkmh(i,kp,l) = dkmh(i,kp,l)/cd(i,l)*cdio(i,l)
          end do
        end do
#endif
      else
c
c
      endif
c
      do 550 k=kv,kx
        s4(k) = dtx2 * gg(k)
  550 continue
c
c   set up loop for variable l (u, v, temperature and mixing ratio)
c
      do 1000 l=1,lmx
#ifdef moist
        if (l .eq. lmx) then
          do 590 i=1,ix
#ifdef wet_ep_exp
            dkmh(i,kp,l) = dkmh(i,kp,l)
#else
            dkmh(i,kp,l) = dw(i) * dkmh(i,kp,l)
#endif
  590     continue
        endif
#endif
        if (l .eq. 2) then
          do 12 k=1,kx
            do 12 i=1,ix
              e(i,k,2) = e(i,k,1)
   12     continue
        else
c
c   compute terms for richtmyer implicit finite difference method
c
c  x, y and z are ak, bk and ck and w is the reciprocal of the denomina-
c    tor of ek and fk in the notes.
c
          do 700 k=kv,kx
            do 700 i=1,ix
              x(i,k)   = s4(k) * dkmh(i,k+1,l)
              z(i,k)   = s4(k) * dkmh(i,k,l)
              y(i,k)   = 1.0 + x(i,k) + z(i,k)
              w(i,k)   = 1.0 / (y(i,k) - z(i,k) * e(i,k-1,l))
              e(i,k,l) = x(i,k) * w(i,k)
  700     continue
        endif
c
c   iterate downward in atmosphere to compute term f for each variable
c
        do 800 k=kv,kx
          do 800 i=1,ix
#ifdef ice_dynamics
c
c put in multiplier in diffused quantity for lead parameterization (cfi)
c
	    f(i,k,l) = (hm(i,k,l)*cdio(i,l)+z(i,k)*f(i,k-1,l))*w(i,k)
#else
            f(i,k,l) = (hm(i,k,l) + z(i,k) * f(i,k-1,l)) * w(i,k)
#endif
  800   continue
        if (l .le. 2) then
          do 1150 i=1,ix
            prox(i,kx,l) = f(i,kx,l)
#ifdef ice_dynamics
c
c divide by additional drag multipler for calculation of ground flux (cfi)
c
            flux(i,l) = dkmh(i,kp,l)/cdio(i,l)*prox(i,kx,l)*gvbrps(i,ih)
#else
            flux(i,l)    = dkmh(i,kp,l) * prox(i,kx,l) * gvbrps(i,ih)
#endif
 1150     continue
#ifdef ice_dynamics
c
c calculate wind stress vectors needed in ice dynamics (cfi)
c
	  if (l .eq. 2) then
            do 19 i = 1, ix
c
c wind stress on ice (cfi)
c
              tauxi = -prox(i,kx,1)*gvbrps(i,ih)*dkmh(i,kp,l)/cdio(i,l)
              tauyi = -prox(i,kx,2)*gvbrps(i,ih)*dkmh(i,kp,l)/cdio(i,l)
c
              tauxic(i,jrow) = tauxic(i,jrow) + tauxi
              tauyic(i,jrow) = tauyic(i,jrow) + tauyi
              windvx(i,jrow) = prox(i,kx,1)/cdio(i,l) + windvx(i,jrow)
              windvy(i,jrow) = prox(i,kx,2)/cdio(i,l) + windvy(i,jrow)
  19       continue
          end if
c
#endif
c
#if defined cpldleg3 || defined cpldleg1
c   compute wind stress components
c    taux & tauy represent the drag of the earth's surface on the
c     atmosphere.
        if (l .eq. 2) then
          do 801 i=1,ix
          taux(i,jrow)=taux(i,jrow) +
     *                 prox(i,kx,1)*gvbrps(i,ih)*dkmh(i,kp,1)
          tauy(i,jrow)=tauy(i,jrow) +
     *                 prox(i,kx,2)*gvbrps(i,ih)*dkmh(i,kp,2)
# ifdef windbc
c   compute scalar wind speed over the sea in cm/sec.
c   prox(i,kx,1) is u-component
c   prox(i,kx,2) is v-component
          swind(i) = sqrt(prox(i,kx,1)**2 + prox(i,kx,2)**2 + eps)
          winds(i,jrow) = winds(i,jrow) + swind(i)
          wind3(i,jrow) = wind3(i,jrow) - s3(i)*gvbrps(i,ih)
     $                   *swind(i)**3
c   s3(i) is minus because of -vrt3 in s3(i); see do 500 loop above,
c    and parmtr
# endif
  801     continue
        endif
#endif
c
c   march up atmosphere to compute implicit winds in boundary layer
c
          do 1200 k=1,kxmv
            do 1200 i=1,ix
              prox(i,kx-k,l) = e(i,kx-k,l) * prox(i,kp-k,l)
     *                       + f(i,kx-k,l)
 1200     continue
c
c   compute wind tendencies from time differences divided by 2*dt
c
          do 1400 k=kv,kx
            do 1400 i=1,ix
#ifdef ice_dynamics
c
c allow for different surface in ice points (cfi)
c
	      vf(i,k,l) = rdtx2*(prox(i,k,l)/cdio(i,l)-hm(i,k,l))
#else
              vf(i,k,l) = rdtx2 * (prox(i,k,l) - hm(i,k,l))
#endif
 1400     continue
        endif
c
c   zero out tendencies at undiffused levels
c
        do 1600 k=1,kvm
          do 1600 i=1,ix
            vf(i,k,l) = 0.0
 1600   continue
        if (l .eq. lmx) then
c
c  after the last variable is done, store e and f at k=kx in common
c    block /surfc/ for temperature and mixing ratio for use by
c    subroutine tsc or tscs.
c
          do 1020 i=1,ix
            ekd (i) = e(i,kx,3)
#ifdef ice_dynamics
c
c correct for different surface type at ice points (cfi)
c
	    fkd(i) = f(i,kx,3) + tfra*cdio(i,l)*(1.0-ekd(i))
#else
            fkd (i) = f(i,kx,3) + tfra * (1.0 - ekd(i))
#endif
#ifdef wet_ep_exp
            ekdr(i) = 1.0/(1.0 - (1.0-dw(i)) * e(i,kx,4))
            fkdr(i) = ekdr(i) * f(i,kx,4)
            ekdr(i) = ekdr(i) * e(i,kx,4) * dw(i)
#else
            ekdr(i) = e(i,kx,4)
            fkdr(i) = f(i,kx,4)
#endif
 1020     continue
        endif
 1000 continue
c
      call tsc
c
c   compute tendencies of temperature and mixing ratio
c
      do 3000 l=3,lmx
        if (l .eq. 3) then
c
c  convert surface temperature t* to celsius
          do 3830 i=1,ix
            tdegc(i) = abs(ts(i,jrow)) - tfra
#ifdef ice_dynamics
	    tdegccd      = abs(tscd(i,jrow)) - tfra*cdio(i,l)
            prox(i,kx,l) = f(i,kx,l) + tdegccd * e(i,kx,l)
            flux(i,l)    = (prox(i,kx,l) - tdegc(i))
     *                   * dkmh(i,kp,l)/cdio(i,l) * gvbrps(i,ih)
#else
            prox(i,kx,l) = f(i,kx,l) + tdegc(i) * e(i,kx,l)
            flux(i,l)    = (prox(i,kx,l) - tdegc(i))
     *                   * dkmh(i,kp,l) * gvbrps(i,ih)
#endif
#ifdef ice_dynamics
c
c compute fluxes over open water and ice for ice points (cfi)
c
	    if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
               proxcdi = prox(i,kx,l)*cdfi(i,l)*area(i,jrow,2)
               proxcdo = prox(i,kx,l)*cdfo(i,l)*(1.-area(i,jrow,2))
               tdegccd      = abs(tscd(i,jrow)) - tfra*cdio(i,l)
               flux(i,l) = (proxcdi+proxcdo-tdegccd)
     >                    *dkmh(i,kp,l)/cdio(i,l)*gvbrps(i,ih)
            end if
#endif
 3830     continue
        else
#ifdef wet_ep_exp
          do 3840 i=1,ix
            tdegc(i) = abs(ts0(i)) - tfra
 3840     continue
#endif
#ifdef ice_dynamics
	  do i = 1, ix
            if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
               tdegc(i) = abs(tsice(i,jrow)) - tfra
            end if
          end do
#endif
c
c  derive saturation vapor pressures for the t* values.
          call esndif(tdegc, ess, dm, -ix)
          do 4100 i=1,ix
            s1(i)       = 0.622 * ess(i) / (gvbrps(i,ih) - ess(i))
#ifdef wet_ep_exp
            s1(i)       = (dw(i) + (1.0-dw(i)) * ekdr(i)) * s1(i)
     *                           + (1.0-dw(i)) * fkdr(i)
#endif
c
#ifdef ice_dynamics
c
c calculate surface mixing ratio for flux computation (cfi)
c
	    if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
               s271 = 0.622*es271/(gvbrps(i,ih)-es271)
               ri = s1(i)
               s1(i) = area(i,jrow,2)*s1(i)*cd(i,l)
     >                 + (1.-area(i,jrow,2))*s271*cdow(i,l)
            else
               s1(i) = s1(i)*cd(i,l)
            end if
#endif
#ifdef ice_dynamics
            prox(i,kx,l) = f(i,kx,l) + s1(i) * e(i,kx,l)
            flux(i,l)    = (prox(i,kx,l) - s1(i))
     *                    * dkmh(i,kp,l)/cdio(i,l) * gvbrps(i,ih)
c
c compute fluxes over open water and ice for ice points (cfi)
c
          if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
c            proxcdi = prox(i,kx,l)*cdfi(i,l)*area(i,jrow,2)
c            proxcdo = prox(i,kx,l)*cdfo(i,l)*(1.-area(i,jrow,2))
c            flux(i,l) = (proxcdi+proxcdo-s1(i))
c    >                  *dkmh(i,kp,l)/cdio(i,l)*gvbrps(i,ih)
c
c separate fluxes over open water and ice for ice points (cfi)
c
             proxcdi = prox(i,kx,l)*cdfi(i,l)
             proxcdo = prox(i,kx,l)*cdfo(i,l)
c
c vapori, vaporo are moisture fluxes over open water and ice (cfi)
c
             vapori(i) = (proxcdi-ri*cd(i,l))*area(i,jrow,2)
     >                  *dkmh(i,kp,l)/cdio(i,l)*gvbrps(i,ih)
             vaporo(i) = (proxcdo-s271*cdow(i,l))*(1.-area(i,jrow,2))
     >                  *dkmh(i,kp,l)/cdio(i,l)*gvbrps(i,ih)
             flux(i,l) = vapori(i) + vaporo(i)
	  end if
#else
	    prox(i,kx,l) = f(i,kx,l) + s1(i) * e(i,kx,l)
            flux(i,l)    = (prox(i,kx,l) - s1(i))
     *                    * dkmh(i,kp,l) * gvbrps(i,ih)
#endif
 4100     continue
        endif
c
        do 4200 k=1,ksm
          do 4200 i=1,ix
            prox(i,kx-k,l) = e(i,kx-k,l) * prox(i,kp-k,l) + f(i,kx-k,l)
 4200   continue
        do 4400 k=kv,kx
          do 4400 i=1,ix
#ifdef ice_dynamics
c
c allow for different surface types in ice points (cfi)
c
	    vf(i,k,l) = rdtx2*(prox(i,k,l)/cdio(i,l)-hm(i,k,l))
#else
            vf(i,k,l) = rdtx2 * (prox(i,k,l) - hm(i,k,l))
#endif
 4400   continue
 3000 continue
c
c the following section computes the mountain stress due to
c gravity wave drag according to a parameterization developed
c by y. hayashi.  the contribution of gravity wave drag is
c then added to the wind tendencies due to vertical diffusion.
c
c  jump around the gravity wave drag code if rk0=0
c
      if (rk0 .eq. 0.0) go to 5556
c
c
c convert potential temperatures to degrees kelvin
      do 4500 k=1,kx
        do 4500 i=1,ix
          hm(i,k,3) = hm(i,k,3) + tfra
 4500 continue
c
c compute low-level mean u, v, t, & theta by taking
c pressure-weighted averages from levels kll to kx
      do 4600 i=1,ix
         ulow(i) = 0.0
         vlow(i) = 0.0
         tlow(i) = 0.0
        thlow(i) = 0.0
 4600 continue
c
      rdqllh = 1./( qmh(kx) - qmh(kll-1) )
      do 4700 k=kll,kx
        do 4700 i=1,ix
          ulow (i) = ulow (i) + dq(k)*umj(i,k)  *rdqllh
          vlow (i) = vlow (i) + dq(k)*vmj(i,k)  *rdqllh
          tlow (i) = tlow (i) + dq(k)*tmj(i,k)  *rdqllh
          thlow(i) = thlow(i) + dq(k)*hm (i,k,3)*rdqllh
 4700 continue
c
c compute reciprocal of low-level wind magnitude
c
      do 4710 i=1,ix
        ruvmag(i) =  ulow(i)**2 + vlow(i)**2
 4710 continue
c
      do 4715 i=1,ix
       if (ruvmag(i) .ne. 0.0) then
         ruvmag(i) = 1.0 / sqrt(ruvmag(i))
       endif
 4715 continue
c
c compute reciprocal of wind magnitude at other levels
      do 4720 k=1,km
        do 4720 i=1,ix
          theta(i,k) = umj(i,k)**2 + vmj(i,k)**2
 4720 continue
c
      do 4725 k=1,km
       do 4725 i=1,ix
        if (theta(i,k) .ne. 0.0) then
          theta(i,k) = 1.0 / sqrt(theta(i,k))
        endif
 4725 continue
c
c compute dot product of wind at each level and low-level wind
      do 5515 k=1,km
        do 5515 i=1,ix
          theta(i,k) = theta(i,k) * ruvmag(i)
     *           * ( ulow(i) * umj(i,k) + vlow(i) * vmj(i,k) )
 5515 continue
c
      do 5525 k=1,km
       do 5525 i=1,ix
        theta(i,k) = min(theta(i,k), 0.99999)
        theta(i,k) = max(theta(i,k),-0.99999)
 5525 continue
c
c determine angle between wind at each level and low-level wind
      do 5530 k=1,km
       do 5530 i=1,ix
        theta(i,k) = acos(theta(i,k))
 5530 continue
c
c in determining critical layer, include only levels
c at which this angle is less than 90 degrees
      do 5535 k=1,km
       do 5535 i=1,ix
        if (theta(i,k) .gt. hpi) then
          test(i,k) = 0.0
        else
          test(i,k) = 1.0
        endif
 5535 continue
c
      do 5100 k=kll,kx
        do 5100 i=1,ix
          test(i,k) = 1.0
 5100 continue
c
c since critical layer begins at the surface and must be continuous,
c eliminate any discontinuous levels above the surface-based layer
      do 5519 i=1,ix
        do 5516 k=kll-1,1,-1
          klast = k
          if ( test(i,k) .ne. 1.0 )  go to 5517
 5516   continue
        klast = 0
 5517   continue
        if ( klast .ne. 0 )  then
          do 5518 k=1,klast
            test(i,k) = 0.0
 5518     continue
        endif
        kcrit(i) = klast + 1
 5519 continue
c
c compute low-level brunt-vaisala frequency
      rdqll = 1./( q(kx) - q(kll) )
      do 5000 i=1,ix
        bvfreq(i) = (hm(i,kx,3) - hm(i,kll,3))*rdqll
 5000 continue
c
      qllbar = ( qmh(kll-1) + qmh(kx) ) / 2.
      do 5050 i=1,ix
        bvfreq(i) = (-1.*grav*grav*qllbar*bvfreq(i)) /
     $     (rgas*tlow(i)*thlow(i))
 5050 continue
c
      do 5541 i=1,ix
       if (bvfreq(i) .gt. 0.0) then
         bvfreq(i) = sqrt(bvfreq(i))
       else
         bvfreq(i) = 0.0
       endif
 5541 continue
c
c add effect of gravity wave drag to wind tendencies
      do 5555 k = 1,kx
        do 5555 i = 1,ix
          vfu(i,k) = vfu(i,k) -
     *       grav*rk0*bvfreq(i)*zsvar(i,jrow)*ulow(i)*test(i,k)
     *             / (rgas * tlow(i) * (1.-q(kcrit(i))))
          vfv(i,k) = vfv(i,k) -
     *       grav*rk0*bvfreq(i)*zsvar(i,jrow)*vlow(i)*test(i,k)
     *             / (rgas * tlow(i) * (1.-q(kcrit(i))))
 5555 continue
c
 5556 continue
c
c  end of gravity wave drag code.
c
c
c  compute zonal mean integrals of torque, heat flux and evaporation
c
       zmdrag(jrow) = 0.0
       zmhtfx(jrow) = 0.0
       zmevap(jrow) = 0.0
      do 81 i=1,ix
       zmdrag(jrow) = zmdrag(jrow) + flux(i,1)
       zmhtfx(jrow) = zmhtfx(jrow) + flux(i,3)
       zmevap(jrow) = zmevap(jrow) + flux(i,4)
   81 continue
      do 82 i=1,ix
        fluxw(i) = dt * flux(i,4)
#ifdef ice_dynamics
c
c evaporation in ice grid points for precipitation in convad (cfi)
c
        if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
           fluxw(i) = dt * vaporo(i)
        end if
c
#endif
   82 continue
      do 83 i=1,ix
#ifdef ice_dynamics
c
c account for different moisture flux types in ice grid points (cfi)
c
	if (pkice(i).gt.0. .and. zin(i).eq.0. .and. .not.swamp) then
           evapor(i,jrow) = evapor(i,jrow) + dt * vaporo(i)
           sublim(i,jrow) = sublim(i,jrow) + dt * vapori(i) * s2(i)
c
c add sublimation to sno here instead of in convad subroutine (cfi)
c
	   sno(i,jrow) = sno(i,jrow) - vapori(i)*dt
c
	else
           evapor(i,jrow) = evapor(i,jrow) + fluxw(i)
           sublim(i,jrow) = sublim(i,jrow) + fluxw(i) * s2(i)
        end if
#else
        evapor(i,jrow) = evapor(i,jrow) + fluxw(i)
        sublim(i,jrow) = sublim(i,jrow) + fluxw(i) * s2(i)
#endif
   83 continue
c
c  compute dissipation of k.e. due to vertical diffusion
c
      do 1042 i=1,ix
        s1(i) = gzs(i,ih) * gpsdx(i,ih) * gvbrps(i,ih)
 1042 continue
      do 1044 k=kv,kx
        do 1044 i=1,ix
          w(i,k) = dq(k) * gvbrps(i,ih) * (gucos(i,k,ih) * vfu(i,k)
     *                                     + gvcos(i,k,ih)*vfv(i,k))
 1044 continue
      zmdrag(jrow) = rix * zmdrag(jrow) * rsekan
      zmhtfx(jrow) = rix * zmhtfx(jrow) * cpp * 60.0
      zmevap(jrow) = rix * zmevap(jrow) * dt
      zmcnvr(jrow) = rix * zmcnvr(jrow)
      zmvdis(jrow) = 0.0
      zmt   (jrow) = 0.0
      do 85 i=1,ix
       zmt(jrow) = zmt(jrow) + s1(i)*rsekan*rix
   85 continue
      do 89 k=kv,kx
       do 89 i=1,ix
        zmvdis(jrow) = zmvdis(jrow) +  w(i,k)*sekan*rix
   89 continue
c
      return
      end
