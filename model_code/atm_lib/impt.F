      subroutine impt(icf)
c
c ****** subroutine to fetch correct latitude data for radiation comp
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
#include "parmrad.h"
      parameter ( nwrds = kp * (kp + 3) * ix )
#ifdef dirnal
      parameter ( kwrds = ix * (kp + 2) )
#endif
c
cSGI
      integer*8 ione
cSGI
c
      dimension pr2i(kx), tranlw(nwrds)
#ifdef dirnal
      dimension srheat(kwrds)
#endif
c
      dimension psq(ix,kx), rsmin(ix,kx)
c
c ****** mean cosine of zenith angle
c ****** and daylight fraction versus longitude
c ****** (x-coordinate) for diurnal variation
c ****** model option dimensioned below
c
#ifdef dirnal
      dimension coszrs(ix),frac(ix)
#endif
#ifdef prdcld
# ifdef covcld
      dimension ifcv(ix,kx), ifcd(ix,kx), nclv(kp), ncld(kp)
c
c  ifcv is the unpacked convective flag bits
c  ifcd is the unpacked unconvective (condensation) flag bits
c  nclv is a counter for the convective flag bits
c  ncld is a counter for the condensation flag bits
c
# endif
#endif
c
#include "albedo.h"
c
#include "albsea.h"
c
#include "anal2.h"
c
c ****** astro is the astronomy common block
c ****** needed for diurnal models for
c ******  communication of astronomical
c ******  parameters from subroutine march
c
#include "astro.h"
c
#include "auxatm.h"
c
#include "cabd.h"
c
#include "cld.h"
c
#include "cld1.h"
c
#include "cldals.h"
c
#include "cldhts.h"
c
#ifdef prdcld
#ifdef covcld
#include "cloudcx.h"
#endif
#endif
c
#include "cloudp.h"
c
#include "ctlblk.h"
c
#include "establ.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "grid.h"
c
#include "iohist.h" 
c
#include "intgrl.h" 
c
#include "latdat.h"
c
#include "lats.h"
c
#include "maimc.h"
c
#include "param.h"
c
#include "rdintl.h"
c
#include "sgwgt1.h"
c
#include "surfc.h"
c
#include "swnuin.h"
c
#include "swnuot.h"
c
#include "trans.h"
c
#ifdef prdcld
      dimension klouds(ix,il), clouds(ix), ncl(kp)
#endif
c
#include "zmflux.h"
c
      equivalence (tranlw(1), transm(1,0,0))
#ifdef dirnal
      equivalence (srheat(1), rheat (1,1  ))
#endif
c
#ifdef prdcld
c  bbcdl is the cirrus blackness coef.
c   bbcdl = 1 => full black
c   bbcdl = x => partially black, where x is the fraction of blackness
# ifdef full_black_cirrus
      data bbcdl / 1.0 /
# endif
# ifdef part_black_cirrus
      data bbcdl / 0.6 /
# endif
# ifndef covcld
      data crz / 1.0 /
c
c  crz is the cloud fraction when there are clouds.
c
# endif
# ifdef covcld
c
c  crx is the cloud fraction for the convective clouds.
c  crz is the cloud fraction for the nonconvective 
c                           (condensation) clouds.
c  crx is a tunable parameter, this number was chosen so that
c  for a R30L9 model the total cloud was ~55% and the net
c  radiation imbalance at the top of the atm was ~0 w/m**2
c
      data crx / 0.30 /
      data crz / 1.0 /
# endif
#endif
c
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
#ifdef dirnal
c
      if (.not. kdisk) go to 2006
c
c ****** fetch hour angle of sunset for appropriate latitude (row)
c
      hr  = ha(jrow)
      dhr = float(intrd) / 3600.0
      dhr = amin1(dhr,24.0)
c
c ****** call zenith to compute coszrs and frac
c ****** from date and astro. parameters
c
      call zenith(fjd,rvct,dlt,slag,xlat(jrow),hr,dhr,ix,coszrs,frac)
c
c ****** integrate cosine of zenith angle and daylight fraction over
c ******  entire latitude circle to obtain zonal means which are
c ******  approximately the daily means for any point on this latitude
c
      ccosz (jrow) = 0.0
      ttauda(jrow) = 0.0
c
      do 2005 i=1,ix
        ttauda(jrow) = ttauda(jrow) + frac(i)
        ccosz (jrow) = ccosz (jrow) + frac(i) * coszrs(i)
 2005 continue
c
      if (ttauda(jrow) .le. 0.0) go to 2006
      ccosz (jrow) = ccosz (jrow) / ttauda(jrow)
      ttauda(jrow) = ttauda(jrow) / float(ix)
 2006 continue
#endif
c
c ****** read in transmission data from disk
c
      if (ldisk) go to 65
#ifdef unitsk
      if (jrow .eq. 1) rewind idsk
      read(idsk) tranlw
#endif
#ifndef unitsk
        do kk=0,kx
         do k=0,kx
          do i=1,ix
           transm(i,k,kk) = glbtrs(i,k,kk,jrow)
          end do
         end do
        end do
        do kk=1,2
         do k=0,kx
          do i=1,ix
           absrbx(i,k,kk) = glbabs(i,k,kk,jrow)
          end do
         end do
        end do
        do k=0,kx
         do i=1,ix
          transz(i,k) = glbtrz(i,k,jrow)
         end do
        end do
#endif
   65 continue
c
c ****** read solar radiation tendencies from disk on days when solar
c ******   radiation is computed
c
#ifdef dirnal
      if (kdisk) go to 75
      read(kdsk) srheat
   75 continue
#endif
c
c * interpolate radiation data. annual mean in both hemis. * * * * * * *
c
      fl = (radian/10.0) * radang(ls)
c
c the constant 5.729578 = (360./(2.*pi))/10.0 = radian/10
c
      if (jrow .le. iy) fl = 18.0 - fl
      ll   = fl
      l    = fl + 0.5
      than = fl - float(ll)
      ll   = ll + 1
      l    = l + 1
c
c  compute the ocean albedo for the current row by interpolating from
c   the table of zonal mean ocean albedoes.
c
      seaalb = oas(ll) + than * ( oas(ll+1) - oas(ll) )
c
c ****** set up appropriate cloud cover for this latitude
c
#ifdef fixcld
      do 200 kc=1,3
        kth  (kc+1) = kkth (l,kc)
        kbh  (kc+1) = kkbh (l,kc)
        cldsw(kc+1) = ccd  (l,kc)
        cldlw(kc+1) = cldsw(kc+1)
  200 continue
#endif
#ifdef prdcld
c  
c  if icf =2 => set clouds = 0
      if (icf .eq. 2) then
        do 206 k=1,kx
         do 206 i=1,ix
          iflg(i,k) = 0
# ifdef covcld
          ifcv(i,k) = 0
# endif
  206   continue
      else
cSGI
        ione = 1
cSGI
        do 201 k=1,kx
         do 201 i=1,ix
cSGI
cc        iflg(i,k) = and(1,kflags(i,jrow)/ipower(k))
          iflg(i,k) = and(ione,kflags(i,jrow)/ipower(k))
cSGI
# ifdef covcld
cSGI
cc        ifcv(i,k) = and(1,kfgcov(i,jrow)/ipower(k))
          ifcv(i,k) = and(ione,kfgcov(i,jrow)/ipower(k))
cSGI
# endif
  201   continue
      endif
# ifdef covcld
      do 30 k=1,kx
      do 30 i=1,ix
      ifcd(i,k) = iflg(i,k)
   30 continue
      do 31 k=1,kx
      do 31 i=1,ix
      if(ifcv(i,k) .ne. 0) ifcd(i,k) = 0
   31 continue
c
# endif
      do 202 k=0,kp
       do 202 i=1,ix
        ccover(i,k) = 0.0
  202 continue
c  Do not allow clouds at k=1 level (top of the atmosphere.
      do 203 k=2,kx
       do 203 i=1,ix
        if (iflg(i,k) .ne. 0) ccover(i,k) = crz
# ifdef covcld
        if (ifcd(i,k) .ne. 0) ccover(i,k) = crz
        if (ifcv(i,k) .ne. 0) ccover(i,k) = crx
# endif
  203 continue
#endif
c
#ifdef prdcld
c
c  compute cloud height boundarys to be used in albcld
c  cldhtb is an entry point in albcld
      call cldhtb(ylat(jrow))
c
c  set cirrus fraction for long wave calculation
c  note ccover is reset only for thick (at least 3 levels)
c       clouds. cldhm is the boundary in sigma units between
c       hi and middle clouds.
c
      do 108 k=1,kx
       if (q(k) .le. cldhm) kcirus = k
  108 continue
      do 109 k=2,kcirus
       do 109 i=1,ix
        if ( ccover(i,k  ).gt.0.0 .and.
     *       ccover(i,k+1).eq.0.0 .and.
     *       ccover(i,k-1).eq.0.0 )
     *       ccover(i,k  ) = bbcdl*crz
  109 continue
#endif
c
c ****** compute mixing ratios for this latitude
c
      do 300 k=1,kx
#ifdef dry
        r(k,1) = rclim(ll,k) + than * (rclim(ll+1,k) - rclim(ll,k))
        r(k,1) = amax1(r(k,1),3.0e-6)
#endif
        r(k,2) = rco2
        r(k,3) = dduo3n(ll,k) + than * ( dduo3n(ll+1,k) - dduo3n(ll,k) )
  300 continue
      r(kp,2) = r(kx,2)
      r(kp,3) = r(kx,3)
#ifdef o3_k1_t80
c  reduce ozone amt by 20 percent at k=1 to improve calc temps there
      r( 1,3) = r( 1,3) * 0.80
#endif
c
c ****** compute tauda and cosz
c
#ifdef annlmn
      tauda = ttauda(ll) + than * (ttauda(ll+1) - ttauda(ll))
      cosz  = ccosz (ll) + than * (ccosz (ll+1) - ccosz (ll))
#endif
#ifdef seanal
      tauda = ttauda(jrow)
      cosz  = ccosz (jrow)
#endif
c
c ****** note: tauda and cosz are redefined below in diurnal models
c
c
c compute solar zenith angle and load appropriate cloud albedos
      zangle = acos(cosz)*57.2957796
c
      call cldalb(zangle, zza)
c
#ifdef fixcld
c
c  set dum1 and dum2 as dummy values, in fixcld models they
c  are only needed to pass as dummy arguments.
      ndum1 = 999
      dum2  = 999
      call albcld(ndum1, dum2)
c
#endif
c
      do 350 ipt = 1,ix
        ralb(ipt) = abs(alb(ipt,jrow))
  350 continue
#ifdef prdcld
      st(1) = 0.11
      st(2) = 0.0
#endif
c
c ----------------------------------------------------------------------
c  call lwcool to compute longwave cooling rates
c ----------------------------------------------------------------------
c
      call lwcool (cgt, gqmix, gvbrps, icf)
c
c ----------------------------------------------------------------------
c  begin main loop on ipt
c ----------------------------------------------------------------------
c
c ****** first compute sqrt(q) for swnew field (used to be in parmtr)
c
      do 799 k=1,kx
        pr2i(k) = sqrt(q(k))
  799 continue
c
c compute pressure at individual levels and use to compute rsmin.
c   rsmin = minimum saturation mixing ratio on es table
c   esmin = minimum saturation vapor pressure on es table
c

      esmin  = etabl(1)
      do 519 k=1,kx
       do 519 i=1,ix
          psq(i,k) = gvbrps(i,ih) * q(k)
        rsmin(i,k) = 0.622*esmin / (psq(i,k)-esmin)
  519 continue
c
c ----------------------------------------------------------------------
c  begin main loop on ipt for sw code
c ----------------------------------------------------------------------
c
      do 800 ipt=1,ix
c
c
c ****** pss    = surface pressure (dynes cm-2)
c ****** zenith = zenith of angle of sun (degrees)
c ****** tauda  = (total daylight hours)/(24 hours)
c
#ifdef dirnal
c ****** set cosz and tauda at each gridpoint,
c ****** ipt, on row in diurnal models
c
        tauda = frac  (ipt)
        cosz  = coszrs(ipt)
#endif
        pss   = gvbrps(ipt,ih)
        pchg  = sqrt(pss / 1013250.0)
        do 243 k=1,kx
          pr2(k) = pr2i(k) * pchg
  243   continue
#ifdef moist
        do 21 k=1,kx
          r(k,1) = max(gqmix(ipt,k,ih),rsmin(ipt,k))
   21   continue
#endif
        r(kp ,1) = r(kx,1)
c
c determine surface albedo 
        call albsfc(ipt,seaalb)
c
        salb(ipt)      = ca(nc)
        calb(ipt,jrow) = ca(nc)
c
#ifdef dirnal
        if (kdisk) go to 84
        do 85 k=1,kp
          heat(k) = rheat(k,ipt)
   85   continue
        sng   = sing(ipt)
        sdtop = song(ipt)
        sdbtm = sng
   84   continue
#endif
#ifdef prdcld
c
c ****** variable cloud computation (no longer used by longwave code)
c
        ncl(kp) = 0
#ifdef covcld
        nclv(kp) = 0
        ncld(kp) = 0
#endif
        do 15 k=1,kx
          ncl(k) = iflg(ipt,k)
#ifdef covcld
          nclv(k) = ifcv(ipt,k)
          ncld(k) = ifcd(ipt,k)
#endif
   15   continue
        ncl(1)   = 0
        do 204 i=1,lx
         cca(i,1) = 0.0
         cca(i,2) = 0.0
         ccb(i,1) = 0.0
         ccb(i,2) = 0.0
          ca(i  ) = 0.0
          cb(i  ) = 0.0
  204   continue
        cldsw(1) = 1.0
        kth  (1) = 0
        kbh  (1) = 0
        ik       = 1
        nd       = 1
#ifndef covcld
    5   if (ncl(ik) .eq. 0) go to 1
        nd      = nd + 1
        kth(nd) = ik
    6   ik      = ik + 1
        if (ncl(ik) .eq. 1) go to 6
        kbh(nd) = ik - 1
    1   ik = ik + 1
        if (ik .le. kx) go to 5
#endif
#ifdef covcld
        do 1 k=2,kx
        if (nclv(k) .eq. 0 .and. ncld(k) .eq. 0) go to 1
        if (ncld(k) .eq. 1) go to 5
        if (nclv(k-1) .eq. 0 .and. nclv(k) .eq. 1) kth(nd+1) = k
        kbh(nd+1) = k
        if (nclv(k) .eq. 1 .and. nclv(k+1) .eq. 0) nd = nd + 1
        go to 1
    5   kth(nd+1) = k
        kbh(nd+1) = kth(nd+1)
        nd = nd + 1
    1   continue
#endif
        kth  (nd+1) = kp
        kbh  (nd+1) = kp
        cldsw(nd+1) = 1.0
        if (nd .eq. 1) go to 9
#ifndef covcld
        do 22 k=2,nd
          cldsw(k) = crz
   22   continue
#endif
#ifdef covcld
        do 2 k=2,nd
          cldsw(k) = crx
    2   continue
#endif
c
c  set the cloud albedos
        call albcld(nd, crz)
c
        do 8 k=2,nd
          ccb(k,1) = 0.0
    8   continue
    9   nccv     = nd
        ncv      = nccv + 1
#endif
#ifdef fixcld
        ncv      = nc
#endif
        ca(ncv)  = calb(ipt,jrow)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  the following section of code used to be the subroutine mn          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     p(k)  = pressure at data-input-level k
c     pp(k) = pressure at flux-level k
c     dp(k) = pressure difference between flux-levels just above and
c             below data-input-level k
c
        p (0) = 0.0
        pp(0) = 0.0
        p(kp) = pss
        do 34 k=1,kx
          p (k) = pss * q(k)
          pp(k) = pss * qmh(k)
   34   continue
c
#ifdef dirnal
        if (.not. kdisk) go to 150
#endif
c
c ****** input for swnew
c
        nc0 = ncv - 2
        nc1 = nc0 + 1
        nc2 = nc0 + 2
        nc3 = nc0 + 3
        nc4 = nc0 + 4
        kbhsw (  1) = 1
        cwca  (  1) = 0.0
        cwcb  (  1) = 0.0
        coca  (  1) = 0.0
        cloudy(  1) = 0.0
        cwca  (nc2) = ca(ncv)
        cwcb  (nc2) = 0.0
        coca  (nc2) = ca(ncv)
        cloudy(nc2) = 1.0
        if (nc0 .ne. 0) then
          do 141 kc=2,nc1
            cwca  (kc) = cca  (kc,2)
            cwcb  (kc) = ccb  (kc,2)
            coca  (kc) = cca  (kc,1)
            cloudy(kc) = cldsw(kc  )
  141     continue
          do 142 kc=2,nc1
            kthsw(kc) = kth(kc)
            kbhsw(kc) = kbh(kc) + 1
  142     continue
          do 143 kc=2,nc1
            kthsw(kc) = max0(kthsw(kc),kbhsw(kc-1))
  143     continue
        endif
#ifdef seanal
        solrsw = solc
#endif
#ifdef annlmn
        solrsw = solarc
#endif
        do 148 k=1,kp
          flx (k) = 0.0
          heat(k) = 0.0
  148   continue
        grdflx = 0.0
        if (cosz .gt. 0.0) then
          call swnew
        endif
c
c  the constant 1.43306e-6 converts erg/(cm**2*sec) to ly/min
c
        sdtop = -flx(1) * 1.43306e-06
        sdbtm = grdflx  * 1.43306e-06
        sng   = sdbtm
#ifdef dirnal
  150   continue
#endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  end of old subroutine mn                                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        if (icf .eq. 1) then
c
c  the constant 1.1574e-5 equals 1.0/86400.
          do 700 k=1,kx
            grad(ipt,k,ih) = grad(ipt,k,ih) + 1.1574e-5 * heat(k)
  700     continue
c
          slwdp(ipt,jrow) = slwdp(ipt,jrow) - sng
c
#ifdef dirnal
          if(kdisk) then
            do 86 k=1,kp
              rheat(k,ipt) = heat(k)
   86       continue
            sing(ipt) = sng
            song(ipt) = sdtop
          endif
#endif
c
          zsnit(ipt) = sdtop
          zsng (ipt) = sdbtm
#ifdef cloud_forcing
          xinsol = solc * tauda * cosz
          if (xinsol .ne. 0.) then
            zplalb(ipt) = (xinsol-sdtop)/xinsol
          else
            zplalb(ipt) = 0.
          endif
        else
c
          zsnitc(ipt) = sdtop
          zsngc (ipt) = sdbtm
          xinsol = solc * tauda * cosz
          if (xinsol .ne. 0.) then
            zplalbc(ipt) = (xinsol-sdtop)/xinsol
          else
            zplalbc(ipt) = 0.
          endif
#endif
        endif
  800 continue
c
c ----------------------------------------------------------------------
c  end main loop on ipt
c ----------------------------------------------------------------------
c
#ifdef prdcld
#ifndef covcld
      do 205 i=1,ix
       if (kflags(i,jrow) .gt. 0) then
         clouds(i) = crz
       else
         clouds(i) = 0.0
       endif
  205 continue
#endif
#ifdef covcld
      do 205 i=1,ix
       clouds(i) = 0.
       do 207 k=1,kx
        if (ifcd(i,k) .gt. 0) clouds(i) = crz
  207 continue
      do 208 k=1,kx
       if (ifcd(i,k) .eq. 0 .and. ifcv(i,k) .gt. 0) clouds(i) = crx
  208 continue
  205 continue
#endif
#endif
c
c ****** write transmission data to disk on full longwave timestep
c
      if (ldisk) then
#ifdef unitsk
        if (jrow .eq. 1) rewind idsk
        write(idsk) tranlw
#endif
#ifndef unitsk
        do kk=0,kx
         do k=0,kx
          do i=1,ix
           glbtrs(i,k,kk,jrow) = transm(i,k,kk)
          end do
         end do
        end do
        do kk=1,2
         do k=0,kx
          do i=1,ix
           glbabs(i,k,kk,jrow) = absrbx(i,k,kk)
          end do
         end do
        end do
        do k=0,kx
         do i=1,ix
          glbtrz(i,k,jrow) = transz(i,k)
         end do
        end do
#endif
      endif
c
c  write solar radiation tendencies on disk at appropriate location
c   defined by gmt index, idtim (in diurnal models).  this is done on
c   days when solar radiation is computed
c
#ifdef dirnal
      if (kdisk) write (kdsk) srheat
c
      do 850 ipt=1,ix
        zshfx (ipt) = shfx(ipt,jrow)
        zaltnt(ipt) = altnt(ipt,jrow)
  850 continue
#endif
c
      tsnit  = 0.0
      ttufxt = 0.0
      tsng   = 0.0
      tufxnt = 0.0
#ifdef prdcld
      tcldz  = 0.0
#endif
      do 851 i=1,ix
       tsnit  = tsnit  + zsnit (i) * rix
       ttufxt = ttufxt + ztufxt(i) * rix
       tsng   = tsng   + zsng  (i) * rix
       tufxnt = tufxnt + zufxnt(i) * rix
#ifdef prdcld
       tcldz  = tcldz  + clouds(i) * rix
#endif
  851 continue
      zmswtp(jrow) = tsnit
      zmswbm(jrow) = tsng
      zmlwtp(jrow) = ttufxt
      zmlwbm(jrow) = tufxnt
#ifdef prdcld
      zmtcld(jrow) = tcldz
#endif
c
c   hemis. mean shfx & altnt are computed in subr. tsc (see /intgrl/)
c
c  compute zonal mean clouds for printout below
c
#ifdef prdcld
      cldjrw(jrow) = zmtcld(jrow)
#endif
c
c  print table of daylight fraction, cosine of zenith angle and solar
c   radiation versus latitude averaged over the whole latitude circle
c   to replace daily mean table skipped in march in diurnal models
c
#ifdef dirnal
      if (kdisk .and. mod(ifcst+idt,86400).eq.0 .and. jrow.eq.1) then
        write(6,2001)
        do 2002 jj=1,il
          xinsol = solc * ttauda(jj) * ccosz(jj)
          write(6,2003) jj,ylat(jj),ttauda(jj),ccosz(jj),xinsol
 2002 continue
      endif
c
 2001 format('  row  latitude    fraction cosz/r**2  insolation')
 2003 format(1x,i3,3x,f8.3,3(2x,f10.6))
#endif
c
      return
      end
