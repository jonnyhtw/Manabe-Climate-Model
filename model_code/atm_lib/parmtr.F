      subroutine parmtr
c     ========== ======
c  this subroutine is called at the beginning of each run from atmos.
c  it calculates  constants, including those in the es table.
c  entry points:       cpoly, lgndre1, lgndre, lgndxy
c  subroutines called : o3init, gaussg
#if defined cpldleg3 || defined cpldleg1
c                       setgrd
#endif
c
#include "parm.h"
#include "parmncc.h"
#include "parmjx.h"
#include "parmkp.h"
#include "clnpol.h"
      parameter ( iym=iy-1 )
c
      real alp,epsi,en,x, zero
c
c
      complex vorm , divm , tm , psm
      complex vordt, divdt, tdt, psdt, ek
c
      dimension repsi(mxp,jxp),en (mxp),ren(mxp),consq(mxp)
      dimension table(5102),esnew(2551),esupc(200)
      dimension cypyx(mx2,jx), cypyy(mx2,jx)
      dimension alp(mxp,jx),paly(mxp,jx)
c
#include "cloudp.h"
c
#include "ctlblk.h"
c
#include "dispat.h"
c
#include "dxordy.h"
c
#include "establ.h"
c
#include "gauss.h"
c
#include "latdat.h"
c
#include "lats.h"
c
#include "lincoef.h"
c
#include "mixlen.h"
c
#include "nxstps.h"
c
#include "param.h"
c
      common / partmp / em   (mxp,jx ),emn(mxp,jx)
c
#include "polyls.h"
c
      common / polyn  / pol(mx,jx,iy)
c
#ifdef prdcld
#include "rhcrit.h"
#endif
c
#include "sgwgt1.h"
c
#include "tsmoth.h"
c
#ifdef ice_dynamics
c
c collect gaussian weights for interpolation (cfi)
c
#include "varice.h"
c
#endif
c
#include "vrtcom.h"
c
c
      equivalence ( cypyx(1,1),colydx(1))
      equivalence ( cypyy(1,1),colydy(1))
c
      data zero/0./
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
#ifdef new_levels
c
c compute qmh using a equation developed by ajb
      qmh(0) = 0.0
      qmh(kx) = 1.0
      do 1 k=1,kx-1
        fk = float(k)/float(kx)
        qmh(k) = 2.0473*fk**5 - 5.4128*fk**4 + 3.0268*fk**3 +
     $           1.0337*fk*fk + 0.3049*fk
    1 continue
#endif
c
c compute sigma values at model full levels
      do 41 k=1,kx
        q(k) = 0.5*(qmh(k)+qmh(k-1))
   41 continue
      write(6,9001) q
      write(6,9002) qmh
c
c  set fundamental atmospheric constants
c
c   compute ozone values for every 10 degrees of latitude and "kx"
c   levels. array dduo3n is filled for annual mean models and arrays
c   xo3sf and xo3ws are filled for seasonal models.
c
      call o3init
c
c  load common block /param/. for additional information,
c  see june 26,1981 notes.
      rgas   = 2.87   e 6
      cp     = 1.00464e 7
      calerg = 2.39   e-8
      ww     = 7.292  e-5
      grav   = 980.6
      sbc    = 8.132e-11
      cappa  = rgas/cp
      gm     = grav/cp
      cpp    =   cp*calerg
      ara    = rgas*calerg
      rqc9   = 1./q(kx)**cappa
      vk0    = 0.4
      z0land = 4.5
      z0sea  = 0.035
      vrt3   = -q(kx)/rgas
#ifdef sponge
c  specify the time constant for damping of winds at top of atm
c  note: dmptop is model dependent
      dmpdays = 0.0
# ifdef dmptop30
      dmpdays = 30.
# endif
# ifdef dmptop20
      dmpdays = 20.
# endif
      if (dmpdays .eq. 0.0) then
        write(6,9006)
        stop 'parmtr'
      endif
      dmptop = dmpdays * 86400.
#endif
#ifdef prdcld
      ipower(1)=1
# ifdef varhy
c
c  compute critical humidity for cloud/no cloud determination
c  in MSTCNV
c  hytop is the crit. RH for cloud formation at the top  of the atm.
c  hybot is the crit. RH for cloud formation at the bot. of the atm.
c
      hybot = 1.00
      hytop = 0.90

      rqkx = 1.0/q(kx)
      do 798 k=1,kx
        hyk(k) = hytop + (hybot - hytop)*q(k)*rqkx
  798 continue
      write (6,787) hyk
# endif
      do 789 j=1,il
       do 789 k=1,kx
c
c  hcrit is the critical relative humidity for cloud formation
c        in non-grid mix models, the magic number is 99%
        hcrit(k,j) = .99
  789 continue
#endif
      factq (1)=0.0
      eradsq = erad * erad
      rasq   = 1.0 / eradsq
      tr2    = 2.0*rasq
      tw     = 2.0*ww
      og     = 1.0/grav
      o2g    = 0.5*og
      absvor = tw/sqrt(1.5)
      rix    = 1.0/float(ix)
      intbal = 86400
      iseg   = 86400
      hmf    = float(3-nhem)
c
c  load dissipation common block /dispat/
      viscxy = 1.0e24
      tconxy = 1.0e24
      qconxy = 1.0e24
      om2al  = 1.0-2.0*alpha
c
#ifdef impcit
c  an implicit timestep is 5 times longer than an explicit timestep.
      if (implct) dt=5.0*dt
#endif
      dtx2 = 2.0*dt
      dtsq = dt**2
c
c  compute various coef and store them in lincoef for use in
c  linear.
c
      do 12 n=1,jx
        do 12 m=1,mx
          elcf(m,n)  = rasq * float((n+isc*(m-1)-1) * (n+isc*(m-1)))
          vxy(m,n)   = -viscxy * (tr2 - elcf(m,n)) ** 2
          vxy(m,n)   = vxy(m,n) / (1.0 - dtx2 * vxy(m,n))
          txy(m,n)   = -tconxy * elcf(m,n)**2
          txy(m,n)   = txy(m,n) / (1.0 - dtx2 * txy(m,n))
          qxy(m,n)     = -qconxy * elcf(m,n)**2
          qxy(m,n)     = qxy(m,n) / (1.0 - dtx2 * qxy(m,n))
   12 continue
c
c Determine the levels and weights to be used for discriminating
c between rain and snow.  The 350m level is used as the critical
c level, and q=0.9575 is used as an approximation to this level
c based on the U. S. standard atmosphere.
      qrs = 0.9575
      do 71 k=1,kx-1
        if (q(k).le.qrs .and. q(k+1).gt.qrs) then
          krs1 = k
          krs2 = k+1
          rswt1 = (q(k+1)-qrs)/(q(k+1)-q(k))
          rswt2 = 1.0 - rswt1
        endif
   71 continue
c
c The following section computes the mixing lengths required by
c subroutine vertdf.
c -------------------------------------------------------------
c
c Define some constants:
c  salaps = standard atmospheric lapse rate
c  tsfc   = idealized global mean surface temperature
c  htcrit = height at which mixing length is a maximum
c  htop   = height at which mixing length vanishes
c
      salaps = -6.5e-5
      tsfc   = 288.16
#ifndef garratt
      htcrit = 7500.0
      htop   = 0.0
c  adjust htop for tuning of cloud properties
c  (default htop is 2500 m)
# ifdef elmix4400
      htop   = 440000.0
# endif
# ifdef elmix4900
      htop   = 490000.0
# endif
      if (htop .eq. 0.0) then
        htop   = 250000.0
      endif
      write (6,9016) htop
#endif
c
c Compute the height of each full level assuming a constant
c standard lapse rate.
c WARNING: These should be used with caution.  They will
c  have large errors above the tropopause.
      do 131 k=1,kx
c
c In order to determine the height, the layer mean temperature
c from the surface to that level is required.  A surface 
c temperature of 15 deg Celsius and a standard lapse rate of 
c -6.5 deg/km will be used to estimate an average temperature 
c profile.  First, determine T at the appropriate full level.
        tcurr = tsfc*(q(k)**(-rgas*salaps/grav))
c
c Find the "layer average" temperature.
        tavg = 0.5*(tsfc+tcurr)
c
c Use this to determine the height of the half-level.
        zstd(k) = -rgas*tavg*alog(q(k))/grav
  131 continue
c
c Compute the height of each half level assuming a constant
c standard lapse rate using the above procedure.
c WARNING: These should be used with caution.  They will
c  have large errors above the tropopause.
      do 132 k=1,kx
        tcurr = tsfc*(qmh(k)**(-rgas*salaps/grav))
        tavg = 0.5*(tsfc+tcurr)
        zhstd(k) = -rgas*tavg*alog(qmh(k))/grav
  132 continue
c
c Compute the mixing length for the top kx-1 model levels.
      do 13 k=1,kx-1
#ifdef garratt
c
c The mixing length is determined from the height by specifying
c an asymptotic value (zasymp).  See Garratt, "The Atmospheric
c Boundary Layer" p. ???
        zasymp = 3000.0
        elmix(k) = (zasymp*vk0*zhstd(k))/(zasymp+vk0*zhstd(k))
#endif
#ifndef garratt
c
c The mixing length is determined from the height by assuming
c that it increases linearly with height to a value of 30m at
c 75m above the ground (htcrit), then decreases linearly from
c there to a value of zero at height (htop).  (See notes, p. IV-6)
        if (zhstd(k).le.htcrit) then
            elmix(k) = vk0*zhstd(k)
          elseif (zhstd(k).gt.htcrit .and. zhstd(k).lt.htop) then
            elmix(k) = vk0*htcrit*(htop-zhstd(k))/(htop-htcrit)
          elseif (zhstd(k).ge.htop) then
            elmix(k) = 0.0
          else
            print *, ' Bad height value.'
        endif
#endif
c
   13 continue
      write(6,9003) elmix
      write(6,9005) (zhstd(k),k=1,kx)
c
c  Compute the number of levels to be diffused (ks)
c
      ks = 0
      do 235 k=1,kx
       if (elmix(k) .ne. 0.0) ks = ks + 1
  235 continue
c  The number of levels is one more than elmix(k) > 0
      ks = ks + 1
      write(6,9007) ks
c
#ifdef garratt
c
c Compute Garratt's f for use in adjusting mixing coefficients
      do 133 k=1,kx-1
        zg = zhstd(k)
        dzg = zstd(k)-zstd(k+1)
        fgarr(k) = ((((zg+dzg)/zg)**(1.0/3.0)-1.0)**1.5) /
     $         ((zg**0.5)*(dzg**1.5))
  133 continue
#endif
c
c  Compute the k level for the top of boundary layer for use in
c  the gravity wave computation in VERTDF. This code tries to
c  find the level closest to crtsig.
c
c  crtsig is the sigma level for the top of the boundary layer for
c         this computation (~850 mb).
      crtsig = 0.85
      do 225 k=kx,1,-1
       if (q(k) .lt. crtsig) then
         ksave = k
         go to 226
       endif
  225 continue
      stop 'parmtr'
  226 continue
      sg1dst = abs(crtsig - q(ksave  ))
      sg2dst = abs(crtsig - q(ksave+1))
      if (sg1dst .le. sg2dst) then
        kll = ksave
      else
        kll = ksave + 1
      endif
      write(6, 425) kll
c
c  Compute other variables used in VERTDF
      ksm  = ks - 1
      kv   = kx - ks + 1
      kvp  = kv + 1
      kvm  = kv - 1
      kxmv = kx - kv
c
c
c
c  compute level dependent constants
c  load common block /sigwgt/.  for additional information,
c  see october 19,1981 notes.
      do 10 k=1,kx
      dq    (k) = qmh(k)-qmh(k-1)
      odq   (k) = 0.5/dq(k)
      qlgc4 (k) = alog(qmh(k))
      akapsg(k) = 0.5*cappa/q(k)
      gg    (k) = -grav/dq(k)
      factr (k) = rgas*alog(qmh(k)/q(k))
      if (k.ne.1) then
#ifdef prdcld
        ipower(k) = 2*ipower(k-1)
#endif
        factq (k) = rgas*alog(q(k)/qmh(k-1))
      endif
 10   continue
      extkx  = factr(kx)
      extkxm = 0.0
      do 20 k=1,km
      dqph  (k  ) = q(k+1)-q(k)
      qlgc3 (k+1) = rgas*(qlgc4(k+1)-qlgc4(k))
      dqmqa5(k  ) = 0.5*(ara/cpp)*dqph(k)/qmh(k)
      den2  (k  ) =  1.0-dqmqa5(k)
      pq    (k  ) = (1.0+dqmqa5(k))/den2(k)
      dqi   (k  ) = 1.0/(dq(k+1)+dq(k  ))
      den   (k  ) = 1.0/(dq(k  )+dq(k+1)*pq(k))
      den1  (k  ) = dqmqa5(k)*546.32
      den2  (k  ) = den1(k)/den2(k)
      vrt4  (k  ) = 2.0*grav*qmh(k)/(rgas*dqph(k))
      vrt2  (k  ) =-((grav*elmix(k)/dqph(k))**2)*((2.*qmh(k)/rgas)**3)
 20   continue
c
c
 65   continue
c
      return
c
      entry cpoly
c     ===== =====
c  this entry point is called from atmos (before parmtr)
c
c  compute pi and friends
         pi  = 4.0*atan(1.)
        hpi  = pi/2.0
      twopi  = 2.0*pi
      pid24  = pi/24.0
      radian = 360.0/twopi
c
c  alphao is a constant used in LWTRAN
c        it is an empirical constant used in computing O3 absorption
      alphao = 0.28
      pialph = pi*alphao
c
c  compute latitude dependent variables in common block /gauss/
c
      call gaussg(iy,coa,weight)
c
#if defined cpldleg3 || defined cpldleg1 || defined ice_dynamics
      call setgrd
#endif
c
      erad  = 6.371e8
      radin = 1.0 / erad
      do 11 l=1,iy
      radang(l) = acos(coa(l))
      cosqr (l) = 1.0-coa(l)**2
      sia   (l) = sqrt(cosqr(l))
#ifdef ice_dynamics
c
c get gaussian weights for use in interpolation (cfi)
c
      weightin(l) = weight(l)
      weightin(il-l+1) = weight(l)
#endif
c
c  compute latitudes
      do 115 j=1,iy
       ylat(   j  ) = radian*(hpi-radang(j))
       ylat(il-j+1) = -ylat(j)
  115 continue
c
      weight(l) = weight(l)/cosqr(l)
 11   continue
c
c
c  ================================================================
c  +                                                              +
c  +             construction of the es table                     +
c  +                                                              +
c  + this table is constructed from es equations from the smithsonian
c  + tables. the es input  is computed from values  (in           +
c  + one-tenth of a degree increments) of es over ice from  -153c +
c  + to 0c  and values of es over water from 0c to 102c.   output +
c  + table contains these data interleaved with their derivatives +
c  + with   respect  to  temperature  except  between -20c and 0c +
c  + where  blended  (over water  and over ice)   es  values  and +
c  + derivatives are calculated.                                  +
c  note: all es computation is done in microbars and not millibars.
c  ================================================================
c
c
c
#ifdef moist
      esbasw = 1013246.0
       tbasw =     373.16
      esbasi =    6107.1
       tbasi =     273.16
c
c  compute es over ice between -153 c and 0 c.
c  see smithsonian meteorological tables page 350.
      do 120 i=1,1530
      tem = 120.16+0.1*float(i-1)
      aa  = -9.09718 *(tbasi/tem-1.0)
      b   = -3.56654 *alog10(tbasi/tem)
      c   =  0.876793*(1.0-tem/tbasi)
      e   = alog10(esbasi)
      table(i)=10**(aa+b+c+e)
 120  continue
c
c  compute es over water between -20c and freezing.
c  see smithsonian meteorological tables page 350.
      do 128 i=1,1221
      tem = 253.16+0.1*float(i-1)
      aa  = -7.90298*(tbasw/tem-1)
      b   =  5.02808*alog10(tbasw/tem)
      c   = -1.3816e-07*(10**((1-tem/tbasw)*11.344)-1)
      d   =  8.1328e-03*(10**((tbasw/tem-1)*(-3.49149))-1)
      e   = alog10(esbasw)
      esh20  = 10**(aa+b+c+d+e)
      if (i.le.200) then
        esupc(i) = esh20
      else
        table(i+1330) = esh20
      endif
  128 continue
c
c  derive blended es over ice and supercooled water between -20c and 0c
      do 130 i=1,200
      tem  = 253.16+0.1*float(i-1)
      wice = 0.05*(273.16-tem)
      wh2o = 0.05*(tem-253.16)
      table(i+1330) = wice*table(i+1330)+wh2o*esupc(i)
  130 continue
c
c  estimate es derivative in microbars per degree celsius by
c    differencing es values in the table
      do 134 i=2,2550
      esnew(i) = 5.0*(table(i+1)-table(i-1))
  134 continue
      esnew(   1) = esnew(   2)
      esnew(2551) = esnew(2550)
c
c  interleave es and its derivative into one table of alternating
c    variables (es entries are odd numbered, derivatives are even)
      do 135 i=1,2551
                  ii = (i-1)*2
      table(5101-ii) = table(2552-i)
  135 continue
      do 136 i=1,2551
            ii    = (i-1)*2
      table(ii+2) = esnew(i)
  136 continue
      do 140 i=1,5102
      etabl(i)=table(i)
  140 continue
#endif
c
c
      return
c
      entry lgndre1
c     ===== ======
c
c  this entry point is called from atmos once per model run to
c  compute legendre polynomials and store them in pol in /polyn/
c  see notes of april 18, 1977 and nov 3, 1981.
c
c  compute epsilon used in deriving the legendre polynomials and store
c  in common block /polyls/. see april 18, 1977 notes.
      epsi(1,1)=zero
      do 2 n=2,jx
      epsi(1,n) = sqrt(float((n-1)**2)/float(4*(n-1)**2-1))
 2    continue
      do 3 n=1,jx
      do 3 m=2,mxp
      epsi(m,n) = sqrt(float((n+m-2)**2-(m-1)**2)/float(4*(n+m-2)**2-1))
 3    continue
c     epsi(*,jxp) = 0.0
      do 96 m=1,mxp
       epsi(m,jxp) = 0.0
   96 continue
      sqrhlf=sqrt(.5)
c
c  consq = sqrt((2m+1)/2m)    where m=lower case m index
      do 266 m=2,mxp
      en   (m) = m-1
      ren  (m) = 1.0/en(m)
      consq(m) = sqrt(.5*(2.*en(m)+1.)*ren(m))
 266  continue
      en(1)=0.
c     repsi=0.
      do 267 n=1,jxp
       do 267 m=1,mxp
        repsi(m,n) = 0.0
  267 continue
      do 268 n=1,jxp
       do 268 m=1,mxp
        if (epsi(m,n) .gt. 0.0) then
          repsi(m,n) = 1.0/epsi(m,n)
        endif
  268 continue
c
c  compute coefficients in common block /dxordy/ for use by uvphrm.
c  see november 4, 1981 notes.
      do 8 mm=1,mx
      m           = isc*(mm-1)+1
      psdxm(mm  ) = radin*float(m-1)
      psdym(mm,1) = radin*float(m-2)*epsi(m,1)
      psdyp(mm,1) = radin*float(m+1)*epsi(m,2)
      uvdxm(mm,1) = -erad/float(m)
      uvdym(mm,1) = 0.0
      uvdyp(mm,1) = -erad*epsi(m,2)/float(m)
      do 8 n=2,jx
      psdym(mm,n) = radin*float(n+m-3)*epsi(m,n)
      psdyp(mm,n) = radin*float(n+m  )*epsi(m,n+1)
      uvdxm(mm,n) = erad*float(m-1  )/float((n+m-2)*(1-n-m))
      uvdym(mm,n) = erad* epsi(m,n  )/float (2-n-m)
      uvdyp(mm,n) = erad* epsi(m,n+1)/float(1-n-m)
  8   continue
c  em=lower case m index, emn=lower case l index
      do 33 n=1,jx
      do 33 m=1,mxp
      em (m,n)=m-1
      emn(m,n)=m+n-2
33    continue
c
c  start latitude loop to load polynomials into pol
      do 44 ls=1,iy
c
c  y = sqrt(1.0-coa**2) = cosine of the latitude
      y = sia(ls)
      x = coa(ls)
      alp(1,1) = sqrhlf
c
c  use recursive formula for n=1 diagonal of rhomboid of polynomials.
      do 4 m=2,mxp
      alp(m,1) = consq(m)*y*alp(m-1,1)
 4    continue
c
c  derive n=2 polynomials from n=1.
      do 55 m=1,mxp
      alp(m,2)=(x*alp(m,1))*repsi(m,2)
 55   continue
c
c  use recursive formula to obtain n=3 to n=jx polynomials.
      do 5 n=3,jx
      do 5 m=1,mxp
      alp(m,n)=(x*alp(m,n-1)-epsi(m,n-1)*alp(m,n-2))*repsi(m,n)
  5   continue
c
c  zero polynomials with absolute values smaller than 10**(-30)
c  and store these in field pol in common block / polyn /.
      do 42 n=1,jx
      do 42 m=1,mxp
      if(abs(alp(m,n)) .ge. 1.0e-30) then
        paly(m,n)=alp(m,n)
      else
        paly(m,n)=0.0
      endif
 42   continue
      do 6 n=1,jx
      do 6 m=1,mx
      pol(m,n,ls)=paly(isc*(m-1)+1,n)
 6    continue
 44   continue
      return
c
      entry lgndre
c     ===== ======
c
c  lgndre loads poly from pol, then the polynomials are cloned
c   and stored in common block /clnpol/.
c
      do 88 n=1,jx
      do 88 m=1,mx
      poly(m,n)=pol(m,n,ls)
 88   continue
c
c  clone polynomials and store into /clnpol/.
      do 89 n=1,jx
       do 89 m=1,mx
        colyl(2*m-1,n) = poly(m,n)
        colyl(2*m  ,n) = poly(m,n)
   89 continue
c
      return
c
      entry lgndxy
c     ===== ======
c
c  lgndxy computes the polydx and polydy, then polydx and dy are
c   cloned and stored in common block /clnpol/.
c
      wght=radin*weight(ls)
c
c  differentiate polynomials.
      do 77 m=1,mx
      polydx(m, 1)= em(isc*(m-1)+1,1)*poly(m,1)*wght
      polydy(m, 1)=-em(isc*(m-1)+1,1)*poly(m,2)*wght*epsi(isc*(m-1)+1,2)
      polydx(m,jx)=0.0
      polydy(m,jx)=0.0
77    continue
      do 7 n=2,jtrun1
      do 7 m=1,mx
      polydx(m,n)=em(isc*(m-1)+1,n)*poly(m,n)*wght
      polydy(m,n)=
     *   ((emn(isc*(m-1)+1,n)+1.)*epsi(isc*(m-1)+1,n  )*poly(m,n-1)
     *   - emn(isc*(m-1)+1,n)    *epsi(isc*(m-1)+1,n+1)*poly(m,n+1))
     *   *wght
 7    continue
      do 90 n=1,jx
       do 90 m=1,mx
        cypyx(2*m-1,n) = polydx(m,n)
        cypyx(2*m  ,n) = polydx(m,n)
        cypyy(2*m-1,n) = polydy(m,n)
        cypyy(2*m  ,n) = polydy(m,n)
   90 continue
c
  425 format(' The top level of the boundary layer for the gravity'
     *       ' wave computation (kll) is: ',i5,/)
  787 format(' hyk (critical humidity for clouds) =',/,(5x,10f8.4))
 9001 format(' q (full levels) =',/,(5x,10f8.4))
 9002 format(' qmh (half levels) =',/,(5x,10f8.4))
 9003 format(' elmix (mixing length [cm] for vertical diffusion) =',/,
     $      (5x,10f8.1))
 9005 format(' zhstd (height [m] of half-levels for elmix) =',/,
     $      (5x,-2p10f8.1))
 9006 format(//,' *** No dmptop specified for sponge option ',
     $          'in parmtr.F ***',//)
 9007 format(' The number of levels to be diffused (ks) is: ',i3,/)
 9016 format(' htop for mixing length is: ',f14.1)
c
      return
      end
