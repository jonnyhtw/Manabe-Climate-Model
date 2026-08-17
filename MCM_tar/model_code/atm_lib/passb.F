      subroutine passb
c
c  passb controls the conversion of spectral variables to grid space,
c   calls subroutines which perform moist and radiative processes on
c   gridpoint data, and controls the conversion of grid varaibles to
c   spectral space.  in addition, the second set of auxiliary variables
c   is written to unit bdsk.  passb is called from atmos each timestep
c   and constitutes the second latitude loop.
c
c  subroutines called:  lgndre(parmtr), getime
c                       convad, mradj, impt, fftm(fftkkg), rhxinv
c                       fftg, scalar
c
c
#include "parm.h"
#include "parmcpu.h"
#include "parmjx.h"
#include "parmkp.h"
#include "parmout.h"
#include "parmdsk.h"
#include "clnpol.h"
c
      complex vor   , div  , t    , psln, ucos , vcos
      complex vorm  , divm , tm   , psm
      complex vordt , divdt, tdt  , psdt
      complex qmix , qmixm, qmxdt
      complex ps    , psdx , psdy , zs, ek
      complex  aift ,  bift, ekift, zaift, zbift
      complex  tift , tuift, tvift, psift, spift
      complex  qeift, quift, qvift
#ifdef dirnal
      complex  qdift, ogift
#endif
#ifdef rad
      complex rdift, radten
#endif
c
#ifdef rad
      dimension  rdift(mx,kx,nhem)
#endif
      dimension anala(awrds),analb(bwrds)
c     dimension anald(bwrds)
      dimension   rtdt(mx2,jx,kx), rrdten(mx2,jx,kx)
      dimension   rrps(mx2,jx, 1)
      dimension rqmxdt(mx2,jx,kx)
c
#include "anal2.h"
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
#include "gridw.h"
c
#include "ift.h"
c
#include "intgrl.h"
c
#include "iohist.h"
c
#include "nxstps.h"
c
#include "param.h"
c
#ifndef unitsk
c     common / pascom  / anlwrd(bwrds,il)
c
c     common / pascom2 /   tdty(mx2,jx,kx,iy), rdteny(mx2,jx,kx,iy)
c     common / pascom2 /    psy(mx2,jx, 1,iy)
c     common / pascom2 / qmxdty(mx2,jx,kx,iy)
#include "psbcom.h"
c
#include "psbcom2.h"
c
#endif
c
#include "polyls.h"
c
#include "pstar.h"
c
#include "radten.h"
c
#include "sfharm.h"
c
#include "surfc.h"
c
#include "surfo.h"
c
#include "timeck.h"
c
#ifdef rad
      dimension hswtop(2), hswbot(2), hlwtop(2)
#endif
      dimension hlwbot(2), hlatht(2), hsenht(2)
      dimension hqao  (2)
#ifdef prdcld
      dimension hcloud(2)
#endif
c
#include "water.h"
c
#include "zmflux.h"
c
      equivalence ( anala  (1), vfu   (1,  1) )
#ifdef moist
      equivalence ( analb  (1), convr (1,  1) )
#endif
#ifdef dry
      equivalence ( analb  (1), convt (1,  1) )
#endif
#ifdef rad
      equivalence ( rdift(1,1,1), tift(1,1,1) )
      equivalence ( rrdten(1,1,1), radten(1,1,1) )
#endif
      equivalence (   rrps(1,1,1),     ps(1,1,1) )
      equivalence (   rtdt(1,1,1),    tdt(1,1,1) )
      equivalence ( rqmxdt(1,1,1),  qmxdt(1,1,1) )
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
#ifdef moist
      gmadj = 0.0
#endif
      if (resetn) then
        tprec = 0.0
        train = 0.0
        tsnow = 0.0
        caht  = 0.0
      endif
c
c  zero qmxdt and tdt for temporary storage of convectively-adjusted
c   spectral temperature and mixing ratio which are accumulated by
c   gaussian quadrature throughout the latitude loop.
c   rqmxdt and rtdt are equivalenced to qmxdt and tdt.
      do 60 k=1,kx
       do 60 n=1,jx
        do 60 m=1,mx2
         rqmxdt(m,n,k) = 0.0
           rtdt(m,n,k) = 0.0
 60   continue
c
c  zero ps for surface pressure (not log p!!) computation
c   rrps is equivalenced to ps.
      do 67 n=1,jx
       do 67 m=1,mx2
        rrps(m,n,1) = 0.0
 67   continue
c
#ifdef rad
c  zero radten for sumation in gaussian quadrature
c   (only on radiation time step).
c   rrdten is equivalenced to radten.
      if (lrad) then
        do 55 k=1,kx
         do 55 n=1,jx
          do 55 m=1,mx2
           rrdten(m,n,k) = 0.0
 55     continue
      endif
#endif
#ifdef dirnal
      if (sgmt) then
        nou  =iabs(nin-2)
        bdski=bdsk+nin
        bdsko=bdsk+nou
      endif
#endif
c
c  begin the latitude loop
      do 62 ls=1,iy
c
      ldum = ls
      call slabb(ldum)
c
 62   continue
c  end of the latitude loop
c
c  if lrad is true, zero radiation integrals
      if (lrad) then
          do 435 ih=1,2
           hswtop(ih) = 0.0
           hswbot(ih) = 0.0
           hlwtop(ih) = 0.0
           hlwbot(ih) = 0.0
           hlatht(ih) = 0.0
           hsenht(ih) = 0.0
           hqao  (ih) = 0.0
#ifdef prdcld
           hcloud(ih) = 0.0
#endif
  435     continue
      endif
c  do loop to compute global integrals
c
      do 68 ls=1,iy
       wght   = weight(ls)
       coswgt = weight(ls)*cosqr(ls)
       rainwt = 0.5*hmf*coswgt
       do 68 ih=1,nhem
        if (ih .eq. 1) then
c         northern hemisphere
          jrow = ls
        else
c         southern hemisphere
          jrow = jlp - ls
        endif
c
c  compute global integrals for printing in atmos and energy
#ifdef moist
        tprec = tprec + rainwt*zmpc (jrow)
        train = train + rainwt*zmrn (jrow)
        tsnow = tsnow + rainwt*zmsn (jrow)
        gmadj = gmadj + rainwt*zmadj(jrow)*og *rix
#endif
        caht  = caht  + rainwt*zmcah(jrow)*cpp*o2g
c
c  on radiation timesteps, compute and print radiation fluxes
c
#ifdef rad
        if (lrad) then
c  compute hemispheric radiation means
          if (ih .eq. 1) then
            hswtop(1) = hswtop(1) + zmswtp(jrow)*coswgt
            hswbot(1) = hswbot(1) + zmswbm(jrow)*coswgt
            hlwtop(1) = hlwtop(1) + zmlwtp(jrow)*coswgt
            hlwbot(1) = hlwbot(1) + zmlwbm(jrow)*coswgt
            hlatht(1) = hlatht(1) + zmltht(jrow)*coswgt
            hsenht(1) = hsenht(1) + zmsnht(jrow)*coswgt
            hqao  (1) = hqao  (1) + zmqao (jrow)*coswgt
#ifdef prdcld
            hcloud(1) = hcloud(1) + zmtcld(jrow)*coswgt
#endif
          else
            hswtop(2) = hswtop(2) + zmswtp(jrow)*coswgt
            hswbot(2) = hswbot(2) + zmswbm(jrow)*coswgt
            hlwtop(2) = hlwtop(2) + zmlwtp(jrow)*coswgt
            hlwbot(2) = hlwbot(2) + zmlwbm(jrow)*coswgt
            hlatht(2) = hlatht(2) + zmltht(jrow)*coswgt
            hsenht(2) = hsenht(2) + zmsnht(jrow)*coswgt
            hqao  (2) = hqao  (2) + zmqao (jrow)*coswgt
#ifdef prdcld
            hcloud(2) = hcloud(2) + zmtcld(jrow)*coswgt
#endif
          endif
        endif
#endif
   68 continue
c
#ifdef rad
c  compute global mean radiation fluxes
c
      if (lrad) then
        tswtop = 0.5 * hmf * (hswtop(1) + hswtop(2))
        tlwtop = 0.5 * hmf * (hlwtop(1) + hlwtop(2))
        tswbot = 0.5 * hmf * (hswbot(1) + hswbot(2))
        tlwbot = 0.5 * hmf * (hlwbot(1) + hlwbot(2))
        tsenht = 0.5 * hmf * (hsenht(1) + hsenht(2))
        tlatht = 0.5 * hmf * (hlatht(1) + hlatht(2))
        tqao   = 0.5 * hmf * (hqao  (1) + hqao  (2))
#ifdef prdcld
        tcloud = 0.5 * hmf * (hcloud(1) + hcloud(2))
#endif
c
        samsw = samsw + tswtop
        samlw = samlw + tlwtop
        kanal = kanal + 1
c
        write(6,900) day
#ifdef dirnal
        if (kdisk) then
#endif
        write(6,1100) hswtop,tswtop,hlwtop,tlwtop
     *               ,hswbot,tswbot,hlwbot,tlwbot
     *               ,hsenht,tsenht,hlatht,tlatht
     *               ,hqao  ,tqao
#ifdef prdcld
     *               ,hcloud,tcloud
#endif
#ifdef dirnal
        endif
#endif
#ifdef prdcld
        write(6,102) zmtcld
#endif
      endif
#endif
c
c  sum complex spectral tendencies which were stored in slabb
#ifndef unitsk
# ifndef use_mpp
      if (lrad) then
cmic$ parallel SHARED(rqmxdt,qmxdty,rtdt,tdty,rrps,psy,rrdten,rdteny)
cmic$ case
      call sumk(rqmxdt, qmxdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(  rtdt,   tdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(  rrps,    psy, mx2, jx,  1, ncpus)
cmic$ case
      call sumk(rrdten, rdteny, mx2, jx, kx, ncpus)
cmic$ end case
cmic$ end parallel
      else
cmic$ parallel SHARED(rqmxdt,qmxdty,rtdt,tdty,rrps,psy)
cmic$ case
      call sumk(rqmxdt, qmxdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(  rtdt,   tdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(  rrps,    psy, mx2, jx,  1, ncpus)
cmic$ end case
cmic$ end parallel
      endif
# endif
# ifdef use_mpp
      if (lrad) then
        call sumk(rrdten, rdteny, mx2, jx, kx, ncpus)
      endif
      call sumk(rqmxdt, qmxdty, mx2, jx, kx, ncpus)
      call sumk(  rtdt,   tdty, mx2, jx, kx, ncpus)
      call sumk(  rrps,    psy, mx2, jx,  1, ncpus)
# endif

c     do 150 ls=1,iy
c      do 150 k=1,kx
c       do 150 n=1,jx
c        do 150 m=1,mx2
c         rqmxdt(m,n,k) = qmxdty(m,n,k,ls) + rqmxdt(m,n,k)
c           rtdt(m,n,k) =   tdty(m,n,k,ls) +   rtdt(m,n,k)
c150  continue
c     if (lrad) then
c       do 151 ls=1,iy
c        do 151 k=1,kx
c         do 151 n=1,jx
c          do 151 m=1,mx2
c           rrdten(m,n,k) = rdteny(m,n,k,ls) + rrdten(m,n,k)
c151    continue
c     endif
c     do 152 ls=1,iy
c      do 152 n=1,jx
c       do 152 m=1,mx2
c         rrps(m,n,1) = psy(m,n,1,ls) + rrps(m,n,1)
c152  continue
#endif
#ifdef triang
c
c  Zero out upper part of the rhombus.
c
      call trunct (tdt,kx)
      call trunct ( ps, 1)
      call trunct (qmxdt,kx)
#ifdef rad
      if (lrad) then
        call trunct(radten,kx)
      endif
#endif
c
#endif
c
c  move convectively-adjusted spectral temperature and mixing ratio
c   to their original locations
      do 71 k=1,kx
       do 71 n=1,jx
        do 71 m=1,mx
         t   (m,n,k) =   tdt(m,n,k)
         qmix(m,n,k) = qmxdt(m,n,k)
 71   continue
c
#ifndef unitsk
c  on restart intervals write data onto unit bdsk
      if (rsrt) then
        do 347 j=1,il
         do 348 n=1,bwrds
          analb(n) = anlwrd(n,j)
  348    continue
         write(bdsk) analb
  347 continue
      endif
#endif
c
   81 format(1x,'shfx for bdsk averaged over ',i5,' timesteps')
#ifdef prdcld
  102 format(/10(2x,f8.4))
#endif
  900 format(' *** radiation timestep at day = ',f11.3)
 1100 format(1x,2(12x,11hnorth. hem.,5x,11hsouth. hem.,7x,6hglobal,13x)
     *                  /8h sntop =,3e16.7,10x,7hlwtop =,3e16.7,
     *                  /8h sngrd =,3e16.7,10x,7hlwgrd =,3e16.7,
     *                  /8h shftx =,3e16.7,10x,7hlatnt =,3e16.7,
     *                  /8h ocnqf =,3e16.7,10x,7hcloud =,3e16.7)
c
      return
      end
