      subroutine passa
c
c  passa controls the conversion of spectral variables to grid space,
c   then calls the sequence of subroutines which evaluate non-linear
c   terms, compute vertical diffusion, and evaluate the surface heat
c   balance.  passa also controls the conversion of grid vaiables to
c   spectral space.  in addition, the first set of auxiliary varaibles
c   is written to unit adsk.  passa represents the first latitude loop
c   and is called from atmos each timestep.
c
c  subroutines called:  slaba
c                       getime
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
      complex vorm  , divm , tm   , psm , ucosm, vcosm
      complex vordt , divdt, tdt  , psdt
      complex qmix , qmixm , qmxdt
      complex ps    , psdx , psdy , zs   , ek
      complex  aift ,  bift, ekift, zaift, zbift
      complex  tift , tuift, tvift, psift, spift
      complex  qeift, quift, qvift
#ifdef dirnal
      complex qdot, somg, qdift, ogift, sqgrd, dqgrd, somgr, domgr
#endif
      complex rzs
c
      dimension anala(awrds),analb(bwrds)
#ifdef dirnal
      dimension analc(awrds)
#endif
      dimension rvordt(mx2,jx,kx), rdivdt(mx2,jx,kx)
      dimension   rtdt(mx2,jx,kx),  rpsdt(mx2,jx, 1)
      dimension   rrps(mx2,jx, 1),    rek(mx2,jx,kx)
      dimension rqmxdt(mx2,jx,kx)
#ifdef dirnal
      dimension  rqdot(mx2,jx,kx),  rsomg(mx2,jx,kx)
#endif
c
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
#include "guvm.h"
c
#include "ift.h"
c
#include "intgrl.h"
c
#include "iohist.h"
c
#include "mrstr.h"
c
#include "nxstps.h"
c
#include "param.h"
c
cifndef unitsk
c     common / pascom  / anlwrd(awrds,il)
c
c     common / pascom2 / vordty(mx2,jx,kx,iy)
c     common / pascom2 / divdty(mx2,jx,kx,iy)
c     common / pascom2 /    eky(mx2,jx,kx,iy)
c     common / pascom2 /   tdty(mx2,jx,kx,iy)
c     common / pascom2 / qmxdty(mx2,jx,kx,iy)
c     common / pascom2 /    psy(mx2,jx, 1,iy)
c     common / pascom2 /  psdty(mx2,jx, 1,iy)
cifdef dirnal
c     common / pascom2 /  qdoty(mx2,jx,kx,iy)
c     common / pascom2 /  somgy(mx2,jx,kx,iy)
cendif
cendif
#ifndef unitsk
# include "psacom.h"
c
# include "psacom2.h"
#endif
c
#include "polyls.h"
c
#include "pstar.h"
c
#include "qa.h"
c
#include "sfharm.h"
c
#include "sftuvm.h"
c
#include "sqdot.h"
c
#include "surfc.h"
c
#include "surfo.h"
c
#include "timeck.h"
c
#include "water.h"
c
#include "zrise.h"
c
      equivalence ( anala(1),  vfu(1,1))
#ifdef moist
      equivalence ( analb(1),convr(1,1))
#endif
#ifdef dry
      equivalence ( analb(1),convt(1,1))
#endif
      equivalence ( rvordt(1,1,1),  vordt(1,1,1))
      equivalence ( rdivdt(1,1,1),  divdt(1,1,1))
      equivalence (   rtdt(1,1,1),    tdt(1,1,1))
      equivalence (  rpsdt(1,1,1),   psdt(1,1,1))
      equivalence (   rrps(1,1,1),     ps(1,1,1))
      equivalence (    rek(1,1,1),     ek(1,1,1))
      equivalence ( rqmxdt(1,1,1),  qmxdt(1,1,1))
#ifdef dirnal
      equivalence (  rqdot(1,1,1),   qdot(1,1,1))
      equivalence (  rsomg(1,1,1),   somg(1,1,1))
#endif
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c  zero complex spectral tendencies to prepare for gaussian quadrature
c   note: rxxxx is equivalenced to xxxx above.
      do 50 k=1,kx
       do 50 n=1,jx
        do 50 m=1,mx2
         rvordt(m,n,k) = 0.0
         rdivdt(m,n,k) = 0.0
          rek  (m,n,k) = 0.0
          rtdt (m,n,k) = 0.0
         rqmxdt(m,n,k) = 0.0
 50   continue
      do 52 n=1,jx
       do 52 m=1,mx2
         rrps(m,n,1) = 0.0
        rpsdt(m,n,1) = 0.0
 52   continue
#ifdef dirnal
      if(sgmt) then
        do 51 k=1,kx
         do 51 n=1,jx
          do 51 m=1,mx2
           rqdot(m,n,k) = 0.0
           rsomg(m,n,k) = 0.0
   51   continue
      endif
#endif
c
c  initialize global integrals to zero
      drvor = 0.0
      cnvr  = 0.0
      tvdis = 0.0
      if(resetn) then
        tevap = 0.0
        gvft  = 0.0
        gvfr  = 0.0
      endif
      zmts  = 0.0
      thtfx = 0.0
#ifdef dirnal
      if(sgmt) then
        nou=iabs(nin-2)
        adski=adsk+nin
        adsko=adsk+nou
      endif
#endif
c
c  begin the latitude loop
      do 46 ls=1,iy
c
      ldum = ls
      call slaba(ldum)
c
 46   continue
c
c  compute latitude weights for use in computing global integrals
c
      do 48 ls=1,iy
       wght   = weight(ls)
       coswgt = weight(ls)*cosqr(ls)
       rainwt = 0.5*hmf*coswgt
       do 48 ih=1,nhem
        if (ih .eq. 1) then
c         northern hemisphere section
          jrow = ls
        else
c         southern hemisphere section
          jrow = jlp-ls
        endif
c
c  compute global integrals for printing in atmos and energy
        drvor = drvor + coswgt*zmdrag(jrow)
        cnvr  =  cnvr + coswgt*zmcnvr(jrow)
        tvdis = tvdis + coswgt*zmvdis(jrow)
        zmts  = zmts  + coswgt*zmt   (jrow)
#ifdef moist
        tevap = tevap + rainwt*zmevap(jrow)
#endif
        thtfx = thtfx + rainwt*zmhtfx(jrow)
#ifdef moist
        gvfr  = gvfr  + rainwt*zmvfr (jrow)*og
#endif
        gvft  = gvft  + rainwt*zmvft (jrow)*cpp*dt*og
   48 continue
c
#ifndef unitsk
c  on restart intervals write data onto unit adsk
      if (rsrt) then
        do 347 j=1,il
         do 348 n=1,awrds
          anala(n) = anlwrd(n,j)
  348    continue
         write(adsk) anala
  347 continue
      endif
c
c  sum complex spectral tendencies which were stored in slaba
      do 150 ls=1,iy
       do 150 k=1,kx
        do 150 n=1,jx
         do 150 m=1,mx2
          rvordt(m,n,k) = vordty(m,n,k,ls) + rvordt(m,n,k)
          rdivdt(m,n,k) = divdty(m,n,k,ls) + rdivdt(m,n,k)
           rek  (m,n,k) =    eky(m,n,k,ls) +    rek(m,n,k)
          rqmxdt(m,n,k) = qmxdty(m,n,k,ls) + rqmxdt(m,n,k)
           rtdt (m,n,k) =   tdty(m,n,k,ls) +   rtdt(m,n,k)
 150  continue
      do 152 ls=1,iy
       do 152 n=1,jx
        do 152 m=1,mx2
          rrps(m,n,1) =   psy(m,n,1,ls) +   rrps(m,n,1)
         rpsdt(m,n,1) = psdty(m,n,1,ls) +  rpsdt(m,n,1)
 152  continue
# ifndef use_mpp
cmic$ parallel SHARED(rvordt,vordty,rdivdt,divdty,rek,eky,rqmxdt,
cmic$*                 qmxdty,rtdt,tdty,rrps,psy,rpsdt,psdty)
cmic$ case
      call sumk(rvordt, vordty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(rdivdt, divdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(   rek,    eky, mx2, jx, kx, ncpus)
#   ifndef grdmx
cmic$ case
      call sumk(rqmxdt, qmxdty, mx2, jx, kx, ncpus)
#   endif
cmic$ case
      call sumk(  rtdt,   tdty, mx2, jx, kx, ncpus)
cmic$ case
      call sumk(  rrps,    psy, mx2, jx,  1, ncpus)
cmic$ case
      call sumk( rpsdt,  psdty, mx2, jx,  1, ncpus)
cmic$ end case
cmic$ end parallel
# endif
# ifdef use_mpp
      call sumk(rvordt, vordty, mx2, jx, kx, ncpus)
      call sumk(rdivdt, divdty, mx2, jx, kx, ncpus)
      call sumk(   rek,    eky, mx2, jx, kx, ncpus)
      call sumk(rqmxdt, qmxdty, mx2, jx, kx, ncpus)
      call sumk(  rtdt,   tdty, mx2, jx, kx, ncpus)
      call sumk(  rrps,    psy, mx2, jx,  1, ncpus)
      call sumk( rpsdt,  psdty, mx2, jx,  1, ncpus)
# endif

c
#endif

#ifdef dirnal
      if(sgmt) then
        call sumk( rqdot,  qdoty, mx2, jx, kx, ncpus)
        call sumk( rsomg,  somgy, mx2, jx, kx, ncpus)
      endif
#endif
cifdef dirnal
c     if(sgmt) then
c       do 151 ls=1,iy
c        do 151 k=1,kx
c         do 151 n=1,jx
c          do 151 m=1,mx2
c           rqdot(m,n,k) = qdoty(m,n,k,ls) + rqdot(m,n,k)
c           rsomg(m,n,k) = somgy(m,n,k,ls) + rsomg(m,n,k)
c 151   continue
c     endif
cendif
c
#ifdef triang
c
c  Zero out upper part of the rhombus.
c
      call trunct (   ek,kx)
      call trunct (vordt,kx)
      call trunct (divdt,kx)
      call trunct (  tdt,kx)
      call trunct ( psdt, 1)
      call trunct (   ps, 1)
      call trunct (qmxdt,kx)
# ifdef dirnal
      if (sgmt) then
        call trunct (qdot,kx)
        call trunct (somg,kx)
      endif
# endif
C
#endif
      call getime(cutime,eltime)
      thex = thex+eltime
c
 4322 format(' iy is not even divisible by ncpus, stop ')
c
c
      return
      end
