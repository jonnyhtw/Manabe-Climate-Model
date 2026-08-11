      subroutine slabb(ldum)
c
c  slabb controls the conversion of spectral variables to grid space,
c   calls subroutines which perform moist and radiative processes on
c   gridpoint data, and controls the conversion of grid varaibles to
c   spectral space.  in addition, the second set of auxiliary variables
c   is written to /pascom/ which are then written to unit 31.  slabb
c   is called from passb for each ls latitude row.  ls runs from the
c   north pole to the equator.
c
c slabb was split out of passb inorder to do multitasking.  ldum is
c  a dummy arguement and should not be used.  use ls in /ctlblk/
c  instead.
c
c  subroutines called:  lgndre(parmtr), hexpft, getime
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
c
      complex vor   , div  , t    , psln, ucos , vcos
      complex vorm  , divm , tm   , psm
      complex vordt , divdt, tdt  , psdt
      complex qmix , qmixm, qmxdt
      complex ps    , psdx , psdy , zs, ek
      complex  aift ,  bift, ekift, zaift, zbift
      complex saift , sbift
      complex daift , dbift
      complex  tift , tuift, tvift, psift, spift
      complex spsift
      complex dpsift
      complex  qeift, quift, qvift
#ifdef dirnal
      complex  qdift, ogift
#endif
#ifdef rad
      complex rdift, srdift, drdift, radten
#endif
      complex ftvor1, ftvor2
c
c
#ifdef rad
      dimension  rdift(mx,kx,nhem)
#endif
      dimension anala(awrds),analb(bwrds)
#ifdef dirnal
      dimension anald(bwrds)
#endif
c
#include "auxatm.h"
c
#include "clnpol.h"
c
#include "ctlblk.h"
c
c
c  *** warning  ****
c  the order of the variables in the common block / fftstr /
c   is linked to the variable order in the common block
c   / grid /.  change the order at your own risk.  adding
c   variables to the end of the common blocks
c   should be safe.
c
      common / fftstr / ftvor1(mx,kx,nhem), ftvor2(mx,nhem,1)
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
      common / pascom  / anlwrd(awrds,il)
c
      common / pascom2 /   tdty(mx2,jx,kx,iy), rdteny(mx2,jx,kx,iy)
      common / pascom2 /    psy(mx2,jx, 1,iy)
      common / pascom2 / qmxdty(mx2,jx,kx,iy)
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
      dimension srdift(mx,kx), drdift(mx,kx)
#endif
      dimension  saift(mx,kx),  sbift(mx,kx)
      dimension  daift(mx,kx),  dbift(mx,kx)
      dimension spsift(mx)   , dpsift(mx)
c
#include "water.h"
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
#endif
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
c  begin the latitude loop
c   note: ls is passed thro /ctlblk/
c
c  zero out spectral tendencies which will be summed in the
c   gaussian quadrature (rhxinv).
#ifndef unitsk
      do 150 k=1,kx
       do 150 n=1,jx
        do 150 m=1,mx2
         qmxdty(m,n,k,ls) = 0.0
           tdty(m,n,k,ls) = 0.0
 150  continue
      do 152 n=1,jx
       do 152 m=1,mx2
          psy(m,n,1,ls) = 0.0
 152  continue
      if(lrad) then
        do 151 k=1,kx
         do 151 n=1,jx
          do 151 m=1,mx2
           rdteny(m,n,k,ls) = 0.0
  151   continue
      endif
#endif
c
      call lgndre
c
c  transform variables from spectral space to grid space
cc    call hexpft (    t,   cgt , kx )
cc    call hexpft ( psln,   gps ,  1 )
      call scalar(ftvor1,    t, kx)
      call scalar(ftvor2, psln,  1)
c  ntrans is the number of parrallel transforms
      ntrans = (kx+1)*nhem
c  do ntrans transforms.  note common block order!!
c
      call fftg  (cgt   , ftvor1, ntrans)
c
      call hexpft(qmix, gqmix, kx)
c
c  compute weight for use in gaussian quadrature (rhxinv).
      coswgt = weight(ls)*cosqr(ls)
c
      call getime(cutime,eltime)
      t2ff = t2ff+eltime
c
c  alternate between hemispheres
      do 622 j=1,jlp,il
c
      if (j.eq.1) then
c       northern hemisphere section
        ih   =  1
        jrow = ls
      else
c       southern hemisphere section
        ih   = 2
        jrow = jlp-ls
      endif
c
#ifdef moist
c  initialize accumulated variables when beginning a restart interval
      if (resetp) then
        if (mod(ktau,intrs).eq.0) then
          do 56 i=1,ix
           precip(i,jrow) = 0.0
           snowfl(i,jrow) = 0.0
#ifdef hydrol
           runoff(i,jrow) = 0.0
           snwmlt(i,jrow) = 0.0
#endif
#if defined cpldleg3 || defined cpldleg1
           snwrnf(i,jrow) = 0.0
#endif
 56       continue
        endif
      endif
#endif
c
      zmpc (jrow)=0.0
      zmrn (jrow)=0.0
      zmsn (jrow)=0.0
      zmcah(jrow)=0.0
c
c  load pkice and the land-sea mask for the first latitude loop
      do 180 i=1,ix
       if (ts(i,jrow) .lt. 0.0) then
         zin(i) = 0.0
       else
         zin(i) = 1.0
       endif
  180 continue
c
c   convert liq eq ice to ice
c     pkice(1;ix) = yce(1,jrow;ix) * 1.06
      do 181 i=1,ix
       pkice(i) = yce(i,jrow) * 1.06
  181 continue
c
#ifdef moist
c
c  perform convective adjustment
c
      call convad
c
      call mradj (gqmix(1,1,ih), gvbrps(1,ih), zmadj)
#endif
c
#ifdef dry
      call passup
#endif
c
#ifdef rad
c  call the radiation package on radiation time steps
c  if icf = 1 => calculate radiation with clouds
c  if icf = 2 => calculate radiation without clouds
c
      if(lrad) then
        icf = 1
        call impt (icf)
# ifdef cloud_forcing
        icf = 2
        call impt (icf)
# endif
      endif
#endif
c
#ifdef dirnal
      if(.not.sgmt) go to 43
      read (bdski) anald
      anald=anald+analb
      write (bdsko) anald
#endif
      if (.not. diurnl) then
        if(rsrt) then
c
c  at each restart interval, compute pfac for division of
c   accumulated s.h. flux by number of timesteps
          pfac=1.0/float(intrs)
          if(jrow.eq.1) write(6,81) intrs
c
c   note that instant. (not time averaged) altnt is stored in analb
          do 80 i=1,ix
           zshfx (i) =  shfx(i,jrow)*pfac
           zaltnt(i) = altnt(i,jrow)
 80       continue
c
c  write auxiliary data to unit 31, for multitasking store into
c   jrow storage and write the data after the jrow do loop.
#ifdef unitsk
          write(bdsk) analb
#endif
#ifndef unitsk
          do 345 n=1,bwrds
           anlwrd(n,jrow) = analb(n)
  345     continue
#endif
        endif
      endif
 43   continue
c
 622  continue
      call getime(cutime,eltime)
      tpas=tpas+eltime
c
c  transform varaibles from grid to spectral domain
c     call fftm (    cgt,  aift, kx*nhem )
c     call fftm ( gvbrps, psift,    nhem )
c
c  note that in the code below:
c        cgt =>  aift
c     gvbrps => spift
c     gqmix  =>  bift
c     grad   => rdift  which is equivalenced into tift
c
c  the loading of gps is done to get longer vectors thro the fft
c
      do 10 i=1,ix
       gps(i,1) = gvbrps(i,1)
       gps(i,2) = gvbrps(i,2)
   10 continue
c  do (kx+1)*nhem transforms.  nore common block order!!
c
      call fftm (    cgt,  aift, (kx+1)*nhem )
c
      call fftm(gqmix, bift, kx*nhem)
c
#ifdef rad
      if (lrad) then
        call fftm ( grad, rdift, kx*nhem )
      endif
#endif
      if (nhem.eq.2) then
        do 91 k=1,kx
         do 91 m=1,mx
         saift (m,k) =  aift (m,k,1) +  aift (m,k,2)
         daift (m,k) =  aift (m,k,1) -  aift (m,k,2)
         sbift (m,k) =  bift (m,k,1) +  bift (m,k,2)
         dbift (m,k) =  bift (m,k,1) -  bift (m,k,2)
 91     continue
        do 92 m=1,mx
         spsift(m  ) = spift (m  ,1) + spift (m  ,2)
         dpsift(m  ) = spift (m  ,1) - spift (m  ,2)
 92     continue
      else
c
c  compute one hemisphere sum and difference fourier coefficients
        do 94 k=1,kx
         do 94 m=1,mx
         saift (m,k) =  aift (m,k,1)
         daift (m,k) =  (0.0,0.0)
         sbift (m,k) =  bift(m,k,1)
         dbift (m,k) =  (0.0, 0.0)
 94     continue
        do 95 m=1,mx
         spsift(m  ) = spift (m  ,1)
         dpsift(m  ) = (0.0,0.0)
 95     continue
      endif
#ifdef rad
      if (lrad) then
        if (nhem.eq.2) then
          do 97 k=1,kx
           do 97 m=1,mx
           srdift (m,k) = rdift (m,k,1) + rdift (m,k,2)
           drdift (m,k) = rdift (m,k,1) - rdift (m,k,2)
 97       continue
        else
          do 98 k=1,kx
           do 98 m=1,mx
           srdift (m,k) = rdift (m,k,1)
           drdift (m,k) = (0.0,0.0)
 98       continue
        endif
#ifndef unitsk
        call rhxinv ( rdteny(1,1,1,ls), srdift, drdift, colyl, coswgt,
     *            kx, rdteny(1,1,1,ls))
#endif
#ifdef unitsk
        call rhxinv ( radten(1,1,1   ), srdift, drdift, colyl, coswgt,
     *            kx, radten(1,1,1   ))
#endif
      endif
#endif
#ifndef unitsk
      call rhxinv ( tdty(1,1,1,ls),  saift,  daift, colyl, coswgt,
     *          kx, tdty(1,1,1,ls))
      call rhxinv ( psy (1,1,1,ls), spsift, dpsift, colyl, coswgt,
     *          1,  psy (1,1,1,ls))
      call rhxinv(qmxdty(1,1,1,ls), sbift, dbift, colyl, coswgt,
     *        kx, qmxdty(1,1,1,ls) )
#endif
#ifdef unitsk
      call rhxinv ( tdt (1,1,1   ),  saift,  daift, colyl, coswgt,
     *          kx, tdt (1,1,1   ))
      call rhxinv ( ps  (1,1,1   ), spsift, dpsift, colyl, coswgt,
     *          1,  ps  (1,1,1   ))
      call rhxinv(qmxdt (1,1,1   ), sbift, dbift, colyl, coswgt,
     *        kx, qmxdt (1,1,1   ) )
#endif
c  in the calls to rhxinv, the spectral harmonics are used twice
c   in order to 'get-a-round' a compiler error on the ymp which
c   slowed down the code in rhxinv.
c
c
   81 format(1x,'shfx for tape31 averaged over ',i5,' timesteps')
c
      return
      end
