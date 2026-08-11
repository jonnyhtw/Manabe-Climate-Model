      subroutine slaba(ldum)
c
c  slaba controls the conversion of spectral variables to grid space,
c   then calls the sequence of subroutines which evaluate non-linear
c   terms, compute vertical diffusion, and evaluate the surface heat
c   balance.  passa also controls the conversion of grid vaiables to
c   spectral space.  in addition, the first set of auxiliary varaibles
c   is written to /pascom / which is then written to unit 30.  slaba
c   represents the latitude loop whcich runs from the north pole to the
c   equator.  it is called from passa for each ls.  this routine
c   was split from passa for multitasking.
c
c  ldum is a passed variable which is not used.  it should be the
c   same value as ls in the cltblk.  use ls, not ldum.
c
c  subroutines called:  lgndre(parmtr), lgndxy(parmtr), hexpft,
c                       getime, grmult,   fftm(fftkkg), rhxinv
c                       fftg(fftkkg), scalar
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
      complex vorm  , divm , tm   , psm , ucosm, vcosm
      complex vordt , divdt, tdt  , psdt
      complex qmix , qmixm, qmxdt
      complex ps    , psdx , psdy , zs, ek, iz
      complex  aift ,  bift, ekift, zaift, zbift
      complex saift , sbift, skift,szaift,szbift
      complex daift , dbift, dkift,dzaift,dzbift
      complex  tift , tuift, tvift, psift, spift
      complex stift ,stuift,stvift,spsift,sspift
      complex dtift ,dtuift,dtvift,dpsift,dspift
      complex  qeift, quift, qvift
      complex sqeift,squift,sqvift
      complex dqeift,dquift,dqvift
#ifdef dirnal
      complex qdot, somg, qdift, ogift, sqgrd, dqgrd, somgr, domgr
#endif
      complex rzs
      complex ftvor1, ftvor2,ftvor3,ftvor4,ftvor5,ftvor6,ftvor7
c
c
      dimension anala(awrds),analb(bwrds)
#ifdef dirnal
      dimension analc(awrds)
#endif
c
c
#include "auxatm.h"
c
#include "clnpol.h"
c
#include "ctlblk.h"
c
c  *** warning  ****
c  the order of the variables in the common block / fftstr /
c   is linked to the variable order in the common block
c   / grid /.  change the order at your own risk.  adding
c   variables to the end of the common blocks
c   should be safe.
c
      common / fftstr / ftvor1(mx,kx,nhem), ftvor2(mx,kx,nhem)
      common / fftstr / ftvor3(mx,kx,nhem), ftvor4(mx,kx,nhem)
      common / fftstr / ftvor5(mx,kx,nhem), ftvor6(mx,nhem, 4)
      common / fftstr / ftvor7(mx,kx,nhem)
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
#ifndef unitsk
c     common /pascom / anlwrd(awrds,il)
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

#include "pascom.h"
c
#include "pascom2.h"

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
      dimension szaift(mx,kx), szbift(mx,kx)
      dimension dzaift(mx,kx), dzbift(mx,kx)
      dimension saift (mx,kx), sbift (mx,kx), skift (mx,kx)
      dimension stift (mx,kx), stuift(mx,kx), stvift(mx,kx)
      dimension squift(mx,kx), sqvift(mx,kx), sqeift(mx,kx)
      dimension spsift(mx   ), sspift(mx   )
      dimension daift (mx,kx), dbift (mx,kx), dkift (mx,kx)
      dimension dtift (mx,kx), dtuift(mx,kx), dtvift(mx,kx)
      dimension dquift(mx,kx), dqvift(mx,kx), dqeift(mx,kx)
      dimension dpsift(mx   ), dspift(mx   )
#ifdef dirnal
      dimension sqgrd (mx,kx), somgr (mx,kx)
      dimension dqgrd (mx,kx), domgr (mx,kx)
#endif
c     dimension anal30(awrds)
c
#include "water.h"
c
#include "zrise.h"
c
#include "zmflux.h"
c
      equivalence ( anala(1),  vfu(1,1))
#ifdef moist
      equivalence ( analb(1),convr(1,1))
#endif
#ifdef dry
      equivalence ( analb(1),convt(1,1))
#endif
c
      data iz / (0.0,1.0)  /
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c  begin the latitude loop
c   note: ls is passed through the ctlblk.
c
c  zero out spectral tendencies which will be summed in the
c   gaussian quadrature (rhxinv).
#ifndef unitsk
      do 150 k=1,kx
       do 150 n=1,jx
        do 150 m=1,mx2
         vordty(m,n,k,ls) = 0.0
         divdty(m,n,k,ls) = 0.0
            eky(m,n,k,ls) = 0.0
         qmxdty(m,n,k,ls) = 0.0
           tdty(m,n,k,ls) = 0.0
 150  continue
      do 152 n=1,jx
       do 152 m=1,mx2
          psy(m,n,1,ls) = 0.0
        psdty(m,n,1,ls) = 0.0
 152  continue
# ifdef dirnal
      if(sgmt) then
        do 151 k=1,kx
         do 151 n=1,jx
          do 151 m=1,mx2
           qdoty(m,n,k,ls) = 0.0
           somgy(m,n,k,ls) = 0.0
  151  continue
      endif
# endif
#endif
c
c  compute latitude weights for use in the gaussian quadrature (rhxinv).
      wght   = weight(ls)
      coswgt = weight(ls)*cosqr(ls)
c
      call lgndre
c
      call lgndxy
c
c  transform variables from spectral to grid space
c     call hexpft (   vor,   gvor, kx )
c     call hexpft (   div,   gdiv, kx )
c     call hexpft ( ucos , gucos , kx )
c     call hexpft ( vcos , gvcos , kx )
c     call hexpft ( ucosm, gucosm, kx )
c     call hexpft ( vcosm, gvcosm, kx )
c     call hexpft (    t ,   cgt , kx )
c     call hexpft (    tm,    gtm, kx )
c     if ( .not. hybrid )  then
c          call hexpft (  qmix,  gqmix, kx )
c          call hexpft ( qmixm, gqmixm, kx )
c     end if
c     call hexpft (  psdx,  gpsdx,  4 )
c
c  the calls to scalar and fftg below convert the spectral vbls
c   into grid space.  in these calls order is important!
c   memory in the common block fftstr is aligned with the
c   memory in the common block grid.  change the order
c   at your own risk.
c
      call scalar(ftvor1       ,   vor, kx)
      call scalar(ftvor2       ,   div, kx)
      call scalar(ftvor3       , ucos , kx)
      call scalar(ftvor4       , vcos , kx)
      call scalar(ftvor5       ,   t  , kx)
      call scalar(ftvor6(1,1,1), psln ,  1)
      call scalar(ftvor6(1,1,2), psdx ,  1)
      call scalar(ftvor6(1,1,3), psdy ,  1)
      if (grow) then
        call scalar(ftvor6(1,1,4),rzs   ,  1)
      else
        call scalar(ftvor6(1,1,4), zs   ,  1)
      endif
c  ntrans is the number of transforms performed in parrallel
      ntrans = (5*kx + 4) * nhem
      call scalar(ftvor7       , qmix, kx)
      ntrans = (6*kx+4) * nhem
c  do ntrans transforms at once.  note common block order!!
c
      call fftg ( gvor, ftvor1, ntrans)
c
      call scalar(ftvor1       , ucosm, kx)
      call scalar(ftvor2       , vcosm, kx)
      call scalar(ftvor3       ,    tm, kx)
      ntrans = 3*kx*nhem
      call scalar(ftvor4       , qmixm, kx)
      ntrans = 4*kx*nhem
c  do ntrans transforms at once.  note common block order!!
c
      call fftg ( gucosm, ftvor1, ntrans)
c
      call getime(cutime,eltime)
      thxp   = thxp+eltime
c
c  alternate between hemispheres
      do 465 j=1,jlp,il
c
      if (j.eq.1) then
c       northern hemisphere section
        ih   = 1
        jrow = ls
        abvr = tw*coa(ls)
      else
c       southern hemisphere section
        ih   = 2
        jrow = jlp-ls
        abvr = -tw*coa(ls)
      endif
c
c  zero out zonal mean fluxes which are loaded in tsc
      zmsnht(jrow) = 0.0
      zmltht(jrow) = 0.0
      zmqao (jrow) = 0.0
c
#ifdef moist
c  initialize accumulated variables when beginning a restart interval
      if (resetp) then
        if (mod(ktau,intrs).eq.0) then
          do 56 i=1,ix
          evapor(i,jrow)=0.0
          sublim(i,jrow)=0.0
          shfx  (i,jrow)=0.0
          ycemlt(i,jrow)=0.0
 56      continue
        endif
      endif
#endif
c
      do 42 k=1,kx
       do 42 i=1,ix
        gvor(i,k,ih) = gvor(i,k,ih) + abvr
 42   continue
c
c  load pkice and the land-sea mask for the first latitude loop
      do 181 i=1,ix
       if (ts(i,jrow) .lt. 0.0) then
         zin(i) = 0.0
       else
         zin(i) = 1.0
       endif
  181 continue
c
c
c   convert liq eq ice to solid sea ice
      do 182 i=1,ix
       pkice(i) = yce(i,jrow) * 1.06
  182 continue
c
c
c  evaluate non-linear terms (grmult), compute vertical diffusion
c   (vertdf), and evaluate surface heat balance (tsc)
c
      call grmult
c
      call getime(cutime,eltime)
      tgrm=tgrm+eltime
c
#ifdef dirnal
      do 36 k=1,kx
      do 36 i=1,ix
      gwdots(i,k,ih)=gwdot(i,k+1)
 36   continue
c
      if (sgmt) then
        read  (adski) analc
        analc=analc+anala
        write (adsko) analc
#endif
c
c  write auxiliary data to unit 30, for multitasking store into
c   jrow storage and write the data after the jrow do loop.
c
        if (.not. diurnl) then
#ifndef unitsk
          if (rsrt) then
            do 345 n=1,awrds
             anlwrd(n,jrow) = anala(n)
  345       continue
          endif
#endif
#ifdef unitsk
          if (rsrt) write(adsk) anala
#endif
        endif
#ifdef dirnal
      endif
#endif
c
#ifdef moist
      do 59 i=1,ix
      altnt(i,jrow)=fluxw(i)
 59   continue
#endif
c
 465  continue
c
c
c  transform variables from grid to spectral space
c     call fftm (  agrid,  aift, kx*nhem)
c     call fftm (  bgrid,  bift, kx*nhem)
c     call fftm (   gteq,  tift, kx*nhem)
c     call fftm (   gtu , tuift, kx*nhem)
c     call fftm (   gtv , tvift, kx*nhem)
c     call fftm ( gekcos, ekift, kx*nhem)
c     call fftm (  gpseq, psift,    nhem)
c     call fftm ( gvbrps, spift,    nhem)
c  ntrans is the number of transforms in parrallel
      ntrans = (6*kx + 2) * nhem
c
cq      call fftm ( gqmxeq, qeift, kx*nhem)
cq      call fftm ( gqmixu, quift, kx*nhem)
cq      call fftm ( gqmixv, qvift, kx*nhem)
        ntrans = (9*kx+2) * nhem
c  do ntrans transforms at once.  note common block order!!
c
      call fftm ( agrid,   aift, ntrans)
c
#ifdef dirnal
      if(sgmt) then
        call fftm ( gwdots, qdift, kx*nhem)
        call fftm ( omg   , ogift, kx*nhem)
      endif
#endif
c
c  differentiate aift, bift, etc. (a,b) in zonal direction by
c   multiplying by square root of -1.  m is incorporated in
c   polydx and polydy. (see march 10, 1982 notes on gaussian.)
c   quadrature.)
      do 47 ihem=1,nhem
       do 47 k=1,kx
        do 47 m=1,mx
         zaift(m,k,ihem) = iz * aift (m,k,ihem)
         zbift(m,k,ihem) = iz * bift (m,k,ihem)
         tuift(m,k,ihem) = iz * tuift(m,k,ihem)
         quift(m,k,ihem) = iz * quift(m,k,ihem)
 47   continue
c
c  perform gaussian quadrature and sum and difference fourier coef.
      if (nhem.eq.2) then
c       global model sum and differences of fourier coef.
        do 80 k=1,kx
         do 80 m=1,mx
         szaift  (m,k) =  zaift(m,k,1) +  zaift(m,k,2)
         dzaift  (m,k) =  zaift(m,k,1) -  zaift(m,k,2)
         szbift  (m,k) =  zbift(m,k,1) +  zbift(m,k,2)
         dzbift  (m,k) =  zbift(m,k,1) -  zbift(m,k,2)
         saift   (m,k) =   aift(m,k,1) +   aift(m,k,2)
         daift   (m,k) =   aift(m,k,1) -   aift(m,k,2)
         sbift   (m,k) =   bift(m,k,1) +   bift(m,k,2)
         dbift   (m,k) =   bift(m,k,1) -   bift(m,k,2)
         skift   (m,k) =  ekift(m,k,1) +  ekift(m,k,2)
         dkift   (m,k) =  ekift(m,k,1) -  ekift(m,k,2)
         stift   (m,k) =   tift(m,k,1) +   tift(m,k,2)
         dtift   (m,k) =   tift(m,k,1) -   tift(m,k,2)
         stuift  (m,k) =  tuift(m,k,1) +  tuift(m,k,2)
         dtuift  (m,k) =  tuift(m,k,1) -  tuift(m,k,2)
         stvift  (m,k) =  tvift(m,k,1) +  tvift(m,k,2)
         dtvift  (m,k) =  tvift(m,k,1) -  tvift(m,k,2)
         sqeift  (m,k) =  qeift(m,k,1) +  qeift(m,k,2)
         dqeift  (m,k) =  qeift(m,k,1) -  qeift(m,k,2)
         squift  (m,k) =  quift(m,k,1) +  quift(m,k,2)
         dquift  (m,k) =  quift(m,k,1) -  quift(m,k,2)
         sqvift  (m,k) =  qvift(m,k,1) +  qvift(m,k,2)
         dqvift  (m,k) =  qvift(m,k,1) -  qvift(m,k,2)
 80     continue
#ifdef dirnal
        if (sgmt) then
          do 801 k=1,kx
           do 801 m=1,mx
           sqgrd   (m,k) =  qdift(m,k,1) +  qdift(m,k,2)
           dqgrd   (m,k) =  qdift(m,k,1) -  qdift(m,k,2)
           somgr   (m,k) =  ogift(m,k,1) +  ogift(m,k,2)
           domgr   (m,k) =  ogift(m,k,1) -  ogift(m,k,2)
 801      continue
        endif
#endif
        do 81 m=1,mx
         spsift  (m  ) =  psift(m,  1) +  psift(m,  2)
         dpsift  (m  ) =  psift(m,  1) -  psift(m,  2)
         sspift  (m  ) =  spift(m,  1) +  spift(m,  2)
         dspift  (m  ) =  spift(m,  1) -  spift(m,  2)
 81     continue
c
      else
c
c  one hemisphere sum and difference fourier coefficients.
c   note that for odd functions (aift,zaift,tvift and qvift) the
c        sums are zero rather than the differences
        do 83 k=1,kx
         do 83 m=1,mx
         szaift (m,k) = (0.0,0.0)
         dzaift (m,k) = zaift(m,k,1)
         szbift (m,k) = zbift(m,k,1)
         dzbift (m,k) = (0.0,0.0)
         saift  (m,k) = (0.0,0.0)
         daift  (m,k) = aift(m,k,1)
         sbift  (m,k) = bift(m,k,1)
         dbift  (m,k) = (0.0,0.0)
         skift  (m,k) = ekift(m,k,1)
         dkift  (m,k) = (0.0,0.0)
         stift  (m,k) = tift(m,k,1)
         dtift  (m,k) = (0.0,0.0)
         stuift (m,k) = tuift(m,k,1)
         dtuift (m,k) = (0.0,0.0)
         stvift (m,k) = (0.0,0.0)
         dtvift (m,k) = tvift(m,k,1)
         sqeift (m,k) = qeift(m,k,1)
         dqeift (m,k) = (0.0, 0.0)
         squift (m,k) = quift(m,k,1)
         dquift (m,k) = (0.0, 0.0)
         sqvift (m,k) = (0.0, 0.0)
         dqvift (m,k) = qvift(m,k,1)
 83     continue
#ifdef dirnal
        if (sgmt) then
          do 802 k=1,kx
           do 802 m=1,mx
           sqgrd   (m,k) =  qdift(m,k,1)
           dqgrd   (m,k) =  (0.0,0.0)
           somgr   (m,k) =  ogift(m,k,1)
           domgr   (m,k) =  (0.0,0.0)
 802      continue
        endif
#endif
        do 84 m=1,mx
         spsift (m  ) = psift(m,  1)
         dpsift (m  ) = (0.0,0.0)
         sspift (m  ) = spift(m,  1)
         dspift (m  ) = (0.0,0.0)
 84     continue
      endif
c
c  in the calls to rhxinv below, the spectral harmonic is passed
c   twice as a 'get-a-round' to a compiler bug on the ymp which
c   slowed down the code.
c
#ifndef unitsk
      call rhxinv (    eky(1,1,1,ls),  skift,  dkift, colyl ,   wght,
     *          kx,    eky(1,1,1,ls) )
      call rhxinv ( vordty(1,1,1,ls), szaift, dzaift, colydx,   -1. ,
     *          kx, vordty(1,1,1,ls) )
      call rhxinv ( vordty(1,1,1,ls),  dbift,  sbift, colydy,    1. ,
     *          kx, vordty(1,1,1,ls) )
      call rhxinv ( divdty(1,1,1,ls), szbift, dzbift, colydx,    1. ,
     *          kx, divdty(1,1,1,ls) )
      call rhxinv ( divdty(1,1,1,ls),  daift,  saift, colydy,    1. ,
     *          kx, divdty(1,1,1,ls) )
      call rhxinv (   tdty(1,1,1,ls),  stift,  dtift, colyl ,   wght,
     *          kx,   tdty(1,1,1,ls) )
      call rhxinv (   tdty(1,1,1,ls), stuift, dtuift, colydx,   -1. ,
     *          kx,   tdty(1,1,1,ls) )
      call rhxinv (   tdty(1,1,1,ls), dtvift, stvift, colydy,    1. ,
     *          kx,   tdty(1,1,1,ls) )
      call rhxinv (  psdty(1,1,1,ls), spsift, dpsift, colyl ,   wght,
     *           1,  psdty(1,1,1,ls) )
      call rhxinv (    psy(1,1,1,ls), sspift, dspift, colyl , coswgt,
     *           1,    psy(1,1,1,ls) )
      call rhxinv ( qmxdty(1,1,1,ls), sqeift, dqeift, colyl , coswgt,
     *          kx, qmxdty(1,1,1,ls) )
      call rhxinv ( qmxdty(1,1,1,ls), squift, dquift, colydx, -1.   ,
     *          kx, qmxdty(1,1,1,ls) )
      call rhxinv ( qmxdty(1,1,1,ls), dqvift, sqvift, colydy,  1.   ,
     *          kx, qmxdty(1,1,1,ls) )
#ifdef dirnal
      if(sgmt) then
        call rhxinv (  qdoty(1,1,1,ls),  sqgrd,  dqgrd, colyl , coswgt,
     *             kx, qdoty(1,1,1,ls) )
        call rhxinv (  somgy(1,1,1,ls),  somgr,  domgr, colyl , coswgt,
     *             kx, somgy(1,1,1,ls) )
      endif
#endif
#endif
#ifdef unitsk
      call rhxinv (    ek (1,1,1   ),  skift,  dkift, colyl ,   wght,
     *          kx,    ek (1,1,1   ) )
      call rhxinv (  psdt (1,1,1   ), spsift, dpsift, colyl ,   wght,
     *           1,  psdt (1,1,1   ) )
      call rhxinv (    ps (1,1,1   ), sspift, dspift, colyl , coswgt,
     *           1,    ps (1,1,1   ) )
# ifndef new_rhxinv
      call rhxinv ( vordt (1,1,1   ), szaift, dzaift, colydx,   -1. ,
     *          kx, vordt (1,1,1   ) )
      call rhxinv ( vordt (1,1,1   ),  dbift,  sbift, colydy,    1. ,
     *          kx, vordt (1,1,1   ) )
      call rhxinv ( divdt (1,1,1   ), szbift, dzbift, colydx,    1. ,
     *          kx, divdt (1,1,1   ) )
      call rhxinv ( divdt (1,1,1   ),  daift,  saift, colydy,    1. ,
     *          kx, divdt (1,1,1   ) )
      call rhxinv (   tdt (1,1,1   ),  stift,  dtift, colyl ,   wght,
     *          kx,   tdt (1,1,1   ) )
      call rhxinv (   tdt (1,1,1   ), stuift, dtuift, colydx,   -1. ,
     *          kx,   tdt (1,1,1   ) )
      call rhxinv (   tdt (1,1,1   ), dtvift, stvift, colydy,    1. ,
     *          kx,   tdt (1,1,1   ) )
      call rhxinv ( qmxdt (1,1,1   ), sqeift, dqeift, colyl , coswgt,
     *          kx, qmxdt (1,1,1   ) )
      call rhxinv ( qmxdt (1,1,1   ), squift, dquift, colydx, -1.   ,
     *          kx, qmxdt (1,1,1   ) )
      call rhxinv ( qmxdt (1,1,1   ), dqvift, sqvift, colydy,  1.   ,
     *          kx, qmxdt (1,1,1   ) )
# endif
# ifdef new_rhxinv
      call rhxin2 ( vordt (1,1,1   ), szaift, dzaift,  dbift,  sbift,
     *colydx,colydy,-1.,1.,kx,vordt (1,1,1   ) )
      call rhxin2 ( divdt (1,1,1   ), szbift, dzbift,  daift,  saift,
     *colydx,colydy, 1.,1.,kx,divdt (1,1,1   ) )
      call rhxin3 (   tdt (1,1,1   ),  stift,  dtift, stuift, dtuift,
     *dtvift, stvift,colyl,colydx,colydy,wght,-1., 1.,kx,tdt(1,1,1))
      call rhxin3 ( qmxdt (1,1,1   ), sqeift, dqeift, squift, dquift,
     *dqvift, sqvift,colyl,colydx,colydy,coswgt,-1., 1.,kx,qmxdt(1,1,1))
# endif
# ifdef dirnal
      if(sgmt) then
        call rhxinv (  qdot (1,1,1   ),  sqgrd,  dqgrd, colyl , coswgt,
     *             kx, qdot (1,1,1   ) )
        call rhxinv (  somg (1,1,1   ),  somgr,  domgr, colyl , coswgt,
     *             kx, somg (1,1,1   ) )
      endif
# endif
#endif
c
c
      return
      end
