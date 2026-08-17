      subroutine convad
c
c    called by: passb
c    calls    : drycnv, mstcnv, smsc
c
c this subroutine controls the computation of convective adjustment.
c the input is pressure, temperature and mixing ratio.
c the output are the new adjusted temps and mixing ratios.
c
#include "parm.h"
#include "parmkp.h"
#include "parmout.h"
#include "parmleng.h"
c
c
      dimension tauran(ix), tausno(ix )
      dimension cmcah (ix), taupre(ix)
      dimension tcc   (ix)
      dimension  hl   (ix)
      dimension r0    (ix,kx)
      dimension t0    (ix,kx)
      dimension temp  (ix,kx)
      dimension es    (ix,kx),     des(ix,kx),      psq(ix,kx)
      dimension ta    (ix,kx,nhem), ra(ix,kx,nhem), psp(ix,nhem)
c
#include "auxatm.h"
c
#include "cloudp.h"
c
#include "ctlblk.h"
c
#include "establ.h"
c
#include "grdfld.h"
c
#include "grid.h"
c
#include "param.h"
c
#include "qa.h"
c
#include "intgrl.h"
c
#include "sgwgt1.h"
c
#include "surfc.h"
c
c
      equivalence (ta  (1,1,1),cgt   (1,1,1)), (ra (1,1,1),gqmix(1,1,1))
      equivalence (psp (1,1)  ,gvbrps(1,1))
      equivalence (r0    (1,1),gdiv  (1,1,1))
      equivalence (tauran  (1),gpsdx  (1,1))
      equivalence (tausno  (1),gpsdy   (1,1))
c
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c  convert log(p*) to p* and store in psp before use
c
      do 1 i=1,ix
        psp(i,ih) = exp(gps(i,ih))
    1 continue
c
c test temperatures for range between -173.16 and +102 degrees c.
c if any temperature is outside of this range, program stops.
c
      fact = 1.0
      if (ktau .eq. 0) fact = 0.5
      do 70 k=1,kx
       do 70 i=1,ix
        temp(i,k) = 102.0 - ta(i,k,ih)
        temp(i,k) = temp(i,k) * (ta(i,k,ih)+173.16)
   70 continue
      ixxx = 0
      do 71 k=1,kx
       do 71 i=1,ix
        if (temp(i,k) .lt. 0.0) go to 771
   71 continue
      ixxx = 0
      go to 772
  771 continue
      ixxx = 1
  772 continue
      if (ixxx .ne. 0)  then
        write(6,9030) jrow, ktau
        do 7770 k=1,kx
          write(6,9031) k, (ta(i,k,ih),i=1,ix)
7770    continue
        call aquit(20)
      endif
c
c save initial mixing ratio and temperature.
c
      do 110 k=1,kx
        do 110 i=1,ix
          r0 (i,k) = ra (i,k,ih)
          t0 (i,k) = ta (i,k,ih)
  110 continue
c
c test to determine whether precipitation occurs as rain or snow.
c heat of condensation is adjusted accordingly.
c
c tcc(i)=+1.0 when the boundary layer temperature .gt.0 deg.c
c tcc(i)= 0.0 when the boundary layer temperature .le.0 deg.c
c
c height of freezing level is taken to be 350 meters based on
c empirical analysis by nmc.  krs1, krs2, rswt1, and rswt2 are
c determined in subroutine parmtr.
c
      do 74 i=1,ix
       tcc(i) = -rswt1*ta(i,krs1,ih) - rswt2*ta(i,krs2,ih)
   74 continue
      do 75 i=1,ix
       if (tcc(i) .lt. 0.0) then
         hl (i) = 600.0
         tcc(i) =   1.0
       else
         hl (i) = 680.0
         tcc(i) =   0.0
       endif
   75 continue
c
#ifndef garratt
c
c  call drycnv which computes the dry adjustment.
c
      call drycnv(t0)
#endif
c
c
c  call mstcnv which computes the moist adjustment.
c
      call mstcnv(hl)
c
c
c compute change of temperature and mixing ratio due to condensation,
c convective and dry adjustment.
c
      do 531 k=1,kx
        do 531 i=1,ix
          convt(i,k) = ta(i,k,ih) - t0(i,k)
          convr(i,k) = ra(i,k,ih) - r0(i,k)
  531 continue
c
c compute total precipitation and heating due to all adjustments
c in column.
c
c convert convt and convr to rates of change.
c
      do 532 i=1,ix
        cmcah (i) =  dq(1) * convt(i,1) * psp(i,ih)
        taupre(i) = -dq(1) * convr(i,1)
  532 continue
      do 534 k=2,kx
        do 534 i=1,ix
          cmcah (i) = cmcah (i) + dq(k) * convt(i,k) * psp(i,ih)
          taupre(i) = taupre(i) - dq(k) * convr(i,k)
  534 continue
      recg = og
      do 550 i=1,ix
        taupre(i) = taupre(i) * psp(i,ih) * recg
  550 continue
      rdt = 1.0 / float ( idt )
      do 535 k=1,kx
        do 535 i=1,ix
          convt(i,k) = rdt * convt(i,k)
          convr(i,k) = rdt * convr(i,k)
  535 continue
c
c  compute relative humidity for analysis and store in field relhum in
c    the analb part of common block /aux/
c
#ifdef dirnal
      if (sgmt) then
        call esndif(ta(1,1,ih),es,des,-leng)
        do 560 k=1,kx
          do 560 i=1,ix
               psq(i,k) = psp(i,ih) * q(k)
            relhum(i,k) = ra(i,k,ih) * (psq(i,k) - es(i,k))
     *                    / (0.622 * es(i,k))
  560   continue
      endif
#endif
      if (rsrt) then
        call esndif(ta(1,1,ih),es,des,-leng)
        do 562 k=1,kx
          do 562 i=1,ix
               psq(i,k) = psp(i,ih) * q(k)
            relhum(i,k) = ra(i,k,ih) * (psq(i,k) - es(i,k))
     *                    / (0.622 * es(i,k))
  562   continue
      endif
c
c set total precipitation as either rain or snow depending upon
c criterion tcc determined near beginning of subroutine.
c
c
c taupre is computed over a 2dt time interval and, therefore,
c must be divided in half.
c
      do 600 i=1,ix
        taupre(i) = taupre(i) * 0.5 * fact
        tauran(i) = taupre(i) *        tcc(i)
        tausno(i) = taupre(i) * (1.0 - tcc(i))
  600 continue
#ifdef moist
      do 630 i=1,ix
        precip(i,jrow) = precip(i,jrow) + taupre(i)
        snowfl(i,jrow) = snowfl(i,jrow) + taupre(i) * (1.0 - tcc(i))
 630  continue
#endif
c
c   for air - sea models
c   tcc=1      pkice>0
c   tcc=0      pkice=0
c   presence of ice determines whether altnt is subtracted from
c    tauran or tausno
c
      do 173 i=1,ix
       if (pkice(i) .gt. 0.0 ) then
         tcc(i) = 1.0
       else
         tcc(i) = 0.0
       endif
  173 continue
      do 625 i=1,ix
#ifdef ice_dynamics
c
c sublimation is included in sno from subroutine vertdf (cfi)
c
      sno(i,jrow) = sno(i,jrow) + tausno(i)*(1.0-zin(i))
#else
      sno(i,jrow)=sno(i,jrow)+(tausno(i)-altnt(i,jrow)*tcc(i))*
     *                        (1.0-zin(i))
#endif
#if defined cpldleg3 || defined cpldleg1
c
c altnt contains evaporation from open water in ice grid points (cfi)
c
#ifdef ice_dynamics
      per(i,jrow)=per(i,jrow)+(tauran(i)-altnt(i,jrow))*
     *                        (1.0-zin(i))
#else
      per(i,jrow)=per(i,jrow)+(tauran(i)-altnt(i,jrow)*
     *                        (1.0-tcc(i)))*
     *                        (1.0-zin(i))
#endif
#endif
  625 continue
c   altnt units are cm/dt, its units are changed in passa.
c   it is calc in vertdf and stored in fluxw
c
c
c compute zonal means of precipitation, rain, snow and total heating
c due to all adjustments.
c
      do 185 i=1,ix
       zmpc (jrow) = zmpc (jrow) + rix * taupre(i)
       zmrn (jrow) = zmrn (jrow) + rix * tauran(i)
       zmsn (jrow) = zmsn (jrow) + rix * tausno(i)
       zmcah(jrow) = zmcah(jrow) + rix * cmcah (i)
  185 continue
c
#ifdef hydrol
c compute ground hydrology.
c
      call smsc
c
#endif
9030  format('1 temperatures bad in convad/passb at row = ',
     *        i3,' ktau = ',i10)
9031  format(1x,'t(i,k), k=',i3,/,8(1x,e14.6))
c
      return
      end
