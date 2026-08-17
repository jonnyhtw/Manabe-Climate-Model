      subroutine mstcnv(hl)
c
c     called by: convad
c     calls    : nothing
c
c     input : ta and ra
c     output: ta and ra
c
c this subroutine computes precipitation due to both large scale
c condensation and small scale moist convection.
c
#include "parm.h"
#include "parmkp.h"
#include "parmjx.h"
      parameter ( iters = 1 )
#include "parmout.h"
#include "parmleng.h"
c
      dimension     hl(ix        )
      dimension tauran(ix        ), tausno(ix     )
#ifdef prdcld
      dimension hh    (ix,kx     ), ihh   (ix,kx     )
#endif
      dimension wk    (ix,kx     ), r0    (ix,kx     ), alrm  (ix,  kx)
      dimension c     (ix,kx     )
      dimension ta    (ix,kx,nhem), ra    (ix,kx,nhem), psp   (ix,nhem)
      dimension rsmin (ix,kx)
c
#include "auxatm.h"
c
#ifdef covcld
#include "cloudcx.h"
#endif
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
#ifdef prdcld
#include "rhcrit.h"
#endif
c
#include "sgwgt1.h"
c
#include "surfc.h"
c
      dimension es    (ix,kx), des   (ix,kx), psq   (ix,kx)
      dimension rs    (ix,kx), mvf   (ix,kx)
      dimension psqm  (ix,kx), it    (ix,kx)
      dimension temnew(ix,kx), ratnew(ix,kx), temavg(ix,kx)
      dimension itest (ix,kx), iwk   (ix,kx), term2 (ix,kx)
      dimension wt1   (ix,kx), wt3   (ix,kx), delta (ix,kx)
      dimension sum0  (ix,kx), term1 (ix,kx), test  (ix,kx)
      dimension sum2  (ix,kp), sum1  (ix,kp)
#ifdef prdcld
      dimension mvd   (ix,kx)
#endif
c
#ifdef prdcld
      equivalence (ihh   (1,1),hh      (1,1))
#endif
      equivalence (ta  (1,1,1),cgt   (1,1,1)), (ra (1,1,1),gqmix(1,1,1))
      equivalence (psp (1,1)  ,gvbrps(1,1)),   (wk   (1,1),gvor (1,1,1))
      equivalence (r0    (1,1),gdiv  (1,1,1)), (alrm (1,1),gucos(1,1,1))
      equivalence (c     (1,1),gvcos (1,1,1)), (tauran (1),gpsdx  (1,1))
      equivalence (tausno  (1),gpsdy   (1,1))
c
c hc = critical relative humidity for condensation.
c
      data hc / 1.00 /
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
c initialize dummy arrays to zero for moist adiabatic adjustment
c computation.
c
      do 72 i=1,ix
       sum0(i,kx) = 0.0
       sum1(i, 1) = 0.0
       sum2(i, 1) = 0.0
   72 continue
c
c compute pressure at middle of layers for use in moist adiabatic
c lapse rate evaluation.
c
      do 130 k=1,km
        do 130 i=1,ix
          psqm(i,k) = psp(i,ih) * qmh(k)
  130 continue
c
c compute pressure at individual levels.
c
      do 140 k=1,kx
        do 140 i=1,ix
          psq(i,k) = psp(i,ih) * q(k)
  140 continue
c
c beginning of code for determining precipitation by both large
c and small scale processes.  compute es, c, and rs at the
c individual levels which are
c   es = saturation vapor pressure.
c    c = derivative of saturation mixing ratio with respect to
c        temperature.
c   rs = saturation mixing ratio.
c   rsmin = minimum saturation mixing ratio on es table
c   esmin = minimum saturation vapor pressure on es table
c   delmin = minimum p-es allowed, this is done to keep p-es
c            positive.
c
c note; repeated iteration through the column is no longer performed
c (i.e. iters=1).
c

      delmin = 1.0e-10
      esmin  = etabl(1)
      do 519 k=1,kx
       do 519 i=1,ix
        rsmin(i,k) = 0.622*esmin / (psq(i,k)-esmin)
  519 continue

      do 520 iter=1,iters
c
c es values are taken from a table(120k-375k) at .1 degree intervals
c
        do 149 k=1,kx
          do 149 i=1,ix
            it(i,k) = 10.0 * (ta(i,k,ih) + 153.05)
  149   continue
        do 80 k=1,kx
         do 80 i=1,ix
          it(i,k) = 2*it(i,k) + 1
          it(i,k) = max(it(i,k), 1)
          it(i,k) = min(it(i,k), 5102-3)
   80   continue
        do 81 k=1,kx
         do 81 i=1,ix
           es(i,k) = etabl(it(i,k)  )
          des(i,k) = etabl(it(i,k)+1)
   81   continue
        do 150 k=1,kx
          do 150 i=1,ix
            wk(i,k)    = 0.622 * hc *  es(i,k)
            des(i,k)   = 0.622 * hc * des(i,k)
            delta(i,k) = psq(i,k) - hc *  es(i,k)
  150   continue
c
c  make sure delta, (press - es press) does not become negative
        do 82 k=1,kx
         do 82 i=1,ix
          delta(i,k) = max(delta(i,k), delmin)
   82   continue
        do 151 k=1,kx
          do 151 i=1,ix
             c(i,k) = psq(i,k) * des(i,k) / (delta(i,k) * delta(i,k))
            rs(i,k) = wk(i,k) / delta(i,k)
  151   continue
        do 83 k=1,kx
         do 83 i=1,ix
          rs(i,k) = max(rs(i,k), rsmin(i,k) )
   83   continue
c
#ifdef prdcld
c determine levels where clouds are to be placed.  used for cloud
c prediction only.
c
        do 180 k=1,kx
         do 180 i=1,ix
# ifndef varhy
          hh(i,k) = r0(i,k) - hcrit(k,jrow)*rs(i,k)
# endif
# ifdef varhy
          hh(i,k) = r0(i,k) - hyk(k)*rs(i,k)
# endif
  180   continue
        do 183 k=1,kx
         do 183 i=1,ix
          if (hh(i,k) .gt. 0.0) then
            ihh(i,k) = 1
          else
            ihh(i,k) = 0
          endif
  183   continue
#endif
c
c compute mean temperature of layer for use in moist adiabatic
c lapse rate evaluation.
c
        do 160 k=1,km
          do 160 i=1,ix
            temavg(i,k) = (ta(i,k,ih) + ta(i,k+1,ih)) * 0.5
  160   continue
c
c compute es and des (derivative of saturation vapor pressure)
c at the mid-point of the layer for use in moist adiabatic lapse
c rate evaluation.
c
        do 169 k=1,km
          do 169 i=1,ix
            it(i,k) = 10.0 * (temavg(i,k) + 153.05)
  169   continue
        do 85 k=1,kx
         do 85 i=1,ix
          it(i,k) = 2*it(i,k) + 1
          it(i,k) = max(it(i,k), 1)
          it(i,k) = min(it(i,k), 5102-3)
   85   continue
        do 86 k=1,kx
         do 86 i=1,ix
           es(i,k) = etabl(it(i,k)  )
          des(i,k) = etabl(it(i,k)+1)
   86   continue
        do 170 k=1,km
          do 170 i=1,ix
            es (i,k) = 0.622 * hc * es (i,k)
            des(i,k) = 0.622 * hc * des(i,k)
  170   continue
c
c compute alrm (moist adiabatic lapse rate) for each layer.
c
        do 171 k=1,km
          do 171 i=1,ix
            alrm(i,k) = (dqph(k) * (ara * psqm(i,k) *
     *                  (temavg(i,k) + tfra) + hl(i) * es(i,k)))
     *                / (qmh(k) * (psqm(i,k) * cpp
     *                               + hl(i) * des(i,k)))
  171   continue
c
c tests are now performed to determine occurrence and type of
c precipitation.  for small scale, moist convection, layer must be
c both super-moist adiabatic and supersaturated.  for large scale
c condensation, just the individual layer must be supersaturated.
c
c here begins the computation of the integer arrays iwk and itest
c which will determine the levels where the moist convective and
c large scale adjustments occur, respectively.  mvf is the same as
c iwk except it applies to the layers rather than the levels
c involved.
c
        if (iter .le. 1) then
c
c         test(i,k) = 1.0 - ra(i,k,ih) / rs(i,k)
c
c set test to negative number for supersaturation.
c
          do 88 k=1,kx
           do 88 i=1,ix
            test(i,k) = wk(i,k) - ra(i,k,ih)*delta(i,k)
   88     continue
          do 89 k=1,kx
           do 89 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   89     continue
c
          do 181 i=1,ix
            mvf(i,kx) = 0
  181     continue
c
          do 90 k=1,km
           do 90 i=1,ix
            mvf(i,k) = itest(i,k) * itest(i,k+1)
   90     continue
c
          do 91 k=1,km
           do 91 i=1,ix
            test(i,k) = alrm(i,k) + 0.01 + ta(i,k,ih) - ta(i,k+1,ih)
   91     continue
c
          do 92 k=1,km
           do 92 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   92     continue
c
          do 93 k=1,km
           do 93 i=1,ix
            mvf(i,k) = mvf(i,k) * itest(i,k)
   93     continue
c
c     mvf(i,k)=1 where moist adiabatic adjustment is required
c     mvf(i,k)=0 where either condensation or nothing is required
c
          do 94 i=1,ix
           iwk(i,1) = mvf(i,1)
   94     continue
c
          do 95 k=2,kx
           do 95 i=1,ix
            iwk(i,k) = max(mvf(i,k),mvf(i,k-1))
   95     continue
c
          do 96 k=1,kx
           do 96 i=1,ix
            test(i,k) = rs(i,k) - ra(i,k,ih)
   96     continue
          do 97 k=1,kx
           do 97 i=1,ix
            if (test(i,k) .ge. 0.0) then
              itest(i,k) = 0
            else
              itest(i,k) = 1
            endif
   97     continue
c
          do 98 k=1,kx
           do 98 i=1,ix
            itest(i,k) = itest(i,k) * (1-iwk(i,k))
   98     continue
c
          do 99 k=1,kx
           do 99 i=1,ix
            wt1(i,k) = mvf(i,k)
            wt3(i,k) = iwk(i,k)
   99     continue
c
c iwk now contains the convection flagbits and itest the
c condensation flagbits: iwk=1 if and only if convection.
c itest=1 if and only if condensation, otherwise they equal zero.
c the value of zero implies no adjustment.
c
#ifdef prdcld
c compute "packed" flagbit array for cloud prediction based on
c previous computation of hh and ihh.
c
          if (lrad) then
            do 182 i=1,ix
             kflags(i,jrow) = 0.0
#ifdef covcld
             kfgcov(i,jrow) = 0.0
#endif
  182       continue
            do 3009 k=1,kx
              do 3009 i=1,ix
                mvd(i,k) = or(iwk(i,k),itest(i,k))
 3009       continue
#ifdef no_kflags_k1
c  Do not allow clouds (kflags) at k=1 level (top of the atm)
            do 3010 i=1,ix
             ihh(i,1) = 0
             iwk(i,1) = 0
 3010       continue
#endif
            do 3005 i=1,ix
              do 3007 k=1,kx
                kflags(i,jrow) = kflags(i,jrow) + ipower(k) * ihh(i,k)
#ifdef covcld
                kfgcov(i,jrow) = kfgcov(i,jrow) + ipower(k) * iwk(i,k)
#endif
 3007         continue
 3005       continue
          endif
#endif
c
c solution of simultaneous equations for the new temperatures and
c mixing ratios as a result of both condensation and convective
c adjustments.  here
c       ta = temperature after dry adjustment.
c       ra = old mixing ratio
c   temnew = new temperature
c   ratnew = new mixing ratio
c
          iflag1 = 0
          iflag2 = 0
          do 188 k=1,kx
           do 188 i=1,ix
            iflag1 = iflag1 + itest(i,k)
            iflag2 = iflag2 +   mvf(i,k)
  188     continue
          if (iflag1 .eq. 0 .and. iflag2 .eq. 0) go to 525
          if (iflag1 .eq. 0) then
            do 3376 k=1,kx
              do 3376 i=1,ix
                temnew(i,k) = ta(i,k,ih)
                ratnew(i,k) = ra(i,k,ih)
 3376       continue
          else
c
c large scale condensation adjustment.
c
            do 3333 k=1,kx
              do 3333 i=1,ix
                test(i,k) = itest(i,k)
 3333       continue
c
c **    itest now has a "1" where condensation is allowed to occur,
c        i.e., where moist adjustment did not occur and where the
c        humidity (ratnew/rs) is greater than 1.0, and "0" otherwise.
c
            do 3004 k=1,kx
              do 3004 i=1,ix
                temnew(i,k) = ta(i,k,ih) +
     *                      test(i,k) * (hl(i) * (ra(i,k,ih) - rs(i,k))
     *                        / (cpp + hl(i) * c(i,k)))
                ratnew(i,k) = test(i,k) * (rs(i,k) + c(i,k)
     *                                     * (temnew(i,k) - ta(i,k,ih))
     *                                     - ra(i,k,ih)) + ra(i,k,ih)
 3004       continue
          endif
c
        endif
c
        if (iter.gt.1) then
          do 3006 k=1,kx
            do 3006 i=1,ix
              temnew(i,k) = ta(i,k,ih)
     *                    + test(i,k) * (hl(i) * (ra(i,k,ih) - rs(i,k))
     *                      / (cpp + hl(i) * c(i,k)))
              ratnew(i,k) = test(i,k) * (rs(i,k) + c(i,k)
     *                                   * (temnew(i,k) - ta(i,k,ih))
     *                                   - ra(i,k,ih)) + ra(i,k,ih)
 3006     continue
        endif
c
c small scale moist convective adjustment.
c
        if (iflag2 .ne. 0) then
          do 700 k=1,km
            do 700 i=1,ix
              sum0(i,kx-k) = wt1(i,kx-k) * (sum0(i,kx+1-k)
     *                                      + alrm(i,kx-k))
  700     continue
          do 701 k=1,kx
            do 701 i=1,ix
              term2(i,k) = dq(k) * (cpp + hl(i) * c(i,k))
              term1(i,k) = term2(i,k) * (ta(i,k,ih) + sum0(i,k))
     *                   + dq(k) * hl(i) * (ra(i,k,ih) - rs(i,k))
  701     continue
          do 702 k=1,kx
            do 706 i=1,ix
              sum1(i,k+1) = wt3(i,k) * (sum1(i,k) + term1(i,k))
              sum2(i,k+1) = wt3(i,k) * (sum2(i,k) + term2(i,k))
  706       continue
            do 172 i=1,ix
             if (wt1(i,k) .ne. wt3(i,k)) then
               temnew(i,k  ) = sum1(i,k+1) / sum2(i,k+1)
                 sum1(i,k+1) = 0.0
                 sum2(i,k+1) = 0.0
             endif
  172       continue
  702     continue
          do 703 k=1,km
            do 703 i=1,ix
              temnew(i,kx-k) = (1.0 - wt1(i,kx-k)) * temnew(i,kx  -k)
     *                       +        wt1(i,kx-k)  * temnew(i,kx+1-k)
  703     continue
          do 704 k=1,km
            do 704 i=1,ix
              temnew(i,k) = temnew(i,k) - sum0(i,k)
  704     continue
          do 705 k=1,kx
            do 705 i=1,ix
              ratnew(i,k) = (rs(i,k) + c(i,k)
     *                       * (temnew(i,k) - ta(i,k,ih))) * wt3(i,k)
     *                    + (1.0 - wt3(i,k)) * ratnew(i,k)
  705     continue
        endif
c
c end of condensation and convective adjustment algorithm.  replace
c previous temperatures and mixing ratios with new ones.
c
        do 510 k=1,kx
          do 510 i=1,ix
            ta(i,k,ih) = temnew(i,k)
            ra(i,k,ih) = ratnew(i,k)
  510   continue
  520 continue
  525 continue
c
      return
      end
