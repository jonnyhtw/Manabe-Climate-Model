      subroutine lwco2 ( tco2,ulog,plog,teff,ln )
c
c  table lookup routine for co2 transmission
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  routine:   lwco2
c  called by: lwtran
c  calls:     none
c  purpose:   compute longwave co2 transmission
c
c  inputs:
c     **passed** ulog           log of co2 amount
c                                  (from top of atmosphere)
c     **passed** plog           log of effective pressure
c     **passed** teff           effective temperature
c     **passed** ln             length of call-argument arrays
c     / co2d   / alogab         absorption(ulog,plog) table
c     / co2d   / alogdt         d(absorption)/d(temperature) table
c
c  output:
c     **passed** tco2           co2 transmission coefficient
c
c  side effects:
c     none
c     (i.e., no input variables or common blocks modified, etc.)
c
c  procedure:
c     a 2d table lookup of co2 absorption is performed.  the axes are
c     log of the co2 amount and log of the effective pressure.  an
c     effective temperature of 300k is assumed for this first step.
c     a correction for temperatures in the range 200k-300k is made by
c     multiplying a temperature correction factor by the difference
c     between the effective temperature and 300k, and adding the
c     result to the previously calculated absorption.  the temperature
c     correction factor is found by 2d table lookup, using the same
c     axes as were used to obtain the 300k absorption.  no temperature
c     correction is applied if the temperature exceeds 300k; if the
c     temperature falls below 200k, the 200k correction is applied.
c     if the co2 concentration is too low, no temperature correction
c     is made.
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c#include "parm.h"
#include "parm.h"
c
      parameter ( lmax = ix * (kx + 1) )
c
      dimension tco2  (  ln), ulog  (  ln), plog  (  ln), teff  (  ln)
      dimension ahi   (lmax), alo   (lmax), cct   (lmax), cdt   (lmax)
      dimension nplog (lmax), iulog (lmax), fulog (lmax), index (lmax)
      dimension plnum (lmax), plden (lmax), plfac (lmax), tcor  (lmax)
      dimension alogoo(lmax), alogpo(lmax), alogop(lmax), alogpp(lmax)
      dimension pltab (  12)
c
#include "co2lw.h"
c
      equivalence (alogoo, plnum), (alogpo, plden, alo, tcor)
cSGI
ccc   equivalence (alogop, iulog), (alogpp, index, ahi, cdt )
      equivalence (alogop, iulog), ( ahi, cdt )
cSGI
c
      data pltab /
     $     -4.50000, -4.00000, -3.40000, -2.88081, -2.27875, -1.88081,
     $     -1.57978, -1.27875, -0.88081, -0.48287,  0.     ,  0.29528 /
c
c ****** check for out-of-range ulog input values
c
      iu = 0
      do 25 l=1,ln
        if ( ulog(l).lt.-4.25 ) iu=iu+1
   25 continue
      if ( iu.ne.0 ) go to 1200
c
c ****** compute pressure range for tables,
c ****** and pressure interpolation factor
c
      do 50 l=1,ln
        nplog(l) = 1
   50 continue
c
      do 100  n=2,11
        do 100  l=1,ln
          if ( plog(l).gt.pltab(n) ) then
            nplog(l) = n
            plnum(l) = plog(l) - pltab(n)
            plden(l) = pltab(n+1) - pltab(n)
          endif
  100 continue
c
      do 200 l=1,ln
        plfac(l) = plnum(l) / plden(l)
  200 continue
c
c ****** compute co2 absorption for an effective temperature of 300k
c
      do 500 l=1,ln
        fulog(l) = 4 * ulog(l) + 18
        iulog(l) = int(fulog(l))
        if ( iulog(l).gt.29 ) iulog(l) = 29
        fulog(l) = fulog(l) - real(iulog(l))
  500 continue
c
      do 600 l = 1, ln
        index(l) = iulog(l) + 30 * (nplog(l) - 1)
  600 continue
c
      do 620 l=1,ln
        alogoo(l) = alogab(index(l),1)
  620 continue
      do 640 l=1,ln
        index (l) = index(l) + 1
        alogpo(l) = alogab(index(l),1)
  640 continue
      do 660 l=1,ln
        index (l) = index(l) + (30 - 1)
        alogop(l) = alogab(index(l),1)
  660 continue
      do 680 l=1,ln
        index (l) = index(l) + 1
        alogpp(l) = alogab(index(l),1)
  680 continue
c
      do 700  l = 1, ln
        ahi(l) = alogop(l) + fulog(l) * (alogpp(l) - alogop(l))
        alo(l) = alogoo(l) + fulog(l) * (alogpo(l) - alogoo(l))
        cct(l) = alo   (l) + plfac(l) * (ahi   (l) - alo   (l))
  700 continue
c
c ****** compute temperature correction coefficient for co2 absorption
c
      do 750 l=1,ln
        fulog(l) = ulog(l) + 4
        iulog(l) = int(fulog(l))
        if ( iulog(l).gt.6 ) iulog(l) = 6
        fulog(l) = fulog(l) - real(iulog(l))
  750 continue
c
      do 800 l = 1, ln
        index(l) = iulog(l) + 7 * (nplog(l) - 1)
  800 continue
c
      do 820 l=1,ln
        alogoo(l) = alogdt(index(l),1)
  820 continue
      do 840 l=1,ln
        index (l) = index(l) + 1
        alogpo(l) = alogdt(index(l),1)
  840 continue
      do 860 l=1,ln
        index (l) = index(l) + (7 - 1)
        alogop(l) = alogdt(index(l),1)
  860 continue
      do 880 l=1,ln
        index (l) = index(l) + 1
        alogpp(l) = alogdt(index(l),1)
  880 continue
c
      do 900 i = 1, ln
        ahi(i) = alogop(i) + fulog(i) * (alogpp(i) - alogop(i))
        alo(i) = alogoo(i) + fulog(i) * (alogpo(i) - alogoo(i))
        cdt(i) = alo   (i) + plfac(i) * (ahi   (i) - alo   (i))
  900 continue
c
c ****** restrict effective temperature input values to proper range
c
      do 950 l=1,ln
        tcor(l) = teff(l)
        if ( tcor(l).gt.320.0 ) tcor(l) = 320
        if ( tcor(l).lt.180.0 ) tcor(l) = 180
  950 continue
c
c ****** apply temperature correction to co2
c ****** absorption if ulog is in range
c
      do 1000  l=1,ln
        cdt (l) = cct(l) + 0.01 * (cdt(l) * (tcor(l) - 300))
        tco2(l) = 1 - 1.3077 * cdt(l)
 1000 continue
c
      do 1100 l=1,ln
        if ( ulog(l).lt.-3 ) tco2(l) = 1 - 1.3077 * cct(l)
 1100 continue
c
      return
c
c ****** fatal error: ulog out-of-range
c
 1200 write(6,1201) ulog
      call panic ( 'program stopped in lwco2.' )
c
 1201 format('1ulog > 3. or < -4.25; execution terminated' /
     $        ' the ulog vector is:' / (8e16.8))
c
      end
