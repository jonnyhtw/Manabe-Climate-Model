      subroutine march
c
c  called by: atmos and indata
c  subroutines called: solar
c
c  subroutine march computes the current month, day, year of the model
c  (based on the forecast time in seconds) every time the subroutine
c  is called. in addition, on a radiation timestep the following are
c  computed:
c   ozone (latitude-height distribution)
c   astronomical data (in subroutine solar)
c   table of solar radiation data (cosz, tauda, etc) vs latitude
c   sea surface temperatures
c
c   this version assumes each year has 365.000 days
c
#include "parm.h"
#include "parmkp.h"
#include "parms.h"
#ifdef fixcld
#include "parmncc.h"
      parameter ( lx = ncc + 1 )
#endif
#ifdef prdcld
      parameter ( lx = kx + 2 )
#endif
c
c
      character*8 months(12)
#ifdef cpldleg1
      character*8 monb
#endif
c
#include "astro.h"
c
#include "calndr.h"
c
#include "cld1.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "latdat.h"
c
#include "lats.h"
c
#include "param.h"
c
#include "sstcof.h"
c
#include "surfo.h"
c
c
cSGI
      integer*8 ndayom
cSGI
      equivalence (ndayom,iday)
c
c  specify phase angles (radians) for ozone & sea surface temperatures
c  based on the middle of the month
c
      data argaug / 4.70916335 /
      data argfeb / 1.56972112 /
c
      data cyear / 365.00 /
c
c
c
      data months / 'jan.    ', 'feb.    ', 'mar.    ', 
     *              'apr.    ', 'may     ', 'june    ',
     *              'july    ', 'aug.    ', 'sep.    ', 
     *              'oct.    ', 'nov.    ', 'dec.    ' /
      data sign / 1h- /, sigb / 1h /
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
c  definitions of arrays and variables for the astronomy:
c         ccosz     cosine of sun zenith angle
c         dayrad    day of the year for radiation
c         dlt       sun declination angle (delta)
c         ha        hour angle of sun at sunset
c         itime     time of day
c         iyear     year
c         month     month (a4)
c         ndayom    day of the month
c         ph        latitude in radians
c         rvct      distance to the sun in astronomical units
c         ttauda    amount of sunshine
c         idyorg    origin of day count.  the julian day corresponds to
c                     january 1, 1970 at noon ut.  (to change day count
c                     change idyorg, preferably in increments of 365)
c
c  compute julian day and fraction from ifcst, the number of seconds
c  since january 1, 1970 at 0 hours ut
c
      idyorg = 2440588
      lday   = ifcst / 86400
      fjd = (float(ifcst) - float(lday)*86400.0)/86400.0 - 0.5
      if (fjd .lt. 0.0) then
        lday = lday - 1
        fjd  = fjd + 1.0
      endif
      jld  = lday + idyorg
c
#ifdef cpldleg1
c  use monb below to set monend=1 if this is the first timestep
c   of a new month
      monb = month
c
#endif
c   compute month,day, and year
      nday = day
      ihour=0
c   date doesnt change from 00z to 12z, just the day for interpolation
      if((day-float(nday)).gt.0.1) ihour=12
      nday = nday+1
      if(nday.gt.iyear*365) iyear = iyear+1
      nday = mod(nday-1,365)+1
      do 56 i=1,12
cSGI
cc    if(nday.gt.dy(i+1)) go to 56
      if(float(nday).gt.dy(i+1)) go to 56
cSGI
      ndayom = nday-dy(i)
      go to 57
   56 continue
   57 month  = months(i)
      astday = amod(day,cyear)
      zmin   = 0.0
c
#ifdef cpldleg1
      monend = 0
      if (monb .ne. month) monend = 1
c
      if (lmon) return
c
#endif
      if (.not.lrad) return
c
c  compute latitudes in radians (xlat) and degrees (alat)
c
      do 30 j=1,iy
        xlat(j) = hpi - acos(coa(j))
        if (nhem .ne. 1) xlat(ilp-j) = -hpi + acos(coa(j))
  30  continue
      do 40 j=1,il
        alat(j) = ylat(j)
  40  continue
c
c  compute ozone by time interpolation of seasonal data from arrays
c  xo3ws and xo3sf in tables above.
c  adjustments to dayrad are phase angles (days)
c
      dayrad = astday
      fday   = amod(dayrad + 46.25, cyear)
      bfday  = twopi * fday / cyear
      arg    = twopi * (dayrad - 12.6875) / cyear
      cosbf  = cos(bfday)
      sinbf  = sin(bfday)
      cos2bf = cos(bfday * 2.0)
      cosfeb = cos(argfeb)
      sinfeb = sin(argfeb)
      cos2fb = cos(argfeb * 2.0)
      cosaug = cos(argaug)
      sinaug = sin(argaug)
      cos2ag = cos(argaug * 2.0)
      cosarg = cos(arg)
      sinarg = sin(arg)
      cos2ar = cos(arg * 2.0)
c
      do 20 k=1,kx
        do 20 l=1,19
          lsc           = 20 - l
          xbar          = 0.25 * (xo3ws(l,k) + xo3ws(lsc,k)
     *                            + xo3sf(l,k) + xo3sf(lsc,k))
          a1            = 0.50 * (xo3sf(l,k) - xo3sf(lsc,k))
          b1            = 0.50 * (xo3ws(l,k) - xo3ws(lsc,k))
          b2            = 0.25 * (xo3ws(l,k) + xo3ws(lsc,k)
     *                            - xo3sf(l,k) - xo3sf(lsc,k))
          dduo3n(lsc,k) = xbar + a1 * sinarg + b1 * cosarg + b2 * cos2ar
   20 continue
c
c
#ifdef dirnal
      if (.not.kdisk) return
      jldsv = jld
c
c  intsr is the interval (days) between solar radiation tendency updates
c
      intsr = intsw / 86400
c
c  compute number of days to advance date to place radiation results in
c  middle of period for which tendencies are used
c
      ladv = (intsr) / 2
      if (resets) ladv = (intsw - mod(ifcst+idt,intsw)) / 172800
      jld = jld + ladv
#endif
c
c  call subroutine solar to compute astronomical data from
c  julian day and fraction
c
#ifdef solar_berger
      call solar (rvct,dlt,alp,slag,il,xlat,ha,ccosz,ttauda)
#else
      call solar (jld,fjd,rvct,dlt,alp,slag,il,xlat,ha,ccosz,ttauda)
#endif
c
#ifdef dirnal
      jld = jldsv
#endif
      write(6,101) ktau, ihour, zmin, ndayom, month, iyear, dayrad
      solc = solarc / ( rvct * rvct )
#ifdef dirnal
      if (.not.kdisk) then
#endif
c
c  print table of cosz, tauda and radiation versus latitude in all but
c  diurnal models
c
#if defined cpldleg3 || defined cpldleg1
      ndayp = day + 0.01
      if (mod(ndayp,5).eq.0) then
#endif
        write(6,102)
        do 19 jj=1,iy
          jr = iy + jj
          xinsol  = solc * ttauda(jj) * ccosz(jj)
          xinsolr = solc * ttauda(jr) * ccosz(jr)
          write(6,103) jj,alat(jj),ttauda(jj),ccosz(jj),xinsol,
     *                 jr,alat(jr),ttauda(jr),ccosz(jr),xinsolr
   19   continue
#if defined cpldleg3 || defined cpldleg1
      endif
#endif
#ifdef dirnal
      endif
#endif
      ddlt = 180.0 * dlt / pi
      idlt = ddlt
      fid  = idlt
      dsig = sigb
      if ((ddlt .lt. 0.0) .and. (idlt .eq. 0)) dsig = sign
      xmin = abs(ddlt-fid) * 60.0
      ixx  = xmin
      fix  = ixx
      xsec = abs(xmin-fix) * 60.0
      halp = 12.0 * alp / pi
      if (halp .lt. 0.0) halp = halp + 24.0
      ihalp = halp
      fih   = ihalp
      ymin  = abs(halp-fih) * 60.0
      iyy   = ymin
      fiy   = iyy
      asec  = abs(ymin-fiy) * 60.0
c
c  print radius vector (distance to sun in astronomical units) and solar
c  right ascension and declination obtained from subroutine solar
c
      write(6,104) rvct, ihalp, iyy, asec, dsig, idlt, ixx, xsec
c
#ifdef dirnal
c  print equation of time (apparent minus mean solar time) in diurnal
c  models, where it is relevant
c
      sslag = 4.0 * radian * slag
      write(6,109) sslag
      if (ladv .gt. 0) write(6,105) ladv
#endif
      return
c
  101 format(' ktau = ',i8,5x,'forecast time',i3,' hr',f5.1,' min z',i4,
     *       1x,a4,i5,5x,'day for sst and ozone =',f10.4)
  102 format(1x,3hrow,4x,8hlatitude,3x,13hmean sunshine,4x,4hcosz,
     *       5x,19hinsolation (ly/min),
     +       4x,3hrow,4x,8hlatitude,3x,13hmean sunshine,4x,4hcosz,
     *       5x,19hinsolation (ly/min))
  103 format(1x,i3,f12.3,3f14.6,10x,i3,f12.3,3f14.6)
  104 format(' radius vector=',f10.7,5x,'r.a. of sun=',i3,'hr',i3,
     *       'min',f6.2,'sec',5x,'declination=',a1,i3,'deg',
     *       i3,'min',f5.1,'sec')
#ifdef dirnal
  105 format('  date for astronomy advanced by',i3,' days')
  109 format('  equation of time=',f9.4,' minutes')
#endif
c
      end
