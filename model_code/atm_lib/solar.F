#ifdef solar_berger
      subroutine solar(r,dlt,alp,slag,n,alat,hang,cosz,frac)
c
c  Note: this solar code was provided by ajb
c
c  subroutine solar is called by march and computes the radius
c    vector, the declination and right ascension of the sun, the
c    apparent sun lag angle (related to the equation of time), the hour
c    angle of the sun at sunset, the mean cosine of the sun zenith
c    angle and daylight fraction for n specified latitudes given the
c    day and fraction.  see notes of august 1989.
c
c    definition of the arguments used by this subroutine -
c
c      r     =        radius vector (distance to sun in a. u.)
c      dlt   =        declination of sun
c      alp   =        right ascension of sun
c      slag  =        apparent sun lag angle (west of mean sun is plus)
c      n     =        number of latitudes
c      alat  =        latitude values dimensioned by n
c      hang  =        hour angle of sun at sunrise or sunset
c      cosz  =        mean cosine of zenith angle
c      frac  =        fraction of daylight
c
c     note - all angles are expressed in radians
c            input variables are day in /ctlblk/, n and alat
c
#include "parm.h"
c
      dimension alat(n), hang(n)
      dimension cosz(n), frac(n)
c
c
#include "ctlblk.h"
c
#include "param.h"
c
c  value of ccr depends on the precision of the machine
c
      data ccr / 1.0e-06 /
c
      piph   = pi  + hpi
ccccc radian = pi  / 180.0
c
c  obliquity in degrees is oblqty
#ifdef solar_6k
      data oblqty  / 24.105 /
#else
#ifdef solar_21k
      data oblqty  / 22.949 /
#else
      data oblqty  / 23.446 /
#endif
#endif
c
c  longitude of earth's perihelion in degrees is peril
#ifdef solar_6k
      data peril  /   0.87 /
#else
#ifdef solar_21k
      data peril  / 114.42 /
#else
      data peril  / 102.04 /
#endif
#endif
c
c  eccentricity of earth's orbit is ec
#ifdef solar_6k
      data ec  / 0.018682 /
#else
#ifdef solar_21k
      data ec  / 0.018994 /
#else
      data ec  / 0.016724 /
#endif
#endif
c
c   set fundamental astronomical parameters
      ecsq = ec**2
ccccc angin  = radian * oblqty
      angin  = oblqty / radian 
      sni    = sin(angin)
c
c   tropical or calendar year in mean solar days is cyear
      cyear = 365.0
c
c   crate is mean daily angular motion of planet in degrees with 
c   respect to equinox
      crate = 360.0 / cyear
c
c   semimajor axis of planetary orbit in astronomical units is sma.
c    one a. u. = 149,600,000 km.  for earth sma=1.
      sma = 1.0
c
c  The reference point for the astronomical calendar is the vernal
c  equinox, which is defined to occur at 0.0 hours on March 21. This
c  occurs 79 days after January 1, which is the reference point for
c  the model calendar.
       en = day - 79.0
c
c   longitude of sun at perihelion in degrees is sunper.  this is
c    exactly 180 degrees from that of earth at perihelion.
c
      sunper = peril - 180.0
      craten = crate * en
c
c  computation of longitude of mean sun follows the procedure
c  described in Berger (JAS, vol. 35, 2362-2367), using eq. 2
c  from annex 2 of Berger's paper "Precession, eccentricity, obliquity,
c  insolation and paleoclimates" from the NATO ASI "Long Term
c  Climatic Variations, Data and Modelling" (J. C. Duplessy, ed.).
c  The equation is a variant of that given on p.2365 of the JAS paper.
      cosphi = sqrt(1.0-ecsq)
      beta = (1.0-cosphi)/ec
ccccc peri = radian*sunper
      peri = sunper / radian
c
c  longitude of mean sun at vernal equinox is velngm
      velngm = -2.0*(beta*(1.0+cosphi)*sin(-peri)) +
     $          beta*beta*(0.5+cosphi)*sin(-2.0*peri) +
     $     beta*beta*beta*(1.0/3.0+cosphi)*sin(-3.0*peri)
c
c  mean longitude of sun in radians is el
ccccc el = velngm + craten*radian
      el = velngm + craten/radian
c
c  constrain el to range of zero to 2 * pi
      el = el - twopi*ifix(el/twopi)
      if (el .lt. 0.0) el= el+twopi
      elt   = el
c
c  ecliptic longitude of sun in radians is rlam
c  Berger (JAS, vol. 35, p. 2365)
      rlam = el + 2.0*ec*sin(el-peri)
     $             + 5.0*ecsq*sin(2.0*(el-peri))/4.0
     $             + 13.0*ec*ecsq*sin(3.0*(el-peri))/12.0
c
c   constrain rlam to range of zero to 2 * pi
      rlam = rlam - twopi*ifix(rlam/twopi)
      if (rlam .lt. 0.0) rlam = rlam + twopi
c
c   compute right ascension of sun (alp) in radians
      alp = atan(cos(angin)*tan(rlam))
c
c   following 6 statements determine proper quadrant for root of
c    equation for alp above and insure that the difference slag
c    below has an absolute magnitude of order pi/10.
c
      if (alp .lt. 0.0) alp = alp + twopi
c
      if ((elt-alp) .gt. piph) elt = elt - twopi
      if ((elt-alp) .lt.-piph) elt = elt + twopi
c
      if (abs(elt-alp) .gt. hpi) then
         if ((elt-alp) .gt. hpi) alp=alp+pi
         if ((elt-alp) .lt.-hpi) alp=alp-pi
      endif
c
c   compute declination of sun (dlt) in radians
      sind = sni*sin(rlam)
      dlt  = asin(sind)
c
c   compute radius vector (distance to sun in a.u.)
      r = (sma*(1.0-ecsq))/(1.0+ec*cos(rlam-peri))
c
c   compute sun lag angle in radians (720*slag/pi = equation
c    of time in minutes)
      slag = elt - alp
c
c  compute hour angle of sunset at all latitudes
c
      do 1 i=1,n
        ph = alat(i)
        if (dlt .ne. 0.0) then
          aap=abs(ph)
c
c  test whether latitude is near either pole
c
          eps = abs(aap-hpi)
          if (eps .le. ccr) then
c
c  hour angle of sunset at pole is either zero or pi depending on dlt
c
             h = hpi*abs(aap/ph+abs(dlt)/dlt)
            ss = sin(ph) * sin(dlt)
            cc = 0.0
          else
            ss = sin(ph) * sin(dlt)
            cc = cos(ph) * cos(dlt)
            ar = -ss / cc
            ac = abs(ar)
            if ((ac - 1.0 + ccr) .eq. 0.0) h = (ac - ar) * hpi
            if ((ac - 1.0 + ccr) .gt. 0.0) then
              if (ar .ge. 0.0) then
                h = 0.0
              else
                h = pi
              endif
            endif
            if ((ac -1.0 + ccr) .lt. 0.0) h = acos(ar)
          endif
        else
          h  = hpi
          ss = 0.0
          cc = cos(ph)
        endif
c
c  below is normal case computation of the sun hour angle
c
        hang(i) = h
c
c  computation of cosz and frac added to solar for models without
c   diurnal variation since subroutine zenith need not be called.
c   cosz is constrained to be positive.
c
        frac(i) = h / pi
        if (h .ne. 0.0) then
          cosz(i) = amax1(ss + cc * sin(h) / h, 0.0)
        else
          cosz(i) = 0.0
        endif
    1 continue
c
      return
      end
#else
      subroutine solar(jld,fjd,r,dlt,alp,slag,n,alat,hang,cosz,frac)
c
c  subroutine solar is called by march and computes the radius
c    vector, the declination and right ascension of the sun, the
c    apparent sun lag angle (related to the equation of time), the hour
c    angle of the sun at sunset, the mean cosine of the sun zenith
c    angle and daylight fraction for n specified latitudes given the
c    day and fraction.  see notes of august 1989.
c
c    definition of the arguments used by this subroutine -
c
c      jld   =        julian day number (an integer)
c      fjd   =        julian date in floating point
c      r     =        radius vector (distance to sun in a. u.)
c      dlt   =        declination of sun
c      alp   =        right ascension of sun
c      slag  =        apparent sun lag angle (west of mean sun is plus)
c      n     =        number of latitudes
c      alat  =        latitude values dimensioned by n
c      hang  =        hour angle of sun at sunrise or sunset
c      cosz  =        mean cosine of zenith angle
c      frac  =        fraction of daylight
c
c     note - all angles are expressed in radians
c            input variables are day in /ctlblk/, n and alat
c
#include "parm.h"
c
      dimension alat(n), hang(n)
      dimension cosz(n), frac(n)
c
c
#include "ctlblk.h"
c
#include "param.h"
c
c  value of ccr depends on the precision of the machine
c
      data ccr / 1.0e-06 /
c
c   epoch for model days is epocm (january 1, 1970 at 0 hours ut
c    or julian date=2440587.5).  this is irrelevant for 365-day
c    years that start at year=0.
c
      data epocm / 2440587.5 /
c
c   epoch for model astronomical calculations is epoca.  this is
c    the standard epoch of january 1, 2000 at 12 hours ut or
c    jd=2451545.0 or j2000.0
c
      data epoca / 2451545.0 /
c
c   obliquity in degrees at epoch (j2000.0) is tilt
c
      data tilt / 23.439 /
c
c   longitude of earth perihelion at epoch (j2000.0) in degrees
c    is perep
c
      data perep / 102.932 /
c
      piph   = pi  + hpi
c
      dat = day
c
c   set fundamental astronomical parameters
c
c   eccentricity of orbit is ec (non-dimensional)
      ec   = 0.01671
      ecsq = ec**2
c
c   i365yr controls which year length is used in solar
c     i365yr = 1 ==> 365 day year, repeats year to year
c     i365yr = 0 ==> 365.25 day year, real time calendar
c
      i365yr = 1
c
c   tropical or calendar year in mean solar days is cyear.
c   anomalistic year in mean solar days (of 86,400 si seconds each)
c    is year.
      if (i365yr .ne. 1) then
c       use statement below for real length tropical
c        and anomalistic  years
        cyear = 365.242191
        year  = 365.259635
      else
c       use statement below for equal lengths for tropical
c        and anomalistic years
        cyear = 365.0
        year  = cyear
      endif
c
c   crate and arate below are mean daily angular motions of planet
c    in degrees with respect to equinox and perihelion respectively
      crate = 360.0 / cyear
      arate = 360.0 /  year
c
c   semimajor axis of planetary orbit in astronomical units is sma.
c    one a. u. = 149,600,000 km.  for earth sma=1.
      sma = 1.0
c
c   mean longitude of sun in degrees at epoch is orsun
      orsun = 280.46
c
c   statements below set julian date and day number
      fjd = dat + epocm
      jld = ifix(fjd)
c
c   days after epoch (j2000.0) is en.  the value of en is adjusted
c    by adding 7.5 days for models with 365-day years to normalize
c    equinox date
      en = fjd - epoca
      if (year .eq. cyear) then
         en  = en + 7.5
         enn = 0.0
       else
         enn = en
      endif
c
c  obliquity of earth axis in degrees is oblqty.  note that secular
c    change in obliquity is suppressed for 365-day years (enn=0).
      oblqty = tilt - 4.0e-7 * enn
      angin  = oblqty / radian
      sni    = sin(angin)
      er     = sqrt((1.0+ec) / (1.0-ec))
c
c   longitude of earth perihelion in degrees is peril.  note that
c    secular change in peril, computed in second term, is suppressed
c    in models with 365-day years because in that case crate=arate.
c
      peril = perep + (crate-arate)*en
c
c   longitude of sun at perihelion in degrees is sunper.  this is
c    exactly 180 degrees from that of earth at perihelion.
c
      sunper = peril - 180.0
      craten = crate * en
c
c   mean longitude of sun in radians is el
      el = (orsun+craten) / radian
c
c   constrain el to range of zero to 2 * pi
      el = el - twopi*ifix(el/twopi)
      if (el .lt. 0.0) el= el+twopi
      elt   = el
      orper = orsun - sunper
      if (orper .gt. 360.0) orper = orper - 360.0
      if (orper .lt. 0.0  ) orper = orper + 360.0
c
c   mean anomaly in radians is em (g in almanac equations)
      em = (orper+craten) / radian
c
c   constrain em to range of zero to 2 * pi
      em = em - twopi*ifix(em/twopi)
      if (em .lt. 0.0) em = em + twopi
      sinem = sin(em)
      cosem = cos(em)
c
c   eccentric anomaly in radians is ecan.  statement below is an
c    approximate solution to Kepler equation that is quite good
c    for small eccentricities (.lt. 0.1).
c
      ecan = em + ec*sinem + ecsq*cosem*sinem
c
c   true anomaly in radians is w
      w = 2.0*atan(er*tan(0.5*ecan))
c
c   ecliptic longitude of sun in radians is rlam
      rlam = w + sunper/radian
c
c   constrain rlam to range of zero to 2 * pi
      rlam = rlam - twopi*ifix(rlam/twopi)
      if (rlam .lt. 0.0) rlam = rlam + twopi
c
c   compute right ascension of sun (alp) in radians
      alp = atan(cos(angin)*tan(rlam))
c
c   following 6 statements determine proper quadrant for root of
c    equation for alp above and insure that the difference slag
c    below has an absolute magnitude of order pi/10.
c
      if (alp .lt. 0.0) alp = alp + twopi
c
      if ((elt-alp) .gt. piph) elt = elt - twopi
      if ((elt-alp) .lt.-piph) elt = elt + twopi
c
      if (abs(elt-alp) .gt. hpi) then
         if ((elt-alp) .gt. hpi) alp=alp+pi
         if ((elt-alp) .lt.-hpi) alp=alp-pi
      endif
c
c   compute declination of sun (dlt) in radians
      sind = sni*sin(rlam)
      dlt  = asin(sind)
c
c   compute radius vector (distance to sun in a.u.)
      r = sma*(1.0-ec*cos(ecan))
c
c   compute sun lag angle in radians (720*slag/pi = equation
c    of time in minutes)
      slag = elt - alp
c
c  compute hour angle of sunset at all latitudes
c
      do 1 i=1,n
        ph = alat(i)
        if (dlt .ne. 0.0) then
          aap=abs(ph)
c
c  test whether latitude is near either pole
c
          eps = abs(aap-hpi)
          if (eps .le. ccr) then
c
c  hour angle of sunset at pole is either zero or pi depending on dlt
c
             h = hpi*abs(aap/ph+abs(dlt)/dlt)
            ss = sin(ph) * sin(dlt)
            cc = 0.0
          else
            ss = sin(ph) * sin(dlt)
            cc = cos(ph) * cos(dlt)
            ar = -ss / cc
            ac = abs(ar)
            if ((ac - 1.0 + ccr) .eq. 0.0) h = (ac - ar) * hpi
            if ((ac - 1.0 + ccr) .gt. 0.0) then
              if (ar .ge. 0.0) then
                h = 0.0
              else
                h = pi
              endif
            endif
            if ((ac -1.0 + ccr) .lt. 0.0) h = acos(ar)
          endif
        else
          h  = hpi
          ss = 0.0
          cc = cos(ph)
        endif
c
c  below is normal case computation of the sun hour angle
c
        hang(i) = h
c
c  computation of cosz and frac added to solar for models without
c   diurnal variation since subroutine zenith need not be called.
c   cosz is constrained to be positive.
c
        frac(i) = h / pi
        if (h .ne. 0.0) then
          cosz(i) = amax1(ss + cc * sin(h) / h, 0.0)
        else
          cosz(i) = 0.0
        endif
    1 continue
c
      return
      end
#endif
