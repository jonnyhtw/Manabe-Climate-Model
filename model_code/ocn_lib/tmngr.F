      subroutine tmngr (dt)
      implicit real*8 (a-h,o-z)
#ifndef gtime
c
c
c=======================================================================
c
c                      t-i-m-e      m-a-n-a-g-e-r
c
c     tmngr measures model time using a julian calendar.
c     date & time are accurate to within one second for arbitrary
c     length time steps. if a calendar with 360 days per year
c     and 30 day months is desired, set "julian" = false (see below).
c     note: the year length is arbitrary. it can be changed below.
c     computer requirements: 64 bit precision is recommended.
c
c     inputs:
c
c     itt    = current model time step (from switch.h).
c     dt     = length of time step in seconds.
c     see "switch.h" for other inputs
c
c     outputs:
c
c     itt    = itt + 1
c     see "switch.h" for other outputs
c
c=======================================================================
c
c
      logical tmnset, eomset
#if defined cpldleg3 || defined cpldleg2
# include "conblk.h"
      character*4 cmonth
      equivalence (cmonth,konblk(22))
#endif
#include "pconst.h"
#include "ctmngr.h"
#include "switch.h"
#ifdef worldocn
# ifdef seasonal
#include "seasnl.h"
# endif
#endif
c
#ifndef worldocn
c     the following variables are put in common "ctmngr" and the
c     data statement in "blkdta" in worldocn cases
c
      save dpm, sumpm, rnorm, rleap, tt0, totda0, clock0, monold
c
c     dpm = number of days per month
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0,
     $          31.0, 30.0, 31.0/
c
#else
c     treat some variables as double precision for added accuracy
cSGI
cc    double precision deltat, dtnorm, dtleap, f1, f2
      real*16 deltat, dtnorm, dtleap, f1, f2
cSGI
#endif
c
c-----------------------------------------------------------------------
c     statement function
c-----------------------------------------------------------------------
c
      frac(r) = r - int(r)
c
c
c-----------------------------------------------------------------------
c     increment time step counter
c-----------------------------------------------------------------------
c
      itt = itt + 1
c
c=======================================================================
c
c     calendar initialization
c
c=======================================================================
c
c
c-----------------------------------------------------------------------
c     set year lengths, partial sums & other misc quantities
c     on the first call to tmngr
c-----------------------------------------------------------------------
c
      if (first) then
#ifdef worldocn
c       store timestep counter for start of this run as "itt1"
        itt1 = itt - 1
#endif
c
c       set initial time in days
c
        clock0 = hour0*r24 + minut0*r1440 + (sec0+p5)*secday
c
c       do bounds checking on initial clock parameters in "ctmngr.h"
c
        if (month0 .lt. 1 .or. month0 .gt. 12) then
          write (stdout,*) ' "month0" is out of bounds'
          stop '=>tmngr'
        endif
        if (day0 .lt. 1 .or. day0 .gt. 31) then
          write (stdout,*) ' "day0" is out of bounds'
          stop '=>tmngr'
        endif
        if (hour0 .lt. 0 .or. hour0 .gt. 23) then
          write (stdout,*) ' "hour0" is out of bounds'
          stop '=>tmngr'
        endif
        if (minut0 .lt. 0 .or. minut0 .gt. 59) then
          write (stdout,*) ' "minut0" is out of bounds'
          stop '=>tmngr'
        endif
        if (sec0 .lt. 0 .or. sec0 .gt. 59) then
          write (stdout,*) ' "sec0" is out of bounds'
          stop '=>tmngr'
        endif
        if (julian) then
c
c         set up a julian calendar
c
          yrnorm = 365.0
#ifdef worldocn
          yrleap = yrnorm
          nleap  = 999999
# if defined cpldleg3 || defined cpldleg2
          if (pintro) then
# endif
          write (stdout,'(/,a48,/)')
     $          ' ==> tmngr: 365 day julian calendar is in effect.'
# if defined cpldleg3 || defined cpldleg2
          endif
# endif
#else
          yrleap = 366.0
          nleap  = 4
          write (stdout,'(/,a41,/)')
     $          ' ==> tmngr: julian calendar is in effect.'
#endif
        else
c
c         set up a simple calendar with 12 equi-length months
c         and no leap years (each year is "yrnorm" days long)
c
          yrnorm = 360.0
          yrleap = yrnorm
          nleap  = 999999
          do 100 m=1,12
            dpm(m) = yrnorm/12
100       continue
          write (stdout,'(/,a47,f9.3,a14/)')
     $     ' ==> tmngr: simple calendar is in effect using ', yrnorm
     $,    ' days per year'
        endif
c
        totday = c0
        tt0    = itt - 1.1
        rleap  = c1/yrleap
        rnorm  = c1/yrnorm
c
        sumpm(1,1) = c0
        do 200 m=2,12
          sumpm(m,1) = sumpm(m-1,1) + dpm(m-1)
          sumpm(m,2) = sumpm(m,1)
          if (m .gt. 2 .and. julian) sumpm(m,2) = sumpm(m,2) + c1
200     continue
c
c-----------------------------------------------------------------------
c       calculate time in years (time1) of model start.
c-----------------------------------------------------------------------
c
        styear = year0 + relyr
        if (mod(int(styear),nleap) .eq. 0) then
          time1 = styear + (sumpm(month0,2) + (day0-1) + clock0)/yrleap
        else
          time1 = styear + (sumpm(month0,1) + (day0-1) + clock0)/yrnorm
        endif
c
        if (julian) then
c
c         1st calc # of days ignoring leap days since 1 jan year 0
c
          totda0 = yrnorm*styear
c
c         add in leap days in years prior to the one listed in styear
c
          lpdays = int(time1 - 29.0/yrleap)/nleap
          if (mod(int(styear),nleap) .eq. 0 .and. (time1-int(time1))
     $        .gt. 29.0/yrleap) lpdays = lpdays - 1
c
c         add in leap day of styear if starting after the end of Feb,
c         year 0
c
          if (time1 .gt. (29.0/yrleap)) lpdays = lpdays + 1
c
c         add in number of days gone in the current year "thisdy"
c
          if (mod(int(styear),nleap) .eq. 0) then
            thisdy = (time1-styear)*yrleap
          else
            thisdy = (time1-styear) * yrnorm
          endif
#ifdef worldocn
          lpdays = 0
#endif
          totda0 = totda0 + lpdays + thisdy
        else
          totda0 = yrnorm*styear
        endif
        totday = totda0 + totday
#ifdef worldocn
c       store off totday1 corresponding to initial conditions
        totday1= totday
#endif
c
c-----------------------------------------------------------------------
c       set time of initial conditions
c-----------------------------------------------------------------------
c
        time0 = time1 - relyr
      endif
c
c
c
c=======================================================================
c
c     calendar calculation
c
c=======================================================================
c
c     increment day counter
c
#ifdef worldocn
      totday = totday1 + (itt-itt1)*dt*secday
#else
      totday = totday + dt*secday
#endif
c
#ifdef testyr
      write(stdout,9026) itt, lpdays, relyr, thisdy, totday1, totday
9026  format(/,1x,'in tmngr: itt lpdays relyr thisdy totday1 totday=',/,
     $         10x,i8,i5,4e19.10)
#endif
c
c-----------------------------------------------------------------------
c     increment relative year according to type of year & time within
c     the year.
c     if clearly within a normal year ... use dtnorm
c     if clearly within a leap year   ... use dtleap
c     if transitioning from a leap year to a normal year   ... use a
c     linear combination of dtleap and dtnorm
c     if transitioning from a normal year to a leap year   ... use a
c     linear combination of dtleap and dtnorm also
c-----------------------------------------------------------------------
c
      dtleap = dt*secday/yrleap
      dtnorm = dt*secday/yrnorm
      curtim = time0 + relyr
      if (mod(int(curtim),nleap) .eq. 0) then
        if (mod(int(curtim+dtleap),nleap) .eq. 0) then
          deltat = dtleap
        else
          f1 =nint(curtim) - curtim
          f2 = dtleap - f1
          deltat = (f1*dtleap + f2*dtnorm)/dtleap
        endif
      else
        if (mod(int(curtim+dtnorm),nleap) .ne. 0) then
          deltat = dtnorm
        else
          f1 = nint(curtim) - curtim
          f2 = dtnorm - f1
          deltat = (f1*dtnorm + f2*dtleap)/dtnorm
        endif
      endif
      relyr = relyr + deltat
#if defined cpldleg3 || defined cpldleg2
c   compute year of the next timestep
      relyrp = relyr + deltat
#endif
c
c-----------------------------------------------------------------------
c     construct julian (or simple) calendar year.
c-----------------------------------------------------------------------
c
      tyear = time0 + relyr
#if defined cpldleg3 || defined cpldleg2
      tyearp = time0 + relyrp
#endif
c
c-----------------------------------------------------------------------
c     decide whether this is a leap year or normal year.
c-----------------------------------------------------------------------
c
      if (mod(int(tyear),nleap) .eq. 0) then
        ntype  = 2
        tday = frac(tyear)*yrleap
#if defined cpldleg3 || defined cpldleg2
        tdayp = frac(tyearp)*yrleap
#endif
      else
        ntype  = 1
        tday = frac(tyear)*yrnorm
#if defined cpldleg3 || defined cpldleg2
        tdayp = frac(tyearp)*yrnorm
#endif
      endif
c
c-----------------------------------------------------------------------
c     construct month, day, hour, minute, & second using fractional
c     part of julian (or simple) year.
c-----------------------------------------------------------------------
c
      month = indp (tday, sumpm(1,ntype), 12)
      if (sumpm(month,ntype) .gt. tday) month = month-1
#if defined cpldleg3 || defined cpldleg2
      monthp = indp (tdayp, sumpm(1,ntype), 12)
      if (sumpm(monthp,ntype) .gt. tdayp) monthp = monthp-1
#endif
      if (first) then
        monold = month
      endif
      tmonth = month
#ifdef worldocn
      curtim = time0 + relyr
# ifdef seasonal
      today  = tday - sumpm(month,ntype)
# endif
#endif
      tday   = tday - sumpm(month,ntype) + 1
      thour  = frac(tday)*c24
      tmin   = frac(thour)*c60
      tsec   = frac(tmin)*c60
      tyear  = int(tyear)
      tday   = int(tday)
      thour  = int(thour)
      tmin   = int(tmin)
      tsec   = int(tsec)
#if defined cfctracers || defined ccl4tracer
c
c     tryr = year in the range 1932.0 to 1990.99
c     tryr0  = "variable "curtim" of model integrations corresponding
c              to 1 Jan 1910 (see "ctmngr.h")
c
      tryr = curtim - tryr0 + 1910.0
#endif
      write (stamp,'(a6,i2,a1,i2,a1,i4,a8,i2,a1,i2,a1,i2)')
     $             'm/d/y=',nint(tmonth),'/',nint(tday),'/'
#if defined cfctracers || defined ccl4tracer
     $,             nint(tryr),', h:m:s:', nint(thour),':'
#else
     $,             nint(tyear),', h:m:s:', nint(thour),':'
#endif
     $,             nint(tmin),':', nint(tsec)
c
c=======================================================================
c
c     set all time dependent logical switches in this section.
c
c=======================================================================
c
c
c-----------------------------------------------------------------------
c     is diagnosis to be done at this time ?
c-----------------------------------------------------------------------
c
c     do an diagnostic timestep when within one half of a timestep
c     of the requested time interval
c
      diagts = tmnset (dgnstc, totday, dt)
c
#if defined cpldleg3 || defined cpldleg2
c-----------------------------------------------------------------------
c     is diagnostic info to be printed this timestep?
c-----------------------------------------------------------------------
c
      prdgts = tmnset (prdiag, totday, dt)
#endif
c-----------------------------------------------------------------------
c     is time to save a snapshot of the data ?
c-----------------------------------------------------------------------
c
c     take a snapshot (save data) when within one half of a timestep
c     of the requested time interval
c
      snapts = tmnset (snaps, totday, dt)
c
c
c-----------------------------------------------------------------------
c     are regional tracer averages to be done at this time ?
c-----------------------------------------------------------------------
c
c     do an average when within one half of a timestep
c     of the requested time interval
c
      tavgts = tmnset (travg, totday, dt)
c
c-----------------------------------------------------------------------
c     is it the end of the month ?
c-----------------------------------------------------------------------
c
      eom = eomset (monold, month)
c
#if defined cpldleg3 || defined cpldleg2
      endmon = .false.                
      if ( monthp .ne. month ) endmon = .true.
#endif
c
#ifdef kdtermb
c-----------------------------------------------------------------------
c     is time averaging of term balances to be done at this time ?
c-----------------------------------------------------------------------
c
      if (avterm .ne. 0.) then
cccc      setting termav to be true at regular intervals...
          termav = tmnset (avterm, totday, dt)
          prterm = termav
c
      else
c       setting termav to be true at the end of each month
        if (eom) then
          termav = .TRUE.
# if defined cpldleg3 || defined cpldleg2
c         print time averaged term balances in Feb and Aug
          prterm = .false.
          if (mod(month,6) .eq. 2) prterm = .true.
# endif
        else
          termav = .FALSE.
        endif
c
      endif
#endif
c
c-----------------------------------------------------------------------
c     is time mixing to be done at this time ?
c-----------------------------------------------------------------------
c
      if (mod(itt,nmix) .eq. 1) then
          mixts = .true.
        else
          mixts = .false.
      endif
c
c-----------------------------------------------------------------------
c     is printing of time step information to be done at this time ?
c-----------------------------------------------------------------------
c
c     print timestep info when within one half of a timestep
c     of the requested time interval
c
      prntsi = tmnset (tsi, totday, dt)
c
#ifdef worldocn
c-----------------------------------------------------------------------
c     is history data for analysis to be written out at this time ?
c-----------------------------------------------------------------------
c
      histts = tmnset (histt, totday, dt)
c
#endif
c-----------------------------------------------------------------------
c     is it time to stop yet ?
c-----------------------------------------------------------------------
c
      if ((itt-tt0)*dt*secday .ge. days) then
          last = .true.
        else
          last = .false.
      endif
c
c-----------------------------------------------------------------------
c     additional logical switches for determining when to bring
c     in time dependent data should go in here.
c     (ie: monthly boundary conditions, assimilations, etc)
c     also logical switches for determining when to write out data.
c     (ie: averaging periods for data, etc)
c-----------------------------------------------------------------------
#ifdef worldocn
# ifdef seasonal
c
c     determine month ahead & month behind & weightings
c     for surface boundary condition seasonal cycle used in "sfcbc.F"
c
#  ifdef oldtime
c     compute weights the same way as was done for Cyber runs
c     "mahd" & "mbhd"   = month ahead & month behind
c     "wmahd" & "wmbhd" = weighting factors for months ahead and behind
c     see "seasnl.h" for details
c
      secmon = (yrnorm * 86400.) / 12.
      ttsec = totday/secday
      mnthfl = (ttsec + .5*secmon) / secmon
      wmahd = (ttsec + .5*secmon - mnthfl*secmon ) / secmon
      wmbhd = 1. - wmahd
      mbhd = mod(mnthfl+11,12) + 1
      mahd = mod(mnthfl,12) + 1
#if defined cpldleg3 || defined cpldleg2
      newmon = 0
      if (mahd .ne. nextm) newmon = 1
      nextm = mahd
#endif
#ifdef january
c     hardwire perpetual January 
      mahd = 1
      mbhd = 1
      wmahd = 1.0
      wmbhd = 0.0
#endif
c
#  else
c     "today" = computed above and is the number of days that have
c     passed in the current month to bring us to the current time
c     "daym"  = number of days in the current month
c     "mahd" & "mbhd"   = month ahead & month behind
c     "wmahd" & "wmbhd" = weighting factors for months ahead and behind
c     see "seasnl.h" for details
c
      daym = dpm(month)
      if (mod(int(tyear),nleap).eq.0 .and. month.eq.2) daym = daym + c1
c
      tplb = dpm(mod(month-1,12) + 1)
      tpla = dpm(mod(month+1,12) + 1)
      if (mod(int(tyear),nleap) .eq. 0) then
        if( (mod(month-1,12)+1) .eq. 2) tplb = tplb + c1
        if( (mod(month+1,12)+1) .eq. 2) tpla = tpla + c1
      endif
      call tmntrp (month, 12, today, tplb,daym,tpla, mahd,mbhd, wmbhd)
c      
cc    call tmntr2 (month, 12, today, daym, mahd, mbhd, wmbhd)
c
      wmahd = c1 - wmbhd
#  endif
# endif
# if defined cpldleg3 || defined cpldleg2
c   set some control block info
      konblk( 3) = itt
      konblk( 9) = totday
      konblk(21) = tday
c   set month and year in fort.60 for restarting, etc (see run script)
      cmonth = months(int(tmonth))
      konblk(23) = tyear
      konblk(65) = month
# endif
#endif
c
      return
      end
      function tmnset (days, totday, dt)
c
c=======================================================================
c     tmnset = true if current time is within 1/2 time step of desired
c              interval
c     days   = time interval between when the logical will be set to
c              true (referenced to jan 1, year 0)
c     totday = accumulated time in days from jan 1, year 0
c     dt     = time step in seconds
c=======================================================================
c
#include "pconst.h"
      logical tmnset
      close  = p5*dt
      i      = nint (totday/days)
      remain = (totday - i*days)*c60*c1440
      if ((remain .ge. c0 .and. remain .lt. close) .or.
     $    (remain .lt. c0 .and. -remain .le. close)) then
        tmnset = .true.
      else
        tmnset = .false.
      endif
      return
      end
      function eomset (oldmon, month)
c
c=======================================================================
c     eomset = true for the first time step when the month changes
c              otherwise  false.
c     oldmon = month when "first" is true
c     month  = current month
c=======================================================================
c
      integer oldmon
      logical eomset
      if (month .ne. oldmon) then
        eomset = .true.
      else
        eomset = .false.
      endif
      oldmon = month
      return
      end
      subroutine tmntrp
     $             (num, maxnum, pos, tplb, tpl, tpla, iahd, ibhd, fbhd)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     time interpolator ... constructs indices & factor needed for
c     linearly interpolating data defined at arbitrary time intervals
c     (midpoints of years, months, days or even random intervals) to
c     the time corresponding to the current time step
c     this routine provides for a smooth, linear interpolation between
c     the midpoints of time intervals of any arbitrary length
c
c=======================================================================
c
c     inputs:
c
c     num    = index of the time interval (containing data for use in
c              interpolating) which includes the current time.
c              (ie: if data is monthly then... num is 2 if the current
c               time is in february)
c     maxnum = max number of time intervals
c              (ie: 12 if using monthly data)
c     pos    = time measured relative to the beginning of interval "num"
c              and in the same units as "tpl".
c              (ie: 14.36 days into the month of february)
c     tplb   = time period length of interval immediately before "num"
c              (ie: 31.0 days if the month is january)
c     tpl    = time period length of interval "num"
c              (ie: 28.0 days if the month is february)
c     tpla   = time period length of interval immediately after "num"
c              (ie: 31.0 days if the month is march)
c
c     outputs:
c
c     iahd   = index for pointing to the time interval with midpoint
c              after (ahead of) the current time.
c              (ie: iahd would be 3 if the current time was beyond the
c                   midpoint of {but still within} february)
c     ibhd   = index for pointing to the time interval with midpoint
c              before (behind) the current time.
c              (ie: ibhd would be 2 if the current time was beyond the
c                   midpoint of {but still within} february)
c     fbhd   = interpolation factor for defining data at the current
c              time (the weight attributed to the data from the ibhd
c              time interval)
c                  schematically the interpolation is defined by:
c              interpolated data = fbhd*data(ibhd) + (1-fbhd)*data(iahd)
c
c     variable:
c
c     rmidpt = the midpoint of the current time interval
c
#include "pconst.h"
c
      rmidpt =  p5*tpl
      if (pos .gt. rmidpt) then
        iahd = num + 1
        fbhd = c1 - c2*(pos-rmidpt)/(tpl+tpla)
      else
        iahd = num
        fbhd = c2*(rmidpt-pos)/(tpl+tplb)
      endif
      iahd = mod(iahd-1,maxnum) + 1
      ibhd = mod(iahd-2+maxnum,maxnum) + 1
c
      return
      end
      subroutine tmntr2
     $             (num, maxnum, pos, tpl, iahd, ibhd, fbhd)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     time interpolator ... constructs indices & factor needed for
c     interpolating data defined at arbitrary time intervals (midpoints
c     of years, months, days or even random intervals) to the time
c     corresponding to the current time step.
c     this routine conducts a linear interpolation between the midpoints
c     of time intervals and the time interval boundaries.
c     at the time interval boundaries the interpolated data values are
c     defined to be the average of the mean values attributed to
c     the two time intervals that border it, regardless of their
c     relative lengths
c
c=======================================================================
c
c     inputs:
c
c     num    = index of the time interval (containing data for use in
c              interpolating) which includes the current time.
c              (ie: if data is monthly then... num is 2 if the current
c               time is in february)
c     maxnum = max number of time intervals
c              (ie: 12 if using monthly data)
c     pos    = time measured relative to the beginning of interval "num"
c              and in the same units as "tpl".
c              (ie: 14.27 days into the month of february would be the
c              designation for 6:28:48 on the morning of february 15th)
c     tpl    = time period length of interval "num"
c              (ie: 28.0 days if the month is february)
c
c     outputs:
c
c     iahd   = index for pointing to the time interval with midpoint
c              after (ahead of) the current time.
c              (ie: iahd would be 3 if the current time was beyond the
c                   midpoint of {but still within} february)
c     ibhd   = index for pointing to the time interval with midpoint
c              before (behind) the current time.
c              (ie: ibhd would be 2 if the current time was beyond the
c                   midpoint of {but still within} february)
c     fbhd   = interpolation factor for defining data at the current
c              time (the weight attributed to the data from the ibhd
c              time interval)
c                  schematically the interpolation is defined by:
c              interpolated data = fbhd*data(ibhd) + (1-fbhd)*data(iahd)
c
c     variable:
c
c     rmidpt = the midpoint of the current time interval
c
c
#include "pconst.h"
c
      rmidpt =  p5*tpl
      if (pos .gt. rmidpt) then
        iahd = mod(num, maxnum) + 1
        ibhd = num
        fbhd = c1 - abs((pos-rmidpt)/tpl)
      else
        iahd = num
        ibhd = mod(num-2+maxnum,maxnum) + 1
        fbhd = abs((pos-rmidpt)/tpl)
      endif
c
#endif
      return
      end
