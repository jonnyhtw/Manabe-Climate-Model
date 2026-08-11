c=======================================================================
c
      subroutine gtmngr (dt)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
#ifdef gtime
c
c                  G-group  t-i-m-e  m-a-n-a-g-e-r
c
c     gtmngr measures model time using a 365 day calendar
c
c     in tests, date, time and setting of logical switches have been
c     shown accurate to within one second for periods of 25,000 years
c
c     the routines and function provided in this file allow for the
c     tracking of time in the model, expressing the time as a character
c     stamp, and the setting of logical switching at times desired by
c     the user.  in addition to setting switches at regular time
c     intervals, logical switches are set at the ends of days, months,
c     years, and decades.
c
c     ASSUMPTIONS:
c  1. Regular repeating 365 day calendar (no leap years).
c  2. Timestep length "dt" is a whole number and an even multiple of 2.
c  3. Logical switches may be set either (A) at "end of a" {day, month,
c     year, decade}", defined as occurring when within one-half
c     timestep of the end of the appropriate calandar indicator, --or--
c     (B) on the "last timestep of a" {day, month, year, decade}", 
c     defined as occurring when the next timestep will have a different
c     calendar indicator (*new feature*), --or-- (C) at a regular time 
c     interval (specified as a number of days: can include fractions of
c     days) as long as there are an integer number of time intervals
c     within a century (i.e., 36500 / (time interval) has no remainder).
c  4. Model is initialized to Jan 1 0:0:0 of user specified year "year0"
c     (i.e., can not initialize to 7:20PM of July 21st).
c  5. Initial time is specified only as a year.  Time stamps are made
c     starting with  1/ 1/XXXXX 0:0:0 (hr:min:sec) where XXXXX is the
c     user specified initial year "year0".
c
c
c     inputs:
c
c     itt    = current model time step (from switch.h).
c     dt     = length of time step in seconds.
c     see "switch.h" for other inputs
c
c     outputs:
c
c     itt    = incremented as itt + 1
c     see "switch.h" for other outputs
c
c   local variables:
c     halfstep = one half of timestep "dt" (integer in units of seconds
c                and the reason "dtts" must be an even multiple of two)
c     iayear, iadayofyr, iasecofday = time level indicators for a point
c                in time ahead of the current time level. 
c     ibyear, ibdayofyr, ibsecofday = time level indicators for a point
c                in time prior to the current time level. 
c     ittseg = counts number of timesteps executed thus far in this
c              segment of the model execution (since "first" was true)
c              (ittseg replaces use of tt0 since precision problems 
c              occur with subtraction on the 32-bit SGI when integers
c              exceed 2^24 or 16777216)
c
      integer halfstep
      logical tmnset
      logical eomset
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     written by k. dixon  (nov 1995)
c
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
      save ittseg, monold
c
c=======================================================================
c-----------------------------------------------------------------------
c     on the first timestep of a job initialize some variables and
c     perform a consistency check
c-----------------------------------------------------------------------
c
      if (first) then
c
c       is the timestep length an even mulitple of two seconds?
c
        halfstep = int(dt/2)
        restime = (dt/2.0) - float(halfstep)
        if (restime .NE. 0) then
          write(stdout,9101) dt
          write(stderr,9101) dt
          stop 'GTIME1'
        endif
c
        write (stdout,'(/,a52,/)')
     $  ' ==> gtmngr: 365 day repeating calendar is in effect.'
c
        sumdayspmon(0) = 0
        do mnth=1,12
          sumdayspmon(mnth) = sumdayspmon(mnth-1) + dayspermon(mnth)
        end do
c
        ittseg = 0
c
        call getmonthnday (icdayofyr, icmonth, icdayofmonth)
        monold = icmonth
c
      endif
c
c-----------------------------------------------------------------------
c     increment and update time step indicators and timestamp
c-----------------------------------------------------------------------
c
c     increment time step counter by one
c
      itt = itt + 1
      ittseg = ittseg + 1
c
c     advance current integer time indicators by "dt"
c
      intdt = int(dt)
      call addtime(intdt, icyear, icdayofyr, icsecofday)
      call getmonthnday (icdayofyr, icmonth, icdayofmonth)
      call gethrmins(icsecofday, ichourofday, icminofhour, icsecofmin)
c
c     the calculation for total number of days completed since
c     initialization ("ctotdays") is done in double precision for
c     increased long term accuracy
c
      ctotdays = dfloat(icyear)*365.0d0 + (dfloat(icdayofyr)-1.0d0)
     $         + (dfloat(icsecofday)/86400.0d0)
      ncentury = (icyear/100)
      nyrofcent= icyear - (100*ncentury)
      cdayofcent=  dfloat(nyrofcent)*365.0d0 + (dfloat(icdayofyr)-1.0d0)
     $         + (dfloat(icsecofday)/86400.0d0)
c
c     create the timestamp for the new time level (the user specified
c     initial year value "year0" is added in when making the stamp)
c
      call mkstmp (stamp, (year0+icyear), icmonth, icdayofmonth,
     $                    ichourofday, icminofhour, icsecofmin)
c
c     for consistency with older time manager, compute "relyr", the
c     floating point number of years since initialization -and-
c     for added precision use the century counters to compute "relyr"
c
cccc  relyr = ctotdays/365.0
      relyr = ncentury*100.0D0 + cdayofcent/365.0D0
c
c=======================================================================
c
c     set all time dependent logical switches in this section.
c
c=======================================================================

c-----------------------------------------------------------------------
c     determine if the current model timestep is the last one of a day,
c     month, year or decade, and set switches accordingly, and
c     determine if the model is within a half a timestep of the end
c     of a day, month, year or decade, and set switches accordingly
c-----------------------------------------------------------------------
c
c     compute what time levels will be a full timestep ahead
c     (the setting of "idtm1" as being one second less than "dt" 
c     ensures that the new month is not considered to begin until after
c     midnight of the last day of the current month)
c
c     set "lo---" (last of...) switches
c
      idtm1 = int(dt-1)
      loday    = .FALSE.
      lomonth  = .FALSE.
      loyear   = .FALSE.
      lodecade = .FALSE.
c
      if (icsecofday .EQ. 0 .OR. (icsecofday+idtm1) .GE. 86400) then
        loday = .TRUE.
c
        iayear     = icyear
        iadayofyr  = icdayofyr
        iasecofday = icsecofday
        iamonth      = icmonth
        iadayofmonth = icdayofmonth
c
        call addtime(idtm1, iayear, iadayofyr, iasecofday)
        call getmonthnday(iadayofyr, iamonth, iadayofmonth)
c
        if ((iayear .NE. icyear) .OR.
     $     (icsecofday .EQ. 0 .AND. icdayofyr .EQ. 1)) then
          loyear = .TRUE.
c
          if ((int(iayear/10) .NE. int(icyear/10)) .OR.
     $        (mod(icyear,10).EQ. 0 .AND. icsecofday .EQ. 0 .AND.
     $         icdayofyr .EQ. 1)) then
            lodecade = .TRUE.
          endif
        endif
c
        if ((iamonth .NE. icmonth) .OR. loyear .OR.
     $      (icsecofday .EQ. 0 .AND. icdayofmonth .EQ. 1)) then
          lomonth = .TRUE.
        endif
c
      endif
c
c-----------------------------------------------------------------------
c
c     set "eo---" (end of...) switches to true if within half
c     a timestep of a calendar change
c
      halfstep = int(dt/2)
      eoday    = .FALSE.
      eomonth  = .FALSE.
      eoyear   = .FALSE.
      eodecade = .FALSE.
c
      if ((icsecofday+halfstep) .GE. 86400 .OR.
     $    (icsecofday-halfstep) .LT. 0) then
        eoday = .TRUE.
c
        iayear     = icyear
        iadayofyr  = icdayofyr
        iasecofday = icsecofday
        iamonth      = icmonth
        iadayofmonth = icdayofmonth
        ibyear     = icyear
        ibdayofyr  = icdayofyr
        ibsecofday = icsecofday
        ibmonth      = icmonth
        ibdayofmonth = icdayofmonth
c
c       compute time levels half a timestep ahead
c
        call addtime( halfstep, iayear, iadayofyr, iasecofday)
        call getmonthnday (iadayofyr, iamonth, iadayofmonth)
c
c       compute time levels half a timestep behind
c
        call addtime(-halfstep, ibyear, ibdayofyr, ibsecofday)
        call getmonthnday (ibdayofyr, ibmonth, ibdayofmonth)
c
        if (iayear .NE. ibyear) then
          eoyear = .TRUE.
          if (int(iayear/10) .NE. int(ibyear/10)) then
            eodecade = .TRUE.
          endif
        endif
        if ((iamonth .NE. ibmonth) .OR. eoyear) then
          eomonth = .TRUE.
        endif
c
      endif
c
#if defined cpldleg3 || defined cpldleg2
c     determine if the next timestep will be in the next month
c     for consistency with previous G group use code, set endmon=lomonth
c
      endmon = lomonth
c
#endif
c-----------------------------------------------------------------------
c     is it the first timestep of a new month?
c-----------------------------------------------------------------------
c
      eom = eomset (monold, icmonth)
c
c-----------------------------------------------------------------------
c     is time mixing to be done at this time ?
c-----------------------------------------------------------------------
c
      if (mod(itt,nmix) .eq. 1) then
          mixts = .TRUE.
        else
          mixts = .FALSE.
      endif
c
c-----------------------------------------------------------------------
c     is it time to stop yet ?
c-----------------------------------------------------------------------
c
      dsofar = (ittseg+0.01)*dt/86400.0
      if (dsofar .ge. days) then
          last = .TRUE.
          write(stdout,9171) itt, days, ittseg
        else
          last = .FALSE.
      endif
c
c-----------------------------------------------------------------------
c     NOTE: calls to tmnset set logical flags at regular intervals.
c           the use of the "cdayofcent" rather than "ctotdays" improves
c           the accuracy, but assumes that the chosen interval divides
c           evenly into the number of days in a century.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     is diagnosis to be done at this time ?
c-----------------------------------------------------------------------
c
c     do an diagnostic timestep when within one half of a timestep
c     of the requested time interval
c
cccc  diagts = tmnset (dgnstc, ctotdays, dt)
      diagts = tmnset (dgnstc, cdayofcent, dt)
c
c-----------------------------------------------------------------------
c     is printing of time step information to be done at this time ?
c-----------------------------------------------------------------------
c
c     print timestep info when within one half of a timestep
c     of the requested time interval
c
cccc  prntsi = tmnset (tsi, ctotdays, dt)
      prntsi = tmnset (tsi, cdayofcent, dt)
c
c-----------------------------------------------------------------------
c     is time to save a snapshot of the data ?
c-----------------------------------------------------------------------
c
c     take a snapshot (save data) when within one half of a timestep
c     of the requested time interval
c
cccc  snapts = tmnset (snaps, ctotdays, dt)
      snapts = tmnset (snaps, cdayofcent, dt)
c
c-----------------------------------------------------------------------
c     are regional tracer averages to be done at this time ?
c-----------------------------------------------------------------------
c
c     do an average when within one half of a timestep
c     of the requested time interval
c
cccc  tavgts = tmnset (travg, ctotdays, dt)
      tavgts = tmnset (travg, cdayofcent, dt)
c
c-----------------------------------------------------------------------
c     additional logical switches for determining when to bring
c     in time dependent data should go in here.
c     (ie: monthly boundary conditions, assimilations, etc)
c     also logical switches for determining when to write out data.
c     (ie: averaging periods for data, etc)
c-----------------------------------------------------------------------
c
#ifdef kdtermb
c-----------------------------------------------------------------------
c     is time averaging of term balances to be done at this time ?
c-----------------------------------------------------------------------
c
      if (avterm .ne. 0.0) then
c         setting termav to be true at regular intervals...
cccc      termav = tmnset (avterm, ctotdays, dt)
          termav = tmnset (avterm, cdayofcent, dt)
          prterm = termav
c
      else
c       setting termav to be true at the end of each month
        if (lomonth) then
          termav = .TRUE.
# if defined cpldleg3 || defined cpldleg2
c   hardwired local modification:
c         print time averaged term balances in Feb and Aug
          prterm = .FALSE.
          if (mod(month,6) .eq. 2) prterm = .TRUE.
# endif
        else
          termav = .FALSE.
        endif
c
      endif
#endif
c
#if defined cpldleg3 || defined cpldleg2
c-----------------------------------------------------------------------
c     is diagnostic info to be printed this timestep?
c-----------------------------------------------------------------------
c
cccc  prdgts = tmnset (prdiag, ctotdays, dt)
      prdgts = tmnset (prdiag, cdayofcent, dt)
c
#endif
#ifdef worldocn
c-----------------------------------------------------------------------
c     is history data for analysis to be written out at this time ?
c-----------------------------------------------------------------------
c
cccc  histts = tmnset (histt, ctotdays, dt)
      histts = tmnset (histt, cdayofcent, dt)
c
# ifdef seasonal
c
c     determine month ahead & month behind & weightings
c     for surface boundary condition seasonal cycle used in "sfcbc.F"
c
c     "today" = is the number of days that have passed in the current
c               month to bring us to the current time
c     "daym"  = number of days in the current month
c     "mahd" & "mbhd"   = month ahead & month behind
c     "wmahd" & "wmbhd" = weighting factors for months ahead and behind
c     see "seasnl.h" for details
c
      today = (icdayofmonth-1.0) + (icsecofday/86400.0)
      daym = dayspermon(icmonth)
c
      tplb = dayspermon(mod(icmonth-1,12) + 1)
      tpla = dayspermon(mod(icmonth+1,12) + 1)
c
      call tmntrp (icmonth, 12, today, tplb,daym,tpla, mahd,mbhd, wmbhd)
c     
      wmahd = c1 - wmbhd
# ifdef oldtime
c     compute weights the same way as was done for Cyber runs
c     "mahd" & "mbhd"   = month ahead & month behind
c     "wmahd" & "wmbhd" = weighting factors for months ahead and behind
c     see "seasnl.h" for details
c
      ayrnorm = 365.0
      secmon = (ayrnorm * 86400.) / 12.
      ttsec = (ctotdays+1.0)/secday
      mnthfl = (ttsec + .5*secmon) / secmon
      wmahd = (ttsec + .5*secmon - mnthfl*secmon ) / secmon
      wmbhd = 1. - wmahd
      mbhd = mod(mnthfl+11,12) + 1
      mahd = mod(mnthfl,12) + 1
#  endif
# endif
# if defined cpldleg3 || defined cpldleg2
c   for consistency with the older version set some control block info
c   and the obsolete variable "month"
c   NOTE: "ib" preface here means time levels from one-half timestep
c   ago (i.e., ibmonth will be 1 on midnight Februay 1st)
c
      month = ibmonth
      konblk( 3) = itt
      totday = ctotdays
      konblk( 9) = totday
      konblk(21) = icdayofmonth
c   set month and year in fort.60 for restarting, etc (see run script)
      cmonth = months(icmonth)
      konblk(23) = icyear + year0
      konblk(65) = icmonth
      newmon = 0
      if (mahd .ne. nextm) newmon = 1
      nextm = mahd
      if (newmon .eq. 1) then
        write(stdout,1919) itt, mahd, mbhd
1919    format('   newmon = 1 at itt=',i12,'  mahd, mbhd=',2i5)
      endif
# endif
#endif
c
c=======================================================================
c

      return
9101  format(//' ##### ERROR ##### ##### ERROR ##### ##### ERROR #####',
     $        /' ##### in gtmngr.F the timestep dt must be a'
     $        /' ##### whole number and an even multiple of two.'
     $        /' ##### here dt is= ',e17.10,
     $        /' ##### ERROR ##### ##### ERROR ##### ##### ERROR #####')
9171  format('   LAST set to TRUE on itt =',i12,
     $ ': MOM ran for the prescribed ',f10.2,' DAYS in',
     $ i10,' timesteps')
      end
c
c - - - - - - - -
c - - - - - - - -
            function tmnset (dayintvl, totday, dt)
c
c=======================================================================
c     tmnset = true if current time is within 1/2 time step of desired
c              interval
c     dayintvl = time interval between when the logical will be set to
c              true (referenced to jan 1 of initial year "year0")
c     totday = accumulated time in days since reference conditions
c              (relative to jan 1 of reference year)
c     dt     = time step in seconds
c=======================================================================
c
      logical tmnset
c
cSGI
cc    double precision totday, xnum, dblintvl, close, remain
      real*16 totday, xnum, dblintvl, close, remain
cSGI
c
cSGI chged dfloat to dble for reals
      close  = 0.5d0*dble(dt)
      dblintvl = dble(dayintvl)
c
      xnum = totday/dblintvl
      i    = nint (xnum)
      remain = (totday - dble(i*dblintvl))*86400.0d0
c
      if ((remain .GE. 0.0 .AND.  remain .LT. close) .or.
     $    (remain .LT. 0.0 .AND. -remain .LE. close)) then
        tmnset = .TRUE.
      else
        tmnset = .FALSE.
      endif
c
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine mkstmp (stamp, year, month, day, hour, min, sec)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     make a 32 character time stamp from integer inputs
c-----------------------------------------------------------------------
c
c     inputs: integer year, month, day, hour, min, sec
c     output: stamp of form m/d/y h:m:s  (works up to year 99,999)
c
      integer year, month, day, hour, min, sec
      character*32 stamp
c
c-----------------------------------------------------------------------
c
      write (stamp,'(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &             'm/d/y=',month,'/',day,'/'
     &,             year,' h:m:s=', hour,':'
     &,             min,':', sec
c
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine rdstmp (stmp, year, month, day, hour, min, sec)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     convert 32 character time stamp into integer year, month, day,
c     hour, minute, second
c-----------------------------------------------------------------------
c
c     input:   32 character time stamp of form m/d/y h:m:s
c     outputs: integer year, month, day, hour, min, sec
c
      character *32 stmp
      character *6  junk6
      character *8  junk8
      character *1  junk1
      integer       year, month, day, hour, min, sec
c
c-----------------------------------------------------------------------
c
      read  (stmp, '(a6,i2,a1,i2,a1,i5,a7,i2,a1,i2,a1,i2)')
     &       junk6, month, junk1, day, junk1, year, junk8, hour,
     &       junk1, min, junk1, sec
c
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine getmonthnday(idayofyr, imonth, idayofmon)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     convert integer day of the year to integer month number and
c     integer day of the month  (reverse of getdayofyr)
c-----------------------------------------------------------------------
#include "ctmngr.h"
c
c     input:   integer day of the year "idayofyr" (1..365) 
c     outputs: integer month (1..12) and day of the month  (1..31)
c
      integer imonth, idayofmon, idayofyr
c
c-----------------------------------------------------------------------
c
      if (idayofyr .LT. 1 .OR. idayofyr .GT. 365) then
        print *, ' ### ### ### ### ### ### ### ### ### '
        print *, ' ### PROBLEM in subroutine getmonthnday'
        print *, ' ### invalid idayofyr = ',idayofyr
        print *, ' ### idayofyr can not be < 1 nor > 365'
        print *, ' ### ### ### ### ### ### ### ### ### '
        stop 'getmon'
      endif
c
      do monnum=1,12
        if (idayofyr .LE. sumdayspmon(monnum)) then
          imonth = monnum
          idayofmon = idayofyr - sumdayspmon(monnum-1) 
          return
         endif
      end do
c     
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine getdayofyr(imonth, idayofmon, idayofyr)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     convert month number and day of the month to day of the year
c     (reverse of getmonthnday)
c-----------------------------------------------------------------------
c
c     inputs: integer month (1..12) and day of the month  (1..31)
c     output: integer day of the year "idayofyr" (1..365) 
c
#include "ctmngr.h"
      integer imonth, idayofmon, idayofyr
c
c-----------------------------------------------------------------------
c
      if (imonth .LT. 1 .OR. imonth .GT. 12) then
        print *, ' ### ### ### ### ### ### ### ### ### '
        print *, ' ### PROBLEM in subroutine getdayofyr'
        print *, ' ### invalid imonth = ',idayofyr
        print *, ' ### imonth can not be < 1 nor > 12'
        print *, ' ### ### ### ### ### ### ### ### ### '
        stop 'getday1'
      endif
      intdpmon = int(dayspermon(imonth))
      if (idayofmon .LT. 1 .OR. idayofmon .GT. intdpmon) then
        print *, ' ### ### ### ### ### ### ### ### ### '
        print *, ' ### PROBLEM in subroutine getdayofyr'
        print *, ' ### invalid idayofmon = ',idayofmon
        print *, ' ### idayofmon can not be < 1 nor > ',intdpmon
        print *, ' ### for month ',imonth
        print *, ' ### ### ### ### ### ### ### ### ### '
        stop 'getday2'
      endif
c
      idayofyr = idayofmon + dayspermon(imonth-1)
c     
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine gethrmins(isecofday, ihourofday, iminofhour, isecofmin)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     convert seconds of the day (seconds since midnight) to
c     hour of the day, minute of the hour and seconds of the minute
c-----------------------------------------------------------------------
c
c     input:   integer seconds of the day isecofday (0..86399)
c     outputs: integers representing hours, minutes and seconds of stamp
c
      integer isecofday, ihourofday, iminofhour, isecofmin
c
c-----------------------------------------------------------------------
c
      if (isecofday .LT. 0 .OR. isecofday .GT. (86400-1)) then
        print *, ' ### ### ### ### ### ### ### ### ### '
        print *, ' ### PROBLEM in subroutine gethrmins'
        print *, ' ### invalid isecofday = ',isecofday
        print *, ' ### isecofday can not be < 0 nor > 86399'
        print *, ' ### ### ### ### ### ### ### ### ### '
        stop 'getmon'
      endif
c
      ihourofday = int(isecofday/3600)
      iminofhour = int((isecofday-(ihourofday*3600))/60)
      isecofmin = (isecofday-ihourofday*3600)-(iminofhour*60)
c
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c=======================================================================
c
      subroutine addtime(ideltat, iyear, idayoyr, isecoday)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     increment integer time keepers of year, day of year and seconds
c     of the day by the timestep "ideltat".  if "ideltat" < 0 then in
c     effect subtract time from the integer time keepers.
c-----------------------------------------------------------------------
c
c     inputs: 
c       ideltat = timestep length in sec (either positive or negative)
c       iyear, idayoyr, isecoday = integer time keepers of year, day of
c                year and seconds of the day
c     outputs: 
c       iyear, idayoyr, isecoday = updated integer time keepers of year,
c                day of year and seconds of the day
c
c-----------------------------------------------------------------------
c
c     1st, add the timestep length to the seconds of the day counter
c
      isecoday = isecoday + ideltat
c
c     2nd, if the updated seconds of the day counter exceeds the number
c     of seconds in a day, figure how many days the day number should
c     be incremented by and then adjust the seconds in a day value.
c     if needed, increment the year counter as well
c
      if (isecoday .GE. 86400) then
        nincday = int(isecoday/86400)
        isecoday = isecoday - nincday*86400
        idayoyr = idayoyr + nincday
c
        if (idayoyr .GT. 365) then
          nincyr = int(idayoyr/365)
          idayoyr = idayoyr - nincyr*365
          iyear = iyear + nincyr
        endif
c
      endif
c
c     3nd, if the updated seconds of the day counter is less than zero
c     figure how many days the day number should be reduced and then
c     adjust the seconds in a day value.  if needed, adjust the year
c     counter also.
c
      if (isecoday .LT. 0) then
        nincday = 1 + int(-isecoday/86400)
        isecoday = nincday*86400 + isecoday
        idayoyr = idayoyr - nincday
c
        if (idayoyr .LT. 1) then
          nincyr = 1 + int(-idayoyr/365)
          idayoyr = 365*nincyr + idayoyr
          iyear = iyear - nincyr
        endif

      endif
c
      return
      end
c
c - - - - - - - -
c - - - - - - - -
c
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
c
c - - - - - - - -
c - - - - - - - -
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
c
      if (month .NE. oldmon) then
        eomset = .TRUE.
      else
        eomset = .FALSE.
      endif
      oldmon = month
c
#endif
      return
      end

