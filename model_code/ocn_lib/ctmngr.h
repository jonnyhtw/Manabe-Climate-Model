c====================== include file "ctmngr.h" ========================
#ifdef gtime
c
c-----------------------------------------------------------------------
c  for use with time manager written by keith dixon for use in 
c  the GFDL Climate Group model integrations.     ( KD:  NOV 1995)
c
c NOTE: KEY ASSUMPTIONS  (1) This code assumes a constant 365-day year
c                        calendar (i.e., no leap years) 
c                        (2) The tracer timestep length "dtts" must be
c                        a whole number and an even multiple of two.
c                        (3) Calculation of logical switches are 
c                        all referenced to initial time specified
c                        by the user in blkdta.F
c                        (4) Initial time is specified only as a year.
c                        Time stamps are made starting with
c                        1/ 1/XXXXX 0:0:0 where XXXXX is the
c                        user specified initial year "year0".
c-----------------------------------------------------------------------

c
c  ---------
c  Variables:
c  ---------
c
c   Integers:
c
c   the following 4 variables are written out on the restart tape
c     itt = current time step counter (since initial cond.)
c     icyear     = years since initial condtns (of current itt)
c     icdayofyr  = days of the year (of current itt)
c     icsecofday = current seconds of the day (of current itt)
c
c   the following 2 variables are set in blkdta.F
c     year0             = year of initial conditions       ( >= 0)
c     dayspermon(1..12) = number of days per month (set in blkdta)
c
c     sumdayspmon(0..12) = sum of the days of the year at the end
c         of the corresponding month
c         (sumdayspmon(0)=0, sumdayspmon(12)=365)
c
c   together with the 4 variables written out on the restart tape
c   the following variables describe the current time of the
c   model simulation (i.e., time passed since initialization)
c     icmonth      =  (1..12)
c     icdayofmonth =  (1..31)
c     ichourofday  =  (0..23)
c     icminofhour  =  (0..59)
c     icsecofmin   =  (0..59)
c
      integer year0
      integer dayspermon, sumdayspmon
c
c     the 4 variables in "tmngrir" are written to the "kontrl" restart
      common /tmngrir/ itt, icyear, icdayofyr, icsecofday
c
      common /tmngri / icmonth, icdayofmonth
      common /tmngri / ichourofday, icminofhour, icsecofmin
      common /tmngri / dayspermon(12), sumdayspmon(0:12), year0
c
#if defined cfctracers || defined ccl4tracer
c     tryr0  = integer variable of model year "icyear" corresponding
c              to 1 Jan 1910 (see "setvbc.F") [in years]
      integer tryr0
c     
      common /tmngri / tryr0
#endif
c
#if defined cpldleg2 || defined cpld_restore
c     icntts = counter for averaging of cumulative "terr" and "serr"
c
      common /tmngrix/ icntts
#endif
#if defined cpldleg3 || defined cpldleg2
c
c     month  = is redundant with "icmonth" (redundancy tolerated to
c              reduce need for code changes)
c     newmon = used as a flag to indicate time to shift monthly surface
c              boundary conditions  for time interpolation (0,1 where
c              1 = shift)
c
      common /tmngrix/ month, newmon, nextm
#endif
c
c   Characters:
c
c     stamp   = date & time stamp corresponding to itt 
c               example:  m/d/y= 7/21/ 1992 h:m:s=19:20: 0
c                         12345678901234567890123456789012
c
      character * 32 stamp
c
      common /tmngrc / stamp
c
#if defined cpldleg3 || defined cpldleg2
c
c     months = 4 character month designations for Leg2 and Leg3 of
c              climate group runs
      character * 4 months
c
      common /tmngrcx/ months(12)
#endif
c
c   Double Precision Real:
c
c     ctotdays = current total number of days since initial conditions
c              (referenced to initial conditions 1 Jan year 0)
c     cdayofcent = current number of days since start of the century
c     relyr   = relative year referenced to time at model
c               initial conditions
c
cSGI
cc    double precision ctotdays, cdayofcent, relyr
      real*16 ctotdays, cdayofcent, relyr
cSGI
c
      common /tmngrdr/ ctotdays, cdayofcent, relyr
c
c-----------------------------------------------------------------------
#else
c-----------------------------------------------------------------------
c
c     time manager clock parameters for setting time of model
c     initial conditions.
c
c     year0  = starting year   ( >= 0)
c     month0 = starting month  (1..12)
c     day0   = starting day    (1..last day of month)
c     hour0  = starting hour   (0..23)
c     minut0 = starting minute (0..59)
c     sec0   = starting second (0..59)
c
      integer year0, month0, day0, hour0, minut0, sec0
c
      common /ctmngr/ year0, month0, day0, hour0, minut0, sec0
c
c     time manager variables corresponding to current model time step
c
c     tyear   = year                     (ie: 1900)
c     tmonth  = month of the year        (ie: 1..12)
c     tday    = day of the month         (ie: 1..31)
c     thour   = hour of the day          (ie: 0..23)
c     tmin    = minute of the hour       (ie: 0..59)
c     tsec    = seconds after the minute (ie: 0..59)
c     relyr   = relative year referenced to time at model
c               initial conditions
c     yrnorm  = number of days in a normal year
c     yrleap  = number of days in a leap year
c     nleap   = number of years between leap years
c     time0   = time (in years) at initial conditions
c     totday  = accumulated time in days from jan 1, year 0
c     julian  = (true, false) = (real, simple) calendar
c     stamp   = date & time stamp corresponding to itt
c
      logical julian
      character*32 stamp
#ifdef worldocn
c     carry relyr in double precision for long runs
cSGI
cc    double precision relyr
      real*16 relyr
cSGI
# if defined cpldleg3 || defined cpldleg2
     $,                relyrp
# endif
      common /ctmngd/ relyr
# if defined cpldleg3 || defined cpldleg2
     $,               relyrp
# endif
c
      common /ctmngr/ tyear, tmonth, tday, thour, tmin, tsec
#else
      common /ctmngr/ tyear, tmonth, tday, thour, tmin, tsec, relyr
#endif
     $,               time0, totday, yrnorm, yrleap, nleap, julian
     $,               stamp
c
c     the following variables are used to update file "kontrl" every 
c     time step and written to a restart file at the end of the run.
c
c     itt     = current time step number
c     years   = number of years of model integration corresponding
c               to "itt".
#ifdef worldocn  
c               years is stored as single precision in file "kontrl"
c               and on the restart tape (relyr is double precision).
#endif
c
      common /ctmngr/ itt, years
c
#if defined cfctracers || defined ccl4tracer
c
c     tryr0  = variable "curtim" of model integrations corresponding 
c             to 1 Jan 1910 (see "setvbc.F") [in years]
c     tryr = tracer year ... year in the range 1910.0 to 1990.99 
c
      common /ctmngr/ tryr0, tryr, curtim
#endif
#ifdef worldocn
c
c     the following variables are put in common rather than kept
c     local in "tmngr" with a save statement (similar to an improvement
c     to be made in the next version of MOM)
c
c     dpm = number of days per month
c     itt1    = timestep number of at startup of the model run
c     totday1 = totday value corresponding to startup of the model run
c
      common /ctmngr/ rleap, tt0, totda0, clock0, monold
     $,               dpm(12), sumpm(12,2)
     $,               itt1, totday1
c
#endif
c
#ifdef cpldleg2
     $,               icntts
#endif
#if defined cpldleg3 || defined cpldleg2
     $,               month, monthp, newmon, nextm
c
      character*4 months(12)
c
      common /ctmngm/ months
c
#endif
c
#endif
