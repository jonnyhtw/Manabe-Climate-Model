      subroutine restio
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     read and write restart file
c
c=======================================================================
c
#include "param.h"
#include "ctmngr.h"
#include "emode.h"
#include "iounit.h"
#include "levind.h"
#ifdef multitasking
#include "cshrbf.h"
#else
#include "slabs.h"
#endif
#ifdef seaice
#include "seayce.h"
#endif
#ifdef gotogtime
      common /gotocmn/ ittblah, yearsblah
#endif
c----------------------------------------------------
      entry rdrest
c----------------------------------------------------
c
#ifdef gotogtime
# ifndef gtime
      write(stdout,9990)
      write(stderr,9990)
      stop  'RDREST1'
9990  format(//60('#')/' ### ERROR ###'/
     $' ### option gotogtime will not work unless option gtime'/
     $' ### is also defined (here gtime was not defined)'/60('#')/)
# endif
#endif
c
#ifdef diskless
c
c-----------------------------------------------------------------------
c     read diskless restart file
c-----------------------------------------------------------------------
c
# ifdef t90_1st
c     when reading in a C90 or YMP restart file, one must read in
c     integers directly into integer variables, and not into floating
c     point variables
c
#  if defined gtime && !defined gotogtime
      read(iorest) itt, icyear, icdayofyr, icsecofday
      rewind iorest
#  else
      read(iorest) itt, years
      write(6,9009) itt, years
9009  format(/,' In restio, itt, years=',i10,e20.10,/)
      ittsav = itt
      savyear = years
      rewind iorest
#  endif
# endif
c
      call ord (iorest)
c
# ifdef t90_1st
#  if defined gtime && defined gotogtime
      itt = ittsav
      years = savyear
#  endif
#  ifdef f90type
#   if defined gtime && !defined gotogtime
      call oputi (kontrl, nvar11, 1, itt)
#   else
      call oputi (kontrl, 1, 1, itt)
      call oput  (kontrl, 1, 2, years)
#   endif
#  else
      call oput (kontrl, nvar11, 1, itt)
#  endif
# endif
# if defined gotogtime
c     read in restart that used old time manager kontrl, but
c     continue integration using the updated G group time manager.
c     also, attempt to accoount for possible previous time truncation
c     by assuming that the restart represents the start of a day,
c     so one will round to the nearest whole day (midnight) and proceed.
c
#  ifdef f90type
      call ogeti (kontrl, 1, 1, ittblah)
      call oget  (kontrl, 1, 2, yearsblah)
#  else
      call oget (kontrl,  2, 1, ittblah)
#  endif
c
      itt = ittblah
      icyear = int(yearsblah)
      daysofyr = (yearsblah-icyear)*365.0
      icdayofyr = int(daysofyr) + 1
      icsecofday = int (86400 * (daysofyr + 1 - icdayofyr))
c
      write(stdout,9018) itt, yearsblah, icyear, icdayofyr, icsecofday
c
c     round to nearest whole day
c
      icdayofyr = nint(daysofyr) + 1
      if (icdayofyr .EQ. 366) then
        icyear = icyear + 1
        icdayofyr = 1
      endif
      icsecofday = 0
c
#  ifdef f90type
#   ifdef gtime
      call oputi (kontrl, nvar11, 1, itt)
#   else
      call oputi (kontrl, 1, 1, itt)
      call oput  (kontrl, 1, 2, years)
#   endif
#  else
      call oput (kontrl, nvar11, 1, itt)
#  endif
c
      write(stdout,9019) icyear, icdayofyr, icsecofday
9018  format(' USING OPTION gotogtime TO SWITCH FROM OLD TIME MANAGER',
     $       ' TO UPDATED G GROUP TIME MANAGER'/' READ IN itt, years=',
     $       i18,2x,e15.9,/' CONVERTS TO icyear, icdayofyr, icsecofday='
     $       ,3i14)
9019  format(' ROUNDING TO THE NEAREST WHOLE DAY',
     $       /' ROUNDS TO.. icyear, icdayofyr, icsecofday=',3i14/)
# else
#  ifdef f90type
#   ifdef gtime
      call ogeti (kontrl, nvar11, 1, itt)
ccc   print *,' itt,icyear,icdayofyr,icsecofday =',itt,icyear,
ccc  &                     icdayofyr,icsecofday
#   else
      call ogeti (kontrl, 1, 1, itt)
      call oget  (kontrl, 1, 2, years)
#   endif
#  else
      call oget (kontrl, nvar11, 1, itt)
#  endif
# endif
c
# ifdef gtime
      years = icyear + (icdayofyr/365.0) + (icsecofday/(365.0*86400.0))
      write(stdout,9090) itt, icyear, icdayofyr, icsecofday
9090  format(1x,'IN RESTIO, ITT=',i14,' icyear=',i10,
     $       ' icdayofyr=',i4,' icsecofday=',i6,/)
#  ifdef f90type
      call oputi (kontrl, nvar11, 1, itt)
#  endif
# else
      write(stdout,9090) itt, years
9090  format(1x,'IN RESTIO, ITT=',i14,' YEARS=',f20.10,/)
#  ifdef f90type
       call oputi (kontrl, 1, 1, itt)
       call oput  (kontrl, 1, 2, years)
#  else
       call oput (kontrl, nvar11, 1, itt)
#  endif
# endif
# if defined cpldleg3 && defined cpldinit
c       reset indicators of ocean time since initialization to zero
c       at start of coupled integration.
c       timestamps will therefore be initialized 1 jan of year "year0"
c       ("year0" is set in blkdta)
c
c ####
c NOTE: USER MUST SET "ittcpldinit" ... the timestep read in from
c ####  the ocean restart tape when first coupling
c
      print*,'cpldinit ITT equals ',itt
cccccccccccccccccccccccccccccccccc
      ittcpldinit = 17003440
cccccccccccccccccccccccccccccccccc
      if (itt .eq. ittcpldinit) then
        ittin = itt
c       reset indicators of ocean time since initialization to zero
c       at start of coupled integration.
c       assuming year0 is 1 in blkdta, then reset the time since
c       initialization counter to zero
c
#  ifdef gtime
        icyearin = icyear
        icdayofyrin = icdayofyr
        icsecofdayin = icsecofday
        icyear = 0
        icdayofyr = 1
        icsecofday = 0
        years = 0.00
        write(stdout,9092) itt, icyearin, icdayofyrin, icsecofdayin
     $                        , icyear,   icdayofyr,   icsecofday
9092    format(/,1x,'IN RESTIO,   AT ITT=',i10,
     $'       icyear=',i10,' icdayofyr=',i4,' icsecofday=',i6,/
     $' RESET icyear=',i10,' icdayofyr=',i4,' icsecofday=',i6,/)
#   ifdef f90type
        call oputi (kontrl, nvar11, 1, itt)
#   else
        call oput (kontrl, nvar11, 1, itt)
#   endif
#  else
        yearin = years
        years = 0.00
        write(stdout,9092) itt, yearin, itt, years
9092    format(/,1x,'IN RESTIO,   AT ITT=',i10,'       YEARS=',f20.10,/,
     $           1x,'                    ',i10,' RESET YEARS=',f20.10,/)
#   ifdef f90type
        call oputi (kontrl, 1, 1, itt)
        call oput  (kontrl, 1, 2, years)
#   else
        call oput (kontrl, nvar11, 1, itt)
#   endif
#  endif
      endif
# endif
#else
c
c-----------------------------------------------------------------------
c     read unformatted i/o restart file
c-----------------------------------------------------------------------
c
c     record 1 is the control block
c     record's 2 & 3 are for stream function time levels
c     record's 4 & 5 are for the stream function guess field
c     record 6 is for the level indicators over "t" points
c
# ifdef gtime
c
#  ifdef gotogtime
c     read in restart that used old time manager kontrl, but
c     continue integration using the updated G group time manager.
c     also, attempt to accoount for possible previous time truncation
c     by assuming that the restart represents the start of a day,
c     so one will round to the nearest whole day (midnight) and proceed.
c
      read (iorest) ittblah, yearsblah
c
      itt = ittblah
      icyear = int(yearsblah)
      daysofyr = (yearsblah-icyear)*365.0
      icdayofyr = int(daysofyr) + 1
      icsecofday = int (86400 * (daysofyr + 1 - icdayofyr))
c
      write(stdout,9018) itt, yearsblah, icyear, icdayofyr, icsecofday
c
c     round to nearest whole day
c
      icdayofyr = nint(daysofyr) + 1
      if (icdayofyr .EQ. 366) then
        icyear = icyear + 1
        icdayofyr = 1
      endif
      icsecofday = 0
c
      write(stdout,9019) icyear, icdayofyr, icsecofday
9018  format(' USING OPTION gotogtime TO SWITCH FROM OLD TIME MANAGER',
     $       ' TO UPDATED G GROUP TIME MANAGER'/' READ IN itt, years=',
     $       i18,2x,f20.10/' CONVERTS TO icyear, icdayofyr, icsecofday='
     $       ,3i14)
9019  format(' ROUNDING TO THE NEAREST WHOLE DAY',
     $       /' ROUNDS TO.. icyear, icdayofyr, icsecofday=',3i14/)
#  else
      read (iorest) itt, icyear, icdayofyr, icsecofday
#  endif
c
      years = icyear + (icdayofyr/365.0) + (icsecofday/(365.0*86400.0))
      write (stdout,9100) iorest, itt, icyear, icdayofyr, icsecofday
#  ifdef f90type
      call oputi (kontrl, nvar11, 1, itt)
#  else
      call oput (kontrl, nvar11, 1, itt)
#  endif
c
# else
      read (iorest) itt, years
      write (stdout,9100) iorest, itt, years
#  ifdef f90type
      call oputi (kontrl, 1, 1, itt)
      call oput  (kontrl, 1, 2, years)
#  else
      call oput (kontrl, nvar11, 1, itt)
#  endif
c
# ifdef cpldleg2
c  "years" decreases from its true value by truncation error, so reset
c   it here. (Do not reset itt, disk pointers are set from it.)
ccccc if (itt .eq. 2376150) then
ccccc   years = 2170.e0
ccccc   write(stdout,9092) itt, years
9092    format(1x,'IN RESTIO, AT ITT=',i10,' RESET YEARS TO',f20.10,/)
ccccc endif
# endif
#endif
c
#if defined cpldleg3 && defined cpldinit
c       reset indicators of ocean time since initialization to zero
c       at start of coupled integration.
c       timestamps will therefore be initialized 1 jan of year "year0"
c       ("year0" is set in blkdta)
c #####
c NOTE: USER MUST SET "ittcpldinit" ... the timestep read in from
c ##### the ocean restart tape when first coupling
c
cccccccccccccccccccccccccccccccccc
      ittcpldinit = 23337141
cccccccccccccccccccccccccccccccccc
      if (itt .eq. ittcpldinit) then
        ittin = itt
c       reset indicators of ocean time since initialization to zero
c       at start of coupled integration.
c       assuming year0 is 1 in blkdta, then reset the time since
c       initialization counter to zero
c
# ifdef gtime
        icyearin = icyear
        icdayofyrin = icdayofyr
        icsecofdayin = icsecofday
        icyear = 0
        icdayofyr = 1
        icsecofday = 0
        years = 0.00
        write(stdout,9092) itt, icyearin, icdayofyrin, icsecofdayin
     $                        , icyear,   icdayofyr,   icsecofday
9092    format(/,1x,'IN RESTIO,   AT ITT=',i10,
     $'       icyear=',i10,' icdayofyr=',i4,' icsecofday=',i6,/
     $' RESET icyear=',i10,' icdayofyr=',i4,' icsecofday=',i6,/)
# else
        yearin = years
        years = 0.00
        write(stdout,9092) itt, yearin, itt, years
9092    format(/,1x,'IN RESTIO,   AT ITT=',i10,'       YEARS=',f20.10,/,
     $           1x,'                    ',i10,' RESET YEARS=',f20.10,/)
# endif
      endif
#endif
c
      relyr = years
#ifdef f90type
# ifdef gtime
      call oputi (kontrl, nvar11, 1, itt)
# else
      call oputi (kontrl, 1, 1, itt)
      call oput  (kontrl, 1, 2, years)
# endif
#else
      call oput (kontrl, nvar11, 1, itt)
#endif
      ndiskb = mod(itt+2,ntlev) + 1
      ndisk  = mod(itt  ,ntlev) + 1
      ndiska = mod(itt+1,ntlev) + 1
c
      read (iorest) ptd
#ifdef twoslabdisks
      call oput (kflds, nwds, mod(itt  ,3)*nwds+1, ptd)
#else
      call oput (kflds, nwds, (ndisk-1)*nwds+1, ptd)
#endif
c
      read (iorest) ptd
#ifdef twoslabdisks
      call oput (kflds, nwds, mod(itt+1,3)*nwds+1, ptd)
#else
      call oput (kflds, nwds, (ndiska-1)*nwds+1, ptd)
#endif
c
c
      read (iorest) ptd
      call oput (kflds, nwds, 3*nwds+1, ptd)
c
      read (iorest) ptd
      call oput (kflds, nwds, 4*nwds+1, ptd)
c
      read (iorest) kmt
# ifdef f90type
      call oputi (kflds, nwds, 5*nwds+1, kmt)
# else
      call oput (kflds, nwds, 5*nwds+1, kmt)
# endif
c
# ifdef seaice
c     restart records 7 & 8 are for the 2 time levels of sea ice
c     first read the tau-1 level, then the tau time level.
c     (using "yca" just for temporary storage, here)
c
      read (iorest) yca
      call oput (kflds, nwds, (5+ndisk )*nwds+1, yca)
c
      read (iorest) yca
      call oput (kflds, nwds, (5+ndiska)*nwds+1, yca)
c
# endif
c     next come the slabs (2 time levels)
c
      do 100 j=1,jmt
       read (iorest) bufsl
       call oput (labs(ndisk), nslab, (j-1)*nslab+1, bufsl)
       if (j .eq. 1 .or. j .eq. jmt) then
         call oput (labs(ndiskb), nslab, (j-1)*nslab+1, bufsl)
       endif
100   continue
c
      do 200 j=1,jmt
       read (iorest) bufsl
       call oput (labs(ndiska), nslab, (j-1)*nslab+1, bufsl)
200   continue
c
# ifdef gtime
9100  format (1x,'===> reading ocean restart data from unit',i3
     $,' on ts=',i10,/,' ===>  icyear=',i10,
     $' icdayofyr=',i4,' icsecofday=',i6,/)
# else
9100  format (1x,'===> reading ocean restart data from unit',i3
     $,      ' on ts=', i10,' years=',f20.10,/)
# endif
c
      rewind iorest
c
#endif
      return
c
c----------------------------------------------------
      entry wrrest
c----------------------------------------------------
# ifdef t90_1st
c     to read in an old YMP or C90 restart file written in the former
c     Cray native form, we have assigned unit iorest (as set in
c     blkdta.F) with the -N cray option for conversion purposes.
c     But, if we wish to have all future restarts written in the
c     T90 64-bit IEEE form, we need to write out the restart to a
c     different file ... one whose assign does not contain the
c     -N cray option.
c     Here the code simply resets the iorest value upon writing out
c     the restart file.  the default is to add 1 to the iorest value
c     set in blkdta ... the user should check and make sure there is
c     no conflict (i.e., that unit number (iorest+1) was not already
c     in use).     -- KD Sept 1996
c
cccc  iorest = iorest + 1
      iorest = 19
# endif
      write (stdout,9200) '==>  Writing  MOM restart on unit ',iorest
     $,                    itt, stamp
#ifdef diskless
cnetcdf
      call bin2nc
      
cnetcdf

      call owrt (iorest)
#else
c
c-----------------------------------------------------------------------
c     write restart file using unformatted fortran i/o
c-----------------------------------------------------------------------
c
c     record 1 is the control block
c     record's 2 & 3 are for stream function time levels
c     record's 4 & 5 are for the stream function guess field
c     record 6 is for the level indicators over "t" points
c
# ifdef gtime
#  ifdef f90type
      call ogeti (kontrl, nvar11, 1, itt)
#  else
      call oget (kontrl, nvar11, 1, itt)
#  endif
      write(stdout,9097) itt, icyear, icdayofyr, icsecofday
9097  format(1x,'IN WRREST, ITT=',i14,' icyear=',i10,
     $       ' icdayofyr=',i4,' icsecofday=',i6,/)
      write (iorest) itt, icyear, icdayofyr, icsecofday
# else
#  ifdef f90type
      call ogeti (kontrl, 1, 1, itt)
      call oget  (kontrl, 1, 2, years)
#  else
      call oget (kontrl, nvar11, 1, itt)
#  endif
      write (iorest) itt, years
# endif
c
#ifdef twoslabdisks
      call oget (kflds, nwds, mod(itt  ,3)*nwds+1, ptd)
#else
      call oget (kflds, nwds, (ndisk-1)*nwds+1, ptd)
#endif
      write (iorest) ptd
c
#ifdef twoslabdisks
      call oget (kflds, nwds, mod(itt+1,3)*nwds+1, ptd)
#else
      call oget (kflds, nwds, (ndiska-1)*nwds+1, ptd)
#endif
      write (iorest) ptd
c
      call oget (kflds, nwds, 3*nwds+1, ptd)
      write (iorest) ptd
c
      call oget (kflds, nwds, 4*nwds+1, ptd)
      write (iorest) ptd
c
# ifdef f90type
      call ogeti (kflds, nwds, 5*nwds+1, kmt)
# else
      call oget (kflds, nwds, 5*nwds+1, kmt)
# endif
      write (iorest) kmt
c
# ifdef seaice
c     record's 7 & 8 are for sea ice
c     first write tau-1 timelevel of sea ice, then the tau time level
c
      call oget (kflds, nwds, (5+ndisk )*nwds+1, yca)
      write (iorest) yca
c
      call oget (kflds, nwds, (5+ndiska)*nwds+1, yca)
      write (iorest) yca
c
# endif
c
c     next come the slabs (2 time levels)
c
      do 400 j=1,jmt
        call oget (labs(ndisk), nslab, (j-1)*nslab+1, bufsl)
        write (iorest) bufsl
400   continue
c
      do 500 j=1,jmt
        call oget (labs(ndiska), nslab, (j-1)*nslab+1, bufsl)
        write (iorest) bufsl
500   continue
c
#endif
      return
9200   format (/1x,a34,i3,', ts = ',i10,1x, a32/)
      end
