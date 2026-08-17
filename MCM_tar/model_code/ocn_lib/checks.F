      subroutine checks
      logical errorc, test
#include "param.h"
#include "accel.h"
#include "chmix.h"
#ifdef isopycmix
# include "cisop.h"
#endif
#include "coord.h"
#ifdef ppvmix
# include "cppmix.h"
#endif
#ifdef tcvmix
# include "ctcmix.h"
#endif
#include "cvmix.h"
#include "ctask.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
c
c-----------------------------------------------------------------------
c     do consistency checks before allowing MOM to continue
c-----------------------------------------------------------------------
c
      errorc = .false.
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
      write (stdout,9000) ' Consistency check messages:'
c
c
      if (ntlev .eq. 2 .and. ntasks .ne. 1) then
        write (stdout,9000)
     $  '==> Error: ntasks > 1 is presently incompatible with the     '
     $, '           diskless option. use ntasks=1 or use another I/O  '
     $, '           option.  (or change "ntlev" to 3 but this will    '
     $, '           increase the memory and was not checked  out)     '
        errorc = .true.
      endif
#if !defined diskless && !defined crayio && !defined fio
        write (stdout,9000)
     $  '==> Error: either "diskless", crayio" or "fio"               '
     $, '           option must be turned on.                         '
        errorc = .true.
#endif
#if defined diskless && defined crayio
        write (stdout,9000)
     $  '==> Error:  "diskless" and "crayio" are incompatible         '
     $, '             choose one or the other option.                 '
        errorc = .true.
#endif
#if defined diskless && defined fio
        write (stdout,9000)
     $  '==> Error:  "diskless" and "fio" are incompatible            '
     $, '             choose one or the other option.                 '
        errorc = .true.
#endif
#if defined crayio && defined fio
        write (stdout,9000)
     $  '==> Error:  "crayio" and "fio" are incompatible              '
     $, '             choose one or the other option.                 '
        errorc = .true.
#endif
c
#ifdef twoslabdisks
# ifdef diskless
        write (stdout,9000)
     $  '==> Warning: since the "diskless" option stores the slabs in '
     $, '             memory, selecting "twoslabdisks" has no effect. '
# endif
# ifdef multitasking
        write (stdout,9000)
     $  '==> Error:   given the unpredictable order in which slabs are'
     $, '             updated when multitasking, 3 slab disks are     '
     $, '             needed, making "twoslabdisks" an incompatible   '
     $, '             option.                                         '
        errorc = .true.
# endif
#endif
c
#ifdef isopycmix
      if (ntasks .gt. 1) then
         write (stdout,9000)
     $  '==> Error:  isopycmix does not yet work with ntasks > 1      '
     $, '            use ntasks=1                                     '
        errorc = .true.
      endif
      if (aidif .lt. p5) then
         write (stdout,9000)
     $  '==> Error:  isopycmix must use aidif >= 0.5                  '
     $, '            (must be fully implicit)                         '
        errorc = .true.
      endif
# ifdef skipland
        write (stdout,9000)
     $  '==> Error:  isopycmix does not yet work with skipland        '
        errorc = .true.
# endif
# ifdef multitasking
      write (stdout,9000)
     $  '==> Error:  isopycmix does not yet work with multitasking    '
        errorc = .true.
# endif
# if !defined ppvmix && !defined tcvmix && !defined constvmix
        write (stdout,9000)
     $  '==> Error: either "ppvmix", "tcvmix"  or "constvmix"         '
     $, '           option must be turned on with isopycmix.          '
        errorc = .true.
# endif
# if !defined consthmix && !defined biharmonic && !defined nlhmix
        write (stdout,9000)
     $  '==> Error: since isopycmix only affects tracer diffusion,    '
     $, '    either "consthmix", "biharmonic" or "nlhmix" must be     '
     $, '    turned on with for use in mixing momentum in clinic      '
     $, '    (only "consthmix" has been checked out)                  '
        errorc = .true.
# endif
# if defined consthmix
        write (stdout,9000)
     $  '==> Note: consthmix will only affect mixing of momentum      '
     $, '    since isopycmix was specified for tracer diffusion       '
# endif
# if defined biharmonic
        write (stdout,9000)
     $  '==> Warning: biharmonic will only affect mixing of momentum  '
     $, '    since isopycmix was specified for tracer diffusion       '
     $, '    (this combination is not fully checked out)              '
# endif
# if defined nlhmix
        write (stdout,9000)
     $  '==> Warning: nlhmix will only affect mixing of momentum      '
     $, '    since isopycmix was specified for tracer diffusion       '
     $, '    (this combination is not fully checked out)              '
# endif
#endif
c
c
#if !defined constvmix && !defined ppvmix && !defined tcvmix
# if !defined isopycmix
        write (stdout,9000)
     $  '==> Error: either "constvmix", ppvmix", "tcvmix"             '
     $, '           or "isopycmix" option must be turned on.          '
        errorc = .true.
# endif
#endif
#ifdef gargettmix
        write (stdout,9000)
     $  '==> Error:  "gargettmix" not implemented yet                 '
        errorc = .true.
#endif
#ifdef windvmix
        write (stdout,9000)
     $  '==> Error:  "windvmix" not implemented yet                   '
        errorc = .true.
#endif
#if defined biharmonic && !defined multitasking
      if (ntasks .gt. 1) then
        write (stdout,9000)
     $  '==> Error: ntasks > 1 is incompatible with the biharmonic    '
     $, '           option. when not multitasking with biharmonic     '
     $, '           mixing use ntasks=1                               '
        errorc = .true.
      endif
#endif
#if defined blvmix && !defined constvmix
        write (stdout,9000)
     $  '==> Error:  "constvmix" must be defined if "blvmix" is used. '
     $, '    when "blvmix" is specified for vertical tracer diffusion '
     $, '    "constvmix" is the only vertical viscosity scheme that   '
     $, '    is checked out.                                          '
        errorc = .true.
#endif
#ifdef blvmix
# ifdef constvmix
        write (stdout,9000)
     $  '==> Warning: blvmix influences vertical tracer diffusion only'
     $, '    vertical mixing of momentum is unaffected                '
# endif
# if defined ppvmix || defined tcvmix
        write (stdout,9000)
     $  '==> Error:   when "blvmix" is specified for vertical tracer  '
     $, '    diffusion, "constvmix" must be defined for the vertical  '
     $, '    mixing of momentum.  Neither "ppvmix" nor "tcvmix" have  '
     $, '    been checked out to work for momentum when "blvmix" is   '
     $, '    set from vertical tracer diffusion.                      '
        errorc = .true.
# endif
#endif
#if defined blhmix
#if !defined consthmix && !defined isopycmix
        write (stdout,9000)
     $  '==> Error:  "consthmix" or "isopycmix" must be defined if    '
     $       "blhmix" is used. '
        errorc = .true.
#endif
# if defined consthmix
        write (stdout,9000)
     $  '==> Warning: blhmix influences horizontal tracer diffusion   '
     $, '    only - vertical mixing of momentum is unaffected         '
# endif
# if defined nlhmix || defined biharmonic
        write (stdout,9000)
     $  '==> Error:   when "blhmix" is specified for horizontal tracer'
     $, '    diffusion, "consthmix"  must be defined for momentum     '
     $, '    neither "nlhmix" nor "biharmonic" have been checked out  '
     $, '    for momentum when "blhmix" is set                        '
        errorc = .true.
# endif
#endif
c
c
#if defined freesfc && defined rigidlid
        write (stdout,9000)
     $  '==> Error:  "freesfc" and "rigidlid" are incompatible        '
     $, '             choose one or the other option.                 '
        errorc = .true.
#endif
#ifdef freesfc
        write (stdout,9000)
     $  '==> Error:  "freesfc" not implemented yet                    '
        errorc = .true.
#endif
#ifndef rigidlid
        write (stdout,9000)
     $  '==> Error: "rigidlid" option is required in this version     '
     $, '            option must be turned on.                        '
        errorc = .true.
#endif
#ifdef rigidlid
# if !defined congrad5pt && !defined congrad9pt && !defined oldrelax
#  if !defined hypergrid
        write (stdout,9000)
     $  '==> Error: either "congrad5pt", congrad9pt", "oldrelax", or  '
     $, '           "hypergrid" option must be turned on.             '
        errorc = .true.
#  endif
# endif
# if defined congrad5pt && defined congrad9pt
        write (stdout,9000)
     $  '==> Error:  "congrad5pt" and "congrad9pt" are incompatible   '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined oldrelax && defined hypergrid
        write (stdout,9000)
     $  '==> Error:  "oldrelax" and "hypergrid" are incompatible      '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined oldrelax && defined congrad5pt
        write (stdout,9000)
     $  '==> Error:  "oldrelax" and "congrad5pt" are incompatible     '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined oldrelax && defined congrad9pt
        write (stdout,9000)
     $  '==> Error:  "oldrelax" and "congrad9pt" are incompatible     '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined hypergrid && defined congrad9pt
        write (stdout,9000)
     $  '==> Error:  "hypergrid" and "congrad9pt" are incompatible    '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
# if defined hypergrid && defined congrad5pt
        write (stdout,9000)
     $  '==> Error:  "hypergrid" and "congrad5pt" are incompatible    '
     $, '             choose one or the other option.                 '
        errorc = .true.
# endif
#endif
c
c
#if !defined biharmonic && !defined nlhmix  && !defined consthmix
# if !defined isopycmix
        write (stdout,9000)
     $  '==> Error: either "biharmonic", nlhmix", "consthmix"         '
     $, '           or "isopycmix" option must be turned on.          '
        errorc = .true.
# endif
#endif
#if defined biharmonic && defined nlhmix
        write (stdout,9000)
     $  '==> Error: "biharmonic" and "nlhmix" are incompatible.       '
     $, '            choose one or the other option.                  '
        errorc = .true.
#endif
#if defined biharmonic && defined consthmix
        write (stdout,9000)
     $  '==> Error: "biharmonic" and "consthmix" are incompatible.    '
     $, '            choose one or the other option.                  '
        errorc = .true.
#endif
#if defined nlhmix && defined consthmix
        write (stdout,9000)
     $  '==> Error: "nlhmix" and "consthmix" are incompatible.        '
     $, '            choose one or the other option.                  '
        errorc = .true.
#endif
#if defined biharmonic
        if (am .lt. 1.e11 .or. ah .lt. 1.e11) then
          write (stdout,9000)
     $  '==> Error: "am" and/or "ah" are too small for the            '
     $, '             "biharmonic"  mixing option                     '
          errorc = .true.
        endif
        if (am .lt. c0 .or. ah .lt. c0) then
          write (stdout,9000)
     $  '==> Error: use positive "am" and "ah" for "biharmonic"       '
     $, '           option.  The (neg) sign is in the equation.       '
          errorc = .true.
        endif
#endif
#if defined consthmix || defined isopycmix
      scdkmx = hdcsck(1)
      scvkmx = hvcsck(1)
      do 90 k=2,km
      if (hdcsck(k) .gt. scdkmx) scdkmx = hdcsck(k)
      if (hvcsck(k) .gt. scvkmx) scvkmx = hvcsck(k)
90    continue
#endif
#if defined consthmix
        if (am*scvkmx .gt. 1.e11) then
          write (stdout,9000)   
     $  '==> Error: "am"*"hvcsck" is too big for the "consthmix"      '
     $, '           mixing option                                     '  
          errorc = .true.
        endif
#endif   
#if defined consthmix && !defined isopycmix
        if (ah*scdkmx .gt. 1.e11) then
          write (stdout,9000)   
     $  '==> Error: "ah"*"hdcsck" is too big for the "consthmix"      '
     $, '           mixing option                                     '  
          errorc = .true.
        endif
#endif
# ifdef isopycmix
        if ((ah+ahisop*scdkmx) .gt. 1.e11) then
          write (stdout,9000)
     $  '==> Error: "ahisop"*"hdcsck" + "ah" is too large for the     '  
     $, '             "isopycmix" mixing option                       '  
          errorc = .true.
        endif
#endif
c
c
#if defined skipland && defined tcmix
        write (stdout,9000)
     $  '==> Warning: "skipland" not yet implemented in "tcmix"       '
#endif
#if defined multitasking && defined nlhmix
        write (stdout,9000)
     $  '==> Error: "nlhmix" not yet implemented for "multitasking"   '
        errorc = .true.
#endif
c
c
c
      if ((dtuv .ne. dtsf) .or. (dtuv .ne. dtts)) then
        write (stdout,9000)
     $  '==> Warning: use of unequal time steps implies the transient '
     $, '             response is unimportant and multiple equilibria '
     $, '             do not exist.                                   '
      endif
c
c     check for mixing coefficients larger than stability permits
c
      dymin  = dyt(2)
      dxmin  = dxt(2)
      dzmin  = dzt(1)
      xlmax  = dtxcel(1)
      do 100 j=2,jmtm1
        dymin  = min(dymin,dyt(j))
100   continue
      do 200 i=2,imtm1
        dxmin  = min(dxmin,dxt(i))
200   continue
      do 300 k=2,km
        xlmax  = max(xlmax,dtxcel(k))
        dzmin  = min(dzmin,dzt(k))
300   continue
c
      if (xlmax .gt. c1) then
        write (stdout,9000)
     $  '==> Warning: use of accelerated time steps implies the       '
     $, '             transient response is unimportant and multiple  '
     $, '             equilibria do not exist. stability tests will   '
     $, '             use "dtts" multiplied by the maximum "dtxcel"   '
      endif
c
c
c
#if defined consthmix || defined isopycmix
      scdkmx = dtts*hdcsck(1)*dtxcel(1)
      scvkmx = dtuv*hvcsck(1)
      do 345 k=2,km
        if (dtts*hdcsck(k)*dtxcel(k) .gt. scdkmx) 
     $                        scdkmx = dtts*hdcsck(k)*dtxcel(k)
        if (dtuv*hvcsck(k) .gt. scvkmx) scvkmx = dtuv*hvcsck(k)
345   continue
      do 350 j=2,jmtm1
        dxdymn = c1/(c1/(dxmin*cst(j))**2 + c1/dymin**2)
# ifdef isopycmix
c  ***not always exactly correct, since xlmax & scdkmx do not have to
c     be greatest for the same k-level**
        if ((ah*xlmax+scdkmx*ahisop)/dxdymn .ge. p25) then
# else
        if ((scdkmx*asubh)/dxdymn .ge. p25) then
# endif
          write (stdout,9000)
# ifdef isopycmix
     $  '==> Warning: lateral diffusive criteria may be exceeded for  '
     $, '             "ah" + "ahisop". use a smaller "dtts", "dtxcel",'
     $, '             "hdcsck", and/or "ah" + "ahisop"                '
# else
     $  '==> Warning: lateral diffusive criteria exceeded for "ah".   '
     $, '             use a smaller "dtts", "dtxcel", "hdcsck",       '
     $, '             and/or "ah"                                     '
# endif
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     $,                                     ',  j=',j
        endif
        if ((scvkmx*asubm)/dxdymn .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: lateral diffusive criteria exceeded for "am".   '
     $, '             use a smaller "dtuv", "hvcsck" and/or "am"      '
          write (stdout,'(a48,f6.2,a5,i3)') ' at latitude ',yt(j)
     $,                                     ',  j=',j
        endif
350   continue
#endif
#ifdef constvmix
      scdimx = vdcsci(1)
      scdjmx = vdcscj(1)
      scvimx = vvcsci(1)
      scvjmx = vvcscj(1)
      do 356 i=1,imt
        if (vdcsci(i) .gt. scdimx) scdimx = vdcsci(i)
        if (vvcsci(i) .gt. scvimx) scvimx = vvcsci(i)
356   continue
      do 358 j=1,jmt
        if (vdcscj(j) .gt. scdjmx) scdjmx = vdcscj(j)
        if (vvcscj(j) .gt. scvjmx) scvjmx = vvcscj(j)
358   continue
      do 360 k=1,km
        if ((dtts*dtxcel(k)*scdimx*scdjmx*vdcsck(k)*fkph/dzt(k)**2)
     $       .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for "fkph".'
     $, '             use a smaller "dtts", "dtxcel", "vdc" scale     '
     $, '             factors, and/or "fkph"                          '
         write (stdout,'(a48,i3)') ' at level =',k
        endif
360   continue
      do 365 k=1,km
        if ((dtuv*scvimx*scvjmx*vvcsck(k)*fkpm/dzt(k)**2) .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for "fkpm".'
     $, '           use a smaller "dtuv", "vvc" scale factors, and/or '
     $, '           "fkpm"                                            '  
         write (stdout,'(a48,i3)') ' at level =',k
      endif
365   continue
#endif
#ifdef implicitvmix
      write (stdout,9000)
     $  '==> Warning: "implicitvmix" turns convective adjustment off. '
     $, '              It will be handled by the vertical mixing      '
#else
      write (stdout,9000)
     $  '==> Warning: convective adjustment is enabled. Verify that   '
     $, '             "ncon" is set appropriately                     '
#endif
#if defined ppvmix && !defined implicitvmix
      do 370 k=1,km
        if ((dtts*dtxcel(k)*fricmx)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '           "fricmx".  use a smaller "dtts", "dtxcel", and/or '
     $, '           "fricmx"                                          '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
        if ((dtts*dtxcel(k)*vdclim)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '            "vdclim". use a smaller "dtts", "dtxcel", and/or '
     $, '            "vdclim"                                         '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
370   continue
      if ((dtuv*fricmx)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '            "fricmx". use a smaller "dtuv" and/or "fricmx"   '
      endif
      if ((dtuv*vvclim)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '           "vvclim".  use a smaller "dtuv" and/or "vvclim"   '
      endif
#endif
#if defined tcvmix && !defined implicitvmix
      do 380 k=1,km
        if ((dtts*dtxcel(k)*vdclim)/dzt(k)**2 .ge. p25) then
          write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '            "vdclim".  use a smaller "dtts", "dtxcel", and/or'
     $, '            "vdclim"                                         '
          write (stdout,'(a48,i3)') ' at level =',k
        endif
380   continue
      if ((dtuv*vvclim)/dzmin**2 .ge. p25) then
        write (stdout,9000)
     $  '==> Warning: vertical diffusive criteria exceeded for        '
     $, '            "vvclim".  use a smaller "dtuv" and/or "vvclim"  '
      endif
#endif
c
      if (cdbot .eq. c0) then
        test = .false.
        do 400 j=2,jmtm1
          do 390 i=2,imtm1
            if (kmt(i,j) .ne. 0 .and. kmt(i,j) .ne. km) test = .true.
390       continue
400     continue
        if (test) then
          write (stdout,9000)
     $  '==> Warning: rough topography may induce checkerboarding     '
     $, '             especially in shelf regions. use of             '
     $, '             "cdbot" > 0.0 may help to supress it.           '
        endif
      endif
c
c     check range of implicit factors
c
#if defined implicitvmix || defined isopycmix
      if (aidif .lt. c0 .or. aidif .gt. c1) then
        write (stdout,9000)
     $  '==> Error: "aidif" must be between 0.0 & 1.0                 '
          errorc = .true.
      endif
#else
      if (aidif .gt. c0) then
        write (stdout,9000)
     $  '==> Warning: "aidif" was set > 0.0  but "implicivmix" is not '
     $, '              defined, so vertical mixing  will be  explicit '
     $, '                                      '
      endif
#endif
c
      if (acor .gt. c0) then
        write (stdout,9000)
     $  '==> Warning: "acor" was set > 0.0  implying that the         '
     $, '              coriolis term will be handled implicitly       '
      endif
      if (acor .lt. c0 .or. acor .gt. c1) then
        write (stdout,9000)
     $  '==> Error: "acor" must be between 0.0 & 1.0                  '
          errorc = .true.
      endif
c
c     check for resolved inertial oscillation
c
      if (acor .eq. c0) then
        fmax = abs(fcor(2))
        do 410 j=3,jmtm2
          fmax = max(fmax,abs(fcor(j)))
410     continue
        if (dtuv .gt. p25*(c2*pi)/fmax) then
          write (stdout,9000)
     $  '==> Warning: the inertial oscillation is not resolved.       '
     $, '             reduce "dtuv" or use "acor" > 0.0               '
        endif
      endif
#ifdef worldocn
# ifdef amazon
#  if defined implicitvmix || defined isopycmix
c     amazon mixing only work thru the implicit treatment of the
c     vertical tracer fluxes
#  else
        write (stdout,9000)
     $  '==> Error: "amazon" will only work when either "isopycmix"   '  
     $, '             or "implicitvmix" are also selected             '
          errorc = .true.
#  endif
# endif
c
c     check for inconsistent ifdef choices for GFDL world ocean runs
c
# if defined ocnr15 && defined ocnr30
        write (stdout,9000)
     $  '==> Error: "ocnr15" & "ocnr30" can not both be defined       '  
     $, '             select just one horizontal resolution           '
          errorc = .true.
# endif
# if defined km12 && defined km18
        write (stdout,9000)
     $  '==> Error: "km12" & "km18" can not both be defined           '  
     $, '             select just one vertical resolution             '
          errorc = .true.
# endif
# if defined restorst
#  if agetracers || defined cfctracers
        write (stdout,9000)
     $  '==> Error: "agetracers" and "cfctracers" have been set up    '  
     $, '           to only work as surface fluxes, and do not work   '
     $, '           when "restorst" is defined.                       '
          errorc = .true.
#  endif
# endif
# if defined agetracers && !defined keepterms
        write (stdout,9000)
     $  '==> Error: "keepterms" must be defined when asking for       '  
     $, '           "agetracers" to exclude diffusion from one tracer '
          errorc = .true.
# endif
# if defined cfctracers
#  if defined cjbcfcs && defined gptcfcs
        write (stdout,9000)
     $  '==> Error: either "cjbcfcs" or "gptcfcs" must be defined     '  
     $, '           when "cfctracers" is defined - not both           '
          errorc = .true.
#  endif
#  if !defined cjbcfcs && !defined gptcfcs
        write (stdout,9000)
     $  '==> Error: either "cjbcfcs" or "gptcfcs" must be defined     '  
     $, '           when "cfctracers" is defined                      '
          errorc = .true.
#  endif
# endif
#endif
#ifdef consthmix
c
c-----------------------------------------------------------------------
c     search for topographic instabilities (based  on the  work of
c     Peter Killworth  ...  eqn 11 from ocean modeling nov 1987)
c-----------------------------------------------------------------------
c
      num   = 50
      do 500 j=2,jmtm1
        dysq = dyt(j)**2
        do 490 i=2,imtm1
          if (kmt(i+1,j-1) .ne. 0 .and. kmt(i+1,j) .ne. 0) then
            dxsq = (dxt(i)*cst(j))**2
            clix = am*dtuv/dxsq
            h1   = zw(kmt(i+1,j-1))
            h2   = zw(kmt(i+1,j))
            hx   = (8.0*h1*h2/(h1+h2)**2 + dxsq/dysq)/(4.0 + dxsq/dysq)
            if (clix .ge. hx .and. num .ge. 0) then
              num = num - 1
              write (stdout,'(1x,a56,a31,i4,a1,i4,a1)')
     $        '==>   Warning: Killworth topographic roughness condition'
     $,       ' exceeded at location (i,j) = (',i+1,',',j,')'
              if (num .eq. 0) then
                write (stdout,9000)
     $         '==> Warning: msgs terminated after 50 cases were found '
              endif
            endif
          endif
490     continue
500   continue
#endif
c
#ifdef sfctopog
# ifndef rigidlid
        write (stdout,9000)
     $  '==> Error: "rigidlid" must be defined along with "sfctopog"   '
          errorc = .true.
# endif
# ifndef selfdoc
        write (stdout,9000)
     $  '==> Error: "selfdoc" must be defined along with "sfctopog"    '
          errorc = .true.
# endif
#endif
c
#if defined cpldleg3
# ifdef seasonal
c   For external mode forcing, annual wind stress is used instead
c    of monthly wind stress when dtts/dtsf > 4.
c   Is annual mean wind stress available for the coupled leg3 model?
        if (dtts/dtsf .gt. 4.0) then
          write (stdout,9000)
     $    '==> Error: dtts/dtsf > 4. Is annual wind stress available? '
          errorc = .true.
        endif
# endif
#endif
c
      write (stdout,9000) ' End of consistency check messages:'
     $,                   ' '
      if (errorc) stop '=>checks'
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
c
c
      return
9000  format (/,(1x,a80))
      end
