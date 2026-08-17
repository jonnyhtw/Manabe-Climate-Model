      subroutine setvbc (j)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
#ifdef worldocn
c=======================================================================
c
c     compute vertical boundary conditions for clinic & tracer
c     for the GFDL world ocean model.
c
c           NOTE: The MOM sample case can NOT be run with this  
c                 version of "setvbc". 
c
c=======================================================================
c
#include "param.h"
c
#if defined cpldleg2 || defined cpldleg3
# include "ascomo.h"
#endif
#include "coord.h"
#if defined cpldleg2 || defined cpldleg3
# include "ctask.h"
#endif
#include "cvbc.h"
#include "grdvar.h"
#include "levind.h"
#include "scalar.h"
#ifdef seaice
# include "seayce.h"
#endif
#include "slabs.h"
#include "ctmngr.h"
#include "iounit.h"
# ifdef seasonal
#include "seasnl.h"
# endif
#include "switch.h"
#ifdef amazon
# include "camazon.h"
#endif
c
# ifdef bctest
c
c     dimension land/masks locally for test of surface boundary
c     conditions that are read in
c
      dimension lsts(imt,jmt), lsuv(imt,jmt)
c
# endif
#ifdef cybertra
      dimension tra(imt)
#endif
c
c=======================================================================
c
c    CAUTION: "j" is in the argument list; do not use it as a loop
c             index. Store j from the argument list to prevent
c             clobbering it; then reload j below.
c
      jarg = j
c
# ifndef seasonal
c
c----------------------------------------------------------------------
c     if not running seasonal, run a perpetual case
c----------------------------------------------------------------------
c
        mon = 1
# endif
c
      if (j.eq.2 .and. .not. mxpas2) then
c
#if defined cpldleg3 && !defined cpld_restore
# ifdef seasonal
c   read surface flux adjustment data
        if (first) then
c   read data for mahd into terr(i,j,1) and mbhd into terr(i,j,2)
          call getfad
          write(stdout,9720) mbhd
          write(stdout,9720) mahd
9720      format(1x,'Read surface flux adjustment data for month=',i3)

        else
c   shift previous "ahead" month to "behind" location for the first
c    timestep needing a new "ahead" month
          if(newmon .eq. 1) then
            do 22 jj=1,jmt
              do 22 i=1,imt
                terr(i,jj,2) = terr(i,jj,1)
                serr(i,jj,2) = serr(i,jj,1)
22          continue
c   read data for new mahd into terr(i,j,1)
            call getfad
            write(stdout,9720) mahd
          endif
        endif
# endif
c
c   For coupled models, boundary data from the atmos model is stored in
c   the bco array in the atmos segment; not necessary to read it here
c
# if defined ocnr30 && defined amazon
      if (loop1) then
c     adjust very large values of p-e+r at Amazon, etc.
        write(6,9991) itt
9991    format(1x,'Adjusting per boundary data at ts=',i8)
        call adjper(bco(1,1,2,1))
      endif
# endif
#endif
c
#if defined cpldleg2
# ifdef seasonal
c   read boundary data from atmos
      if (first) then

c   read data for mahd into bco(i,j,n,1) and mbhd into bco(i,j,n,2)
        call getvbc(iobc, mahd, bcfld)
        write(stdout,9992) iobc, mahd
9992    format(1x,'Read ocean surface boundary data from unit',i3,
     $            ' for month=',i3)
        write(stdout,9992) iobc, mbhd
c
# if defined ocnr30 && defined amazon
c   adjust very large values of p-e+r
        write(6,9991) mahd
9991    format(1x,'Adjusting per boundary data for month=',i3)
        call adjper(bco(1,1,2,1))
        write(6,9991) mbhd
        call adjper(bco(1,1,2,2))
# endif
c
      else
c   shift previous "ahead" month to "behind" location for the first
c    timestep needing a new "ahead" month
        if(newmon .eq. 1) then
          do 22 n=1,nbo
            do 22 jj=1,jmt
              do 22 i=1,imt
                bco(i,jj,n,2) = bco(i,jj,n,1)
22        continue
c   read data for new month "ahead" into bco(1,1,n,1)
          call getvbc(iobc, mahd, bcfld)
          write(stdout,9992) iobc, mahd
c
# if defined ocnr30 && defined amazon
c   adjust very large values of p-e+r
          write(6,9991) mahd
          call adjper(bco(1,1,2,1))
# endif
c
        endif
      endif
# endif
#endif
        if (first) then
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c       read in GFDL world ocean surface boundary conditions
c       (done only for the first timestep of each run)
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
#ifdef cybertra
c   set tracera = 50 everywhere at the surface:
        do 5 i=1,imt
          tra(i) = 50.0
5       continue
#endif
c
#if defined cpldleg2 || defined cpld_restore
c   Initialize terr and serr arrays
        icntts = 0
        do 745 jj=1,jmt
          do 745 i=1,imt
            terr(i,jj) = c0
            serr(i,jj) = c0
745     continue
c
c    Read observed surface temperature, salinity, ice
        rewind ioobs
        do 760 mon=1,12
          read(ioobs) ((sstobs(i1,j1,mon),i1=1,imt),j1=1,jmt)
          read(ioobs) ((salobs(i2,j2,mon),i2=1,imt),j2=1,jmt)
          read(ioobs) ((yceobs(i3,j3,mon),i3=1,imt),j3=1,jmt)
          do 765 j2=1,jmt
            do 762 i2=1,imt
              salobs(i2,j2,mon) = (salobs(i2,j2,mon)-35.0) / 1000.0
762         continue
765       continue
760     continue
        write(stdout,9760) itt, mahd, sstobs(14,14,mahd),
     $                                salobs(14,14,mahd),
     $                                yceobs(14,14,mahd)
9760    format(1x,
     $   'Read observed surface temperature, salinity, ice at itt=',i10,
     $   /,'  sstobs, salobs, yceobs for mahd=',i3,' at (14,14,mahd):',
     $   3e13.5)
#endif
c
#ifdef cybervbc
c----------------------------------------------------------------------
c       read in surface bc used in cyber r15 coupled model runs
c       stored in /ymp/o/kd/models/model_input/R15/seasonal_bc
c       (same as used in cyber r15 coupled model experiments;
c       taken from cyber205 tape=8791,part=13)
c----------------------------------------------------------------------
c
        open (iobc,form='UNFORMATTED',
     $             file='seasonalbc')
        read(iobc) sstobs
        read(iobc) salobs
        read(iobc) wsxobs
        read(iobc) wsyobs
        read(iobc) yceobs
        close(iobc,status='KEEP')
#endif
#ifdef r15vbc
c
c----------------------------------------------------------------------
c       read in surface bc used in Cray r15 model runs
c       (one type of boundary condition is kept in each file)
c----------------------------------------------------------------------
c      
#   ifdef bctest
c
c       to check boundary conditions, the surface land/sea masks
c       are needed
c
        write(stdout,8111)
8111    format(/' checking on Cray boundary conditions')
c
        call Setlsm (lsts, lsuv)
c
#   endif
c
c       get sea surface temperatures, one month at a time
c
        open (iobc,form='UNFORMATTED',
     $       file = './sst')
c
        do 10 mon=1,12
          call Getvbc (iobc, mon, sstobs(1,1,mon) )
10      continue
c
        close(iobc,status='KEEP')
#   ifdef bctest
c
c       call subrountine to compute zonal averages for a quick check
c       and set flag used for land values
c
        flagln = -99.99
c
        call testbc (sstobs, lsts, flag)
#   endif
c
c       get sea surface salinities, one month at a time
c
        open (iobc,form='UNFORMATTED',
     $       file = './sfcsalt')
c
        do 20 mon=1,12
          call Getvbc (iobc, mon, salobs(1,1,mon) )
          do 15 jrow=1,jmt
            do 12 i=1,imt
              salobs(i,jrow,mon) = (salobs(i,jrow,mon)-35.0) / 1000.0
12          continue
15        continue
20      continue
c
        close(iobc,status='KEEP')
#   ifdef bctest
c
c       call subrountine to compute zonal averages for a quick check
c       and set flag used for land values
c
        flagln = -99.99
c
        call testbc (salobs, lsts, flag)
#   endif
c
c       get sea x-component of wind stress, one month at a time
c
        open (iobc,form='UNFORMATTED',
     $       file = './taux')
c
        do 30 mon=1,12
          call Getvbc (iobc, mon, wsxobs(1,1,mon) )
30      continue
c
        close(iobc,status='KEEP')
#   ifdef bctest
c
c       call subrountine to compute zonal averages for a quick check
c       and set flag used for land values
c
        flagln = -99.99
c
        call testbc (wsxobs, lsuv, flag)
#   endif
c
c       get sea y-component of wind stress, one month at a time
c
        open (iobc,form='UNFORMATTED',
     $       file = './tauy')
c
        do 40 mon=1,12
          call Getvbc (iobc, mon, wsyobs(1,1,mon) )
40      continue
c
        close(iobc,status='KEEP')
#   ifdef bctest
c
c       call subrountine to compute zonal averages for a quick check
c       and set flag used for land values
c
        flagln = -99.99
c
        call testbc (wsyobs, lsuv, flag)
#   endif
        write(stdout,9101) j
9101    format( '  === read in obs since first j = ',i8)
#  endif
#  ifdef seasonal
c
c       compute annual mean windstresses (are sometimes used for
c       external mode forcing when DTTS/DTUV is large)
c
        r12 = c1 / 12.0
c
        do 55 jj=1,jmt
          do 50 i=1,imt
            wsxavg(i,jj) = c0
            wsyavg(i,jj) = c0
50        continue
55      continue
c
#  ifdef cpldleg3 
cc   Annual mean wind stress is not used in current versions of the
cc    coupled model. Store the latest bco data instead.
c
        do 155 jj=1,jmt
          do 150 i=1,imt
            wsxavg(i,jj) = bco(i,jj,4,1)
            wsyavg(i,jj) = bco(i,jj,5,1)
150       continue
155     continue
c
#  endif
#  ifdef cybervbc
        do 80 mon=1,12
          do 70 jj=1,jmt
            do 60 i=1,imt
              wsxavg(i,jj) = wsxavg(i,jj) + wsxobs(i,jj,mon)*r12
              wsyavg(i,jj) = wsyavg(i,jj) + wsyobs(i,jj,mon)*r12
60          continue
70        continue
80      continue
#  endif
#  ifdef cpldleg2
c   compute annual mean wind stress for use in subr clinic
        call getwsa
#  endif
c
        write(6,9102)  mbhd, mahd, wmbhd, wmahd
9102    format(/' === running with seasonal surface bc'/
     $     ' mbhd, mahd, wmbhd, wmahd =',2i4,2(1x,f9.6))
#  else
        write(6,9102) mon
9102    format(/' === running with perpetual month ',i2,' sfcbc'/)
#  endif
        endif
      endif
c
c   Reset "j" from the argument list in case it was clobbered above.
      j = jarg
c
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     begin the code that is executed for each row each timestep
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
# if !defined cpldleg3 && !defined cpldleg2
c      realize a heat flux of 100 langleys per day per deg C
c                      ( -or- 100 cal per sq cm per day per deg C)
c      of air-sea temperature difference
c
c      dT/dt = {100.0/[86400.*dzt(1)]} * (Tair - Tsea)
c      ...here, replace Tair with obs SST...
c      after Haney, JPO, (1971)
# endif
c      use "ramp" to start the model run more gently for heat, salt
c      and momentum fluxes
c      ittrmp = number of timesteps that will pass before ramp = 1.0
c----------------------------------------------------------------------
c
c   Note that dzt(1) acts as a time constant here: 50 days for 50m depth,
c    32 days for 32m depth, etc.
      gamma = 100.0/(8.64e4 * dzt(1))
      ittrmp = 720
#ifdef ocnr30
c   For "r30" model, increase ittrmp 
      ittrmp = ittrmp * 3
#endif
      ramp = c1
      if (itt .lt. ittrmp) then
#ifdef noramp
c     for testing purposes, ignore the ramp-up for the sfc boundary
        ramp = c1  
#else
cc 
        fitt = float(itt)
cc 
        ramp = p5 *(c1 + cos(pi * (c1 + (fitt / float(ittrmp)))))
#endif
        gamma = gamma * ramp
      endif
c
#ifdef cpldleg2
      if (first .or. newmon.eq.1) then
#endif
#ifdef cpldleg3
      if (loop1 .or. newmon.eq.1) then
#endif
#if defined cpldleg2 || defined cpldleg3
       if (.not. mxpas2) then
c   compute area integrals for consistency check:
        if (j .eq. jsrow+1) then
          do 82 nri=1,nbo
            tbco(nri) = c0
82        continue
# ifdef cpldleg2
          do 84 nri=1,3
            tobs(nri) = c0
84        continue
# endif
# if defined cpldleg3 && !defined cpld_restore
          tterr = c0
          tserr = c0
# endif
          sumt = c0
          sumu = c0
        endif
c
        fxt = cst(j) * dxt(2) * dyt(j)
        fxu = csu(j) * dxu(2) * dyu(j)
        do 86 i=2,imt-1
          tbco(1) = tbco(1) + bco(i,j,1,1) * fm(i,1,jc) * fxt
          tbco(2) = tbco(2) + bco(i,j,2,1) * fm(i,1,jc) * fxt
          tbco(3) = tbco(3) + bco(i,j,3,1) * fm(i,1,jc) * fxt
          tbco(4) = tbco(4) + bco(i,j,4,1) * gm(i,1,jc) * fxu
          tbco(5) = tbco(5) + bco(i,j,5,1) * gm(i,1,jc) * fxu
# ifdef cpldleg2
          tobs(1) = tobs(1) + sstobs(i,j,mahd) * fm(i,1,jc) * fxt
          tobs(2) = tobs(2) + salobs(i,j,mahd) * fm(i,1,jc) * fxt
          tobs(3) = tobs(3) + yceobs(i,j,mahd) * fm(i,1,jc) * fxt
# endif
# if defined cpldleg3 && !defined cpld_restore
          tterr = tterr + terr(i,j,1) * fm(i,1,jc) * fxt
          tserr = tserr + serr(i,j,1) * fm(i,1,jc) * fxt
# endif
          sumt = sumt + fm(i,1,jc) * fxt
          sumu = sumu + gm(i,1,jc) * fxu
86      continue
c
        if ( j .eq. jerow ) then
          tbco(1) = tbco(1) / sumt
          tbco(2) = tbco(2) / sumt
          tbco(3) = tbco(3) / sumt
          tbco(4) = tbco(4) / sumu
          tbco(5) = tbco(5) / sumu
# ifdef cpldleg2
          tobs(1) = tobs(1) / sumt
          tobs(2) = tobs(2) / sumt
          tobs(3) = tobs(3) / sumt
# endif
# if defined cpldleg3 && !defined cpld_restore
          tterr = tterr / sumt
          tserr = tserr / sumt
# endif
          write(stdout,9694) itt, mahd, (tbco(itb),itb=1,5)
9694      format(1x,'Surface boundary data stored at itt=',i10,
     $              '   Global integrals for mahd=',i3,' are:',/,
     $           5x,'heat p-e+r s+i tx ty=',2x,5e15.7)
# ifdef cpldleg2
          write(stdout,9695) tobs
9695      format(5x,'sstobs salobs yceobs=',2x,3e15.7)
# endif
# if defined cpldleg3 && !defined cpld_restore
          write(stdout,9695) tterr, tserr
9695      format(5x,'           terr serr=',2x,2e15.7)
# endif
        endif
c
       endif
      endif
#endif
c
# if !defined restorst      
c----------------------------------------------------------------------
c     set world ocean surface tracer fluxes off of restoring terms
c     set no flux condition for tracer at bottom.
c----------------------------------------------------------------------
c
#ifdef cybervbc
c     convert gamma so that a flux will result (as oposed to a time
c     rate of change for the entire top layer)
c
      rdays = gamma
      gamma = gamma * dzt(1)
c
cc    do 100 m=1,2
        do 90 i=2,imtm1
#  ifdef seasonal
          stf(i,1) = gamma * ((sstobs(i,j,mahd)*wmahd + 
     $                 sstobs(i,j,mbhd)*wmbhd) - t(i,1,jc,nm,1))
          stf(i,2) = gamma * ((salobs(i,j,mahd)*wmahd + 
     $                 salobs(i,j,mbhd)*wmbhd) - t(i,1,jc,nm,2))
#    ifdef notracer
c     for testing purposes, set surface tracer flux to zero
          stf(i,1) = c0
          stf(i,2) = c0
#    endif
#   else
          stf(i,1) = gamma * ((sstobs(i,j,mon) - t(i,1,jc,nm,1))
          stf(i,2) = gamma * ((salobs(i,j,mon) - t(i,1,jc,nm,2))
#   endif
cc        btf(i,m) = c0
          btf(i,1) = c0
          btf(i,2) = c0
90      continue       
100   continue       
#endif
c
#if defined cpldleg2 || defined cpldleg3
      do 400 i=1,imt
#ifdef cybertra
c  Specify tracera as in cyber code:
        stf(i,3) = gamma * (tra(i)-t(i,1,jc,nm,3)) * dzt(1)
        btf(i,3) = c0
#endif
c
c  Set bottom tracer flux conditions:
        btf(i,1) = c0
        btf(i,2) = c0
400   continue       
#endif
c
#if defined cpldleg3 
      do 500 i=1,imt
        stf(i,1) = bco(i,j,1,1) 
        stf(i,2) = bco(i,j,2,1) * 0.035 
# ifdef seaice
       snowfx(i) = bco(i,j,3,1)
# endif
# ifndef cpld_restore
# ifdef seasonal
c   (terr & serr are in flux form; dzt(1) is not needed here)
        stfadj(i,1) = (terr(i,j,1)*wmahd + 
     $                 terr(i,j,2)*wmbhd)
        stfadj(i,2) = (serr(i,j,1)*wmahd + 
     $                 serr(i,j,2)*wmbhd)
# else
        stfadj(i,1) = terr(i,j,1)
        stfadj(i,2) = serr(i,j,1)
# endif
# ifdef agetracers
        stfadj(i,3) = c0
# endif
# endif
500   continue       
#endif
c
#if defined cpldleg2 || defined cpld_restore
      do 600 i=1,imt
# ifdef seasonal
# ifndef cpld_restore
        stf(i,1) = (bco(i,j,1,1)*wmahd + 
     $              bco(i,j,1,2)*wmbhd) 
        stf(i,2) = (bco(i,j,2,1)*wmahd + 
     $              bco(i,j,2,2)*wmbhd) * 0.035 
# endif
        restr(i,1) = gamma * (sstobs(i,j,mahd)*wmahd + 
     $                        sstobs(i,j,mbhd)*wmbhd 
     $                           - t(i,1,jc,nm,1)) * dzt(1)
        restr(i,2) = gamma * (salobs(i,j,mahd)*wmahd + 
     $                        salobs(i,j,mbhd)*wmbhd 
     $                           - t(i,1,jc,nm,2)) * dzt(1)
#  ifdef seaice
#  ifndef cpld_restore
        snowfx(i) = (bco(i,j,3,1)*wmahd + 
     $               bco(i,j,3,2)*wmbhd) 
#  endif
        adjyc(i)  = - gamma * (yceobs(i,j,mahd)*wmahd + 
     $                         yceobs(i,j,mbhd)*wmbhd 
     $                           - yc(i,j,nm))
#  endif
# else
        stf(i,1) = bco(i,j,1,1) 
        stf(i,2) = bco(i,j,2,1) * 0.035 
        restr(i,1) = gamma * (sstobs(i,j,1)
     $                           - t(i,1,jc,nm,1)) * dzt(1)
        restr(i,2) = gamma * (salobs(i,j,1)  
     $                           - t(i,1,jc,nm,2)) * dzt(1)
#  ifdef seaice
        snowfx(i) = bco(i,j,3,1)
        adjyc(i)  = - gamma * (yceobs(i,j,1)
     $                           - yc(i,j,nm))
#  endif
# endif
600   continue       
c
      if (.not. mxpas2) then
c   Store FLUXES of terr (cal/cm2/sec) and serr (gm/cm2/sec)
        do 610 i=1,imt
          terr(i,j) = terr(i,j) + restr(i,1)
          serr(i,j) = serr(i,j) + restr(i,2)            
c   Store terrj and serrj for history tape, current timestep only
c    (also updated in ice.F)
          terrj(i)  = restr(i,1)
          serrj(i)  = restr(i,2)            
610     continue
c
        if (j .eq. jerow) then
          icntts = icntts + 1
        endif
c   save terr and serr on the last timestep of the month in subr step.
c
# ifdef testterr
c   also print terr, serr in ice.F after adding ice effects.
c   note that terr and serr are summed arrays here.
        if (j .eq. 12) then
          write(stdout,9604) itt, icntts, 
     $               terr(12,12), serr(12,12), terrj(12), serrj(12)
     $              ,adjyc(12)
9604      format(1x,'in setvbc: at (12,12), itt,icntts,terr,serr=',
     $              2i10,5e13.5)
        endif
# endif
      endif
c
# endif
#   ifdef agetracers
c
c     for the idealized age tracer, set the surface flux equal to that
c     needed to instantaneously reset the uppermost layer values to zero
c     Note: this means the flux is dependent upon the depth of the top
c           layer and the timestep length
c     also, set the bottom flux to zero
c
      m = ntage1
      do 110 i=2,imtm1
        stf(i,m) = -t(i,1,jc,nm,m)*dzt(1)/c2dtts
        btf(i,m) = c0
110   continue
#   endif
#   if defined cfctracers || defined ccl4tracer
cdemo(
c (Note: the do 133 loop is for demostration purposes only...
c  in the absense of wind speeds being carried in this test case,
c  a simple latitude dependent wind speed field will be used.
c  for the equations of liss and merlivat are to be applied, scalar
c  wind speed data must be available, from either observations or
c  the atmospheric component of a coupled model.)
      do 133 i=1,imt
        wspd(i) = abs(5.0*cst(j)) + abs(3.5*sine(j))
133   continue
cdemo)
c
#   endif
#   ifdef ccl4tracer
c
c     simulate the air-sea flux of carbon tetrachloride
c
c     1st, determine the current atmospheric mixing ratio
c     (varies with time, but not latitude)
c     tryr = year in the range 1910.0 to 1990.99 (see "tmngr.F")
c     tryr0  = variable "curtim" of model integrations corresponding
c             to 1 Jan 1910 (see "ctmngr.h")
c     cl4air = atmospheric mixing ratios of carbon tetrachloride
c              [parts per trillion volume]
c
      call ctcair (tryr , cl4air)
c
c     if (j.eq. 12) write(stdout,1959) j, tryr, cl4air(44)
c1959  format(' in setvbc: j, tryr, cl4air(44)=',i7,2x,2e14.6)
c
c     the solubility of carbon tetrachloride in seawater for the
c     appropriate temperatures and salinities (after Hunter-Smith, 
c     Balls and Liss)   (here, model predicted level one t,s values 
c     are used ... one could also consider using observed ssts and
c     surface salinities instead)
c
      call CTCsol (t(1,1,jc,nc,1), t(1,1,jc,nc,2), fm(1,1,jc),
     $             CL4sol)
c
c      if (j.eq.12) write(stdout,1960) cl4sol(45) 
c1960  format(' cl4sol(45,12)=',e14.6)
c
c     next, approximate the transfer velocity for CCl4
c     (here, the relations of Liss and Merlivat are used, requiring
c     Schmidt numbers to be estimated and then used in the L & M
c     equations with scalar wind speeds.  one could also consider other
c     functional relationships or simply set the transfer velocities to
c     represent a restoring timescale for the uppermost level)
c
c     get temperature dependent schmidt numbers for CCl4
c
      call CTCsc (t(1,1,jc,nc,1), fm(1,1,jc), CCL4sc)
c
c      if (j.eq.12) write(stdout,1961) ccl4sc(45) 
c1961  format(' ccl4sc(45,12)=',e14.6)
c
c     get wind speed and schmidt number dependent transfer velocities
c     for CCl4
c
      call lismer (wspd, CCL4sc, fm(1,1,jc), CL4tvl)
c
c      if (j.eq.12) write(stdout,1962)  wspd(45), cl4tvl(45) 
c1962  format(' in setvbc: wspd(45), cl4tvl(45,12)=',2e14.6)
c
c
      do 170 i=2,imtm1
        stf(i,ntccl4) = CL4tvl(i)*
     $                 (CL4sol(i)*CL4air(i) - t(i,1,jc,nc,ntccl4))
        btf(i,ntccl4) = c0
170   continue
c
c
#   endif
#   ifdef cfctracers
c
c     simulate the air-sea flux of cfc-11 and cfc-12
c
c     1st, determine the current atmospheric mixing ratios of the two
c     compounds (vary in time and with latitude)
c     tryr = year in the range 1910.0 to 1990.99 (see "tmngr.F")
c     tryr0  = variable "curtim" of model integrations corresponding
c              to 1 Jan 1910 (see "ctmngr.h")
c     f11air, f12air = atmospheric mixing ratios of cfc-11 and cfc-12
c
      call cfcair (yt(j), tryr, acfc11, acfc12)
c
      do 130 i=1,imt
        f11air(i) = acfc11
        f12air(i) = acfc12
130   continue
c
c     2nd, get the solubility of the cfc's in seawater ("f11sol" &
c     "f12sol")for the appropriate temperatures and salinities 
c     (after Warner & Weiss)
c     (here, model predicted level one t,s values are used ...
c     one could also consider using observed ssts and surface 
c     salinities instead)
c
      call cfcsol (t(1,1,jc,nc,1), t(1,1,jc,nc,2), fm(1,1,jc),
     $             f11sol, f12sol)
c
c     3rd, approximate the transfer velocities for cfc-11 and cfc-12
c     (here, the relations of Liss and Merlivat are used, requiring
c     Schmidt numbers to be estimated and then used in the l & M
c     equations with scalar wind speeds.  one could also consider other
c     functional relationships or simply set the transfer velocities to
c     represent a restoring timescale for the uppermost level)
c
c     get temperature dependent schmidt numbers ("f11sc" & "f12sc")
c     for cfc-11 and cfc-12
c
      call cfcsc (t(1,1,jc,nc,1), fm(1,1,jc), f11sc, f12sc)
c 
c  get wind speed and scmidt number dependent transfer velocities
c  for cfc-11 and cfc-12
c
c
      call lismer (wspd, f11sc, fm(1,1,jc), f11tvl)
      call lismer (wspd, f12sc, fm(1,1,jc), f12tvl)
c
      do 140 i=2,imtm1
        stf(i,ntcfc1) = f11tvl(i)*
     $                 (f11sol(i)*f11air(i) - t(i,1,jc,nc,ntcfc1))
        stf(i,ntcfc2) = f12tvl(i)*
     $                 (f12sol(i)*f12air(i) - t(i,1,jc,nc,ntcfc2))
        btf(i,ntcfc1) = c0
        btf(i,ntcfc2) = c0
c      if ( i.eq.44 .and. j.eq.11) then
c        write(6,1319) ntcfc2, f12air(i), f12sol(i), f12tvl(i)
c1319    format('  ntcfc2, f12air(i), f12sol(i), f12tvl(i)=',i5,3e14.6)
c      endif
140   continue
c
#   endif
# endif
c
#ifdef restorst
c
c----------------------------------------------------------------------
c     restore surface tracers to prescribed values via
c     newtonian damping term with restoring time scale of "rests" days
c----------------------------------------------------------------------
c
cmom     gamma = secday/rests
cc    do 200 m=1,2  
        do 190 i=2,imtm1
# ifdef seasonal
          restr(i,1) = gamma * ((sstobs(i,j,mahd)*wmahd +
     $                 sstobs(i,j,mbhd)*wmbhd) - t(i,1,jc,nm,1))
          restr(i,2) = gamma * ((salobs(i,j,mahd)*wmahd +
     $                 salobs(i,j,mbhd)*wmbhd) - t(i,1,jc,nm,2))
# else
          restr(i,1) = gamma * ((sstobs(i,j,mon) - t(i,1,jc,nm,1))
          restr(i,2) = gamma * ((salobs(i,j,mon) - t(i,1,jc,nm,2))
# endif
cmom          restr(i,m) = gamma * (star(m)-t(i,1,jc,nm,m))
190     continue
200   continue
#endif
c
c----------------------------------------------------------------------
c     set surface & bottom stresses
c----------------------------------------------------------------------
c
      do 300 i=1,imt
        kz = kmu(i,j)
#ifdef seasonal
# ifdef noclinic
        smf(i,1) = c0   
        smf(i,2) = c0   
        smfavg(i,1) = c0
        smfavg(i,2) = c0
cccc    if (j.eq.10 .and. i.eq. 30) smf(i,1) = 1.0e-5
# else
#  ifdef cybervbc
        smf(i,1) = (wsxobs(i,j,mahd)*wmahd + wsxobs(i,j,mbhd)*wmbhd)
     $             *ramp
        smf(i,2) = (wsyobs(i,j,mahd)*wmahd + wsyobs(i,j,mbhd)*wmbhd)
     $             *ramp
#  endif
#  ifdef cpldleg2
        smf(i,1) = (bco(i,j,4,1)*wmahd + bco(i,j,4,2)*wmbhd)
     $             *ramp
        smf(i,2) = (bco(i,j,5,1)*wmahd + bco(i,j,5,2)*wmbhd)
     $             *ramp
#  endif
#  ifdef cpldleg3
        smf(i,1) =  bco(i,j,4,1) * ramp
        smf(i,2) =  bco(i,j,5,1) * ramp
#  endif
        smfavg(i,1) = wsxavg(i,j)*ramp
        smfavg(i,2) = wsyavg(i,j)*ramp
# endif
#else
# ifdef noclinic
        smf(i,1) = c0
        smf(i,2) = c0
# else
#  ifdef cybervbc
        smf(i,1) = ramp*wsxobs(i,j,mon)
        smf(i,2) = ramp*wsyobs(i,j,mon)
#  endif
#  if defined cpldleg2 || defined cpldleg3
        smf(i,1) = ramp*bco(i,j,4,1)
        smf(i,2) = ramp*bco(i,j,5,1)
#  endif
# endif
#endif
        if (kz .ne. 0) then
          uvmag    = sqrt(u(i,kz,jc,nm)**2 +v(i,kz,jc,nm)**2)
          bmf(i,1) = cdbot*u(i,kz,jc,nm)*uvmag
          bmf(i,2) = cdbot*v(i,kz,jc,nm)*uvmag
        else
          bmf(i,1) = c0
          bmf(i,2) = c0
        endif
300   continue
c
#endif
c
      return
      end
