      subroutine histw (ktau)
c
#if defined worldocn && defined selfdoc
c     NOTE: this code is specifically designed for GFDL global ocean
c           climate studies and produces self documenting output for
c           use with the smmaps / latxht family of graphical analysis
c           programs
c           this version writes out the history info in j-row slabs to
c           reduce the number of i/o operations    ( Keith Dixon - 1993 )
c
#ifdef cyclic
c     NOTE: this code assumes cyclic boundary conditions
c           a few minor changes would be needed to i loops to remove
c           this restriction
c
c=======================================================================
c
c     printout self documenting history data for analysis
c
c     "histw" is called from the 3800 timestep loop of subroutine ocean
c
c     Data written to the "intrnl" unit in "histwj" is read in, arranged
c     into a form used by SMMAPS and other off-line analysis routines,
c     and written in a slab-wise manner to the "iohist" unit.
c
c     before history data from the first j-row is written out,
c     "self-documenting" info is written to the "iohist" unit via a
c     call to "histw0"
c
c     a "cheat sheet" relating variable names and their units of
c     measure is created through calls to routine "setcht and written
c     to unit "iochet", using unit "iodumy" along the way for formatting
c
c=======================================================================
c
#include "param.h"
#include "accel.h"
#include "convt.h"
#include "coord.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "docnam.h"
#include "emode.h"
#include "grdvar.h"
#include "histjr.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef seasonal
#include "seasnl.h"
#endif
#include "slabs.h"
#include "switch.h"
# if defined cpldleg2 || defined cpldleg3
c   this code assumes that the ifdef seaice option is seleted only
c   during leg2 or leg3 of an air sea project
# include "seayce.h"
# endif
c
c-----------------------------------------------------------------------
c
c     NOTE: modifications of this routine should be made only after
c           reviewing the way in which it works, to ensure that
c           understanding that the modified code will work for the
c           various pre-exiting options.
c     explanations of the "iblk" value, how it is computed and the
c     point-wise data arrangement for each i,j point can be found in
c     file "histjr.h"
c
c-----------------------------------------------------------------------
c
c     "uext" and "vext" are the external mode velocity components
c          calculated in subroutine step
c     "twr" contains the tracers to be written
c     "uwr" and "vwr" contain the full velocity components to be written
c     "wwr" contains the vertical velocity as computed in tracer
c     "var" is the history slab for a given j-row, containing info in
c           a form compatible with version of smmaps and other analysis
c           programs
c     "lcheat" is a logical variable used to have the cheat sheet 
c           initialized before handling the first i,j point.
c
      dimension uext(imt), vext(imt), fkmt(imt), fkmu(imt)
      dimension twr(imt,km,nt), uwr(imt,km), vwr(imt,km), wwr(imt,km)
c
      dimension var (iblk,imt-2)
cAltix
      real*4    var4(iblk,imt-2)
cAltix
c
      logical lcheat
c
      character * 16 vname, units
# ifdef sfctopog
      character * 16 chsupp
# endif
c
c=======================================================================
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     generate the variable names and units of measure for this case
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
      rewind (intrnl)
c
c---------------------------------------------------------------------
c     read thru "intrnl" unit row by row, arrange data point by point
c     and within the history slab array named "var", and then write
c     the history info  to unit "iohist" in a slab-wise form 
c     compatible with versions of smmaps, latxht, etc.
c---------------------------------------------------------------------
c
      nvtot = iblk
c
      write(stdout,239) nvtot, jhists, jhiste
239   format(' writing history tape using: iblk, jhists, jhiste =',3i8)
c
      if ((jhiste-jhists+1) .NE. (jmt-2)) then
        print *, jhiste, jhists, (jmt-2), ' PROBLEM'
        stop ' histwr'
      endif
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     read in history data from internal i/o unit one row at a time,
c     then load the data into slab-wise data vectors point-by-point ,
c     and finally write a history slab to the history file for each
c     j-row.
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     read in history data from internal i/o unit one row at a time,
c---------------------------------------------------------------------
c
      do 200 jj=jhists,jhiste
        read (intrnl) fkmt, fkmu, twr, uwr, vwr, wwr
     $             ,uext, vext, stf, smf, flgcnv
# ifdef ccl4tracer
     $             ,CL4air, CL4sol, CCL4sc, CL4tvl, Cl4sat, CL4inv
# endif
# ifdef cfctracers
     $             ,f11air, f12air, f11sol, f12sol, f11sc,  f12sc
     $             ,f11tvl, f12tvl, f11sat, f12sat, f11inv, f12inv
# endif
# if defined cfctracers || defined  ccl4tracer
#  if defined cpldleg2 || defined cpldleg3
     $             ,wspd
#  else
     $             ,wspnow, ycenow
#  endif
# endif
# if defined cpldleg2 || defined cpldleg3
# ifdef seaice
     $             ,snowfx
# endif
# endif
# if defined cpldleg2 || defined cpld_restore
     $             ,restr, adjyc
     $             ,terrj, serrj
# endif
# if defined cpldleg3 && !defined cpld_restore
     $             ,stfadj
# endif
c
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     load the data into slab-wise data vector "var", generating
c     cheat sheet info during loop thru the first i,j point data
c---------------------------------------------------------------------
c
        if (jj.ge.2 .and. jj.le.jmtm1) then
c
          do 190 i=2,imtm1
c
            ii = i-1
c
            if (jj .eq. 2 .and. i .eq. 2) then
c     initialize cheat sheet variables
              lcheat = .TRUE.
              nvchk = 0
              do 110 nv=1,iblk
                varnam(nv) = '    MISSING     '
                varunt(nv) = '    MISSING     '
                varhgr(nv) = 'MM'
                varvgr(nv) = 'MM'
                klevar(nv) =  -9e6
110           continue
            else
              lcheat = .FALSE.
            endif
c
c---------------------------------------------------------------------
c     begin with base variables associated with "iblk0"
c---------------------------------------------------------------------
c
c     the first two locations are always given to the number of ocean
c     levels at that gridpoint on the ts & uv grids (.le. "km")
c
            nv = 1
            if (lcheat) call setcht('fkmt            ',
     $                              'k-levels ts grid',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = fkmt(i)
c
            nv = 2
            if (lcheat) call setcht('fkmu            ',
     $                              'k-levels uv grid',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = fkmu(i)
c
c     the next km*nt positions are filled with the tracers,
c     at each of the "km" levels (including land points).
c
            do 130 m=1,nt
              if (lcheat) then
              if (m.eq.nttemp) then
                  units = 'degrees C       '
                elseif (m.eq.ntsalt) then
                  units = 'diff from 0.0350'
# ifdef agetracers
                elseif (m.eq.ntage1) then
                  units = 'age in years    '
# endif
# ifdef ccl4tracer
                elseif (m.eq.ntccl4) then
                  units = 'picomoles/liter '
# endif
# ifdef cfctracers
                elseif (m.eq.ntcfc1 .or. m.eq.ntcfc2) then
                  units = 'picomoles/liter '
# endif
                else
                  units = ' unknown units  '
                endif
              endif
c
              do 120 k=1,km
                nv = 2+(m-1)*km+k
                if (lcheat) then
                  rewind(iodumy)
                  write(iodumy,9051) trname(m), k
                  rewind(iodumy)
                  read (iodumy,9008) vname 
                  call setcht(vname, units,
     $                        'ts', 'zt', k, nv)
                endif
                var(nv,ii) = twr(i,k,m)
120           continue
130         continue
c
c
c     the next km*3 positions are filled with ocean velocity data,
c     u & v (total velocity) & w (from tracer) at each of the "km"
c     levels (including land points).
c           u component of velocity (total velocity)
            do 140 k=1,km
              nv = 2 + km*nt  + k
              if (lcheat) then
                rewind(iodumy)
                write(iodumy,9101) k
                rewind(iodumy)
                read (iodumy,9008) vname
                call setcht(vname, 'cm per second   ',
     $                      'uv', 'zt', k, nv)
              endif
              var(nv,ii) = uwr(i,k)
140         continue
c
c           v component of velocity (total velocity)
            do 150 k=1,km
              nv = 2 + km*(1+nt) + k
              if (lcheat) then
                rewind(iodumy)
                write(iodumy,9102) k
                rewind(iodumy)
                read (iodumy,9008) vname
                call setcht(vname, 'cm per second   ',
     $                      'uv', 'zt', k, nv)
              endif
              var(nv,ii) = vwr(i,k)
150         continue
c
c           w component of velocity
            do 160 k=1,km
              nv = 2 + km*(2+nt) + k
              if (lcheat) then
                rewind(iodumy)
                write(iodumy,9103) k
                rewind(iodumy)
                read (iodumy,9008) vname
                call setcht(vname, 'cm per second   ',
     $                      'ts', 'zw', k, nv)
              endif
              var(nv,ii) = wwr(i,k)
160         continue
c
c     the next 2 positions are filled with external mode velocities
c
            nv = 2 + km*(nt+3) + 1
            if (lcheat) call setcht('uext extrnl mode',
     $                              'cm per second   ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = uext(i)
c
            nv = 2 + km*(nt+3) + 2
            if (lcheat) call setcht('vext extrnl mode',
     $                              'cm per second   ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = vext(i)
c
c     the next nt positions are filled with the surface tracer fluxes
c
            do 170 m=1,nt
              nv = 2 + km*(nt+3) + 2 + m
              if (lcheat) then
                rewind(iodumy)
                write(iodumy,9151) trname(m)
                rewind(iodumy)
                read (iodumy,9008) vname
                if (m .eq. nttemp) then
                  units = 'cal/(cm2 sec)   '
                elseif (m .eq. ntsalt) then
                  units = 'grams/(cm2 sec) '
                else
                  units = 'surface flux    '
                endif
                call setcht(vname, units,
     $                      'ts', 'zw', 0, nv)
              endif
              var(nv,ii) = stf(i,m)
170         continue
c
c     the next 2 positions are filled with the surface momentum fluxes
c     otherwise referred to as the wind stress components
c
            nv = 2 + km*(nt+3) + 2 + nt + 1
            if (lcheat) call setcht('smf wind stressX',
     $                              'dynes per cm**2 ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = smf(i,1)
c
            nv = 2 + km*(nt+3) + 2 + nt + 2
            if (lcheat) call setcht('smf wind stressY',
     $                              'dynes per cm**2 ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = smf(i,2)
c
c     the next km positions are filled with the convection flags
c
            do 180 k=1,km
              nv = 2 + km*(nt+3) + 2 + nt + 2 + k
              if (lcheat) then
                rewind(iodumy)
                write(iodumy,9161) k
                rewind(iodumy)
                read (iodumy,9008) vname
                call setcht(vname, '0=no conv, 1=con',
     $                      'ts', 'zt', k, nv)
              endif
              var(nv,ii) = flgcnv(i,k)
180         continue
c
c     convert the barotropic streamfunction to units of sverdrups and
c     load it into its spot at the end of the base of the data vector
c     (position iblk0)
c
            nv = iblk0
            if (lcheat) call setcht('streamfunction  ',
     $                              'Sv=10**6 m**3/s ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = p(i,jj,1) * 1.0e-12
c
c---------------------------------------------------------------------
c     done with base variables associated with "iblk0"
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     consider gaseous tracer variables associated with "iblktr"
c---------------------------------------------------------------------
c
# ifdef ccl4tracer
c
c           add in variables associated with carbon tetrachloride
c
            nv = iblk0 + 1
            if (lcheat) call setcht('CCl4 air concent',
     $                              'pp trillion vol ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = CL4air(i)
c
            nv = iblk0 + 2
            if (lcheat) call setcht('CCl4 solubility ',
     $                              'moles/kg/atmos  ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = CL4sol(i)
c
            nv = iblk0 + 3
            if (lcheat) call setcht('CCl4 Schmidt num',
     $                              'dimensionless   ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = CCL4sc(i)
c
            nv = iblk0 + 4
            if (lcheat) call setcht('CCl4 transfer vl',
     $                              'cm per sec      ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = CL4tvl(i)
c
            nv = iblk0 + 5
            if (lcheat) call setcht('CCl4 pct saturtn',
     $                              'percent         ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = CL4sat(i)
c
            nv = iblk0 + 6
            if (lcheat) call setcht('CCl4 cumlt input',
     $                              'nanomoles/m**2  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = cinpct(i,jj)
c
            nv = iblk0 + 7
            if (lcheat) call setcht('CCl4 col inventr',
     $                              'nanomoles/m**2  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = CL4inv(i)
# endif
# ifdef cfctracers   
c
c           add in variables associated with cfc-11 and cfc-12
c
            nv = iblk0 + ibccl4 + 1
            if (lcheat) call setcht('cfc-11 air conc ',
     $                              'pp trillion vol ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = f11air(i)
c
            nv = iblk0 + ibccl4 + 2
            if (lcheat) call setcht('cfc-12 air conc ',
     $                              'pp trillion vol ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = f12air(i)
c
#  ifdef cfcperkg
            units = 'moles/kg/atmos  '
#  else
            units = 'mole/liter/atmos'
#  endif
            nv = iblk0 + ibccl4 + 3
            if (lcheat) call setcht('cfc-11 solubilit',
     $                              units,
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = f11sol(i)
c
            nv = iblk0 + ibccl4 + 4
            if (lcheat) call setcht('cfc-12 solubilit',
     $                              units,
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = f12sol(i)
c
            nv = iblk0 + ibccl4 + 5
            if (lcheat) call setcht('cfc-11 Schmidt N',
     $                              'dimensionless   ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = f11sc (i)
c
            nv = iblk0 + ibccl4 + 6
            if (lcheat) call setcht('cfc-12 Schmidt N',
     $                              'dimensionless   ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = f12sc (i)
c
            nv = iblk0 + ibccl4 + 7
            if (lcheat) call setcht('cfc-11 trans vel',
     $                              'cm per sec      ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = f11tvl(i)
c
            nv = iblk0 + ibccl4 + 8
            if (lcheat) call setcht('cfc-12 trans vel',
     $                              'cm per sec      ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = f12tvl(i)
c
            nv = iblk0 + ibccl4 + 9
            if (lcheat) call setcht('cfc-11 satrtn eq',
     $                              'picomoles/liter ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = f11sat(i)
c
            nv = iblk0 + ibccl4 + 10
            if (lcheat) call setcht('cfc-12 satrtn eq',
     $                              'picomoles/liter ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii)= f12sat(i)
c
            nv = iblk0 + ibccl4 + 11
            if (lcheat) call setcht('cfc-11 cumlt inp',
     $                              'nanomoles/m**2  '
     $                              'ts', 'zw', 0, nv)
            var(nv,ii)= cinp11(i,jj)
c
            nv = iblk0 + ibccl4 + 12
            if (lcheat) call setcht('cfc-12 cumlt inp',
     $                              'nanomoles/m**2  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii)= cinp12(i,jj)
c
            nv = iblk0 + ibccl4 + 13
            if (lcheat) call setcht('cfc-11 col invnt',
     $                              'nanomoles/m**2  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii)= f11inv(i)
c
            nv = iblk0 + ibccl4 + 14
            if (lcheat) call setcht('cfc-12 col invnt',
     $                              'nanomoles/m**2  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii)= f12inv(i)
# endif
# if defined cfctracers || defined ccl4tracer
c
c     variable associated with cfcs or carbon tet tracers
c
            nv = iblk0 + ibccl4 + ibcfcs + 1
            if (lcheat) call setcht('wind speed      ',
     $                              'meters per sec  ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = windsp(i) 
#  if !defined cpldleg2 && !defined cpldleg3
c
c     variable associated with cfcs or carbon tet tracers in versions
c     without prognostic sea ice
c
            nv = iblk0 + ibccl4 + ibcfcs + 2
            if (lcheat) call setcht('sea ice coverage',
     $                              'fractional cover',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = ycenow(i)
#  endif
# endif
c
c---------------------------------------------------------------------
c     done with gaseous tracer variables associated with "iblktr"
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     consider variables associated with "iblkas" ... the air-sea
c     model project
c---------------------------------------------------------------------
c     next, load air-sea variables associated with "iblkas"
c
c
# if defined cpldleg3 || defined cpldleg2
c
c     load in variables associated with sea ice that are found in both
c     the leg 2 and coupled model runs
c
            nv = iblk0 + iblktr + 1
# ifdef seaice
            if (lcheat) call setcht('snow flux       ',
     $                              'cm per sec      ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = snowfx(i)
c
            nv = iblk0 + iblktr + 2
            if (lcheat) call setcht('sea ice         ',
     $                              'cm              ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = yc(i,jj,nc)
c
            nv = iblk0 + iblktr + 3
            if (lcheat) call setcht('stfyc (latnt ht)',
     $                              'cal/(cm2 sec)   ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = stfyc(i,jj,1)
c
            nv = iblk0 + iblktr + 4
            if (lcheat) call setcht('stfyc (salt)    ',
     $                              'grams/(cm2 sec) ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = stfyc(i,jj,2)
c
            nv = iblk0 + iblktr + 5
# endif
            if (lcheat) call setcht('flux adjust (T) ',
     $                              'cal/(cm2 sec)   ',
     $                              'ts', 'zw', 0, nv)
# if defined cpldleg2 || defined cpld_restore
            var(nv,ii) = terrj(i)
# endif
# if defined cpldleg3 && !defined cpld_restore
            var(nv,ii) = stfadj(i,1)
# endif
c
            nv = nv + 1
            if (lcheat) call setcht('flux adjust (S) ',
     $                              'grams/(cm2 sec) ',
     $                              'ts', 'zw', 0, nv)
# if defined cpldleg2 || defined cpld_restore
            var(nv,ii) = serrj(i)
# endif
# if defined cpldleg3 && !defined cpld_restore
            var(nv,ii) = stfadj(i,2)
# endif
# endif
# if defined cpldleg2 || defined cpld_restore
c
c     load in variables associated with the leg2 ocean only run
c
            nv = iblk0 + iblktr + iblegs + 1
            if (lcheat) call setcht('restor heat flux',
     $                              'cal/(cm2 sec)   ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = restr(i,1)
c
            nv = iblk0 + iblktr + iblegs + 2
            if (lcheat) call setcht('restor salt flux',
     $                              'grams/(cm2 sec) ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = restr(i,2)
c
            nv = iblk0 + iblktr + iblegs + 3
            if (lcheat) call setcht('restor ice rate ',
     $                              'cm per sec      ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = adjyc(i)
c
            nv = iblk0 + iblktr + iblegs + 4
            if (lcheat) call setcht('dtyce (t-restor)',
     $                              'cal/(cm2 sec)   ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = dtyce(i,jj,1)
c
            nv = iblk0 + iblktr + iblegs + 5
            if (lcheat) call setcht('dtyce (s-restor)',
     $                              'grams/(cm2 sec) ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = dtyce(i,jj,2)
c
            nv = iblk0 + iblktr + iblegs + 6
            if (lcheat) call setcht('observed sst    ',
     $                              'degrees C       ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = wmahd*sstobs(i,jj,mahd) +
     $                   wmbhd*sstobs(i,jj,mbhd)
c
            nv = iblk0 + iblktr + iblegs + 7
            if (lcheat) call setcht('obs sfc salinity',
     $                              'parts/1000      ',
     $                              'ts', 'zt', 1, nv)
            var(nv,ii) = wmahd*salobs(i,jj,mahd) +
     $                   wmbhd*salobs(i,jj,mbhd)
c
            nv = iblk0 + iblktr + iblegs + 8
            if (lcheat) call setcht('observed sea ice',
     $                              'cm              ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = wmahd*yceobs(i,jj,mahd) +
     $                   wmbhd*yceobs(i,jj,mbhd)
c
# ifdef seasonal
            nv = iblk0 + iblktr + iblegs + 9
            if (lcheat) call setcht('smfavg wind strX',
     $                              'dynes per cm**2 ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = wsxavg(i,jj)
c
            nv = iblk0 + iblktr + iblegs + 10
            if (lcheat) call setcht('smfavg wind strY',
     $                              'dynes per cm**2 ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = wsyavg(i,jj)
c
# endif
# endif
c
# if defined cpldleg3 && !defined cpld_restore
c
c     load in variables associated with the leg3 coupled model run
c
# endif
c
c---------------------------------------------------------------------
c     done with variables associated with "iblkas" ... the air-sea
c     model project
c---------------------------------------------------------------------
c
c---------------------------------------------------------------------
c     consider extra history variables associated with "iblkxt"
c     (those that do not fall in the above categories)
c---------------------------------------------------------------------
c
# ifdef sfctopog
            nv = iblk0 + iblktr + iblkas + 1
            if (lcheat) call setcht('zu vert avg frcg',
     $                              'cm per second**2',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = zu(i,jj)
c
            nv = iblk0 + iblktr + iblkas + 2
            if (lcheat) call setcht('zv vert avg frcg',
     $                              'cm per second**2',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = zv(i,jj)
c
            nv = iblk0 + iblktr + iblkas + 3
            if (lcheat) call setcht('ptd  d psi / d t',
     $                              'Sv / second     ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = ptd(i,jj) * 1.0e-12 / c2dtsf 
c
            nv = iblk0 + iblktr + iblkas + 4
            if (lcheat) call setcht('psx (sfc pr grd)',
     $                              'cgs units       ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = psx(i,jj)
c
            nv = iblk0 + iblktr + iblkas + 5
            if (lcheat) call setcht('psy (sfc pr grd)',
     $                              'cgs units       ',
     $                              'uv', 'zw', 0, nv)
            var(nv,ii) = psy(i,jj)
c
            nv = iblk0 + iblktr + iblkas + 6
            if (lcheat) call setcht('sealev divxdepth',
     $                              'cgs             ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = divslh(i,jj)
c
            nv = iblk0 + iblktr + iblkas + 7
            if (lcheat) call setcht('sea level topog ',
     $                              'cm              ',
     $                              'ts', 'zw', 0, nv)
            var(nv,ii) = sealev(i,jj)
c
# endif
c
c---------------------------------------------------------------------
c     done with extra history variables associated with "iblkxt"
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
            if (lcheat) then
c---------------------------------------------------------------------
c     Writing out self documenting info to "iohist" via "histw0"
c     and print out list of history variables for "cheat sheet"
c---------------------------------------------------------------------
c
#ifdef gtime
              cmonth = months(icmonth)
#else
              cmonth = months(int(tmonth))
#endif
c
c     set info about 1st and last timestamps used in averaging
c     the history data, and the number of instantaneous snapshots
c     that were averaged to compute the mean field
c     (here, no averaging has been done yet, so "knthst" is set to 1
c     and the first and last timestamps are the same.
c
              stampe = stamp
              stampb = stamp
              knthst = 1
c
c
              call histw0
c
              if (lcheat1st .AND. cmonth .EQ. 'jan.') then
                open (iochet, file='cheat_sheet', access='SEQUENTIAL',
     $               form='FORMATTED', position='REWIND')
                write(stdout,9189)
                write(iochet,9189)
                do 189 nv=1,iblk
                  write(iochet,9199) nv, varnam(nv), varunt(nv),
     $                               varhgr(nv), varvgr(nv), klevar(nv)
189             continue
                lcheat1st = .FALSE.
                close (iochet, STATUS='KEEP')
              endif
c
            endif
c
190       continue
c
c
c---------------------------------------------------------------------
c     check that data to be written out can be converted to IEEE without
c     underflow or overflow conversion problems
c---------------------------------------------------------------------
c
          nums = iblk*(imt-2)
c
          call ieeeok (var, nums, .true.)
c
c---------------------------------------------------------------------
c     write the history data slab for the current j-row to unit "iohist"
c---------------------------------------------------------------------
c
cAltix
          var4 = var
          write(iohist) var4
ccc       write(iohist) var
cAltix
c
        endif
c
200   continue
c
c---------------------------------------------------------------------
c
      write(stdout,9002) iohist, ktau, nvtot
     $                  ,month
c
c---------------------------------------------------------------------
c
#endif
#endif
c
      return
c
9002  format(1x,'*** history data written for analysis to unit',i3,
     $          ' on itt=',i8,'   iblk=',i4,'  month=',i3,/)
9008  format(a16)
9051  format(a12,'(',i2,')')
9101  format('u(',i2,') total vel ',1x,a16)
9102  format('v(',i2,') total vel ',1x,a16)
9103  format('w(',i2,') box bottom',1x,a16)
9151  format('stf ',a12)
9161  format('convectn flag',i3)
9189  format(/,1x,'List of history tape variables')
9199  format(1x,i4,2x,a16,2x,a16,2x,a2,2x,a2,2x,i3)
c
      end
