      subroutine histwj (j, uext, vext)
      implicit real*8 (a-h,o-z)
c
#if defined worldocn && defined selfdoc
c     NOTE: this code is specifically designed for GFDL global ocean
c           climate studies and produces self documenting output for
c           use with the smmaps / latxht family of graphical analysis
c           programs    ( Keith Dixon - 1993 )
#ifdef cyclic
c     NOTE: this code assumes cyclic boundary conditions
c           a few minor changes would be needed to i loops to remove
c           this restriction
c
c=======================================================================
c
c     prepare and write history data to an internal i/o unit
c
c     "histwj" is called from the 700 j row loop of subroutine step
c
c     for all rows that call "histwj", the history information is
c     written to a sequential "intrnl" unit for later processing
c     and writing to the "iohist" unit.  The reading and processing of
c     the data written to the "intrnl" unit is done in subroutine
c     "histw"
c     (time averaging is a post-processing step).
c
c=======================================================================
c
#include "param.h"
#include "convt.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "emode.h"
#include "grdvar.h"
#include "histjr.h"
#include "iounit.h"
#include "levind.h"
# if defined cpldleg2 || defined cpldleg3
#include "seayce.h"
# endif
#include "slabs.h"
#include "switch.h"
c 
c   variables:
c     "uext" and "vext" are the external mode velocity components
c          calculated in subroutine step
c     "twr" contains the tracers to be written
c     "uwr" and "vwr" contain the full velocity components to be written
c     "wwr" contains the vertical velocity as computed in tracer
c
      dimension uext(imt), vext(imt), fkmt(imt), fkmu(imt)
      dimension twr(imt,km,nt), uwr(imt,km), vwr(imt,km), wwr(imt,km)
c
c
c=======================================================================
c
c---------------------------------------------------------------------
c     open needed unit(s) the first time "histwj" is called each run
c---------------------------------------------------------------------
c
      if (hist0) then
c
        open (UNIT=intrnl,FILE='internal', FORM='UNFORMATTED',ERR=998)
        write(6,9001) intrnl
c
        write(6,9001) iohist
c
        jhists = j
        hist0 = .false.
        goto 50
c
998     write(stdout,9990) intrnl
        stop ' histwj'
9990    format(/' error => problem opening unit ',i3,' in histwj')
50      continue
c
      endif
c
c---------------------------------------------------------------------
c     prepare data for write out to internal i/o unit
c---------------------------------------------------------------------
c
      do 100 i=1,imt
        fkmt(i) = kmt(i,j)
        fkmu(i) = kmu(i,j)
100   continue
c
      do 200 m=1,nt
        do 190 k=1,km
          do 180 i=1,imt
            twr(i,k,m) = t(i,k,jc,nc,m)
180       continue
190     continue
200   continue
c
      do 300 k=1,km
        do 290 i=1,imt
          uwr(i,k) = u(i,k,jc,nc)
          vwr(i,k) = v(i,k,jc,nc)
          wwr(i,k) = w(i,k)
290     continue
300   continue
c
# ifdef cfctracers
c     calculate column inventories for cfc-11 and cfc-12
c     and load cumulative inputs into 1-D vector
c
      do 400 i=1,imt
        f11inv(i) = c0
        f12inv(i) = c0
        do 390 k=1,km
          f11inv(i) = f11inv(i) + fm(i,k,jc)*t(i,k,jc,nc,ntcfc1)
     $                            *dzt(k)*0.01
          f12inv(i) = f12inv(i) + fm(i,k,jc)*t(i,k,jc,nc,ntcfc2)
     $                            *dzt(k)*0.01
390     continue
400   continue
# endif
c
# ifdef ccl4tracer
c     calculate column inventories for carbon tetrachloride
c     and load cumulative inputs into 1-D vector
c
      do 500 i=1,imt
        CL4inv(i) = c0
        do 490 k=1,km
          CL4inv(i) = CL4inv(i) + fm(i,k,jc)*t(i,k,jc,nc,ntccl4)
     $                            *dzt(k)*0.01
490     continue
500   continue
# endif
c
c---------------------------------------------------------------------
c     write out analysis data to internal i/o unit
c---------------------------------------------------------------------
c
      if (j .eq. jhists) rewind(intrnl)
c
      write(intrnl) fkmt, fkmu, twr, uwr, vwr, wwr
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
      jhiste = j
c
#endif
#endif
      return
9001  format(/,'  *** unit',i3,' opened for analysis data ***',/)
      end
