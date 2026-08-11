      subroutine step (ntask)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     step is called once per timestep. it initializes various
c     quantities and cycles data from tau & tau-1 discs into the memory
c     slab window a row at a time. after prognostic variables are
c     predicted for each row, it writes them out to the tau+1 disc.
c
c=======================================================================
c
c
#include "param.h"
#include "cdiag.h"
#ifdef isopycmix
#include "cisop.h"
#endif
#include "crelax.h"
#include "ctask.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "cvmix.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#include "slabs.h"
#ifdef saltnorm
# include "conblk.h"
# include "sltdat.h"
      equivalence (sltsum, konblk(100))
#endif
#include "switch.h"
#include "temp.h"
c
      dimension ubar(imt), vbar(imt)
c
c-----------------------------------------------------------------------
c     set up the starting & ending j rows for this task
c-----------------------------------------------------------------------
c
      if (ntask .ne. 0) then
        jsrow =  jsrowb(ntask)
        jerow =  jerowb(ntask)
      else
        jsrow =  1
        jerow =  jmt-1
      endif
c
c-----------------------------------------------------------------------
c     queue up disk reads
c-----------------------------------------------------------------------
c
      do 100 j=jsrow,jerow
        if (mixts) then
          call ofind(labs(ndisk),nslab,(j-1)*nslab+1)
        else
          call ofind(labs(ndiskb),nslab,(j-1)*nslab+1)
        endif
        call ofind(labs(ndisk),nslab,(j-1)*nslab+1)
100   continue
c
#ifdef biharmonic
c
c-----------------------------------------------------------------------
c     read  row  "jsrow-1", "jsrow" & "jsrow+1" variables from disk into
c     the memory slab window. save internal modes, construct masks &
c     total velocity for tau & tau-1
c-----------------------------------------------------------------------
c
      jm  = mod(jsrow-1,nslabs) + 1
      jc  = mod(jsrow,  nslabs) + 1
      jp  = mod(jsrow+1,nslabs) + 1
c
# ifdef worldocn
      call getvar (max(1,jsrow-1),jm,ubar,vbar)
      call getvar (jsrow,  jc,ubar,vbar)
      call getvar (jsrow+1,jp,ubar,vbar)
# else
      call getvar (max(1,jsrow-1),jm)
      call getvar (jsrow,  jc)
      call getvar (jsrow+1,jp)
# endif
c
c-----------------------------------------------------------------------
c     compute del**2 quantities for "jsrow"
c-----------------------------------------------------------------------
c
      call delsqr (jsrow, jm, jc, jp)
#else
c
c-----------------------------------------------------------------------
c     read  row "jsrow" variables from disk into the memory slab window.
c     save internal modes, construct masks & total velocity for tau &
c     tau-1
c-----------------------------------------------------------------------
c
      jc = mod(jsrow,nslabs) + 1
c
# ifdef worldocn
      call getvar (jsrow, jc, ubar, vbar)
# else
      call getvar (jsrow, jc)
# endif
#endif
c
#ifdef saltnorm
c-----------------------------------------------------------------------
c     set salt normalization variables
c-----------------------------------------------------------------------
c
      sltsmj = c0
cc    totslv = totslt * volume
cc    sltnrm = (totslv - sltsum) / volume
      sltnrm =  totslt - (sltsum / volume)
      if (sltsum .eq. c0) sltnrm = c0
cc    write(stdout,1028) itt, sltnrm, totslt, sltsum, volume
1028  format(1x,' at itt',i7,
     $       2x,' sltnrm, totslt, sltsum, volume=',4e18.8)
#endif
c
c=======================================================================
c=======================================================================
c
c     begin row-by-row computation of prognostic variables within task
c
c=======================================================================
c=======================================================================
c
      do 700 j=jsrow,jerow
c
c-----------------------------------------------------------------------
c       set the memory slab window pointers (indices)
c       jm = southern slab (row j-1, m stands for minus one)
c       jc = central slab  (row j)
c       jp = northern slab (row j+1, p stands for plus one)
c-----------------------------------------------------------------------
c
        jm = mod(j-1,nslabs) + 1
        jc = mod(j  ,nslabs) + 1
        jp = mod(j+1,nslabs) + 1
#ifdef biharmonic
c
c-----------------------------------------------------------------------
c       set additional memory slab window pointer (index)
c       jpp = northern most slab (row j+2, p stands for plus one)
c       read row j+2 variables from disk into the memory slab window.
c       save internal modes, construct masks & total velocity for tau &
c       tau-1 then ... compute del**2 quantities for row "j+1"
c-----------------------------------------------------------------------
c
        jpp = mod(j+2,nslabs) + 1
c
# ifdef worldocn
        call getvar (j+2, jpp, ubar, vbar)
# else
        call getvar (j+2, jpp)
# endif
        call delsqr (j+1, jc, jp, jpp)
#else
c
c-----------------------------------------------------------------------
c       read row j+1 variables from disk into the memory slab window.
c       save internal modes, construct masks & total velocity for tau &
c       tau-1)
c-----------------------------------------------------------------------
c
# ifdef worldocn
        call getvar (j+1, jp, ubar, vbar)
# else
        call getvar (j+1, jp)
# endif
#endif
c
c-----------------------------------------------------------------------
c       initiate writeout of newly computed data from previous row
c       should a copy of the data be saved for off line analysis also?
c-----------------------------------------------------------------------
c
        if (j .gt. jsrow + 1) then
c
          call oput (labs(ndiska), nslab, (j-2)*nslab+1, ta)
c
          if (eots .and. snapts) then
            if (j .eq. jsrow + 2) then
              write (stdout,920)
     $         '==>  Writing  a snapshot of MOM"s data on unit '
     $,                    iosnap, itt, stamp
            endif
            write (iosnap) smf, stf, ta, ua, va
#ifdef tcvmix
     $,                    vvca, vdca, vdqca, q2a
# ifdef leq
     $,                    q2la
# endif
#endif
          endif
        endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c       set symmetry boundary conditions on last row
c-----------------------------------------------------------------------
c
        if (j .eq. jmtm1) then
          do 500 m=1,nt
            do 490 k=1,km
              do 480 i=1,imt
                t(i,k,jp,nm,m) = t(i,k,jc,nm,m)
                t(i,k,jp,nc,m) = t(i,k,jc,nc,m)
480           continue
490         continue
500       continue
#ifdef seaice
          do 510 i=1,imt
            yc(i,jp,nm) = yc(i,jc,nm)
            yc(i,jp,nc) = yc(i,jc,nc)
510       continue
#endif
        endif
#endif
c
c-----------------------------------------------------------------------
c       set all source terms
c-----------------------------------------------------------------------
c
        if (j .ne. jsrow) then
          do 660 k=1,km
            do 650 i=1,imt
              sourcu(i,k) = c0
              sourcv(i,k) = c0
650         continue
660       continue
          do 690 m=1,nt
            do 680 k=1,km
              do 670 i=1,imt
                sourct(i,k,m) = c0
670           continue
680         continue
690       continue
        endif
c
c-----------------------------------------------------------------------
c       set vertical boundary conditions (surface & bottom)
c-----------------------------------------------------------------------
c
        if (j .ne. jsrow) call setvbc (j)
c
c-----------------------------------------------------------------------
c       set vertical mixing coefficients
c-----------------------------------------------------------------------
c
#ifdef ppvmix
        call ppmix (j)
#endif
#if defined constvmix || defined blvmix
        if (j .ne. jsrow) call cnvmix (j)
#endif
#ifdef tcvmix
        call tcmix(j)
#endif
#ifdef isopycmix
c
c       calculate slopes and set mixing coefficients
c       to orient diffusion of tracers along isopycnal surfaces
c
        call isop0 (j)
#endif
#ifdef testcfl
c
c-----------------------------------------------------------------------
c       monitor velocities at current time level "nc"
c-----------------------------------------------------------------------
c
# ifdef worldocn
        if (eots .and. diagts) call cfl (j)
# else
        call cfl (j)
# endif
#endif
c
c-----------------------------------------------------------------------
c       calculate internal mode velocities
c-----------------------------------------------------------------------
c
#ifdef symmetry
        call clinic (j)
#else
        if (j .ne. jmtm1) call clinic (j)
#endif
c
c-----------------------------------------------------------------------
c       calculate tracers
c-----------------------------------------------------------------------
c
        call tracer (j)
c
c-----------------------------------------------------------------------
c       calculate & print diagnostics on diagnostic timesteps
c-----------------------------------------------------------------------
c
        if (diagts .and. eots) call diag (j)
c
#ifdef worldocn
c-----------------------------------------------------------------------
c     write history for analysis when histts is 'true'.
c-----------------------------------------------------------------------
c
      if (histts .and. eots) then
        if (j .gt. 1) call histwj (j, ubar, vbar)
      endif
#endif
c
700   continue
c
#ifdef saltnorm
c-----------------------------------------------------------------------
c     set the total salt variable
c-----------------------------------------------------------------------
c
      sltsum = sltsmj
c
#endif
c
#if defined cpldleg2 || defined cpld_restore
c-----------------------------------------------------------------------
c   save terr and serr on the last timestep of the month
c-----------------------------------------------------------------------
      if (endmon .and. .not. mxpas2) then
ccccc   rewind iobca
c    write only the current month of terr & serr
          if (icntts .gt. 0) then
            rcnt = c1/float(icntts)
            do 625 jj=1,jmt
              do 625 i=1,imt
                terr(i,jj) = terr(i,jj) * rcnt
                serr(i,jj) = serr(i,jj) * rcnt
625         continue
            write (iobca) terr
            write (iobca) serr
          else
            write(6,9609)
9609        format(///,1x,'icntts = 0. Stop the run in step.')
            stop
          endif
        write(stdout,9620) itt, month, icntts
9620    format(1x,'Write surface flux adjustment terms (terr, serr)'
     $            ' at itt=',i10, '  for month=',i3,' icntts=',i10)
        icntts = 0
        do 645 jj=1,jmt
          do 645 i=1,imt
            terr(i,jj) = c0
            serr(i,jj) = c0
645     continue
      endif
#endif
c
c-----------------------------------------------------------------------
c     initiate writeout of newly computed data from row "jerow"
c     should a copy of the data be saved for off line analysis also?
c-----------------------------------------------------------------------
c
      call oput (labs(ndiska), nslab, (jerow-1)*nslab+1, ta)
c
      if (eots .and. snapts) then
        write (iosnap) smf, stf, ta, ua, va
#ifdef tcvmix
     $,                vvca, vdca, vdqca, q2a
# ifdef leq
     $,                q2la
# endif
#endif
      endif
      return
920   format (/1x,a47,i3,', ts = ',i7,1x, a32/)
      end
#ifdef worldocn
      subroutine getvar (jrowt, jptr, ubar, vbar)
#else
      subroutine getvar (jrowt, jptr)
#endif
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     read prognostic variables from row "jrowt" on disk into the
c     memory slab window at position "jptr". set the masks for "jrowt" &
c     save a copy of the internal modes before constructing the full
c     velocity.
c
c=======================================================================
c
#include "param.h"
#include "ctask.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "slabs.h"
#include "switch.h"
c
      dimension ubar(imt), vbar(imt)
c
c     limit size of jrow
c
      jrow = min(jrowt,jmt)
c
c-----------------------------------------------------------------------
c     read slabs from row "jrow" on disk into the memory slab window
c     at row "jptr"
c-----------------------------------------------------------------------
c
      if (jrow .lt. jmt) then
        if (mixts) then
          call oget (labs(ndisk), nslab, (jrow-1)*nslab+1, bufsl)
        else
          call oget (labs(ndiskb), nslab, (jrow-1)*nslab+1, bufsl)
        endif
        call xfer (bufsl, t(1,1,jptr,nm,1))
        call oget (labs(ndisk ), nslab, (jrow-1)*nslab+1, bufsl)
        call xfer (bufsl, t(1,1,jptr,nc,1))
      endif
c
c-----------------------------------------------------------------------
c     set masks for row "jrow"
c-----------------------------------------------------------------------
c
      do 100 k=1,km
        do 90 i=1,imt
          if (kmt(i,jrow) .ge. k) then
            fm(i,k,jptr) = c1
          else
            fm(i,k,jptr) = c0
          endif
          if (kmu(i,jrow) .ge. k) then
            gm(i,k,jptr) = c1
          else
            gm(i,k,jptr) = c0
          endif
90      continue
100   continue
c
c-----------------------------------------------------------------------
c     set pointers (indices) to cycle internal modes
c     (also used to cycle del**2 quantities for biharmonic option)
c-----------------------------------------------------------------------
c
      jpt1 = mod (jrowt+1,numjpt) + 1
      jpt2 = mod (jrowt+2,numjpt) + 1
      jpt  = jpt2
#ifdef biharmonic
      jpt3 = mod (jrowt+3,numjpt) + 1
      jpt  = jpt3
#endif
c
c-----------------------------------------------------------------------
c     save a copy of the internal mode velocity from row "jrow"
c     (row "jrow + 1"  if the biharmonic option is enabled) for use
c     in constructing "fvsu" in clinic and diagnostic caluclations
c-----------------------------------------------------------------------
c
      do 200 k=1,km
        do 190 i=1,imt
          uclin(i,k,jpt) = u(i,k,jptr,nc)
          vclin(i,k,jpt) = v(i,k,jptr,nc)
190     continue
200   continue
c
c-----------------------------------------------------------------------
c     add external mode velocity to internal mode velocity for row
c     row "jsrow" (tau & tau-1)
c-----------------------------------------------------------------------
c
      if (jrow .le. jmtm1) then
        do 300 m=1,2
          if (m .eq. 1) then
            n = nc
          else
            n = nm
          endif
          do 270 i=1,imtm1
            diag1   = p(i+1,jrow+1,m) - p(i  ,jrow,m)
            diag0   = p(i  ,jrow+1,m) - p(i+1,jrow,m)
            ubar(i) = -(diag1+diag0)*dyu2r(jrow)*hr(i,jrow)
            vbar(i) =  (diag1-diag0)*dxu2r(i)*hr(i,jrow)*csur(jrow)
270       continue
#ifdef cyclic
c
c         set cyclic boundary conditions
c
          ubar(imt) = ubar(2)
          vbar(imt) = vbar(2)
#endif
          do 290 k=1,km
            do 280 i=1,imu
              u(i,k,jptr,n) = (u(i,k,jptr,n) + ubar(i))*gm(i,k,jptr)
              v(i,k,jptr,n) = (v(i,k,jptr,n) + vbar(i))*gm(i,k,jptr)
280         continue
290       continue
300     continue
      endif
c
      return
      end
