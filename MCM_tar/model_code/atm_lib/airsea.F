      program airsea
      implicit real*8 (a-h,o-z)
c
#if !defined cpldleg3 && !defined cpldleg2
#include "parm.h"
#endif
c
c
      character*80 ocndata, atmdata, albdata
#ifdef cpldleg2
c   This param.h is from the ocean library; it contains imt & jmt:
#include "param.h"
      parameter ( ix=(imt-2)/2, il=jmt-2 )
#endif
#ifndef unitsk
c
#include "parmcpu.h"
c
#endif
c
#if !defined cpldleg3 && !defined cpldleg2
      external maim1
      external s1gwgt
      external swtbls
#endif
c
#include "asblck.h"
c
#ifndef cpldleg2
#include "ctlblk.h"
#endif
c
#if defined cpldleg3 || defined cpldleg2
#include "conblk.h"

cAltix
      real*4 totdyo
      character*4 omonth
      equivalence (omonth, konblk(22))
cAltix
#endif
c
#if !defined cpldleg2
#include "iohist.h"
#endif
c
#include "cexpnam.h"
c
      character*16 expern
      data cyear /365.0/
c
c  Define all fortran units used in model
c    See iohist.h for details
c
cSGI
cif !defined cpldleg3 && !defined cpldleg1
c     data itap / 11 /
c     data irsf / 12 /
cendif
cif defined cpldleg3 || defined cpldleg1
c     data itap   / 71 /
c     data irsf   / 72 /
c     data irunof / 38 /
cendif
c     data irs    / 20 /
c     data kdsk   / 26 /
c     data adsk   / 30 /
c     data bdsk   / 31 /
c     data iunt   / 55 /
c     data ibd    / 56 /
c     data iunitr / 41 /
c     data iunitw / 42 /
c     data ialb   / 39 /
c     data idsk   / 25 /
c     data izsvar / 85 /
c     data intgl  / 77 /
c     data ictocn / 60 /
c     data isstyc / 36 /
c     data iflux  / 35 /
c     data iflxan / 37 /
c     data iflxac / 48 /
c     data iycthk / 58 /
c     data iycrst /  1 /
cSGI
c
c
c   <<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
cSGI
#if !defined cpldleg3 && !defined cpldleg1
      itap = 11
      irsf = 12
#endif
#if defined cpldleg3 || defined cpldleg1
      itap   = 71
      irsf   = 72
      irunof = 38 
#endif
      irs    = 20 
      kdsk   = 26 
      adsk   = 30 
      bdsk   = 31 
      iunt   = 55 
      ibd    = 56 
      iunitr = 41 
      iunitw = 42 
      ialb   = 39 
      idsk   = 25 
      izsvar = 85 
      intgl  = 77 
      ictocn = 60 
      isstyc = 36 
      iflux  = 35 
      iflxan = 37 
      iflxac = 48 
      iycthk = 58 
      iycrst =  1 
cSGI
c
#ifndef unitsk
cc    ncpus = max(1,ncpu())
cc    if (ncpus .gt. maxcpu) then
cc      print *,'maxcpu must be at least ',ncpus,' in parm.h.',
cc   &          '  Run terminated !!'
cc      stop
cc    end if
      read(5,4320) ncpus
      write(io_out, 4321) ncpus
      if (ncpus .gt. maxcpu) then
        write(6, 4319) maxcpu
        call aquit(4319)
        stop
      endif
#endif

c
c   Read and store the experiment name
c
      read (5,51) expern
c  store the experiment name
      expn1 = expern(1:8)
      expn2 = expern(9:16)
      write (6,61) expn1, expn2
c
#if defined cpldleg3 || defined cpldleg1
c  open fortran units for use by the model.
c
#ifdef ice_dynamics
      open(iycrst,form='unformatted')
#endif
      open(itap ,form='unformatted')
      open(iunt ,form='unformatted')
      open(ialb ,form='unformatted')
      open(irsf ,form='unformatted')
      open(irs  ,form='unformatted')
      open(adsk ,form='unformatted')
      open(bdsk ,form='unformatted')
      open(ibd  ,form='unformatted')
      open(idsk ,form='unformatted')
#endif
c
      la=.true.
      lo=.true.
c
c   Read boundary data at beginning of the run
      call bcio
c
#if !defined cpldleg3 && !defined cpldleg2 && !defined cpldleg1
      call rdctl(itap)
      rewind itap
c
      call indato
#endif
c
#ifndef cpldleg2
c
      call iniprm
c
c   atmosphere leg
 100  continue
c
#if defined fixed_sst || defined q_flux || defined restore
c
c get sea surface boundary conditions when these are prescribed
c rather than predicted
      astday = amod(day,cyear)
      call bciol(astday)
#endif
c
      call atmos
c
      la=.false.
#endif
c
#ifndef cpldleg1
 200  continue
# ifdef ice_dynamics
c
c sea-ice model (cfi)
c
      call sphercav
# endif
c
c   ocean leg
      call ocean
c
      lo=.false.
c
# if defined cpldleg3 || defined cpldleg2
      rewind ictocn
cAltix
ccc   write (ictocn) konblk
      ittocn  = konblk(3)
      totdyo  = konblk(9)
      idaymno = konblk(21)
      iyearo  = konblk(23)
      imontho = konblk(65)
      sltsmo  = konblk(100)
      print *,' in airsea '
      print *,' ittocn,totdyo,idaymno,iyearo,imontho,sltsmo,omonth='
      print *,  ittocn,totdyo,idaymno,iyearo,imontho,sltsmo,omonth
      write (ictocn) ittocn,totdyo,idaymno,iyearo,imontho,sltsmo,omonth
cAltix
# endif
#endif
c
#ifndef cpldleg2
      if(ktau .ge. ktaust) go to 300
      go to 100
c
c  stop when ktau  .ge. ktaust 
 300  continue
c
c
c  Write control block at end of run
      iunit  = 6
      ifrmat = 1
      call wrtctl(iunit,ifrmat)
#endif
c
# if defined cpldleg3 || defined cpldleg2
c
c   Write boundary data at end of the run
      call bcout
c
c  Write ocean control block at end of run
      write (6,302) konblk
# endif
c
#ifndef cpldleg2
c  Close FORTRAN units at end of run
      close(itap)
      close(iunt)
      close(ialb)
      close(irsf)
      close(irs )
      close(adsk)
      close(bdsk)
      close(ibd )
      close(idsk,status='delete')
#endif
c
c
  51  format(a16)
  61  format(/,1x,'The experiment name is ',2a8,/)
 302  format(//,' konblk',/,(1x,10i12))
 4319 format(' ncpus gt maxcpu, maxcpu = ',i5)
 4320 format(i5)
 4321 format(' for this run, ncpus = ',i5)
c
      stop
      end
