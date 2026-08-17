      subroutine ocn1st
      implicit real*8 (a-h,o-z)
c
c
c=======================================================================
c
c     initialize quantites at initial condition time
c
c=======================================================================
c
c
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef multitasking
      parameter (ntp1=nt+1, ntp2=nt+2, ntp3=nt+3, ntp4=nt+4
     $,          ntp5=nt+5, ntp6=nt+6, ntp7=nt+7)
#include "cshrbf.h"
      dimension ua(imt,km), va(imt,km)
      equivalence (ta(1,1,ntp1), ua(1,1))
      equivalence (ta(1,1,ntp2), va(1,1))
# ifdef tcvmix
      dimension vvca(imt,km), vdca(imt,km), vdqca(imt,km), q2a(imt,km)
      equivalence (ta(1,1,ntp3), vvca(1,1))
      equivalence (ta(1,1,ntp4), vdca(1,1))
      equivalence (ta(1,1,ntp5), vdqca(1,1))
      equivalence (ta(1,1,ntp6), q2a(1,1))
#  ifdef leq
      dimension q2la(imt,km)
      equivalence (ta(1,1,ntp7), q2la(1,1))
#  endif
# endif
#else
#include "slabs.h"
#endif
#ifdef seaice
#include "seayce.h"
#endif
#ifdef saltnorm
# include "conblk.h"
# include "sltdat.h"
      equivalence (sltsum,konblk(100))
#endif
#include "switch.h"
c
      dimension tinit(km,nt)
c
c-----------------------------------------------------------------------
c     generate topography over "t" points
c-----------------------------------------------------------------------
c
      call topog
c
c-----------------------------------------------------------------------
c     initialize slab data on disk
c-----------------------------------------------------------------------
c
      do 200 j=1,jmt
c
c      construct initial density profile
c
       do 110 k=1,km
#ifdef worldocn
c       match initial conditions of previous cyber r15 model runs
c
        tinit(k,1) = 3.9281
        tinit(k,2) = 0.0346786
# ifdef agetracers
c
c       initialize age tracer to zero
c
        tinit(k,ntage1) = c0
# endif
# ifdef cfctracers
c
c       initialize chlorofluorcarbon tracers to zero
c
        tinit(k,ntcfc1) = c0
        tinit(k,ntcfc2) = c0
# endif
#else
         tinit(k,1) = theta0 (yt(j), zt(k))
         tinit(k,2) = 0.0349
#endif
110    continue
c
c       set internal mode velocities to zero
c
        do 130 k=1,km
          do 120 i=1,imt
            ua(i,k) = c0
            va(i,k) = c0
#ifdef tcvmix
c
c      set turbulent kinetic energy
c
            q2a(i,k)  = 1.e-6
#  ifdef leq
            q2la(i,k) = 1.e-4
#  endif
#endif
120       continue
130     continue
c
        do 160 m=1,nt
          do 150 k=1,km
            do 140 i=1,imt
              ta(i,k,m) = c0
140         continue
150       continue
160     continue
c
c       set tracers over ocean points to specified values
c
        do 190 i=1,imt
          kz = kmt(i,j)
          if (kz .ne. 0) then
            do 180 k=1,kz
              do 170 m=1,nt
                ta(i,k,m) = tinit(k,m)
cc              if(m .eq. 1) ta(i,k,m) = tinit(k,m)
                if(m .eq. 2) ta(i,k,m) = tinit(k,m) - 0.035
170           continue
180         continue
          endif
190     continue
c
c       send the initial slabs to disk. use "oput" since data is
c       contiguous in buffer slab "ta"
c
        call oput (labs(1), nslab, (j-1)*nslab+1, ta)
        call oput (labs(2), nslab, (j-1)*nslab+1, ta)
        if (ntlev .eq. 3 .and. (j .eq. 1 .or. j .eq. jmt))
     $     call oput (labs(3), nslab,(j-1)*nslab+1, ta)
200   continue
c
#ifdef saltnorm
c
c-----------------------------------------------------------------------
c     store the total salt content of the ocean for use by the
c     salt normalization code.  Note this code assumes tinit
c     is NOT a fuction of i or j or k!
c-----------------------------------------------------------------------
c
      totslt = (tinit(1,2) - 0.035)
      sltsum = c0
c
#endif
c
c-----------------------------------------------------------------------
c     end slab initialization on disk
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c     initialize two dimensional fields on disk
c-----------------------------------------------------------------------
c
      do 300 n=1,ntau
        do 290 j=1,jmt
          do 280 i=1,imt
            p (i,j,n) = c0
#ifdef seaice
            yc(i,j,n) = c0
#endif
280       continue
290     continue
300   continue
c
c     block's 1 through 3 are for stream function time levels
c     block's 4 & 5 are for the stream function guess field
c     block 6 is for the level indicators over "t" points
#ifdef seaice
c     block 7-8 are the two time levels (tau, tau-1...order depends
c     on the timestep number) of sea ice thickness
#endif
c
      call oput (kflds, nwds, 1, p)
      call oput (kflds, nwds, nwds+1, p)
      call oput (kflds, nwds, 3*nwds+1, p)
      call oput (kflds, nwds, 4*nwds+1, p)
# ifdef f90type
      call oputi (kflds, nwds, 5*nwds+1, kmt)
# else
      call oput (kflds, nwds, 5*nwds+1, kmt)
# endif
#ifdef seaice
      call oput (kflds, nwds, 6*nwds+1, yc)
      call oput (kflds, nwds, 7*nwds+1, yc)
      if (ntlev .EQ. 3) call oput (kflds, nwds, 8*nwds+1, yc)
#endif
c
c-----------------------------------------------------------------------
c     initialize control block on disk
c-----------------------------------------------------------------------
c
      itt   = 0
#ifdef gtime
c     initialize current time at initialization to be equivalent to
c     Jan 1 0:0:0 of year number 0
c     (elsewhere, the timestamp is formed by adding year0 to icyear)
c
      icyear     = 0
      icdayofyr  = 1
      icsecofday = 0
#else
      years = c0
#endif
c
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
c
      return
      end
