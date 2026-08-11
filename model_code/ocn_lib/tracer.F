      subroutine tracer(j)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     tracer computes tracer quantities at time level n+1 for row j
c
c=======================================================================
c
c
#include "param.h"
#include "accel.h"
#ifdef biharmonic
# include "cbihar.h"
#endif
#include "cdiag.h"
#include "chmix.h"
#ifdef inland
#include "cinlnd.h"
#endif
#ifdef isopycmix
#include "cisop.h"
#endif
#ifdef nlhmix
# include "cnlmix.h"
#endif
#include "cregin.h"
#include "ctask.h"
#if defined saltnorm || defined agetracers
#include "ctmngr.h"
#endif
#ifdef worldocn
#include "convt.h"
#endif
#include "cvbc.h"
#include "cvmix.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#include "seasnl.h"
#include "slabs.h"
#ifdef saltnorm
# include "conblk.h"
# include "sltdat.h"
      equivalence (sltsum,konblk(100))
#endif
#include "switch.h"
#include "temp.h"
      dimension twodt(km)
#ifdef kdtermb
c     dimension arrays to keep track of convective tracer changes
      dimension tchang(imt,km,nt), tbefor(imt,km,nt)
#endif
#ifdef saltnorm
# include "thick.h"
#endif
#ifdef amazon 
# include "camazon.h"
c
c     variable "vdc2" is used to force complete mixing of levels    
c     1 thru kamzn within the range of i,j points given by   
c     jamzns, jamzne and iamzns, iamzne. this is accomplished by 
c     setting the vertical diffusion coeff to vdclim, which simulates
c     free convection via the implicit vertical mixing.
c
      dimension vdc2(imt,km)
#endif
#ifndef keepterms
c
c---------------------------------------------------------------------
c     define operators using statement functions ...
c     otherwise they will be dimensioned in "slabs.h"
c---------------------------------------------------------------------
c
#include "fdift.h"
#endif
c
      if (j .eq. jsrow) then
c
c-----------------------------------------------------------------------
c       on the first row of a task:
c
c       1) compute the advective coefficients at the south face of the
c          "t" grid box for the next row (jsrow+1)
c
c       2) return
c-----------------------------------------------------------------------
c
c       calculate "fvst" this way to insure that answers do not change
c       regardless of how "ntasks" is set
c
        fxb = cstr(j)*dytr(j)*csu(j)
        do 100 k=1,km
          do 90 i=2,imt
            fvst(i,k) = (v(i,k,jc,nc)*dxu(i) + v(i-1,k,jc,nc)*dxu(i-1))
     $               *fxb*dxt4r(i)
#ifdef noadvection
c        for test purposes, ignore advection of tracers
            fvst(i,k) = c0
#endif
90        continue
100     continue
        fx = cst(j)*dyt(j)*cstr(j+1)*dytr(j+1)
        do 170 k=1,km
          do 160 i=2,imt
            fvst(i,k) = fvst(i,k)*fx
160       continue
170     continue
c
        return
      else
c
c-----------------------------------------------------------------------
c       on all but the first row of a task:
c
c       1) compute the advective coefficients "fuw" at the west face of
c         the "t" grid box & "fvn" at the north face of the "t" grid box
c
c       2) initialize "w" & "tracers"
c-----------------------------------------------------------------------
c
        fxa = cstr(j)*dytr(j)
        fxb = cstr(j)*dytr(j)*csu(j)
        do 200 k=1,km
          do 180 i=2,imt
            fuw(i,k) = (u(i-1,k,jc,nc)*dyu(j) + u(i-1,k,jm,nc)*dyu(j-1))
     $               *fxa
            fvn(i,k) = (v(i  ,k,jc,nc)*dxu(i) + v(i-1,k,jc,nc)*dxu(i-1))
     $               *fxb*dxt4r(i)
#ifdef noadvection
c        for test purposes, ignore advection of tracers
            fuw(i,k) = c0
            fvn(i,k) = c0
#endif
180       continue
#ifdef cyclic
          fuw(imt,k) = fuw(2,k)
#else
          fuw(imt,k) = c0
#endif
c
#ifdef skipland
          do 190 i=1,imt
            w(i,k) = c0
190       continue
#endif
200     continue
c
#ifdef skipland
        do 300 m=1,nt
          do 290 k=1,km
            do 280 i=1,imt
              ta(i,k,m) = c0
280         continue
290       continue
300     continue
#endif
      endif
#ifdef nlhmix
c
c-----------------------------------------------------------------------
c     set diffusion factors for the row
c-----------------------------------------------------------------------
c
c     use pure deformation for nonlinear diffusion coeff
c
      hfact = c1
      bah   = c0
      do 250 k=1,km
        do 240 i=2,imt
          dpt1 = p25*(dpure(i,k)  + dpure(i-1,k)
     $              + dpures(i,k) + dpures(i-1,k))
          dpt2 = p25*(dpures(i,k) + dpures(i-1,k)
     $              + dpuret(i,k) + dpuret(i-1,k))
          bbtj(i,k) = (effk(i)*hfact*dpures(i,k) + bah)*bbt(j)
          cctj(i,k) = (effky(j  )*hfact*dpt1 + bah)*cct(j)
          ddtj(i,k) = (effky(j-1)*hfact*dpt2 + bah)*ddt(j)
240     continue
250   continue
#endif
#ifdef skipland
c
c-----------------------------------------------------------------------
c     solve "tracers" over longitudinal strips of ocean points
c     if skipland is not set them "istart" & "iend" are parameters
c-----------------------------------------------------------------------
c
      do 1300 lstrip=1,lseg
        istart = ist(j,lstrip)
        iend   = iet(j,lstrip)
        if (istart .ne. 0) then
#endif
c
c-----------------------------------------------------------------------
c         compute vertical velocity in t columns
c-----------------------------------------------------------------------
c
c         compute divergence of mass at each level
c
          do 400 k=1,kmm1
            do 390 i=istart,iend
              w(i,k) = c2dzt(k)*((fuw(i+1,k)-fuw (i,k))*dxt4r(i)
     $                           +fvn(i  ,k)-fvst(i,k))
#ifdef noadvection
c        for test purposes, ignore advection of tracers
              w(i,k) = c0
#endif
390         continue
400       continue
c
c         integrate downward from the surface so vertical velocity is
c         defined at the bottom of each level
c
          do 500 k=2,kmm1
            do 490 i=istart,iend
              w(i,k) = w(i,k) + w(i,k-1)
490         continue
500       continue
c
c         set vertical velocity at the surface to zero (rigid-lid)
c         set vertical velocity at the bottom to zero also
c
          do 600 i=istart,iend
            w(i,0)        = c0
            w(i,kmt(i,j)) = c0
 600      continue
c
#ifdef isopycmix
c
c         add fk3(i,k+1,3) (the zz term from isopycnal mixing) to
c         vdc(i,k) for use in implicit mixing
c
          do 610 i=1,imt
            vdc(i,km) = c0
# ifdef amazon
            vdc2(i,km) = c0
# endif
610       continue
c
          do 630 k=1,kmm1
            do 620 i=1,imt
              vdc(i,k) = vdc(i,k) + fk3(i,k+1,3)
# ifdef amazon
              vdc2(i,k) = vdc(i,k)
# endif
620         continue
630       continue
c
#endif
# ifdef worldocn
c
c        initialize diagnostic convection flag to zero
c
         do 670 i=1,imt
           do 660 k=1,km
             flgcnv(i,k) = c0
660        continue
670      continue
# endif
          do 1200 m=1,nt
c
c-----------------------------------------------------------------------
c           calculate quantities for the computation of
c           vertical diffusion and advection of tracers
c-----------------------------------------------------------------------
c
            do 700 k=1,kmm1
              do 690 i=istart,iend
#ifdef tcvmix
                vtf(i,k) = vdc(i,k,jc,nm)*(t(i,k,jc,nm,m)
     $                                 - t(i,k+1,jc,nm,m))*dzwr(k)
#else
# ifdef amazon
c     force complete mixing of salinities from the surface thru level
c     kamzn (see "camazon.h") within the range of i,j points given by
c     jamzns, jamzne and iamzns, iamzne. this is accomplished by
c     setting the vertical diffusion coeff to a big number,
c     which simulates free convection via the implicit vertical mixing.
                if ( m .eq. 2
     $               .and.  k.lt.kamzn  .and. (k+1).le.kmt(i,j)
     $               .and. (j.ge.jamzns .and. j.le.jamzne)
     $               .and. (i.ge.iamzns .and. i.le.iamzne)) then
                  vdc2(i,k) = 1.0e6 
                  vtf(i,k) = 1.0e6 *(t(i,k,jc,nm,m)-t(i,k+1,jc,nm,m))
     $                     *dzwr(k)
ccccccccccccccccc
ccccc write(stdout,6513) i,j,k,m, vdc2(i,k)
6513  format(' AMAZON MIXING AT I,J,K,M,vdc2= ',4i5,e16.5)
ccccccccccccccccc
                else
                  vtf(i,k) = vdc(i,k)*(t(i,k,jc,nm,m)-t(i,k+1,jc,nm,m))
     $                     *dzwr(k)
                endif
# else
                vtf(i,k) = vdc(i,k)*(t(i,k,jc,nm,m) - t(i,k+1,jc,nm,m))
     $                     *dzwr(k)
# endif
#endif
                fwb1(i,k) = w(i,k)*(t(i,k,jc,nc,m) + t(i,k+1,jc,nc,m))
690           continue
700         continue
c
c           for vertical diffusion of tracers:
c           set k=0 elements to reflect surface tracer flux.
c           set bottom level "kz" to reflect insulation condition.
c           for vertical advection of tracers:
c           set k=0 elements of vertical advection array to the rigid
c           lid condition (w(0) = 0). also assume w(km)=0. If "kz" <
c           "km" and there is a bottom slope, w(kz) may be non zero.
c
c
            do 800 i=istart,iend
              fwb1(i,0)        = c0
              fwb1(i,km)       = c0
              vtf(i,0)         = stf(i,m)
              vtf(i,kmt(i,j))  = btf(i,m)
800         continue
c
c-----------------------------------------------------------------------
c           construct auxiliary arrays for statement functions to
c           prevent unnecessary computation.
c-----------------------------------------------------------------------
c
            do 900 k=1,km
              do 890 i=istart,iend+1
                aux1(i,k) = fuw(i,k)*(t(i,k,jc,nc,m)+t(i-1,k,jc,nc,m))
#if defined biharmonic && !defined isopycmix
                aux2(i,k) = dxu2r(i-1)*(del2(i,k,m,jpt2)
     $                               -  del2(i-1,k,m,jpt2))
#else
# ifndef isopycmix
                aux2(i,k) = dxu2r(i-1)*(t(i,k,jc,nm,m)-t(i-1,k,jc,nm,m))
# endif
#endif
890           continue
900         continue
c
c-----------------------------------------------------------------------
c           construct tracer source terms
c-----------------------------------------------------------------------
c
#ifdef agetracers
c
            if (m .eq. ntage1 ) then
c
c             put in place the sub-surface sources for the idealized
c             age tracers.  the "sourct" variable is in terms of a time
c             rate of change, and the age tracer is to be reported
c             in units of years, so the growth rate "growrt" is set to 
c             1/[#sec in a year] per second
# ifdef gtime
c             when using gtime option, year length is hardwired to 365
              yrnorm = 365.0
# endif
c
            growrt = c1 / (yrnorm * 86400.0)
c
            do 940 k=2,km
              do 935 i=istart,iend
                sourct(i,k,m) = sourct(i,k,m) + growrt
935           continue
940         continue
c
          endif
c
#endif
#ifdef ccl4tracer
c
          if (m .eq. ntccl4 ) then
c
c           put in place the hydrolysis sink for carbon tetrachloride
c           from Jeffers, P.M. and N.L. Wolfe (1989) Hydrolysis of 
c           Carbon Tetrachloride.  Science, 246:1638.
c           rate of change, and the age tracer is to be reported
c           in units of years, so the growth rate "growrt" is set to 
c           1/[#sec in a year] per second
c
            degk = 273.16
c
            do 955 k=2,km
              do 953 i=istart,iend
                sourct(i,k,m) = sourct(i,k,m)
     $            - t(i,k,jc,nm,ntccl4)*
     $            (4.07e12*exp(-13790.0/(degk+t(i,k,jc,nm,nttemp))))
     $            / 60.0
953           continue
955         continue
            endif
#endif
#if defined restorst || defined cpldleg2 || defined cpld_restore
c
c-----------------------------------------------------------------------
c           restore surface tracers to prescribed values by adding
c           newtonian damping to source term. (see "setvbc.F")
c-----------------------------------------------------------------------
c
            do 960 i=istart,iend
              sourct(i,1,m) = sourct(i,1,m) + restr(i,m) / dzt(1)
960         continue
#endif
#if defined cpldleg3 && !defined cpld_restore
c----------------------------------------------------------------------
c           add the surface flux correction terms for t and s
c           (see setvbc.F)
c----------------------------------------------------------------------
            do 965 i=istart,iend
              sourct(i,1,m) = sourct(i,1,m) + stfadj(i,m) / dzt(1)
965         continue
#endif
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c           construct arrays for use in isopycnal diffusion
c           calculation
c-----------------------------------------------------------------------
c
            call isop1 (j, m)
#endif
#ifdef keepterms
c
c-----------------------------------------------------------------------
c           if terms are defined as arrays, then calculate them here
c           otherwise statement functions will be used
c-----------------------------------------------------------------------
c
            do 980 k=1,km
              do 970 i=istart,iend
#include "fdift.h"
970           continue
980         continue
#endif
c
c-----------------------------------------------------------------------
c           calculate the new tracer quantities allowing for implicit
c           treatment of vertical diffusion
c-----------------------------------------------------------------------
c
            do 1000 k=1,km
              twodt(k) = c2dtts*dtxcel(k)
              do 990 i=istart,iend
                ta(i,k,m) = (t(i,k,jc,nm,m) + twodt(k)*(
     $                       UTx(i,k) + VTy(i,k) + WTz(i,k)
#ifdef isopycmix
     $                     + Tisox(i,k) + Tisoy(i,k) + Tisoz(i,k)
#else
     $                     + Txx(i,k) +Tyy(i,k)
#endif
#if defined implicitvmix || defined isopycmix
     $                     + Tzz(i,k)*(c1-aidif)
#else
     $                     + Tzz(i,k)
#endif
     $                     + sourct(i,k,m))
     $                                   )*fm(i,k,jc)
990           continue
1000        continue
c
#ifdef inland
c-----------------------------------------------------------------------
c     establish communication with inland seas
c-----------------------------------------------------------------------
c
            do 1020 nin=1,ninlnd
              if (jinlnd(nin) .eq. j) then
                i1 = iinlnd(nin,1)
                i2 = iinlnd(nin,2)
                if (i2 .ge. istart .and. i1 .le. iend) then
                  do 1010 k=kinlnd(nin,1),kinlnd(nin,2)
                    ta(i1,k,m) = ta(i1,k,m) + c2dtts*dtxcel(k)*
     $                          (t(i2,k,jc,nm,m)-t(i1,k,jc,nm,m)) *
     $                          binlnd(nin,1)*fm(i1,k,jc)*fm(i2,k,jc)
                    ta(i2,k,m) = ta(i2,k,m) + c2dtts*dtxcel(k)*
     $                          (t(i1,k,jc,nm,m)-t(i2,k,jc,nm,m)) *
     $                          binlnd(nin,2)*fm(i1,k,jc)*fm(i2,k,jc)
1010              continue
                endif
              endif
1020        continue
#endif  
#if defined implicitvmix || defined isopycmix
c
c----------------------------------------------------------------------
c           solve vertical diffusion implicitly
c----------------------------------------------------------------------
c
c           store terms to compute implicit vertical mixing on
c           diagnostic time steps
c
            if (diagts .and. eots) then
              do 1050 k=1,km
                do 1040 i=istart,iend
                  zzi(i,k,1) = ta(i,k,m)
1040            continue
1050          continue
            endif
c
# ifdef tcvmix
            call invtri (ta(1,1,m),  stf(1,m), btf(1,m)
     $,             vdca, twodt, kmt, fm(1,1,jc), istart, iend, j)
# else
#  ifdef amazon
            if ( m .eq. 2 ) then
ccccccccccccc
ccccc           if ( j.ge.jamzns .and. j.le.jamzne) then
ccccc             do 1634 i=iamzns,iamzne
ccccc             write(stdout,6575) i,j,m, vdc2(i,1), vdc(i,1),
ccccc$                       ta(i,1,m), ta(i,2,m)
1634              continue
ccccc           endif
6575  format(3i4,' AMAZON: before invtri vdc2,vdc,s(1),s(2) =',4e12.5)
ccccccccccccc
              call invtri (ta(1,1,m),  stf(1,m), btf(1,m)
     $,             vdc2, twodt, kmt, fm(1,1,jc), istart, iend, j)
ccccccccccccc
ccccc           if ( j.ge.jamzns .and. j.le.jamzne) then
ccccc             do 1635 i=iamzns,iamzne
ccccc               write(stdout,6576) i,j,vdc2(i,1), vdc(i,1),
ccccc$                       ta(i,1,m), ta(i,2,m)
1635              continue
ccccc           endif
6576  format(2i4,' AMAZON: after invtri vdc2,vdc,s(1),s(2) =',4e12.5)
ccccccccccccc
            else
              call invtri (ta(1,1,m),  stf(1,m), btf(1,m)
     $,              vdc, twodt, kmt, fm(1,1,jc), istart, iend, j)
            endif
#  else
            call invtri (ta(1,1,m),  stf(1,m), btf(1,m)
     $,              vdc, twodt, kmt, fm(1,1,jc), istart, iend, j)
#  endif
# endif
c
c           compute residual implicit vertical mixing
c
            if (diagts .and. eots) then
              do 1070 k=1,km
                rc2dt = c1/twodt(k)
                do 1060 i=istart,iend
                  zzi(i,k,1) = rc2dt*(ta(i,k,m) - zzi(i,k,1))
1060            continue
1070          continue
            endif
#endif
#ifdef saltnorm
c
c-----------------------------------------------------------------------
c           perform the salt normalization before every 5th mixing
c           timestep
c-----------------------------------------------------------------------
c
            if (m .eq. 2) then
              n5mix = nmix * 5
              if (mod(itt,n5mix) .eq. 0) then
                if (j.eq.2) write(stdout,888) sltnrm
                do 1081 k=1,km
                 do 1081 i=istart,iend
                  ta(i,k,m) = ta(i,k,m) + sltnrm*fm(i,k,jc)
 1081           continue
              endif
            endif
c
#endif       
c
c-----------------------------------------------------------------------
c           do analysis of tracers on diagnostic timesteps over
c           longitudinal strips of ocean points
c-----------------------------------------------------------------------
c
            if (diagts .and. eots) then
              f1 = cst(j)*dyt(j)
              do 1100 i=istart,iend
                kz = kmt(i,j)
                if (kz .ne. 0) then
                  area       = f1*dxt(i)
                  sstf       = area*stf(i,m)
                  ssst       = area*t(i,1,jc,nc,m)
                  asst(m,0)  = asst(m,0)  + ssst
                  stflx(m,0) = stflx(m,0) + sstf
                  n          = mskhr(i,j)
                  if (n .gt. 0) then
                    stflx(m,n) = stflx(m,n) + sstf
                    asst(m,n)  = asst(m,n)  + ssst
                  endif
                  do 1080 k=1,kz
                    boxvol = area*dzt(k)
                    satwv  = boxvol*t(i,k,jc,nc,m)
                    sutx   = boxvol*UTx(i,k)
                    svty   = boxvol*VTy(i,k)
                    swtz   = boxvol*WTz(i,k)
#ifdef isopycmix
                    stxx   = boxvol*Tisox(i,k)
                    styy   = boxvol*Tisoy(i,k)
#else
                    stxx   = boxvol*Txx(i,k)
                    styy   = boxvol*Tyy(i,k)
#endif
#if defined implicitvmix || defined isopycmix
                    stzz   = boxvol*(Tzz(i,k)*(c1-aidif) + zzi(i,k,1))
# ifdef isopycmix
     $                     + boxvol*Tisoz(i,k)
# endif
#else
                    stzz   = boxvol*Tzz(i,k)
#endif
                    ssourc = boxvol*sourct(i,k,m)
c
                    n = 0
                    termbt(k,2,m,n) = termbt(k,2,m,n) + sutx
                    termbt(k,3,m,n) = termbt(k,3,m,n) + svty
                    termbt(k,4,m,n) = termbt(k,4,m,n) + swtz
                    termbt(k,5,m,n) = termbt(k,5,m,n) + stxx
                    termbt(k,6,m,n) = termbt(k,6,m,n) + styy
                    termbt(k,7,m,n) = termbt(k,7,m,n) + stzz
                    termbt(k,8,m,n) = termbt(k,8,m,n) + ssourc
                    atwv(m,n)       = atwv(m,n)       + satwv
                    n      = nhreg*(mskvr(k)-1) + mskhr(i,j)
                    if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                      termbt(k,2,m,n) = termbt(k,2,m,n) + sutx
                      termbt(k,3,m,n) = termbt(k,3,m,n) + svty
                      termbt(k,4,m,n) = termbt(k,4,m,n) + swtz
                      termbt(k,5,m,n) = termbt(k,5,m,n) + stxx
                      termbt(k,6,m,n) = termbt(k,6,m,n) + styy
                      termbt(k,7,m,n) = termbt(k,7,m,n) + stzz
                      termbt(k,8,m,n) = termbt(k,8,m,n) + ssourc
                      atwv(m,n)       = atwv(m,n)       + satwv
                    endif
1080              continue
c
                  if (m .eq. 1 .and. kz .ge. 2) then
                    fx = area*grav*p5
                    do 1090 k=2,kz
                      buoy(k) = buoy(k) - fx*dzw(k-1)*w(i,k-1)*
     $                                (rhos(i,k-1)+rhos(i,k))
1090                continue
                  endif
                endif
1100          continue
            endif
#ifdef isopycmix
c
c     save the north - central isopycnal diffusion of "e(i,k,2)"
c     (yy + zy terms ... xy is neglected) into "esav(i,k,m)" for
c     use during the next j row calculation
c
      do 1190 k=1,km
        do 1180 i=1,imt
          esav(i,k,m) = e(i,k,2)
1180    continue
1190  continue
#endif
1200      continue
#ifdef skipland
        endif
1300  continue
#endif
#ifdef seaice
c-----------------------------------------------------------------------
      call icedyn (j)
      call icethm (j, .false.)
c
c-----------------------------------------------------------------------
#endif
#ifndef implicitvmix
c
c-----------------------------------------------------------------------
c         convectively adjust water column if gravitationally unstable
c         over longitudinal strips of ocean points
c-----------------------------------------------------------------------
c
c         ncon determines number of passes through convection loop
c
c         ks=1: compare lev. 1 to 2; 3 to 4; etc.
c         ks=2: compare lev. 2 to 3; 4 to 5; etc.
# ifdef kdtermb
c
c    for tracers, keep track of change due to convective adjustment
c
          if (diagts .and. eots) then
            do 1368 m=1,nt
              do 1367 kq=1,km
                do 1366 iq=1,imt
                  tchang(iq,kq,m) = c0
1366            continue
1367          continue
1368        continue
          endif
# endif
c
          do 1400 n=1,ncon
# ifdef kdtermb
c
c    for tracers, keep track of change due to convective adjustment
c
              if (diagts .and. eots) then
                do 1308 m=1,nt
                  do 1307 kq=1,km
                    do 1306 iq=1,imt
                      tbefor(iq,kq,m) = ta(iq,kq,m)*fm(iq,kq,jc)
1306                continue
1307              continue
1308            continue
              endif
c
# endif
            do 1390 ks=1,2
c
c             find density for entire slab for stability determination
c
              call statec (ta(1,1,1), ta(1,1,2), tempik(1,1,1), ks)
# ifdef skipland
              do 1380 lstrip=1,lseg
                istart = ist(j,lstrip)
                iend   = iet(j,lstrip)
                if (istart .ne. 0) then
# endif
c
c                 set "heavy water" in land to stop convection
c
                  dense = 1.e15
                  do 1340 i=istart,iend
                    k = kmt(i,j) + 1
                    if (k .le. km) then
# ifdef worldocn
                      do 1330 kk=k,km
                        tempik(i,kk,1) = dense
1330                  continue
# else
                      tempik(i,k,1) = dense
# endif
                    endif
1340              continue
c
c                 if unstable,  mix tracers on adjoining levels
c
                  do 1370 m=1,nt
                    do 1360 k=ks,kmm1,2
                      do 1350 i=istart,iend
                        if (tempik(i,k,1) .gt. tempik(i,k+1,1)) then
                          ta(i,k,m) = (dztxcl(k)  *ta(i,k,m) +
     $                                 dztxcl(k+1)*ta(i,k+1,m))
     $                                 *dzwxcl(k)
                          ta(i,k+1,m) = ta(i,k,m)
# ifdef worldocn
c                         set diagnostic convection flags to one
                          flgcnv(i,k) = c1
                          flgcnv(i,k+1) = c1
# endif
                        endif
1350                  continue
1360                continue
1370              continue
# ifdef skipland
                endif
1380          continue
# endif
c
1390        continue
# ifdef kdtermb
c    for tracers, keep track of change due to convection
c
          if (diagts .and. eots) then
            do 1398 m=1,nt
              do 1397 kq=1,km
                do 1396 iq=1,imt
                  tchang(iq,kq,m) = tchang(iq,kq,m) + 
     $                  (ta(iq,kq,m)-tbefor(iq,kq,m))*fm(iq,kq,jc)
1396            continue
1397          continue
1398        continue
          endif
# endif
c
1400      continue
c
# ifdef seaice
c
c-----------------------------------------------------------------------
c     adjust sea ice values after conv adjustment loop is completed
c-----------------------------------------------------------------------
c
            call icethm (j, .true.)
# endif
c
#endif
c
#ifndef kdtermb
c-----------------------------------------------------------------------
c         integrate total changes in t,s and squared t,s on diagnostic
c         timesteps over longitudinal strips of ocean points
c-----------------------------------------------------------------------
c
          if (diagts .and. eots) then
            do 1500 m=1,nt
              do 1490 k=1,km
                fx = cst(j)*dyt(j)*dzt(k)/(c2dtts*dtxcel(k))
                do 1480 i=2,imtm1
                  boxfx           = fx*dxt(i)*fm(i,k,jc)
                  sddt            = (ta(i,k,m)-t(i,k,jc,nm,m))*boxfx
                  svar            = (ta(i,k,m)**2-t(i,k,jc,nm,m)**2)
     $                              *boxfx
                  n               = 0
                  termbt(k,1,m,n) = termbt(k,1,m,n) + sddt
                  tvar(k,m,n)     = tvar(k,m,n)     + svar
                  n     = nhreg*(mskvr(k)-1) + mskhr(i,j)
                  if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                    termbt(k,1,m,n) = termbt(k,1,m,n) + sddt
                    tvar(k,m,n)     = tvar(k,m,n)     + svar
                  endif
1480            continue
1490          continue
1500        continue
          endif
#endif
c
#if defined cpldleg3 || defined cpldleg2
c
c     save unfiltered sst and sea ice for the atmos segment
c
#ifdef cpldleg2
      if(eots) then
#endif
#ifdef cpldleg3
c  For coupled model, store sst and sea ice on last timestep of the day
      if(endseg .and. eots) then
#endif
        call bcatm(j)
      endif  
c
c
#endif
c
#if defined fourfil || defined firfil
c
c-----------------------------------------------------------------------
c     fourier filter tracers at high latitudes
c-----------------------------------------------------------------------
c
#ifdef cyclic
c
c     set cyclic boundary conditions on newly computed tracers
c
      do 1600 m=1,nt
        do 1590 k=1,km
          ta(1  ,k,m) = ta(imtm1,k,m)
          ta(imt,k,m) = ta(2    ,k,m)
1590    continue
1600  continue
# ifdef seaice
      yca(  1,j) = yca(imtm1,j)
      yca(imt,j) = yca(    2,j)
# endif
#endif
      call filt (j)
c
#ifdef seaice
c-----------------------------------------------------------------------
c     fourier filter sea ice at high latitudes
c-----------------------------------------------------------------------
      call filtyc (j)
#endif
#endif
c
#ifdef kdtermb
c-----------------------------------------------------------------------
c         integrate total changes in t,s and squared t,s on diagnostic
c         timesteps over longitudinal strips of ocean points
c         (do so after filtering to account for possible zonal
c         transfers due to filtering)
c         keep track of change due to convection in termbt(k,9,m,n)
c-----------------------------------------------------------------------
c
          if (diagts .and. eots) then
            do 1650 m=1,nt
              do 1640 k=1,km
                fx = cst(j)*dyt(j)*dzt(k)/(c2dtts*dtxcel(k))
                do 1630 i=2,imtm1
                  boxfx           = fx*dxt(i)*fm(i,k,jc)
                  sddt            = (ta(i,k,m)-t(i,k,jc,nm,m))*boxfx
                  svar            = (ta(i,k,m)**2-t(i,k,jc,nm,m)**2)
     $                              *boxfx
                  n               = 0
                  termbt(k,1,m,n) = termbt(k,1,m,n) + sddt
                  tvar(k,m,n)     = tvar(k,m,n)     + svar
                  n     = nhreg*(mskvr(k)-1) + mskhr(i,j)
                  if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                    termbt(k,1,m,n) = termbt(k,1,m,n) + sddt
                    tvar(k,m,n)     = tvar(k,m,n)     + svar
                    termbt(k,9,m,n) = termbt(k,9,m,n) +  
     $                                boxfx*tchang(i,k,m)
                  endif
1630            continue
1640          continue
1650        continue
          endif
#endif
      if (tavgts .and. eots) then
#ifdef kdtravg
        call regnj (j)
c
c       write formatted tracer average data to iotavg file
        if (j .eq. jmtm1) call region (iotavg, .true., .true., .true.)
c
c       write unformatted tracer average data to unit iotavg2
        if (j .eq. jmtm1) call region (iotavg2, .false., .true., .true.)
#else
        call region (j, iotavg, .true., .true., .true.)
#endif
      endif
c
#ifdef saltnorm
c-----------------------------------------------------------------------
c     find the total salt content of the ocean for normalization
c-----------------------------------------------------------------------
c
      facj = cst(j)*dyt(j)
      do 1610 k=1,km
       do 1610 i=istart,iend
        tempik(i,k,1) = (ta(i,k,2)*fm(i,k,jc))*facj*dzt(k)*dxt(i)
 1610 continue
      do 1620 k=1,km
       do 1620 i=istart,iend
        sltsmj = sltsmj + tempik(i,k,1)
 1620 continue
c
#endif
c-----------------------------------------------------------------------
c     accumulate integrated absolute changes in tracers over
c     longitudinal strips of ocean points when printout of timestep
c     info is desired
c-----------------------------------------------------------------------
c
      if (prntsi) then
        fx1 = p5*cst(j)*dyt(j)/c2dtts
#ifdef skipland
        do 1700 lstrip=1,lseg
          istart = ist(j,lstrip)
          iend   = iet(j,lstrip)
          if (istart .ne. 0) then
#endif
            do 1690 m=1,nt
              do 1680 k=1,km
                fx = fx1/dtxcel(k)
                do 1660 i=istart,iend
                  tempik(i,k,1) = abs(ta(i,k,m)-t(i,k,jc,nm,m))*
     $                            c2dzt(k)*fx*dxt(i)*fm(i,k,jc)
1660            continue
                do 1670 i=istart,iend
                  dtabs(k,m,j) = dtabs(k,m,j) + tempik(i,k,1)
1670            continue
1680          continue
1690        continue
#ifdef skipland
          endif
1700    continue
#endif
      endif
c
c-----------------------------------------------------------------------
c     transfer quantities computed to the north of the present row
c     to be defined to the south in the computation of the next row
c-----------------------------------------------------------------------
c
      if (j .ne. jerow) then
        fx = cst(j)*dyt(j)*cstr(j+1)*dytr(j+1)
        do 1800 k=1,km
          do 1790 i=1,imt
            fvst(i,k) = fvn(i,k)*fx
            rhos(i,k) = rhon(i,k)
#ifdef nlhmix
            dtens(i,k)  = dten(i,k)
            dshers(i,k) = dshear(i,k)
            dpuret(i,k) = dpures(i,k)
            dpures(i,k) = dpure(i,k)
#endif
1790      continue
1800    continue
      endif
#ifdef cyclic
c
c     set cyclic boundary conditions on newly computed tracers
c
      do 1900 m=1,nt
        do 1890 k=1,km
          ta(1  ,k,m) = ta(imtm1,k,m)
          ta(imt,k,m) = ta(2    ,k,m)
1890    continue
1900  continue
# ifdef seaice
      yca(  1,j) = yca(imtm1,j)
      yca(imt,j) = yca(    2,j)
# endif
#endif
#ifdef icetest
      if (j.eq.74) then
        do 1974 i=11,15
        write(stdout,9074) j, i, yca(i,j), yc(i,j,nc), yc(i,j,nm)
9074    format(1x,'in subr tracer: j, i, yca, ycnc, ycnm=',/,
     $         2i5,3e13.5)
1974    continue
      endif
#endif
#ifndef symmetry
c
c-----------------------------------------------------------------------
c     set new velocities at northern wall to zero since no pass through
c     clinic is made for this row
c-----------------------------------------------------------------------
c
      if (j .eq. jmtm1) then
        do 2000 k=1,km
          do 1990 i=1,imt
            ua(i,k) = c0
            va(i,k) = c0
1990      continue
2000    continue
      endif
#endif
      return
#ifdef saltnorm
 888  format(' salt normalization done - adjustment factor = ',3pf10.6)
#endif
      end
