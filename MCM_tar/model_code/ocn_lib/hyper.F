      subroutine hyper
      implicit real*8 (a-h,o-z)
#ifdef hypergrid
c
c=======================================================================
c
c     hyper is similar to relax except it solves the equation along
c     diagonals to allow for complete vectorization & parallelization
c
c=======================================================================
c
c---------------------------------------------------------------------
c  define global data
c---------------------------------------------------------------------
c
#include "param.h"
#include "crelax.h"
#include "ctmngr.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "iounit.h"
#include "scalar.h"
#include "switch.h"
#include "temp.h"
c
      dimension ptdb(imt,jmt), resmi(jmt)
      equivalence (res,ptdb)
c
c=======================================================================
c     begin introductory section to prepare for the relaxation
c=======================================================================
c
c-----------------------------------------------------------------------
c    initiate readin of relaxation solution of 2 timesteps back
c    for the purpose of computing an initial guess for the present pass
c    (input unit alternates on tsteps between 4 & 5)
c-----------------------------------------------------------------------
c
      luptdb = 5 - mod(itt,2)
      luptd  = 4 + mod(itt,2)
      call ofind (kflds, nwds, (luptdb-1)*nwds+1)
      call ofind (kflds, nwds, (luptd-1)*nwds+1)
c
c
      if (first .or. mixts .or. mod(itt,nmix) .eq. 2) then
c
c-----------------------------------------------------------------------
c       initialize the work area
c-----------------------------------------------------------------------
c
        do 100 j=1,jmt
          do 90 i=1,imt
            cfn(i,j)  = c0
            cfs(i,j)  = c0
            cfe(i,j)  = c0
            cfw(i,j)  = c0
            rmsk(i,j) = c0
#ifdef islands
            cof(i,j)  = c1
#endif
90        continue
100     continue
c
c-----------------------------------------------------------------------
c       generate arrays of coefficients for relaxation
c-----------------------------------------------------------------------
c
c       1st, compute coefficients of the laplacian star
c       (hold non-interior points to zero using start and end indices)
c
        do 600 j=3,jmtm1
          do 590 l=1,lseg
            if (isz(j,l) .ne. 0) then
              ie = iez(j,l)
              fxa = c2*cstr(j)*cstr(j)
              fxb = c2*csu(j  )*cstr(j)*dytr(j)*dyur(j  )
              fxc = c2*csu(j-1)*cstr(j)*dytr(j)*dyur(j-1)
              do 560 i=isz(j,l),iez(j,l)
                rmsk(i,j) = c1
                cfn(i,j)  = fxb/(h(i-1,j)+h(i,j))
                cfs(i,j)  = fxc/(h(i-1,j-1)+h(i,j-1))
                cfe(i,j)  = fxa*dxur(i)*dxtr(i)/(h(i,j)+h(i,j-1))
                cfw(i,j)  = fxa*dxur(i-1)*dxtr(i)/(h(i-1,j)+h(i-1,j-1))
                cpf(i,j) = c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
560          continue
c
c             augment coeffs for implicit treatment of coriolis term
c
              if (acor .ne. c0) then
                fx = -c2dtsf*acor*cstr(j)*dytr(j)*omega
                do 570 i=isz(j,l),iez(j,l)
                  cfn(i,j) = cfn(i,j)+(hr(i,j  )-hr(i-1,j  ))*sine(j  )
     $                        *fx*dxtr(i)
                  cfs(i,j) = cfs(i,j)-(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
     $                        *fx*dxtr(i)
                  cfe(i,j) = cfe(i,j)-(hr(i  ,j)*sine(j)-hr(i  ,j-1)*
     $                        sine(j-1))*fx*dxtr(i)
                  cfw(i,j) = cfw(i,j)+(hr(i-1,j)*sine(j)-hr(i-1,j-1)*
     $                        sine(j-1))*fx*dxtr(i)
570            continue
              endif
c
c             normalize coefficients
c
              do 580 i=isz(j,l),iez(j,l)
                cfn(i,j) = cfn(i,j)*cpf(i,j)
                cfs(i,j) = cfs(i,j)*cpf(i,j)
                cfe(i,j) = cfe(i,j)*cpf(i,j)
                cfw(i,j) = cfw(i,j)*cpf(i,j)
580           continue
            endif
590       continue
600     continue
#ifdef islands
c
c       modify coefficients to include island perimeter points
c
        do 700 isle=1,nisle
          cofis(isle) = c0
          do 690 n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            fxa = c2*cstr(j)*cstr(j)
            fxb = c2*csu(j  )*dyur(j  )*dytr(j)*cstr(j)
            fxc = c2*csu(j-1)*dyur(j-1)*dytr(j)*cstr(j)
            fxd = -c2dtsf*acor*cstr(j)*dytr(j)*omega
            if(hr(i-1,j  ) .ne. c0 .or. hr(i  ,j  ) .ne. c0)
     $        cfn(i,j) = fxb/(h(i-1,j)+h(i,j))
            if(hr(i-1,j-1) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $        cfs(i,j) = fxc/(h(i-1,j-1)+h(i,j-1))
            if(hr(i  ,j  ) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $        cfe(i,j) = fxa*dxtr(i)*dxur(i)/(h(i,j)+h(i,j-1))
            if(hr(i-1,j  ) .ne. c0 .or. hr(i-1,j-1) .ne. c0)
     $        cfw(i,j) = fxa*dxtr(i)*dxur(i-1)/(h(i-1,j)+h(i-1,j-1))
c
            cpf(i,j) = c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
c
            if(hr(i-1,j  ) .ne. c0 .or. hr(i  ,j  ) .ne. c0)
     $        cfn(i,j) = cfn(i,j)
     $        +fxd*dxtr(i)*(hr(i,j)-hr(i-1,j))*sine(j)
            if(hr(i-1,j-1) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $        cfs(i,j) = cfs(i,j)
     $        -fxd*dxtr(i)*(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
            if(hr(i  ,j  ) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $        cfe(i,j) = cfe(i,j)
     $        -fxd*dxtr(i)*(hr(i,j)*sine(j)-hr(i,j-1)*sine(j-1))
            if(hr(i-1,j  ) .ne. c0 .or. hr(i-1,j-1) .ne. c0)
     $        cfw(i,j) = cfw(i,j)
     $        +fxd*dxtr(i)*(hr(i-1,j)*sine(j)-hr(i-1,j-1)*sine(j-1))
            cfn(i,j) = cfn(i,j)*cpf(i,j)
            cfs(i,j) = cfs(i,j)*cpf(i,j)
            cfe(i,j) = cfe(i,j)*cpf(i,j)
            cfw(i,j) = cfw(i,j)*cpf(i,j)
            cof(i,j) = cst(j)*dxt(i)*dyt(j)/cpf(i,j)
            cofis(isle) = cofis(isle)+cof(i,j)
690       continue
          cofis(isle) = c1/cofis(isle)
700     continue
#endif
      endif
c
c-----------------------------------------------------------------------
c     normalize vorticity "ztd"
c-----------------------------------------------------------------------
c
      do 800 j=1,jmt
        do 790 i=1,imt
          ztd(i,j) = ztd(i,j)*cpf(i,j)
790     continue
800   continue
c
c-----------------------------------------------------------------------
c     complete readin of relaxation solution of 2 timesteps prior and
c     of previous timestep.
c-----------------------------------------------------------------------
c
      call oget (kflds, nwds, (luptdb-1)*nwds+1, ptdb)
      call oget (kflds, nwds, (luptd-1)*nwds+1, ptd)
c
c-----------------------------------------------------------------------
c     erstimate a first guess by extrapolating the two
c     previous solutions forward in time.
c-----------------------------------------------------------------------
c
      fxa=c1
      if (mixts .or. mxpas2) fxa=p5
      do 900 j=1,jmt
        do 890 i=1,imt
          ptd(i,j) = fxa*(c2*ptd(i,j)-ptdb(i,j))
890    continue
900   continue
c
c     compute criterion for convergence and set residuals to zero
c
      crtp = crit*fxa
      do 1000 j=1,jmt
        do 990 i=1,imt
          res(i,j) = c0
990     continue
1000  continue
c
c=======================================================================
c     begin section to do the relaxation
c=======================================================================
c
      nwcnt = 0
      nwmax = 300
      mscan = 0
#ifdef testrelx
c    initialize array for storing resmax
      do 999 l=1,mxwres
        wresmx(l) = c0
999   continue
#endif
1001  continue
      mscan = mscan + 1
c
c-----------------------------------------------------------------------
c     consider the arrays as being defined on the squares of an
c     "imt by jmt" checkerboard. take two passes: first solve the
c     equation on the black squares, then on the red squares.
c-----------------------------------------------------------------------
c
      do 1200 icolor=1,2
        do 1190 j=3,jscan
          resmi(j) = c0
          ishypr   = 1 + mod(icolor+j,2)
          if (ishypr .eq. 1) ishypr = 3
          do 1160 i=ishypr,imtm1,2
            res(i,j) = rmsk(i,j)*(cfn(i,j)*ptd(i,j+1)
     $                 +cfs(i,j)*ptd(i,j-1)
     $                 +cfe(i,j)*ptd(i+1,j)
     $                 +cfw(i,j)*ptd(i-1,j)
     $                 -ptd(i,j)-ztd(i,j))
1160      continue
          if (ishypr .eq. 3) ishypr = 1
#ifdef cyclic
          res(1,j)   = res(imum1,j)
          res(imt,j) = res(2,j)
#else
          res(1,j)   = c0
          res(imt,j) = c0
#endif
c
c-----------------------------------------------------------------------
c         make a correction to ptd based on the residuals
c-----------------------------------------------------------------------
c
          do 1170 i=ishypr,imt,2
            ptd(i,j) = ptd(i,j) + sor*res(i,j)
1170      continue
c
c-----------------------------------------------------------------------
c         find the maximum absolute residual to determine convergence
c-----------------------------------------------------------------------
c
#ifdef testrelx
          ires = 0
          jres = 0
#endif
          do 1180 i=ishypr,imt,2
c          do 1180 i=2,imtm1
#ifdef testrelx
            resi = resmi(j)  
#endif
            resmi(j) = max(abs(res(i,j)),resmi(j))
#ifdef testrelx
            if(resmi(j) .ne. resi) ires = i
#endif
1180      continue
1190    continue
1200  continue
c
      resmax = c0
      do 1300 j=3,jscan
#ifdef testrelx
        resj = resmax
#endif
        resmax = max(resmi(j),resmax)
#ifdef testrelx
        if(resmax .ne. resj) jres = j
#endif
1300  continue
#ifdef islands
c
c-----------------------------------------------------------------------
c     do hole relaxation for each island
c-----------------------------------------------------------------------
c
      do 1400 isle=1,nisle
        resis = c0
        do 1360 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          resis = resis + (cfn(i,j)*ptd(i,j+1)
     $                  +  cfs(i,j)*ptd(i,j-1)
     $                  +  cfe(i,j)*ptd(i+1,j)
     $                  +  cfw(i,j)*ptd(i-1,j)
     $                  -  ptd(i,j) - ztd(i,j))*cof(i,j)
1360    continue
c
c       normalize the island residual and update the maximum
c       absolute residual of the relaxation if necessary
c
        resis  = resis*cofis(isle)
#ifdef testrelx
        resmis = resmax
#endif
        resmax = max(abs(resis),resmax)
#ifdef testrelx
        if(resmis .ne. resmax) then
          ires = 99
          jres = 99
        endif
#endif
        resis  = resis*sor
c
c       make a correction to ptd over the island and its perimeter
c
        do 1370 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          ptd(i,j) = ptd(i,j) + resis
1370    continue
1400  continue
#endif
#ifdef testrelx
       wresmx(mscan) = resmax
#endif
#ifdef symmetry
c
c-----------------------------------------------------------------------
c  set symmetry boundary condition
c-----------------------------------------------------------------------
c
      do 1500 i=1,imt
        ptd(i,jmt) = -ptd(i,jmtm1)
1500  continue
#endif
#ifdef cyclic
c
c-----------------------------------------------------------------------
c  set cyclic boundary condition
c-----------------------------------------------------------------------
c
      do 1600 j=1,jmt
        ptd(1,j)   = ptd(imum1,j)
        ptd(imu,j) = ptd(2,j)
1600  continue
#endif
c
c-----------------------------------------------------------------------
c     test maximum residual for convergence of the relaxation.
c     if not converged, proceed with another scan.
c     (if the number of scans reaches mxscan, leave the loop)
c-----------------------------------------------------------------------
c
      if (resmax.ge.crtp .and. mscan.lt.mxscan) go to 1001
c
c=======================================================================
c  end of the relaxation
c=======================================================================
c
c-----------------------------------------------------------------------
c  update the stream function based on the relaxation solution
c-----------------------------------------------------------------------
c
      if (.not. mxpas2) then
        do 1700 j=1,jmt
          do 1690 i=1,imt
            ptdb(i,j) = p(i,j,2) + ptd(i,j)
            p(i,j,2)  = p(i,j,1)
            p(i,j,1)  = ptdb(i,j)
1690      continue
1700    continue
      else
        do 1800 j=1,jmt
          do 1790 i=1,imt
            p(i,j,1) = p(i,j,2) + ptd(i,j)
1790      continue
1800    continue
      endif
c
c-----------------------------------------------------------------------
c  save ptd to compute 1st guess for relaxation next timestep
c  (..note.. on 1st pass of euler backward timestep, bypass this
c            save, since it will be done on the 2nd pass)
c  (..note.. on a mixing timestep, alter ptd to be consistent with
c            normal, leap-frog stepping)
c-----------------------------------------------------------------------
c
      if (mixts .and. eb) return
c
      if (mxpas2 .or. mixts) then
        do 1900 j=1,jmt
          do 1890 i=1,imt
            ptd(i,j) = c2*ptd(i,j)
1890      continue
1900    continue
      endif
c
      call oput (kflds, nwds, (luptdb-1)*nwds+1, ptd)
#endif
      return
      end
