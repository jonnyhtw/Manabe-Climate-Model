      subroutine congr5
      implicit real*8 (a-h,o-z)
#ifdef congrad5pt
c
c=======================================================================
c
c     congr5 solves for the time derivitive of the mass transport
c     stream function using a preconditioned conjugate gradient
c     algorithm with 5 point numerics.
c
c     note: abrupt topographic changes may lead to non convergence. In
c     practice this may only be a problem for coarse resolution
c     model with non smooth topography.
c
c     this algorithm is described in gfdl tech report by j. derber, 1987
c     (see also: gill, p.e., w. murray & m.h. wright, 1981,
c      practical optimizations. academic press)
c
c=======================================================================
c
#include "param.h"
#include "crelax.h"
#include "ctmngr.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "temp.h"
c
c=======================================================================
c     begin introductory section to prepare for the relaxation
c=======================================================================
c
c-----------------------------------------------------------------------
c     initiate readin of solution of 2 timesteps before
c     for the purpose of computing an initial guess for the present
c     solution. (input unit alternates on time steps between 4 & 5)
c-----------------------------------------------------------------------
c
      luptdb = 5 - mod(itt,2)
      luptd  = 4 + mod(itt,2)
      call ofind (kflds, nwds, (luptdb-1)*nwds+1)
      call ofind (kflds, nwds, (luptd-1)*nwds+1)
c
      if (first .or. mixts .or. mod(itt,nmix) .eq. 2) then
c
c-----------------------------------------------------------------------
c       initialize the work area on the first entry. also during
c       a mixing time step and the following time step
c-----------------------------------------------------------------------
c
        do 100 j=1,jmt
          do 90 i=1,imt
            cfn(i,j) = c0
            cfs(i,j) = c0
            cfe(i,j) = c0
            cfw(i,j) = c0
            cpf(i,j) = c0
#ifdef islands
            cof(i,j) = c1
# ifdef oldcongrad
congrad
            ismask(i,j) = 0
# endif
#endif
90        continue
100     continue
#ifdef islands
# ifdef oldcongrad
c
c-----------------------------------------------------------------------
c       form island mask  (ismask=0 over interior ocean points
c                          ismask=1 over perimeter ocean points
c                          ismask=2 over land points)
c-----------------------------------------------------------------------
c
        do 200 j=2,jmtm1
          do 190 i=2,imu
            test1 = hr(i,j)*hr(i-1,j)*hr(i,j-1)*hr(i-1,j-1)
            test2 = hr(i,j)+hr(i-1,j)+hr(i,j-1)+hr(i-1,j-1)
            if (test1. eq. c0) ismask(i,j) = 1
            if (test2 .eq. c0) ismask(i,j) = 2
190       continue
200     continue
# endif
#endif
c
c-----------------------------------------------------------------------
c       generate arrays of coefficients for relaxation
c-----------------------------------------------------------------------
c
c       compute coefficients for the five point del square operator
c       coefficients are zero in land & on perimeter of land masses
c
        do 300 j=3,jmtm1
          do 290 l=1,lseg
            if (isz(j,l) .ne. 0) then
              fxa=c2*cstr(j)*cstr(j)
              fxb=c2*csu(j  )*cstr(j)*dytr(j)*dyur(j  )
              fxc=c2*csu(j-1)*cstr(j)*dytr(j)*dyur(j-1)
              do 260 i=isz(j,l),iez(j,l)
                cfn(i,j)=fxb/(h(i-1,j)+h(i,j))
                cfs(i,j)=fxc/(h(i-1,j-1)+h(i,j-1))
                cfe(i,j)=fxa*dxur(i)*dxtr(i)/(h(i,j)+h(i,j-1))
                cfw(i,j)=fxa*dxur(i-1)*dxtr(i)/(h(i-1,j)+h(i-1,j-1))
                cpf(i,j)=c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
260           continue
c
c             augment coefficients for implicit coriolis treatment
c
              if (acor .ne. c0) then
                fx = -c2dtsf*acor*cstr(j)*dytr(j)*omega
                do 270 i=isz(j,l),iez(j,l)
                  cfn(i,j) = cfn(i,j)+(hr(i,j  )-hr(i-1,j  ))*sine(j  )
     $                          *fx*dxtr(i)
                  cfs(i,j) = cfs(i,j)-(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
     $                          *fx*dxtr(i)
                  cfe(i,j) = cfe(i,j)-(hr(i  ,j)*sine(j)-hr(i  ,j-1)
     $                          *sine(j-1))*fx*dxtr(i)
                  cfw(i,j) = cfw(i,j)+(hr(i-1,j)*sine(j)-hr(i-1,j-1)
     $                          *sine(j-1))*fx*dxtr(i)
270             continue
              endif
c
c             normalize coefficients
c
              do 280 i=isz(j,l),iez(j,l)
                cfn(i,j) = cfn(i,j)*cpf(i,j)
                cfs(i,j) = cfs(i,j)*cpf(i,j)
                cfe(i,j) = cfe(i,j)*cpf(i,j)
                cfw(i,j) = cfw(i,j)*cpf(i,j)
280           continue
            endif
290       continue
300     continue
#ifdef islands
c
c       modify coefficients to include island perimeter points
c
        do 400 isle=1,nisle
          cofis(isle) = c0
# ifdef oldcongrad
congrad
          do 390 j=jsis(isle),jeis(isle)
            fxa = c2*cstr(j)*cstr(j)
            fxb = c2*csu(j  )*dyur(j  )*dytr(j)*cstr(j)
            fxc = c2*csu(j-1)*dyur(j-1)*dytr(j)*cstr(j)
            fxd = -c2dtsf*acor*cstr(j)*dytr(j)*omega
            do 380 i=isis(isle),ieis(isle)
              if(ismask(i,j).eq.1) then
                if(hr(i-1,j  ) .ne. c0 .or. hr(i  ,j  ) .ne. c0)
     $            cfn(i,j) = fxb/(h(i-1,j)+h(i,j))
     $            +fxd*dxtr(i)*(hr(i,j)-hr(i-1,j))*sine(j)
                if(hr(i-1,j-1) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $            cfs(i,j) = fxc/(h(i-1,j-1)+h(i,j-1))
     $            -fxd*dxtr(i)*(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
                if(hr(i  ,j  ) .ne. c0 .or. hr(i  ,j-1) .ne. c0)
     $            cfe(i,j) = fxa*dxtr(i)*dxur(i)/(h(i,j)+h(i,j-1))
     $            -fxd*dxtr(i)*(hr(i,j)*sine(j)-hr(i,j-1)*sine(j-1))
                if(hr(i-1,j  ) .ne. c0 .or. hr(i-1,j-1) .ne. c0)
     $            cfw(i,j) = fxa*dxtr(i)*dxur(i-1)/(h(i-1,j)+h(i-1,j-1))
     $            +fxd*dxtr(i)*(hr(i-1,j)*sine(j)-hr(i-1,j-1)*sine(j-1))
                cpf(i,j) = c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
                cfn(i,j) = cfn(i,j)*cpf(i,j)
                cfs(i,j) = cfs(i,j)*cpf(i,j)
                cfe(i,j) = cfe(i,j)*cpf(i,j)
                cfw(i,j) = cfw(i,j)*cpf(i,j)
                cof(i,j) = cst(j)*dxt(i)*dyt(j)/cpf(i,j)
                cofis(isle) = cofis(isle)+cof(i,j)
              endif
380         continue
390       continue
          cofis(isle) = c1/cofis(isle)
400     continue
c     
        do 500 j=1,jmt
          do 490 i=1,imt
            hdir(i,j) = cof(i,j)
490       continue
500     continue
c
        do 600 j=1,jmt
          do 590 i=1,imt
            if (ismask(i,j) .ne. 0) then
              cof(i,j) = c0
            else
              cof(i,j) = c1
            endif
590       continue
600     continue
c
# else
          do 390 n=1,nippts(isle)
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
390       continue
          cofis(isle) = c1/cofis(isle)
400     continue
# endif
c
# ifdef oldcongrad
congrad
        do 700 isle=1,nisle
          do 690 j=jsis(isle),jeis(isle)
            do 680  i=isis(isle),ieis(isle)
              if (ismask(i,j) .eq. 1) cof(i,j) = hdir(i,j)*cofis(isle)
680          continue
690       continue
700     continue
# else
        do 700 isle=1,nisle
          do 690 n=1,nippts(isle)
            i = iperm(iofs(isle)+n)
            j = jperm(iofs(isle)+n)
            cof(i,j) = cof(i,j)*cofis(isle)
690       continue
700     continue
# endif
# ifdef cyclic
c
        do 800 j=1,jmt
          cof(1,j)   = c0
          cof(imt,j) = c0
800     continue
# endif
#endif
      endif
c
c-----------------------------------------------------------------------
c     normalize vorticity "ztd"
c-----------------------------------------------------------------------
c
      do 900 j=1,jmt
        do 890 i=1,imt
          ztd(i,j) = ztd(i,j)*cpf(i,j)
890     continue
900   continue
c
c-----------------------------------------------------------------------
c     complete readin of relaxation solution of 2 timesteps prior and
c     of previous timestep.
c-----------------------------------------------------------------------
c
      call oget (kflds, nwds, (luptdb-1)*nwds+1, res)
      call oget (kflds, nwds, (luptd-1)*nwds+1, ptd)
c
c-----------------------------------------------------------------------
c     erstimate a first guess by extrapolating the two
c     previous solutions forward in time.
c-----------------------------------------------------------------------
c
      fxa=c1
      if (mixts .or. mxpas2) fxa=p5
      do 1000 j=1,jmt
        do 990 i=1,imt
          ptd(i,j) = fxa*(c2*ptd(i,j)-res(i,j))
990     continue
1000  continue
#ifdef cyclic
c
c     set cyclic bc
c
      do 1110 j=2,jmtm1
        ptd(1,j)   = ptd(imtm1,j)
        ptd(imt,j) = ptd(2,j)
1110  continue
#endif
c
c-----------------------------------------------------------------------
c     compute convergence criterion & find initial residuals
c-----------------------------------------------------------------------
c
      crtp = crit*fxa
      do 1200 j=1,jmt
        do 1190 i=1,imt
          res(i,j) = c0
1190    continue
1200  continue
c
      do 1300 j=2,jmtm1
        do 1290 i=2,imtm1
          res(i,j) = cfn(i,j)*ptd(i,j+1) +
     $               cfs(i,j)*ptd(i,j-1) +
     $               cfe(i,j)*ptd(i+1,j) +
     $               cfw(i,j)*ptd(i-1,j) -
     $               ptd(i,j) - ztd(i,j)
#ifdef islands
          res(i,j) = res(i,j)*cof(i,j)
#endif
1290    continue
1300  continue
#ifdef islands
# ifdef oldcongrad
congrad
c
c     find initial island residual using stokes theorem for the
c     line integral. update res(i,j) on perimeter & interior of island
c
      do 1400 isle=1,nisle
        resis = c0
        do 1370 j=jsis(isle),jeis(isle)
          do 1360 i=isis(isle),ieis(isle)
            if (ismask(i,j) .eq. 1) resis = resis + res(i,j)
1360      continue
1370    continue
        do 1390 j=jsis(isle),jeis(isle)
          do 1380 i=isis(isle),ieis(isle)
            if (ismask(i,j) .ge. 1) res(i,j) = resis
1380      continue
1390    continue
1400  continue
#else
c
c     find initial island residual using a line integral.
c      update res(i,j) on perimeter & interior of island
c
      do 1400 isle=1,nisle
        resis = c0
        do 1360 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          resis = resis + res(i,j)
1360    continue
        do 1370 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          res(i,j) = resis
1370    continue
1400  continue
# endif
#endif
#ifdef cyclic
c
c     set cyclic bc
c
      do 1500 j=2,jmtm1
        res(1,j)   = res(imtm1,j)
        res(imt,j) = res(2,j)
1500  continue
#endif
#ifdef symmetry
c
c     set symmetry bc
c
      do 1600 i=1,imt
        res(i,jmt) = - res(i,jmtm1)
1600  continue
#endif
c
c     initialize beta & gdir
c
      beta = c0
      do 1700 j=1,jmt
        do 1690 i=1,imt
          gdir(i,j) = c0
1690    continue
1700  continue
c
c-----------------------------------------------------------------------
c     begin iteration loop using conjugate gradient technique
c-----------------------------------------------------------------------
c
      mscan=0
#ifdef testrelx
c    initialize array for storing resmax
      do 999 l=1,mxwres
        wresmx(l) = c0
999   continue
#endif
1701  continue
      mscan = mscan + 1
c
c-----------------------------------------------------------------------
c     compute field of residuals
c-----------------------------------------------------------------------
c
      do 1800 j=2,jmtm1
        do 1790 i=2,imtm1
          gdir(i,j) = cfn(i,j)*res(i,j+1) +
     $                cfs(i,j)*res(i,j-1) +
     $                cfe(i,j)*res(i+1,j) +
     $                cfw(i,j)*res(i-1,j) -
     $                res(i,j)
#ifdef islands
          gdir(i,j) = gdir(i,j)*cof(i,j)
#endif
1790    continue
1800  continue
#ifdef islands
# ifdef oldcongrad
congrad
c     find island residuals using stokes theorem for the
c     line integral. update gdir on perimeter & interior of island
c
      do 1900 isle=1,nisle
        resis = c0
        js    = jsis(isle)
        je    = jeis(isle)
        is    = isis(isle)
        ie    = ieis(isle)
        do 1870 j=js,je
          do 1860 i=is,ie
            if (ismask(i,j) .eq. 1) resis = resis + gdir(i,j)
1860      continue
1870    continue
        do 1890 j=js,je
          do 1880 i=is,ie
            if (ismask(i,j) .ge. 1) gdir(i,j) = resis
1880      continue
1890    continue
1900  continue
# else
c
c     find island residuals using a line integral.
c     update gdir on perimeter & interior of island
c
      do 1900 isle=1,nisle
        resis = c0
        do 1860 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          resis = resis + gdir(i,j)
1860    continue
        do 1870 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          gdir(i,j) = resis
1870    continue
1900  continue
# endif
#endif
#ifdef cyclic
c
c     set cyclic condition
c
      do 2000 j=2,jmtm1
        gdir(1,j)   = gdir(imtm1,j)
        gdir(imt,j) = gdir(2,j)
2000  continue
#endif
c
      do 2100 j=2,jmtm1
        asum(j) = c0
        do 2090 i=2,imtm1
          asum(j) = asum(j) + gdir(i,j)*res(i,j)
2090    continue
2100  continue
c
      apr4 = c0
      do 2200 j=2,jmtm1
        apr4 = apr4 + asum(j)
2200  continue
c
      if (mscan .gt. 1) then
        if (oapr .ne. c0) then
          beta = apr4/oapr
        else
          beta = c0
        endif
      endif
c
c-----------------------------------------------------------------------
c     calculate new values for hdir and fdir
c-----------------------------------------------------------------------
c
      do 2300 j=2,jmtm1
        do 2290 i=1,imt
          hdir(i,j)   = res(i,j)  + beta*hdir(i,j)
          fdir(i,j)   = gdir(i,j) + beta*fdir(i,j)
2290    continue
2300  continue
c
      oapr = apr4
      apr2 = c0
      do 2400 j=2,jmtm1
        asum(j) = c0
        do 2390 i=2,imtm1
          asum(j) = asum(j) + fdir(i,j)*fdir(i,j)
2390    continue
2400  continue
c
      do 2500 j=2,jmtm1
        apr2 = apr2 + asum(j)
2500  continue
c
c-----------------------------------------------------------------------
c     calculate stepsize (asor) & make a correction to ptd res
c-----------------------------------------------------------------------
c
      if (apr2 .ne. c0) then
        asor = - apr4/apr2
      else
        asor = c0
      endif
c
      do 2600 j=2,jmtm1
        do 2590 i=1,imt
          ptd(i,j) = ptd(i,j) + asor*hdir(i,j)
          res(i,j) = res(i,j) + asor*fdir(i,j)
2590    continue
2600  continue
#ifdef cyclic
c
c     set cyclic condition
c
      do 2700 j=2,jmtm1
        res(1,j)   = res(imtm1,j)
        res(imt,j) = res(2,j)
2700  continue
#endif
c
c-----------------------------------------------------------------------
c     find the maximum absolute residual to determine convergence
c-----------------------------------------------------------------------
c
#ifdef testrelx
      ires = 0
      jres = 0
#endif
      resmax = c0
      do 2800 j=3,jscan
        do 2790 i=2,imtm1
#ifdef testrelx
          resmij = resmax
#endif
          resmax = max(abs(res(i,j)), resmax)
#ifdef testrelx
          if(resmij .ne. resmax) then
            ires = i
            jres = j
          endif
#endif
2790    continue
2800  continue
#ifdef testrelx
      wresmx(mscan) = resmax
#endif
#ifdef symmetry
c
c     set symmetry on northern wall
c
      do 2900 i=1,imt
        res(i,jmt) = - res(i,jmtm1)
2900  continue
#endif
c
c-----------------------------------------------------------------------
c      take another iteration unless the max residual < crtp or
c      the number of scans reaches mxscan
c-----------------------------------------------------------------------
c
      if (resmax .ge. crtp .and. mscan .lt. mxscan) go to 1701
c
c=======================================================================
c     end of the iteration loop
c=======================================================================
c
#ifdef symmetry
c
c     set symmetry on northern wall
c
      do 3000 i=1,imt
        ptd(i,jmt) = - ptd(i,jmtm1)
3000  continue
#endif
c-----------------------------------------------------------------------
c     update the stream function based upon the relaxation solution
c-----------------------------------------------------------------------
c
      if (.not. mxpas2) then
        do 3100 j=1,jmt
          do 3090 i=1,imt
            res(i,j) = p(i,j,2) + ptd(i,j)
            p(i,j,2)  = p(i,j,1)
            p(i,j,1)  = res(i,j)
3090      continue
3100    continue
      else
        do 3200 j=1,jmt
          do 3190 i=1,imt
            p(i,j,1) = p(i,j,2) + ptd(i,j)
3190      continue
3200    continue
      endif
c
c-----------------------------------------------------------------------
c     save ptd to compute 1st guess for relaxation next timestep
c     (..note.. on 1st pass of euler backward timestep, bypass this
c            save, since it will be done on the 2nd pass)
c     (..note.. on a mixing timestep, alter ptd to be consistent with
c            normal, leap-frog stepping)
c-----------------------------------------------------------------------
c
      if (mixts .and. eb) return
c
      if (mxpas2. or. mixts) then
        do 3300 j=1,jmt
          do 3290 i=1,imt
            ptd(i,j)=c2*ptd(i,j)
3290      continue
3300    continue
      endif
c
      call oput (kflds, nwds, (luptdb-1)*nwds+1, ptd)
c
#endif
      return
      end
