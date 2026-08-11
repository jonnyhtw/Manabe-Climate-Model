      subroutine relax
      implicit real*8 (a-h,o-z)
#ifdef oldrelax
c
c=======================================================================
c                                                                    ===
c  relax takes as input the vorticity driving function computed in   ===
c        "clinic" (ztd) and, using sequential over-relaxation,       ===
c        solves the laplacian equation for the external mode of      ===
c        velocity in terms of a mass transport stream function (p).  ===
c                                                                    ===
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
c---------------------------------------------------------------------
c  demension and equivalence local data
c---------------------------------------------------------------------
c
      dimension ptdb(imt,jmt)
      equivalence (res,ptdb)
c
c---------------------------------------------------------------------
c  begin executable code
c---------------------------------------------------------------------
c
c=======================================================================
c  begin introductory section to prepare for the relaxation  ===========
c=======================================================================
c
c---------------------------------------------------------------------
c  initiate readin of relaxation solution of 2 timesteps back
c  for the purpose of computing an initial guess for the present pass
c  (input unit alternates on tsteps between 4 & 5)
c---------------------------------------------------------------------
c
      luptdb=5-mod(itt,2)
      luptd =4+mod(itt,2)
      call ofind(kflds,nwds,(luptdb-1)*nwds+1)
c
c---------------------------------------------------------------------
c  initialize the work area
c---------------------------------------------------------------------
c
      fx=c0
      do 100 j=1,jmt
        do 90 i=1,imt
          cfn(i,j)=fx
          cfs(i,j)=fx
          cfe(i,j)=fx
          cfw(i,j)=fx
#ifdef islands
          cof(i,j)=c1
#endif
90      continue
100   continue
c
#ifdef kdtest
c---------------------------------------------------------------------
c  calculate depth field from its reciprocal
c  (..note.. h is equivalenced with ptd)
c   !!!actually, this is not needed since "h" is computed in ocean.
c   however this code is kept to duplicate answers.
c---------------------------------------------------------------------
c
c  compute reciprocals with an epsilon added to avoid zero divide
c
      fxa=1.e-20
      do 300 j=2,jmtm1
        do 290 i=1,imt
          h(i,j)=c1/(hr(i,j)+fxa)
290     continue
300   continue
# ifdef islands
c
c  then reset h over land to zero
c
      do 400 j=1,jmt
        do 390 i=1,imt
          if(hr(i,j).eq.c0)h(i,j)=c0
390     continue
400   continue
# endif
# ifdef cyclic
c
c  then set cyclic boundary conditions on h
c
      do 500 j=1,jmt
        h(1,j)=h(imum1,j)
500   continue
# endif
c
#endif
c---------------------------------------------------------------------
c  complete readin of relaxation solution of 2 timesteps back and
c  initiate readin of relaxation solution of previous timestep
c  (input unit alternates on timesteps between 5 & 4)
c---------------------------------------------------------------------
c
      call oget(kflds,nwds,(luptdb-1)*nwds+1,ptdb)
      call ofind(kflds,nwds,(luptd-1)*nwds+1)
c
c---------------------------------------------------------------------
c  generate arrays of coefficients for relaxation
c---------------------------------------------------------------------
c
c  1st, compute coefficients of the laplacian star
c  (hold non-interior points to zero using start and end indices)
c
      do 600 j=3,jmtm1
      do 600 l=1,lseg
        is=isz(j,l)
        if(is.eq.0) go to 600
        ie=iez(j,l)
        fxa=c2*cstr(j)*cstr(j)
        fxb=c2*csu(j  )*cstr(j)*dytr(j)*dyur(j  )
        fxc=c2*csu(j-1)*cstr(j)*dytr(j)*dyur(j-1)
        fx=c1
        do 570 i=is,ie
          cfn(i,j)=fxb/(h(i-1,j)+h(i,j))
          cfs(i,j)=fxc/(h(i-1,j-1)+h(i,j-1))
          cfe(i,j)=fxa*dxur(i)*dxtr(i)/(h(i,j)+h(i,j-1))
          cfw(i,j)=fxa*dxur(i-1)*dxtr(i)/(h(i-1,j)+h(i-1,j-1))
          cpf(i,j)=fx/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
570     continue
c
c  2nd, augment coefficients for implicit treatment of coriolis term
c
        if(acor.ne.c0) then
          fx=-c2dtsf*acor*cstr(j)*dytr(j)*omega
          do 580 i=is,ie
            cfn(i,j)=cfn(i,j)+(hr(i,j  )-hr(i-1,j  ))*sine(j  )
     *                        *fx*dxtr(i)
            cfs(i,j)=cfs(i,j)-(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
     *                        *fx*dxtr(i)
            cfe(i,j)=cfe(i,j)-(hr(i  ,j)*sine(j)-hr(i  ,j-1)*sine(j-1))
     *                        *fx*dxtr(i)
            cfw(i,j)=cfw(i,j)+(hr(i-1,j)*sine(j)-hr(i-1,j-1)*sine(j-1))
     *                        *fx*dxtr(i)
580       continue
        endif
c
c  3rd, normalize coefficients and forcing term for efficiency
c
        do 590 i=is,ie
          cfn(i,j)=cfn(i,j)*cpf(i,j)
          cfs(i,j)=cfs(i,j)*cpf(i,j)
          cfe(i,j)=cfe(i,j)*cpf(i,j)
          cfw(i,j)=cfw(i,j)*cpf(i,j)
          ztd(i,j)=ztd(i,j)*cpf(i,j)
590     continue
600   continue
#ifdef islands
c
c  4th, compute coefficients on island perimeter points
c
      do 700 isle=1,nisle
        cofis(isle)=c0
        do 690 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          fxa=c2*cstr(j)*cstr(j)
          fxb=c2*csu(j  )*dyur(j  )*dytr(j)*cstr(j)
          fxc=c2*csu(j-1)*dyur(j-1)*dytr(j)*cstr(j)
          fxd=-c2dtsf*acor*cstr(j)*dytr(j)*omega
          if(hr(i-1,j  ).ne.c0 .or. hr(i  ,j  ).ne.c0)
     *      cfn(i,j)=fxb/(h(i-1,j)+h(i,j))
     *         +fxd*dxtr(i)*(hr(i,j)-hr(i-1,j))*sine(j)
          if(hr(i-1,j-1).ne.c0 .or. hr(i  ,j-1).ne.c0)
     *      cfs(i,j)=fxc/(h(i-1,j-1)+h(i,j-1))
     *         -fxd*dxtr(i)*(hr(i,j-1)-hr(i-1,j-1))*sine(j-1)
          if(hr(i  ,j  ).ne.c0 .or. hr(i  ,j-1).ne.c0)
     *      cfe(i,j)=fxa*dxtr(i)*dxur(i)/(h(i,j)+h(i,j-1))
     *         -fxd*dxtr(i)*(hr(i,j)*sine(j)-hr(i,j-1)*sine(j-1))
          if(hr(i-1,j  ).ne.c0 .or. hr(i-1,j-1).ne.c0)
     *      cfw(i,j)=fxa*dxtr(i)*dxur(i-1)/(h(i-1,j)+h(i-1,j-1))
     *         +fxd*dxtr(i)*(hr(i-1,j)*sine(j)-hr(i-1,j-1)*sine(j-1))
          cof(i,j)=c1/(cfn(i,j)+cfs(i,j)+cfe(i,j)+cfw(i,j))
          cfn(i,j)=cfn(i,j)*cof(i,j)
          cfs(i,j)=cfs(i,j)*cof(i,j)
          cfe(i,j)=cfe(i,j)*cof(i,j)
          cfw(i,j)=cfw(i,j)*cof(i,j)
          ztd(i,j)=ztd(i,j)*cof(i,j)
          cof(i,j)=cst(j)*dxt(i)*dyt(j)/cof(i,j)
          cofis(isle)=cofis(isle)+cof(i,j)
690     continue
        cofis(isle)=c1/cofis(isle)
700   continue
#endif
c
c  finally, multiply the coefficients by the overrelaxation factor
c
      do 800 j=1,jmt
        do 790 i=1,imt
          cfn(i,j)=cfn(i,j)*sor
          cfs(i,j)=cfs(i,j)*sor
          cfe(i,j)=cfe(i,j)*sor
          cfw(i,j)=cfw(i,j)*sor
790     continue
800   continue
c
c---------------------------------------------------------------------
c  compute a first guess for the relaxation by extrapolating the two
c  previous solutions forward in time.
c---------------------------------------------------------------------
c
c  1st, complete readin of relaxation solution of previous timestep
c
      call oget(kflds,nwds,(luptd-1)*nwds+1,ptd)
c
c  2nd, perform time extrapolation (accounting for mixing timestep)
c
      fxa=c1
      fxb=c2
      if(mixts.or.mxpas2) fxa=p5
      do 900 j=1,jmt
        do 890 i=1,imt
          ptd(i,j)=fxa*(fxb*ptd(i,j)-ptdb(i,j))
890     continue
900   continue
c
c  compute criterion for convergence of relaxation and set residuals 0
c
      crtp=crit*fxa*sor
      do 1000 j=1,jmt
        do 990 i=1,imt
          res(i,j)=c0
990     continue
1000  continue
c
c=======================================================================
c  end introductory section  ===========================================
c=======================================================================
c
c=======================================================================
c  begin section to do the relaxation  =================================
c=======================================================================
c
      mscan=0
10000 mscan=mscan+1
c
c---------------------------------------------------------------------
c  compute entire field of residuals as in simultaneous relaxation
c---------------------------------------------------------------------
c
      do 1100 j=3,jscan
#ifdef kdtest
        do 1090 i=1,imt
#else
        do 1090 i=2,imtm1
#endif
        res(i,j)=cfn(i,j)*ptd(i,j+1)
     *          +cfs(i,j)*ptd(i,j-1)
     *          +cfe(i,j)*ptd(i+1,j)
     *          +cfw(i,j)*ptd(i-1,j)
     *          -sor*(ptd(i,j)+ztd(i,j))
1090    continue
1100  continue
#ifdef islands
c
c---------------------------------------------------------------------
c  reset residuals to zero over land
c---------------------------------------------------------------------
c
      do 1200 j=3,jscan
        do 1170 i=1,imt
          tempik(i,1,1)=c0
1170    continue
        do 1180 l=1,lseg
          is=isz(j,l)
          if(is.eq.0) go to 1181
          ie=iez(j,l)
        do 1180 i=is,ie
          tempik(i,1,1)=res(i,j)
1180     continue
1181   continue
        do 1190 i=1,imt
          res(i,j)=tempik(i,1,1)
1190    continue
1200  continue
#endif
#ifdef cyclic
c
c---------------------------------------------------------------------
c  set cyclic boundary conditions on the residuals
c---------------------------------------------------------------------
c
      do 1300 j=3,jscan
        res(1,j)   = res(imum1,j)
        res(imt,j) = res(2,j)
1300  continue
#else
      do 1300 j=3,jscan
        res(1,j)   = c0
        res(imt,j) = c0
1300  continue
#endif
c
c---------------------------------------------------------------------
c  perform correction on southern point to yield sequential relaxation
c---------------------------------------------------------------------
c
      do 1400 j=3,jscan
#ifdef islands
        do 1360 l=1,lseg
          is=isz(j,l)
          if(is.eq.0) go to 1361
          ie=iez(j,l)
          do 1350 i=is,ie
            res(i,j)=res(i,j)+cfs(i,j)*res(i,j-1)
 1350     continue
 1360   continue
 1361   continue
#else
        do 1370 i=2,imtm1
          res(i,j)=res(i,j)+cfs(i,j)*res(i,j-1)
 1370   continue
#endif
c
c---------------------------------------------------------------------
c  perform correction on western point to yield sequential relaxation
c---------------------------------------------------------------------
c
        do 1390 l=1,lseg
          is=isz(j,l)
          if(is.eq.0) go to 1400
          ie=iez(j,l)
          do 1380 i=is,ie
            res(i,j)=res(i,j)+cfw(i,j)*res(i-1,j)
1380      continue
1390    continue
1400  continue
c
c---------------------------------------------------------------------
c  make a correction to ptd based on the residuals
c---------------------------------------------------------------------
c
      do 1500 j=3,jscan
        do 1490 i=1,imt
          ptd(i,j)=ptd(i,j)+res(i,j)
1490    continue
1500  continue
c
c---------------------------------------------------------------------
c  find the maximum absolute residual to determine convergence
c---------------------------------------------------------------------
c
      resmax=c0
      do 1600 j=3,jscan
        do 1590 i=2,imtm1
          resmax=max(abs(res(i,j)),resmax)
1590    continue
1600  continue
#ifdef islands
c
c---------------------------------------------------------------------
c  do hole relaxation for each island
c---------------------------------------------------------------------
c
      fx=c0
      do 1700 isle=1,nisle
        resis=fx
        do 1660 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
              resis=resis+(cfn(i,j)*ptd(i  ,j+1)
     *                    +cfs(i,j)*ptd(i  ,j-1)
     *                    +cfe(i,j)*ptd(i+1,j  )
     *                    +cfw(i,j)*ptd(i-1,j  )
     *                    -sor*(ptd(i,j)+ztd(i,j)))*cof(i,j)
1660    continue
c
c  normalize the island residual and update the maximum
c  absolute residual of the relaxation if necessary
c
        resis=resis*cofis(isle)
        resmax=max(abs(resis),resmax)
c
c  make a correction to ptd over the island and its perimeter points
c
        do 1690 n=1,nippts(isle)
          i = iperm(iofs(isle)+n)
          j = jperm(iofs(isle)+n)
          ptd(i,j)=ptd(i,j)+resis
1690    continue
1700  continue
#endif
#ifdef symmetry
c
c---------------------------------------------------------------------
c  set symmetry boundary condition
c---------------------------------------------------------------------
c
      do 1800 i=1,imt
        ptd(i,jmt)=-ptd(i,jmtm1)
1800  continue
#endif
#ifdef cyclic
c
c---------------------------------------------------------------------
c  set cyclic boundary condition
c---------------------------------------------------------------------
c
      do 1900 j=1,jmt
        ptd(1,j)=ptd(imum1,j)
        ptd(imu,j)=ptd(2,j)
1900  continue
#endif
c
c---------------------------------------------------------------------
c  test maximum residual for convergence of the relaxation.
c  if not converged, proceed with another scan.
c  (..note.. if the number of scans reaches mxscan, leave the loop)
c---------------------------------------------------------------------
c
      if(resmax.ge.crtp .and. mscan.lt.mxscan) go to 10000
#ifdef kdtest
      write(stdout,1313) crtp, resmax
1313  format(' crtp, resmax =',2e14.6)
#endif
c
c=======================================================================
c  end of the relaxation  ==============================================
c=======================================================================
c
c---------------------------------------------------------------------
c  update the stream function based upon the relaxation solution
c---------------------------------------------------------------------
c
      if(.not. mxpas2) then
        do 2000 j=1,jmt
          do 1990 i=1,imt
            ptdb(i,j) = p(i,j,2)+ptd(i,j)
            p(i,j,2) = p(i,j,1)
            p(i,j,1) = ptdb(i,j)
1990      continue
2000     continue
      else
        do 2100 j=1,jmt
          do 2090 i=1,imt
            p(i,j,1) = p(i,j,2)+ptd(i,j)
2090      continue
2100    continue
      endif
c
c---------------------------------------------------------------------
c  save ptd to compute 1st guess for relaxation next timestep
c  (..note.. on 1st pass of euler backward timestep, bypass this
c            save, since it will be done on the 2nd pass)
c  (..note.. on a mixing timestep, alter ptd to be consistent with
c            normal, leap-frog stepping)
c---------------------------------------------------------------------
c
      if(mixts .and. eb) return
      if(mxpas2 .or. mixts) then
        do 2200 j=1,jmt
          do 2190 i=1,imt
            ptd(i,j)=c2*ptd(i,j)
2190      continue
2200    continue
      endif
      call oput(kflds,nwds,(luptdb-1)*nwds+1,ptd)
#endif
      return
      end
