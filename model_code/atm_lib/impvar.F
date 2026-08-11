      subroutine impvar
c
c  impvar is called by linear for implicit models.
c    entry points:       kernal
c    subroutines called: invert
c    see semi-implicit notes in supersource manual.
c    note  in this subroutine complex variables are represented by
c          real arrays whose first dimension is doubled.
c                    this subroutine computes implicit (bar tau) type
c                    variables and adds the contribution to the total
c                    tendencies due to implicit terms.
c                    then subroutine tderiv can march the integration
c                    ahead in time, as in the purely explicit case.
c                    the assumption that sigma dot=0 at sigma=0,1 is
c                    used in the semi-implicit algorithm.
c
#include "parm.h"
#include "parmkp.h"
      parameter (mtrun1=2*(mtrun+1))
      parameter (mx=mtrun1, jx=jtrun +2)
      parameter (lpmax=isc*mtrun+jtrun+2,mxc=mtrun+1)
      parameter (kxm2=kx-2)
      parameter (kx3=3, kx11a=11, lpmx=lpmax-1)
c
      dimension eigjm(mx,jx)
      dimension vect(mx,jx,kx),gp(mx,jx,kx)
      dimension acoef(kx,kx11a),bcoef(kx,kx11a),tcoef(kx,kx3)
c
      common / imptmp / eigsng(lpmx)
      common / imptmp / dttbar(kx)
c
#include "matrix.h"
c
#include "nxstps.h"
c
#include "param.h"
c
#include "sfharm.h"
c
#include "sgwgt1.h"
c
#include "wimplt.h"
c
#include "zrise.h"
c
      equivalence (gp(1,1,1),vcos(1,1,1)),(vect(1,1,1),ek(1,1,1))
      equivalence (acoef(1,1),acof1(1)),(bcoef(1,1),bcof1(1))
      equivalence (tcoef(1,1),tcof1(1))
      equivalence (dtstep,dt),(dtxdt,dtsq)
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c   compute (phi+rgas*ts*ln(ps)) (at t-dt) + dt*d(phi+rgas*ts*ln(ps))/dt
c   store result temporarily in gp(n,m,k) for later use.
c   first zero gp field to ensure jx index values are zero.
c
      do 10 k=1,kx
      do 10 n=1,jx
      do 10 m=1,mx
 10   gp(m,n,k)=0.0
c
      if (grow) then
c   change topography gradually at the start of a new experiment.
c   also see atmos, grmult, indata, passa to raise topography
        do 7 n=1,jx
         do 7 m=1,mx
          gp(m,n,kx)=grav*rzs(m,n,1)+extkx*(tm(m,n,kx)+dtstep*
     *                 tdt(m,n,kx))+extkxm*(tm(m,n,km)+dtstep*
     *                 tdt(m,n,km))
   7    continue
      else
        do 28 n=1,jx
         do 28 m=1,mx
          gp(m,n,kx)=grav*zs(m,n,1)+extkx*(tm(m,n,kx)+dtstep*
     *                tdt(m,n,kx))+extkxm*(tm(m,n,km)+dtstep*
     *                tdt(m,n,km))
   28   continue
      endif
      gp(1,1,kx)=gp(1,1,kx)+sqrt(2.0)*tfra*(extkx+extkxm)
      do 23 k=1,km
      do 29 n=1,jx
      do 29 m=1,mx
 29   gp(m,n,kx-k)=gp(m,n,kp-k)+factq(kp-k)*(tm(m,n,kp-k)+dtstep*
     1           tdt(m,n,kp-k))+factr(kx-k)*(tm(m,n,kx-k)+dtstep*
     2           tdt(m,n,kx-k))
 23   gp(1,1,kx-k)=gp(1,1,kx-k)+sqrt(2.0)*tfra*(factq(kp-k)+factr(kx-k))
c
c   previous time step mean temperature (tmean) should be used below
c   for consistency. however, tmean is stored for only one timestep.
      do 129 k=1,kx
      do 129 n=1,jx
      do 129 m=1,mx
      vect(m,n,k)=tmean(k)*psm(m,n,1)
      vect(m,n,k)=vect(m,n,k)+dttbar(k)*psdt(m,n,1)
      gp(m,n,k)=gp(m,n,k)+rgas*vect(m,n,k)
  129 continue
c
c   fetch appropriate values of -l(l+1)/a**2 for del square operator
c   and store in eigjm field with first dimension doubled for paired
c   values of real aand imaginary parts
      do 132 n=1,jx
      do 132 m=1,mxc
      lp=n+isc*(m-1)
      if(lp.le.lpmx) then
      eigjm(2*m-1,n)=eigsng(lp)
      eigjm(2*m  ,n)=eigsng(lp)
      else
      eigjm(2*m-1,n) = 0.0
      eigjm(2*m  ,n) = 0.0
      endif
 132  continue
      do 30 k=1,km
      do 30 n=1,jx
      do 30 m=1,mx
      vect(m,n,k)=divm(m,n,k)-divm(m,n,k+1)+
     1             dtstep*(divdt(m,n,k)-divdt(m,n,k+1)
     2             -eigjm(m,n)*(gp(m,n,k)-gp(m,n,k+1)))
   30 continue
      do 31 n=1,jx
      do 31 m=1,mx
      vect(m,n,kx)=(divm(m,n,kx)+dtstep*divdt(m,n,kx))-
     1              dtstep*eigjm(m,n)*gp(m,n,kx)
   31 continue
c
c   now compute the div bar tau spectral coefficients.
c   first zero gp field, which is reused for div bar tau.
      do 32 kprm=1,kx
      do 32 n=1,jx
      do 32 m=1,mx
      gp(m,n,kprm)=0.0
   32 continue
      do 33 kprm=1,kx
      do 33 kdprm=1,kx
      do 33 n=1,jtrun1
      do 33 m=1,mxc
      mr=2*m-1
      mi=2*m
      lp=n+isc*(m-1)
      gp(mr,n,kprm)=gp(mr,n,kprm)+divmtx(lp,kdprm,kprm)*vect(mr,n,kdprm)
      gp(mi,n,kprm)=gp(mi,n,kprm)+divmtx(lp,kdprm,kprm)*vect(mi,n,kdprm)
   33 continue
c
c   compute implicit part of spectral sigma dot.
c   temporarily store sigma dot implicit at level k+1/2 in
c   vect(n,m,k+1), for 1 .le. k .le. kx-1.
c   also, since sigma dot implicit at level k=1/2 is equal to 0,
c   store div wiggle bar tau there.
c   sigma dot implicit at level kx+1/2 = 0, and will not be stored.
c
c   re-initialize vect(n,m,1)
      do 34 n=1,jx
      do 34 m=1,mx
      vect(m,n,1)=0.0
   34 continue
c
c   compute div wiggle bar tau spectral coeff.
      do 36 k=1,kx
      do 36 n=1,jx
      do 36 m=1,mx
      vect(m,n,1)=vect(m,n,1)+dq(k)*gp(m,n,k)
   36 continue
c
c   compute sigma dot implicit at level 3/2
      do 37 n=1,jx
      do 37 m=1,mx
      vect(m,n,2)=dq(1)*(vect(m,n,1)-gp(m,n,1))
   37 continue
c
c   compute sigma dot implicit at levels 5/2,7/2,...kx-1/2
      do 38 k=2,km
      do 38 n=1,jx
      do 38 m=1,mx
      vect(m,n,k+1)=vect(m,n,k)+dq(k)*(vect(m,n,1)-gp(m,n,k))
   38 continue
c
c   add on implicit contributions to psdt and tdt tendencies,
c   then store the true total divergence tendency in divdt.
      do 40 n=1,jx
      do 40 m=1,mx
      tdt(m,n,kx)=tdt(m,n,kx)+tcof2(kx)*vect(m,n,kx)+
     1                        tcof3(kx)*vect(m,n,1)
      psdt(m,n,1)=psdt(m,n,1)-vect(m,n,1)
   40 continue
      do 45 k=2,km
      do 45 n=1,jx
      do 45 m=1,mx
      tdt(m,n,k)=tdt(m,n,k)+tcof2(k)*vect(m,n,k)
   45 continue
      do 145 k=1,km
      do 145 n=1,jx
      do 145 m=1,mx
      tdt(m,n,k)=tdt(m,n,k)+tcof1(k)*vect(m,n,k+1)
      tdt(m,n,k)=tdt(m,n,k)+tcof3(k)*vect(m,n,1)
  145 continue
      do 46 k=1,kx
      do 46 n=1,jx
      do 46 m=1,mx
      divdt(m,n,k)=dstpin*(gp(m,n,k)-divm(m,n,k))
   46 continue
c
c   now subroutine tderiv may be called from atmos.
c   the algorithm there will be the same both
c   for implicit and explicit differencing.
c
      return
c
      entry kernal
c     ===== ======
c
c  kernal computes elements of the a-matrix and then of the p-matrix
c    which is then inverted.  kernal is called from linear at the
c    beginning of every run and whenever tmean is reset (resetm.eq.
c    true).
c
      dstpin=1.0/dtstep
      do 9 k=1,kx
      dttbar(k)=dtstep*tmean(k)
    9 continue
c
c   compute the elements of the matrix. see appendix a of
c   september 1978 notes: semi-implicit differencing method for gfdl
c   spectral models
      do 303 i=1,kx11a
      do 303 k=1,kx
      bcoef(k,i)=0.0
  303 continue
      do 6 k=1,km
      bcof1(k)=-1.0
      bcof6(k)=-1.0
      bcof11(k)=-1.0
      kprm=0
    8 kprm=kprm+1
      bcof2(k)=bcof2(k)-dq(kprm)
      if (k .gt. 1 .and. kprm .lt. k) bcof3(k)=bcof3(k)-dq(kprm)
      if (kprm.lt.k) go to 8
      bcof4(k)=-bcof2(k)
      bcof5(k)=-bcof3(k)
      bcof7(k)=+bcof2(k)-dq(k+1)
      bcof8(k)= bcof2(k)
      bcof9(k)=-bcof7(k)
      bcof10(k)=bcof4(k)
    6 continue
      bcof1(kx)=-1.0
      bcof6(kx)=bcof1(kx)
      bcof11(kx)=bcof1(kx)
c     bcof8(kx)=sum(dq(i)),i=1,kx-2 -- collapses to assignment below
      bcof8(kx)=-qmh(kx-2)
      bcof10(kx)=-bcof8(kx)
      bcof3(kx)=bcof8(kx)-dq(km)
      bcof5(kx)=-bcof3(kx)
      bcof7(kx)=bcof3(kx)
      bcof9(kx)=-bcof3(kx)
c
c   re-initialize tcoef, acoef, dbrcof, and pmatrx.
      do 2 i=1,kx3
      do 2 k=1,kx
      tcoef(k,i)=0.0
    2 continue
      do 3 i=1,kx11a
      do 3 k=1,kx
      acoef(k,i)=0.0
    3 continue
      do 4 k=1,kx
      dbrcof(k)=0.0
      do 4 kprm=1,kx
      pmatrx(kprm,k)=0.0
    4 continue
      do 5 k=1,km
      acof1(k)=rgas*(tmean(k)-tmean(k+1))
      acof2(k)= odq (k)*factr(k)*(tmean(k+1)-tmean(k) )
      if (k .gt. 1) acof3(k)= odq (k)*factr(k)*(tmean(k)-tmean(k-1))
      acof4(k)=akapsg(k)*factr(k)*tmean(k)
      if (k .gt. 1) acof5(k)=acof4(k)
      acof6(k)= cappa*factr(k)*tmean(k)
      if (k .lt. km) acof7(k)=odq(k+1)*factq(k+1)*(tmean(k+2)
     1                       -tmean(k+1))
      acof8(k)= odq (k+1)*factq(k+1)*(tmean(k+1)-tmean(k) )
      acof9(k)=akapsg(k+1)*factq(k+1)*tmean(k+1)
      acof10(k)=acof9(k)
      acof11(k)= cappa*factq(k+1)*tmean(k+1)
    5 continue
      acof1(kx)=rgas*tmean(kx)
      acof3(kx)=extkx* odq (kx)*(tmean(kx)-tmean(kx-1))
      acof5(kx)=extkx*akapsg(kx)*tmean(kx)
      acof6(kx)=extkx* cappa*tmean(kx)
      acof7(kx)=extkxm* odq (kx-1)*(tmean(kx)-tmean(kx-1))
      acof8(kx)=extkxm* odq (kx-1)*(tmean(kx-1)-tmean(kx-2))
      acof9(kx)=extkxm*akapsg(kx-1)*tmean(kx-1)
      acof10(kx)=acof9(kx)
      acof11(kx)=extkxm* cappa*tmean(kx-1)
      do 13 k=1,kx
      if (k .lt. kx) tcof1(k)=akapsg(k) *tmean(k)
     1    -odq(k)*(tmean(k+1)-tmean(k) )
      if (k .gt. 1) tcof2(k)=akapsg(k) *tmean(k)
     1    -odq(k)*(tmean(k)-tmean(k-1))
      tcof3(k) =-cappa*tmean(k)
   13 continue
      do 14 i=1,kx11a
      do 14 k=1,kx
      dbrcof(k)=dbrcof(k)+acoef(k,i)*bcoef(k,i)
   14 continue
c
c
c   compute elements of p-matrix, initially.
c   recompute them whenever tmean(k) is recomputed.
c   the element of row k and column kprm is stored in
c   array pmatrx(kprm,k) = p-matrix sub k sub kprm.
c   the p-matrix is further manipulated and then
c   inverted to solve a matrix equation for the
c   div bar tau (n,m,kprm) spectral coefficients alone.
c   then the matrix product delt*(pmatrx(kprm,k))
c   *(div bar tau(n,m,kprm)) yields the solution to
c   p(kp,m,k)-p(n,m,k+1) for 1 .le. k .le. kx+1, and to
c   p(n,m,kx), where summation over kprm is implied, and
c   p(n,m,k) = gp(n,m,k)+rgas*tmean(k)*(ln(ps(n,m,k))).
c   note that ln(ps), d(ln(ps))/dx, and d(ln(ps))/dy are
c   stored respectively in the arrays ps, psdx, and psdy.
c
c   compute elements of p-matrix.
      do 15 k=1,km
      if (k .le. 1) go to 17
      km1=k-1
      kprm=0
   16 kprm=kprm+1
      pmatrx(kprm,k)=pmatrx(kprm,k)+(acof3(k)-acof5(k))*dq(kprm)
      if (kprm.lt.km1) go to 16
   17 continue
      kprm=0
   18 kprm=kprm+1
      pmatrx(kprm,k)=pmatrx(kprm,k)+(acof2(k)-acof4(k)+acof7(k)+acof8(k)
     1  -acof9(k)-acof10(k))*dq(kprm)
      if (kprm.lt.k) go to 18
      pmatrx(k+1,k)=pmatrx(k+1,k)+(acof7(k)-acof9(k))*dq(k+1)
   15 continue
      do 20 kprm=1,kxm2
      pmatrx(kprm,kx)=pmatrx(kprm,kx)+(acof8(kx)-acof10(kx))*dq(kprm)
   20 continue
      do 21 kprm=1,km
      pmatrx(kprm,kx)=pmatrx(kprm,kx)+(acof3(kx)-acof5(kx)+acof7(kx)
     1  -acof9(kx) )*dq(kprm)
   21 continue
c
c   add on contribution of div wiggle bar tau to pmatrx.
      do 22 k=1,kx
       do 22 kprm=1,kx
      pmatrx(kprm,k)=pmatrx(kprm,k)+dbrcof(k)*dq(kprm)
   22 continue
c
c   compute elements of div bar tau-matrix, then call subroutine invert
c   to obtain the inverse matrix.  then the spectral coefficients
c   div bar tau(jp,mp,k) can be determined.
c   there is a separate kx*kx div bar tau matrix for each value of l.
c   the elements (and later, the inverse elements) are stored in the
c   array divmtr(lp,kprm,k).
c   also, if tmean(k) is a function of time, then divmtr must be
c   recomputed whenever tmean(k) is.
c
c   construct div bar tau matrix elements
      do 125 lp=1,lpmx
 125  eigsng(lp)=-rasq*float(lp*(lp-1))
      do 25 k=1,kx
      do 25 kprm=1,kx
      do 25 lp=1,lpmx
      divmtx(lp,kprm,k)=dtxdt*eigsng(lp)*pmatrx(kprm,k)
   25 continue
      do 26 k=1,km
      kr=k
      do 26 lp=1,lpmx
      divmtx(lp,k,kr)=divmtx(lp,k,kr)+1.0
      divmtx(lp,k+1,kr)=divmtx(lp,k+1,kr)-1.0
   26 continue
      do 27 lp=1,lpmx
      divmtx(lp,kx,kx)=divmtx(lp,kx,kx)+1.0
   27 continue
c
c   invert the matrix divmtx
c
      call invert(divmtx,kx,lpmx)
c
      return
c
      end
