      subroutine isop0 (j)
      implicit real*8 (a-h,o-z)
c
#ifdef isopycmix
c=======================================================================
c
c     an adaptation of M. Cox's isopycnal diffusion code, after M. Redi
c     (this portion is called from subroutine step)
c
c     isopycmix does not yet work with microtasking, ntasks > 1, or
c     skipland
c     (because of need to have proper "esav" for all jsrows)
c
c   References:
c     Cox, M.D., 1987, Isopycnal diffusion in a z-coordinate model,
c       Ocean Modelling, issue 74, 1-5.
c     Redi, M.H., 1982, Oceanic isopycnal mixing by coordinate rotation,
c       J. Phys. Oceanogr., vol 12, 1154-1158, 1982.
c
c=======================================================================
c
#include "param.h"
#include "accel.h"
#include "chmix.h"
#include "cisop.h"
#include "ctask.h"
#include "cvmix.h"
#include "grdvar.h"
#include "slabs.h"
#include "temp.h"
c
      parameter ( c1e10 = 1.0e10 )
      parameter ( eps = 1.0e-25 )
c
c    "rx", "ry" & "rz" terms = store density differences used to compute
c               isopycnal diffusion tensor components "fk#"
c    "e" = density gradients (1=zonal; 2=meridional; 3=vertical)
c    "slmxr" = reciprocal of maximum allowable slope (limit mixing to
c        avoid numerical instabilities)
c    see "isop.h" for more details on variables listed above
c
c---------------------------------------------------------------------
c
c     compute five of the nine components of isopycnal diffusion tensor
c     and load coefficients into "fk..." arrays
c     components (1,1) & (2,2) of tensor are assumed = 1.0 (xx & yy)
c     components (1,2) & (2,1) of tensor are assumed = 0.0 (yx & xy)
c
c     scale density gradients by c1e10 to keep exponents in range later
c
c     numerical constraint: slmxr > sqrt(4*ahisop*hdcsck(k)*dt/dx*dz)
c     ("slmxr" can generally be 100 or less)
c
c---------------------------------------------------------------------
c
      if (j .eq. jsrow) then
c
c     for just the first jrow of a task, initialize the three "fk3"
c     coefficients (defined at the top of level one) to zero
c     (no diffusion of tracers through sea surface)
c
        do 100 i=1,imt
          fk3(i,1,1) = c0
          fk3(i,1,2) = c0
          fk3(i,1,3) = c0
100     continue
c
      else
c
c     for all tracer rows but the first row of a task, shift slabwise
c     rows of density differences
c
        do 200 k=1,km
          do 190 i=1,imt
            rx (i,k) = rxp(i,k)
            rys(i,k) = ryn(i,k)
            rz (i,k) = rzp(i,k)
190       continue
200     continue
c
        do 300 i=1,imt
          rz (i,kmp1) = rzp(i,kmp1)
300     continue
c
      endif
c
c     for all jrows...
c     load densities from row j+1 into tempik(i,k,1) & tempik(i,k,2)
c
      call statec ( t(1,1,jp,nm,1), t(1,1,jp,nm,2), tempik(1,1,1), 1)
      call statec ( t(1,1,jp,nm,1), t(1,1,jp,nm,2), tempik(1,1,2), 2)
c
c     calculate vertical density differences for j+1 row and load
c     into "rzp" (mult by phantom 1.e10)
c
c     first, subtract even numbered levels from the level above them
c
      do 400 k=2,km,2
        do 390 i=1,imt
          rzp(i,k) = fm(i,k-1,jp)*fm(i,k,jp)*c1e10*
     $               (tempik(i,k-1,1) - tempik(i,k,1))
390     continue
400   continue
c
c     second, subtract odd numbered levels from the level above them
c
      do 500 k=3,km,2
        do 490 i=1,imt
          rzp(i,k) = fm(i,k-1,jp)*fm(i,k,jp)*c1e10*
     $               (tempik(i,k-1,2) - tempik(i,k,2))
490     continue
500   continue
c
c     set vertical density differences to zero at top and bottom
c     interfaces
c
      do 600 i=1,imt
        rzp(i,1)    = c0
        rzp(i,kmp1) = c0
600   continue
c
c     load densities at center row j into tempik(i,k,1)
c     (compatible densities at row j+1 are still in tempik(i,k,2))
c
      call statec ( t(1,1,jc,nm,1), t(1,1,jc,nm,2), tempik(1,1,1), 2)
c
c     load east-center density differences at row j+1 into "rxp"
c     load north-center density differences (row j+1 - row j) into "ryn"
c     (multply "rxp" & "ryn"  by phantom 1.e10)
c
      do 700 k=1,km
        do 690 i=1,imtm1
          rxp(i,k) = fm(i+1,k,jp)*fm(i,k,jp)*c1e10*
     $               (tempik(i+1,k,2) - tempik(i,k,2))
          ryn(i,k) = fm(i  ,k,jp)*fm(i,k,jc)*c1e10*
     $               (tempik(i  ,k,2) - tempik(i,k,1))
690     continue
700   continue
c
      do 800 k=1,km
#ifdef cyclic
        rxp(imt,k) = rxp(2,k)
        ryn(imt,k) = ryn(2,k)
#else
        rxp(imt,k) = c0
        ryn(imt,k) = c0
#endif
800   continue
c
c     if this is the first row of a task, set up is finished
c     otherwise, continue with calculation of mixing tensors
c
      if (j .eq. jsrow .and. j .eq. 1) then
        do 900 m=1,nt
          do 890 k=1,km
            do 880 i=1,imt
              esav(i,k,m) = c0
880         continue
890       continue
900     continue
      endif

c NOTE: isopycmix does not yet work with microtasking yet
c       or for code using more than 1 tasks, or skipland option
c       (because of need to have proper "esav" for all jsrows)
c
c---------------------------------------------------------------------
c
       if (j .eq. jsrow) return
c
c---------------------------------------------------------------------
c
c     convert density differences into density gradients and load into
c     "e" to allow calculation of "fk1(3)" mixing component (xz).
c     resolved at center of eastern face of t,s box; at grid point
c     level "zt"; as indicated by "X" below.
c            lat/lon view               depth/lon view
c        u,v j   -> +-------+           +-------+ <- zw(k-1)
c                   |       |           |       |	
c        t,s j   -> |   i   X           |   i   X <- zt(k)
c                   |       |           |       |
c        u,v j-1 -> +-------+           +-------+ <- zw(k)
c
      do 1000 k=1,km
        do 990 i=1,imtm1
          e(i,k,1) = rx(i,k)*dxur(i)*cstr(j)
          e(i,k,2) = dyt4r(j)*
     $              (rys(i,k) + rys(i+1,k) + ryn(i,k) + ryn(i+1,k))
          e(i,k,3) = dzt2r(k)*p5*
     $              (rz(i,k) + rz(i+1,k) + rz(i,k+1) + rz(i+1,k+1))
990     continue
1000  continue
c
      do 1100 k=1,km
#ifdef cyclic
        e(imt,k,1) = e(2,k,1)
        e(imt,k,2) = e(2,k,2)
        e(imt,k,3) = e(2,k,3)
#else
        e(imt,k,1) = c0
        e(imt,k,2) = c0
        e(imt,k,3) = c0
#endif
1100  continue
c
c     use slmxr to limit slope of mixing in the vertical to guard
c     against numerical instabilities and compute the "fk1(3)"
c     mixing coefficient (zx).
c     "chkslp" = criterion used to check slope steepness
c
      do 1200 k=1,km
        do 1190 i=1,imt
          chkslp = -sqrt (e(i,k,1)**2 + e(i,k,2)**2) * slmxr * dtxsqr(k)
          if (e(i,k,3) .gt. chkslp) e(i,k,3) = chkslp
          fk1(i,k,3) = (-e(i,k,1)*e(i,k,3)*ahisop*hdcsck(k)) / 
     $                 (e(i,k,3)**2 + eps)
1190    continue
1200  continue
c
c     convert density differences into density gradients and load into
c     "e" to allow calculation of "fk2(3)" mixing component. (zy)
c     resolved at center of northern face of t,s box; at grid point
c     level "zt"; as indicated by "X" below.
c            lat/lon view               depth/lon view
c        u,v j   -> +---X---+           +-------+ <- zw(k-1)
c                   |       |           |       |
c        t,s j   -> |   i   |           |   X   | <- zt(k)
c                   |       |           |       |
c        u,v j-1 -> +-------+           +-------+ <- zw(k)
c
      do 1300 k=1,km
        do 1290 i=2,imt
          e(i,k,1) = dxt4r(i)*csur(j)*
     $              (rx(i,k) + rx(i-1,k) + rxp(i,k) + rxp(i-1,k))
          e(i,k,2) = ryn(i,k)*dyur(j)
          e(i,k,3) = dzt2r(k)*p5*
     $              (rz(i,k) + rz(i,k+1) + rzp(i,k) + rzp(i,k+1))
1290    continue
1300  continue
c
      do 1400 k=1,km
#ifdef cyclic
        e(1,k,1) = e(imtm1,k,1)
        e(1,k,2) = e(imtm1,k,2)
        e(1,k,3) = e(imtm1,k,3)
#else
        e(1,k,1) = c0
        e(1,k,2) = c0
        e(1,k,3) = c0
#endif
1400  continue
c
c     use slmxr to limit slope of mixing in the vertical to guard
c     against numerical instabilities and compute the "fk2(3)"
c     mixing coefficient (zy).
c     "chkslp" = criterion used to check slope steepness
c
      do 1500 k=1,km
        do 1490 i=1,imt
          chkslp = -sqrt (e(i,k,1)**2 + e(i,k,2)**2) * slmxr * dtxsqr(k)
          if (e(i,k,3) .gt. chkslp) e(i,k,3) = chkslp
          fk2(i,k,3) = (-e(i,k,2)*e(i,k,3)*ahisop*hdcsck(k)) /
     $                 (e(i,k,3)**2 + eps)
1490    continue
1500  continue
c
c     convert density differences into density gradients and load into
c     "e" to allow calculation of the three "fk3(1:3)" tensor
c     components resolved at the center of upper face of t,s box; at
c     depth of "zw(k-1)" as indicated by "X" below.
c            lat/lon view               depth/lon view
c        u,v j   -> +-------+           +---X---+ <- zw(k-1)
c                   |       |           |       |
c        t,s j   -> |   X   |           |   i   | <- zt(k)
c                   |       |           |       |
c        u,v j-1 -> +-------+           +-------+ <- zw(k)
c
c     calculate density gradients at the top of all levels except k=1
c
      do 1600 k=2,km
        do 1590 i=2,imt
          e(i,k,1) = dxt4r(i)*cstr(j)*
     $              (rx(i,k-1) + rx(i-1,k-1) + rx(i,k) + rx(i-1,k))
          e(i,k,2) = dyt4r(j)*
     $              (rys(i,k) + rys(i,k-1) + ryn(i,k) + ryn(i,k-1))
          e(i,k,3) = dzwr(k-1)*rz(i,k)
1590    continue
1600  continue
c
      do 1700 k=2,km
#ifdef cyclic
        e(1,k,1) = e(imtm1,k,1)
        e(1,k,2) = e(imtm1,k,2)
        e(1,k,3) = e(imtm1,k,3)
#else
        e(1,k,1) = c0
        e(1,k,2) = c0
        e(1,k,3) = c0
#endif
1700  continue
c
c     use slmxr to limit slope of mixing in the vertical to guard
c     against numerical instabilities and compute the three "fk3"
c     mixing coefficients ("fk3(1)"=xz, "fk3(2)"=yz, "fk3(3)"=zz)
c     "chkslp" = criterion used to check slope steepness
c
      do 1800 k=2,km
        do 1790 i=1,imt
          chkslp = -sqrt (e(i,k,1)**2 + e(i,k,2)**2) * slmxr * dtxsqr(k)
          if (e(i,k,3) .gt. chkslp) e(i,k,3) = chkslp
1790    continue
1800  continue
c
c       "fk3(i,k,1)" and "fk3(i,k,2)" go to zero for uniform horizontal
c     densities.  otherwise, "fk3(i,k,1)" represents mixing across the
c     eastern face of the t,s box due to sloping isopycnals (xz) and
c     "fk3(i,k,2)" represents mixing across the northern face of the
c     t,s box due to sloping isopycnals (yz).
c       "fk3(i,k,3)" includes diffusion in the vertical due to sloping
c     isopycnals (zz).
c     later, in tracer, "fk3(i,k,3)" will be added into the vertical
c     (or almost diapycnic) mixing represented by the "vdc" variable
c     to get the full zz vertical mixing contribution.
c     (values for fk3 at k=1 were set at the top of this routine)
c
      do 1900 k=2,km
        do 1890 i=1,imt
          ahfctr = ahisop*hdcsck(k) / (e(i,k,3)**2 + eps)
          fk3(i,k,1) = -e(i,k,3)*e(i,k,1)*ahfctr
          fk3(i,k,2) = -e(i,k,3)*e(i,k,2)*ahfctr
          fk3(i,k,3) = (e(i,k,1)**2 + e(i,k,2)**2)*ahfctr
1890    continue
1900  continue
c
#endif
      return
      end
      subroutine isop1 (j, mtrace)
      implicit real*8 (a-h,o-z)
c
#ifdef isopycmix
c=======================================================================
c
c     an adaptation of M. Cox's isopycnal diffusion code, after M. Redi
c     (this portion is called from subroutine tracer)
c
c=======================================================================
c
#include "param.h"
#include "chmix.h"
#include "cisop.h"
#include "cvmix.h"
#include "grdvar.h"
#include "slabs.h"
#include "temp.h"
c
c----------------------------------------------------------------------
c  compute isopycnal diffusion of tracers (except zz term)
c----------------------------------------------------------------------
c
      m = mtrace
c
c     use tracer gradients in the east - central (e(1)) and in the
c     north - central (e(2)) directions at all levels to calculate
c     effects of purely horizontal diffusion (uses the sum of the
c     along isopycnal and background horizontal diffusion coefficients
c     as the total diffusion coefficient)
c
      do 100 k=1,km
        do 90 i=1,imtm1
          e(i,k,1) = (ah + ahisop*hdcsck(k))*dxur(i)*cstr(j)*
     $               (t(i+1,k,jc,nm,m) - t(i,k,jc,nm,m))
          e(i,k,2) = (ah + ahisop*hdcsck(k))*dyur(j)*
     $               (t(i  ,k,jp,nm,m) - t(i,k,jc,nm,m))
90      continue
100   continue
c
c     for all but the top and bottom levels, add in the effects of
c     the along isopycnal diffusion computed using the "fk1" and "fk2"
c     mixing tensors (see subroutine isop0 for details)
c
      do 200 k=2,kmm1
        do 190 i=1,imtm1
          e(i,k,1) = e(i,k,1) + fk1(i,k,3)*dzt2r(k)*p5*
     $       ( fm(i  ,k-1,jc)*(t(i  ,k-1,jc,nm,m) - t(i  ,k  ,jc,nm,m))
     $       + fm(i+1,k-1,jc)*(t(i+1,k-1,jc,nm,m) - t(i+1,k  ,jc,nm,m))
     $       + fm(i  ,k+1,jc)*(t(i  ,k  ,jc,nm,m) - t(i  ,k+1,jc,nm,m))
     $       + fm(i+1,k+1,jc)*(t(i+1,k  ,jc,nm,m) - t(i+1,k+1,jc,nm,m)))
          e(i,k,2) = e(i,k,2) + fk2(i,k,3)*dzt2r(k)*p5*
     $       ( fm(i  ,k-1,jc)*(t(i  ,k-1,jc,nm,m) - t(i  ,k  ,jc,nm,m))
     $       + fm(i  ,k+1,jc)*(t(i  ,k  ,jc,nm,m) - t(i  ,k+1,jc,nm,m))
     $       + fm(i  ,k-1,jp)*(t(i  ,k-1,jp,nm,m) - t(i  ,k  ,jp,nm,m))
     $       + fm(i  ,k+1,jp)*(t(i  ,k  ,jp,nm,m) - t(i  ,k+1,jp,nm,m)))
190     continue
200   continue
c
c     do the special case of the top and bottom levels for e(1) & e(2).
c
      do 300 i=1,imtm1
        e(i, 1,1) = e(i, 1,1) + fk1(i, 1,3)*dzt2r( 1)*p5*
     $      ( fm(i  ,   2,jc)*(t(i  ,   1,jc,nm,m) - t(i  , 2,jc,nm,m))
     $      + fm(i+1,   2,jc)*(t(i+1,   1,jc,nm,m) - t(i+1, 2,jc,nm,m)))
        e(i,km,1) = e(i,km,1) + fk1(i,km,3)*dzt2r(km)*p5*
     $      ( fm(i  ,kmm1,jc)*(t(i  ,kmm1,jc,nm,m) - t(i  ,km,jc,nm,m))
     $      + fm(i+1,kmm1,jc)*(t(i+1,kmm1,jc,nm,m) - t(i+1,km,jc,nm,m)))
        e(i, 1,2) = e(i, 1,2) + fk2(i, 1,3)*dzt2r( 1)*p5*
     $      ( fm(i  ,   2,jc)*(t(i  ,  1,jc,nm,m) - t(i  ,  2,jc,nm,m))
     $      + fm(i  ,   2,jp)*(t(i  ,  1,jp,nm,m) - t(i  ,  2,jp,nm,m)))
        e(i,km,2) = e(i,km,2) + fk2(i,km,3)*dzt2r(km)*p5*
     $      ( fm(i  ,kmm1,jc)*(t(i  ,kmm1,jc,nm,m) - t(i  ,km,jc,nm,m))
     $      + fm(i  ,kmm1,jp)*(t(i  ,kmm1,jp,nm,m) - t(i  ,km,jp,nm,m)))
300   continue
c
c     apply land/sea masks: set mixing to zero if central point or
c     adjacent one is land.
c     "e(i,k,1)" represents diffusion at the eastern face arising from
c     the sloping isopycnal surfaces and from purely horizontal mixing
c     "e(i,k,2)" represents diffusion at the northern face arising from
c     the sloping isopycnal surfaces and from purely horizontal mixing
c
      do 400 k=1,km
        do 390 i=1,imtm1
          e(i,k,1) = fm(i+1,k,jc)*fm(i,k,jc)*e(i,k,1)
          e(i,k,2) = fm(i  ,k,jp)*fm(i,k,jc)*e(i,k,2)*csu(j)
390     continue
400   continue
c
c     e(i,k,3) represents diffusion at the top face arising from the
c     sloping isopycnal surfaces (xz + yz terms)
c
      do 500 k=2,km
        do 490 i=2,imtm1
          e(i,k,3) = fm(i  ,k-1,jc)*fm(i  ,k,jc)*
     $             ( fk3(i,k,1)*dxt4r(i)*cstr(j)*
     $       ( fm(i-1,k  ,jc)*(t(i  ,k  ,jc,nm,m) - t(i-1,k  ,jc,nm,m))
     $       + fm(i-1,k-1,jc)*(t(i  ,k-1,jc,nm,m) - t(i-1,k-1,jc,nm,m))
     $       + fm(i+1,k  ,jc)*(t(i+1,k  ,jc,nm,m) - t(i  ,k  ,jc,nm,m))
     $       + fm(i+1,k-1,jc)*(t(i+1,k-1,jc,nm,m) - t(i  ,k-1,jc,nm,m)))
     $             + fk3(i,k,2)*dyt4r(j)*
     $       ( fm(i  ,k  ,jm) * (t(i  ,k  ,jc,nm,m)-t(i  ,k  ,jm,nm,m))
     $       + fm(i  ,k-1,jm) * (t(i  ,k-1,jc,nm,m)-t(i  ,k-1,jm,nm,m))
     $       + fm(i  ,k  ,jp) * (t(i  ,k  ,jp,nm,m)-t(i  ,k  ,jc,nm,m))
     $       + fm(i  ,k-1,jp) * (t(i  ,k-1,jp,nm,m)-t(i  ,k-1,jc,nm,m)))
     $              )
490     continue
500   continue
c
c     set the values of e(3) to zero for level 1 and kmp1
c     for zero flux conditions through the surface and bottom
c
      do 600 i=1,imt
        e(i,   1,3) = c0
        e(i,kmp1,3) = c0
600   continue
c
      do 700 k=1,km
#ifdef cyclic
        e(imt,k,1) = e(2,k,1)
        e(imt,k,2) = e(2,k,2)
        e(imt,k,3) = e(2,k,3)
        e(1,k,3) = e(imtm1,k,3)
#else
        e(imt,k,1) = c0
        e(imt,k,2) = c0
        e(imt,k,3) = c0
        e(  1,k,3) = c0
#endif
700   continue
c
#endif
      return
      end
