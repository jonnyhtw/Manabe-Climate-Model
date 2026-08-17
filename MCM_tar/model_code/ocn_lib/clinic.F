      subroutine clinic(j)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     clinic computes internal mode velocity at time level n+1 for jow j
c
c=======================================================================
c
#include "param.h"
#ifdef biharmonic
# include "cbihar.h"
#endif
#include "cdiag.h"
#include "chmix.h"
#ifdef nlhmix
# include "cnlmix.h"
#endif
#include "cregin.h"
#include "ctask.h"
#include "cvbc.h"
#include "cvmix.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "slabs.h"
#include "temp.h"
c
      dimension ubar(imt), vbar(imt), ueng(imt), veng(imt), twodt(km)
c
#ifndef keepterms
c
c-----------------------------------------------------------------------
c     define operators using statement functions ...
c     otherwise they will be dimensioned in "slabs.h"
c-----------------------------------------------------------------------
c
#include "fdifm.h"
#endif
c
c-----------------------------------------------------------------------
c       calculate external mode velocities
c-----------------------------------------------------------------------
c
        do 100 i=2,imtm1
          ubar(i) = -(p(i,j+1,1)-p(i,j,1))*dyur(j)*
     $                   min(hr(i-1,j),hr(i,j))
          vbar(i) =  (p(i+1,j+1,1)-p(i,j+1,1))*dxur(i)*
     $                   min(hr(i,j+1),hr(i,j))*cstr(j+1)
100     continue
c
      if (j .eq. jsrow) then
c
c-----------------------------------------------------------------------
c       on the first row of a task:
c
c       1) compute advective coefficients at the south face of the
c          "u,v" box for the next row (jsrow+1)
c
c       2) compute density for the next row (jsrow+1) for ocean points
c          on the "t" grid
c
c       3) return
c-----------------------------------------------------------------------
c
#ifdef nlhmix
c
c       initialize deformation fields for nonlinear horizontal mixing
c
        do 180 k=1,km
          do 170 i=1,imt
            dtens(i,k)  = c0
            dshers(i,k) = c0
            dpuret(i,k) = c0
            dpures(i,k) = c0
170       continue
180     continue
#endif
c
c       calculate "fvsu" this way to insure that answers do not change
c       regardless of how "ntasks" is set
c
        fx = dyu2r(j)*csur(j)*cst(j+1)
        do 200 k=1,km
          do 190 i=2,imtm1
            fvsu(i,k)  = ((vclin(i,k,jpt2) + vclin(i,k,jpt1))*p5
     $                   + vbar(i))*fx
190       continue
200     continue
        fxb = csu(j)*dyu(j)*csur(j+1)*dyur(j+1)
        do 260 k=1,km
          do 250 i=2,imtm1
            fvsu(i,k) = fvsu(i,k)*fxb
250       continue
260     continue
c
        call state (t(1,1,jp,nc,1), t(1,1,jp,nc,2), rhos)
c
        return
c
      else
c
c-----------------------------------------------------------------------
c       on all but the first row of a task:
c
c       1) calculate advective coefficients "fuw" at the west face of
c          the "u,v" box & "fvn" at the north face of the "u,v" box
c          by combining external & internal modes with a grid factor
c
c       2) initialize velocities
c-----------------------------------------------------------------------
c
        fx = dyu2r(j)*csur(j)*cst(j+1)
        do 300 k=1,km
          do 270 i=2,imtm1
            fuw(i,k) = ((uclin(i,k,jpt1) + uclin(i-1,k,jpt1))*p5
     $                  + ubar(i))*csur(j)
            fvn(i,k) = ((vclin(i,k,jpt2)+vclin(i,k,jpt1))*p5
     $                  + vbar(i))*fx
270       continue
#ifdef cyclic
          fuw(imt,k) = fuw(2,k)
#else
          fuw(imt,k) = c0
#endif
#ifdef skipland
          do 290 i=1,imt
            ua(i,k) = c0
            va(i,k) = c0
290       continue
#endif
300     continue
      endif
c
#ifdef symmetry
c
c-----------------------------------------------------------------------
c     set symmetry conditions on the last row
c-----------------------------------------------------------------------
c
      if (j .eq. jmtm1) then
        do 600 k=1,km
          do 590 i=1,imt
            fvn(i,k)     = -fvsu(i,k)
            u(i,k,jp,nm) = u(i,k,jm,nm)
            u(i,k,jp,nc) = u(i,k,jm,nc)
590       continue
600     continue
      endif
#endif
c
c-----------------------------------------------------------------------
c     compute density on ocean points of "t" grid row j+1
c-----------------------------------------------------------------------
c
      call state (t(1,1,jp,nc,1), t(1,1,jp,nc,2), rhon)
c
#ifdef nlhmix
c
c-----------------------------------------------------------------------
c     set up diffusion factors for this row
c-----------------------------------------------------------------------
c
      call nlmix(j)
      bbu(j) = p5*csur(j)
      ccu(j) = effky(j)*dyur(j)*csur(j)*csur(j)
#endif
c
c-----------------------------------------------------------------------
c     initialize integral arrays for external mode diagnostics
c-----------------------------------------------------------------------
c
      if (diagts .and. eots) then
        do 700 ll=1,8
          do 690 i=1,imt
            zuneng(i,ll) = c0
            zvneng(i,ll) = c0
690       continue
700     continue
      endif
#ifdef skipland
c
c-----------------------------------------------------------------------
c     solve for velocity on longitudinal strips of ocean points.
c     if skipland is not set them "istart" & "iend" are parameters
c-----------------------------------------------------------------------
c
      do 3500 lstrip=1,lseg
        istart = isu(j,lstrip)
        iend   = ieu(j,lstrip)
        if (istart .ne. 0) then
#endif
c
c-----------------------------------------------------------------------
c         accumulate kinetic energy on row "j" every "ntsi" timesteps
c-----------------------------------------------------------------------
c
          if (prntsi) then
            fx = p25*csu(j)*dyu(j)
#ifdef symmetry
c
c           weight symmetry row by one half
c
            if (j .eq. jmtm1) fx = fx*p5
#endif
            do 800 k=1,km
              do 780 i=istart,iend
                tempik(i,k,1) = (fx*(u(i,k,jc,nc)*u(i,k,jc,nc) +
     $                      v(i,k,jc,nc)*v(i,k,jc,nc)))*c2dzt(k)*dxu(i)
780           continue
              do 790 i=istart,iend
                ektot(k,j) = ektot(k,j) + tempik(i,k,1)
790           continue
800         continue
          endif
c
c-----------------------------------------------------------------------
c         compute vertical velocity at the bottoms of u,v boxes
c-----------------------------------------------------------------------
c
c         compute divergence of mass at each level
c
          do 900 k=1,kmm1
            do 890 i=istart,iend
              w(i,k) = c2dzt(k)*((fuw(i+1,k)-fuw (i,k))*dxu2r(i) +
     $                            fvn(i  ,k)-fvsu(i,k))
890         continue
900       continue
c
c         integrate downward from the surface to define vertical
c         velocity at the bottom of each level
c
          do 1000 k=2,kmm1
            do 990 i=istart,iend
              w(i,k) = w(i,k) + w(i,k-1)
990         continue
1000      continue
c
c         set vertical velocity at the surface to zero (rigid-lid)
c         set vertical velocity at level k=km to zero
c         note: vertical velocity may be non zero at physical bottom
c         (if kmu(i,j) < km & there is a bottom slope).
c
          do 1100 i=istart,iend
            w(i,0)  = c0
            w(i,km) = c0
1100      continue
c
c-----------------------------------------------------------------------
c         compute hydrostatic pressure gradient
c-----------------------------------------------------------------------
c
c         compute it at the first level
c
          fxa = grav*dzw(0)*csur(j)
          fxb = grav*dzw(0)*dyu2r(j)
          do 1200 i=istart-1,iend
            t1        = rhon(i+1,1) - rhos(i  ,1)
            t2        = rhon(i  ,1) - rhos(i+1,1)
            dpdx(i,1) = (t1-t2)*fxa*dxu2r(i)
            dpdy(i,1) = (t1+t2)*fxb
1200      continue
c
c         compute the change in pressure gradient between levels
c
          fxa = grav*csur(j)*p5
          fxb = grav*dyu4r(j)
          do 1300 k=2,km
            do 1280 i=istart-1,iend+1
              tempik(i,k,1) = rhon(i,k-1) + rhon(i,k)
              tempik(i,k,2) = rhos(i,k-1) + rhos(i,k)
1280        continue
            do 1290 i=istart-1,iend
              t1        = tempik(i+1,k,1) - tempik(i  ,k,2)
              t2        = tempik(i  ,k,1) - tempik(i+1,k,2)
              dpdx(i,k) = fxa*(t1-t2)*dzw(k-1)*dxu2r(i)
              dpdy(i,k) = fxb*(t1+t2)*dzw(k-1)
1290        continue
1300      continue
c
c         integrate downward from the first level
c
          do 1400 k=1,kmm1
            do 1390 i=istart-1,iend
              dpdx(i,k+1) = dpdx(i,k) + dpdx(i,k+1)
              dpdy(i,k+1) = dpdy(i,k) + dpdy(i,k+1)
1390        continue
1400      continue
c
#ifdef skipland
          if (lstrip .eq. 1) then
#endif
            do 1500 k=1,km
#ifdef cyclic
              dpdx(imt,k) = dpdx(2,k)
              dpdy(imt,k) = dpdy(2,k)
#else
              dpdx(imt,k) = c0
              dpdy(imt,k) = c0
#endif
1500        continue
#ifdef skipland
          endif
#endif
c
c-----------------------------------------------------------------------
c         calculate quantities for the computation of
c         vertical diffusion and advection of momentum
c-----------------------------------------------------------------------
c
          do 1600 k=1,kmm1
            do 1590 i=istart,iend
#ifdef tcvmix
              vmf(i,k,1) = vvc(i,k,jc,nm)*(u(i,k,jc,nm)
     $                                   - u(i,k+1,jc,nm))*dzwr(k)
              vmf(i,k,2) = vvc(i,k,jc,nm)*(v(i,k,jc,nm)
     $                                   - v(i,k+1,jc,nm))*dzwr(k)
#else
              vmf(i,k,1) = vvc(i,k)*(u(i,k,jc,nm) - u(i,k+1,jc,nm))
     $                     *dzwr(k)
              vmf(i,k,2) = vvc(i,k)*(v(i,k,jc,nm) - v(i,k+1,jc,nm))
     $                     *dzwr(k)
#endif
              fwb1(i,k) = w(i,k)*(u(i,k,jc,nc) + u(i,k+1,jc,nc))
              fwb2(i,k) = w(i,k)*(v(i,k,jc,nc) + v(i,k+1,jc,nc))
1590        continue
1600      continue
c
c         set k=0 elements to reflect wind stress
c         set level "kz" for a no bottom drag condition.
c         set k=0 elements of vertical advection array to the rigid lid
c         condition (w(0) = 0). also assume w(km)=0. If "kz" < "km"
c         and there is a bottom slope, w(kz) may be non zero.
c
          do 1700 i=istart,iend
            vmf(i,0,1) = smf(i,1)
            vmf(i,0,2) = smf(i,2)
            kz         = kmu(i,j)
            vmf(i,kz,1)= bmf(i,1)
            vmf(i,kz,2)= bmf(i,2)
            fwb1(i,0)  = c0
            fwb2(i,0)  = c0
            fwb1(i,km) = c0
            fwb2(i,km) = c0
1700      continue
c
c
c-----------------------------------------------------------------------
c         construct auxiliary arrays for to prevent unnecessary
c         computation.
c-----------------------------------------------------------------------
c
#ifdef biharmonic
c
c         set up the momentum pointers (indices)
c         m1 = u component
c         m2 = v component
c
          m1 = nt + 1
          m2 = nt + 2
#endif
          do 1800 k=1,km
            do 1790 i=istart,iend+1
              aux1(i,k) = fuw(i,k)*(u(i-1,k,jc,nc) + u(i,k,jc,nc))
              aux2(i,k) = fuw(i,k)*(v(i-1,k,jc,nc) + v(i,k,jc,nc))
#ifdef biharmonic
              aux3(i,k) = dxt4r(i)*(del2(i,k,m1,jpt2)
     $                            - del2(i-1,k,m1,jpt2))
              aux4(i,k) = dxt4r(i)*(del2(i,k,m2,jpt2)
     $                            - del2(i-1,k,m2,jpt2))
#endif
#ifdef consthmix
              aux3(i,k) = dxt4r(i)*(u(i,k,jc,nm) - u(i-1,k,jc,nm))
              aux4(i,k) = dxt4r(i)*(v(i,k,jc,nm) - v(i-1,k,jc,nm))
#endif
1790        continue
1800      continue
c
c-----------------------------------------------------------------------
c         set up for the coriolis term:
c         use tau timestep for explicit trtmnt ... or
c         use tau-1 timestep for implicit treatment with remainder of
c         term to be added later
c-----------------------------------------------------------------------
c
          if (acor .eq. c0) then
            lcor = nc
          else
            lcor = nm
          endif
c
c-----------------------------------------------------------------------
c         if momentum terms are defined as arrays, then calculate
c         them here otherwise statement functions will be used
c-----------------------------------------------------------------------
c
#ifdef keepterms
          do 1880 k=1,km
            do 1870 i=istart,iend
#include "fdifm.h"
1870        continue
1880      continue
#endif
c
c-----------------------------------------------------------------------
c         construct the forcing terms for momentum
c-----------------------------------------------------------------------
c
          do 1900 k=1,km
            do 1890 i=istart,iend
              ua(i,k) = (UUx(i,k) + VUy(i,k) + WUz(i,k)
     $                 + Uxx(i,k) + Uyy(i,k) + Umet(i,k)
#ifndef nohilats
     $                 + UVtan(i,k)
#endif
#ifdef implicitvmix
     $                 + Uzz(i,k)*(c1-aidif)
#else
     $                 + Uzz(i,k)
#endif
     $                 + fV(i,k) - dpdx(i,k) + sourcu(i,k))*gm(i,k,jc)
c
              va(i,k) = (UVx(i,k) + VVy(i,k) + WVz(i,k)
     $                 + Vxx(i,k) + Vyy(i,k) + Vmet(i,k)
#ifndef nohilats
     $                 - UUtan(i,k)
#endif
#ifdef implicitvmix
     $                 + Vzz(i,k)*(c1-aidif)
#else
     $                 + Vzz(i,k)
#endif
     $                 - fU(i,k) - dpdy(i,k) + sourcv(i,k))*gm(i,k,jc)
1890        continue
1900      continue
#ifdef implicitvmix
c
c-----------------------------------------------------------------------
c         solve vertical diffusion implicitly
c-----------------------------------------------------------------------
c
c        update the velocity without implicit vertical diffusion
c
          do 1930 k=1,km
            twodt(k) = c2dtuv
            do 1920 i=istart,iend
              ua(i,k) = u(i,k,jc,nm) + c2dtuv*ua(i,k)
              va(i,k) = v(i,k,jc,nm) + c2dtuv*va(i,k)
1920        continue
1930      continue
c
c         store terms to compute implicit vertical diffusion on
c         diagnostic time steps
c
          if (diagts .and. eots) then
            do 1050 k=1,km
              do 1040 i=istart,iend
                zzi(i,k,1) = ua(i,k)
                zzi(i,k,2) = va(i,k)
1040          continue
1050        continue
          endif
c
c         add in the implicit vertical diffusion
c
#ifdef tcvmix
          call invtri (ua, smf(1,1), bmf(1,1), vvca, twodt, kmu
     $,                gm(1,1,jc), istart, iend, j)
          call invtri (va, smf(1,2), bmf(1,2), vvca, twodt, kmu
     $,                gm(1,1,jc), istart, iend, j)
#else
          call invtri (ua, smf(1,1), bmf(1,1), vvc, twodt, kmu
     $,                gm(1,1,jc), istart, iend, j)
          call invtri (va, smf(1,2), bmf(1,2), vvc, twodt, kmu
     $,                gm(1,1,jc), istart, iend, j)
#endif
c
c         compute residual implicit vertical diffusion
c
          r2dtuv = c1/c2dtuv
          if (diagts .and. eots) then
            do 1070 k=1,km
              do 1060 i=istart,iend
                zzi(i,k,1) = r2dtuv*(ua(i,k) - zzi(i,k,1))
                zzi(i,k,2) = r2dtuv*(va(i,k) - zzi(i,k,2))
1060          continue
1070        continue
          endif
c
c         convert back to time change of velocity
c
          do 1990 k=1,km
            do 1980 i=istart,iend
              ua(i,k) = r2dtuv*(ua(i,k) - u(i,k,jc,nm))
              va(i,k) = r2dtuv*(va(i,k) - v(i,k,jc,nm))
1980        continue
1990      continue
#endif
c
c-----------------------------------------------------------------------
c         form time change of vertically averaged forcing
c-----------------------------------------------------------------------
c
c         integrate time change vertically
c
          do 2000 i=istart,iend
            zu(i,j) = c0
            zv(i,j) = c0
2000      continue
          do 2100 k=1,km
            fx = c2dtsf*dzt(k)
            do 2090 i=istart,iend
              zu(i,j) = zu(i,j) + ua(i,k)*fx
              zv(i,j) = zv(i,j) + va(i,k)*fx
2090        continue
2100      continue
c
#ifdef worldocn
# ifdef seasonal
c   For external mode forcing, subtract term derived from monthly
c    wind stress and add in term derived from annual wind stress
          if (dtts/dtsf .gt. 4.0) then
            do 2120 i=istart,iend
              zu(i,j) = zu(i,j) + (smfavg(i,1) - smf(i,1))*c2dtsf
              zv(i,j) = zv(i,j) + (smfavg(i,2) - smf(i,2))*c2dtsf
2120        continue
          endif
# endif
#endif
#ifdef symmetry
c
c         set symmetry row to zero
c
          if (j .eq. jmtm1) then
            do 2200 i=istart,iend
              zv(i,j) = c0
2200        continue
          endif
#endif
c
c         form vertical average
c
          do 2300 i=istart,iend
            zu(i,j) = zu(i,j)*hr(i,j)
            zv(i,j) = zv(i,j)*hr(i,j)
2300      continue
c
c-----------------------------------------------------------------------
c         do analysis of internal mode forcing on diagnostic timesteps
c         also, form vert ave. for use later in ext. mode analysis
c         along with term balances
c-----------------------------------------------------------------------
c
          if (diagts .and. eots) then
c
            fx = csu(j)*dyu(j)
#ifdef symmetry
c
c           weight symmetry row by one half
c
            if (j .eq. jmtm1) fx = fx*p5
#endif
            do 2400 k=1,km
c
c             compute ke change due to pressure term globally
c             also integrate pressure term over the regional volume
c
              n0 = nhreg*(mskvr(k)-1)
              do 2320 i=istart,iend
                ueng(i)   = -gm(i,k,jc)*dpdx(i,k)
                veng(i)   = -gm(i,k,jc)*dpdy(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,6) = engint(k,6) + (uclin(i,k,jpt1)*ueng(i) +
     $                        vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,6) = zuneng(i,6) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,6) = zvneng(i,6) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,2,1,n) = termbm(k,2,1,n) + ueng(i)*boxvol
                  termbm(k,2,2,n) = termbm(k,2,2,n) + veng(i)*boxvol
                endif
2320          continue
c
c             compute ke change due to zonal advection of momentum
c             globally. also integrate zonal advection term over the
c             regional volumes
c
              do 2330 i=istart,iend
                ueng(i)     = gm(i,k,jc)*UUx(i,k)
                veng(i)     = gm(i,k,jc)*UVx(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,2) = engint(k,2) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,2) = zuneng(i,2) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,2) = zvneng(i,2) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,3,1,n) = termbm(k,3,1,n) + ueng(i)*boxvol
                  termbm(k,3,2,n) = termbm(k,3,2,n) + veng(i)*boxvol
#ifndef nohilats
                  termbm(k,13,1,n) = termbm(k,13,1,n)
     $                             + UVtan(i,k)*boxvol
                  termbm(k,13,2,n) = termbm(k,13,2,n)
     $                             - UUtan(i,k)*boxvol
#endif
                endif
2330          continue
c
c             compute ke change due to meridional advection of momentum
c             globally. also integrate meridional advection term over
c             the regional volumes
c
              do 2340 i=istart,iend
                ueng(i)     = gm(i,k,jc)*VUy(i,k)
                veng(i)     = gm(i,k,jc)*VVy(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,2) = engint(k,2) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,2) = zuneng(i,2) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,2) = zvneng(i,2) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,4,1,n) = termbm(k,4,1,n) + ueng(i)*boxvol
                  termbm(k,4,2,n) = termbm(k,4,2,n) + veng(i)*boxvol
                endif
2340          continue
c
c             compute ke change due to vertical advection of momentum
c             globally. also integrate vertical advection term over the
c             regional volumes
c
              do 2350 i=istart,iend
                ueng(i)     = gm(i,k,jc)*WUz(i,k)
                veng(i)     = gm(i,k,jc)*WVz(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,3) = engint(k,3) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,3) = zuneng(i,3) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,3) = zvneng(i,3) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,5,1,n) = termbm(k,5,1,n) + ueng(i)*boxvol
                  termbm(k,5,2,n) = termbm(k,5,2,n) + veng(i)*boxvol
                endif
2350          continue
c
c             compute ke change due to zonal diffusion of momentum
c             globally. also integrate zonal diffusion term over the
c             regional volumes
c
              do 2360 i=istart,iend
                ueng(i)     = gm(i,k,jc)*Uxx(i,k)
                veng(i)     = gm(i,k,jc)*Vxx(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,4) = engint(k,4) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,4) = zuneng(i,4) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,4) = zvneng(i,4) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,6,1,n) = termbm(k,6,1,n) + ueng(i)*boxvol
                  termbm(k,6,2,n) = termbm(k,6,2,n) + veng(i)*boxvol
                endif
2360          continue
c
c             compute ke change due to meridional diffusion of momentum
c             lobally. also integrate meridional diffusion term over
c             the regional volumes
c
              do 2370 i=istart,iend
                ueng(i)     = gm(i,k,jc)*Uyy(i,k)
                veng(i)     = gm(i,k,jc)*Vyy(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,4) = engint(k,4) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,4) = zuneng(i,4) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,4) = zvneng(i,4) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,7,1,n) = termbm(k,7,1,n) + ueng(i)*boxvol
                  termbm(k,7,2,n) = termbm(k,7,2,n) + veng(i)*boxvol
                endif
2370          continue
c
c             compute ke change due to vertical diffusion of momentum
c             globally. also integrate vertical diffusion term over the
c             regional volumes
c
              do 2380 i=istart,iend
#ifdef implicitvmix
                ueng(i) = gm(i,k,jc)*(Uzz(i,k)*(c1-aidif) + zzi(i,k,1))
                veng(i) = gm(i,k,jc)*(Vzz(i,k)*(c1-aidif) + zzi(i,k,2))
#else
                ueng(i) = gm(i,k,jc)*Uzz(i,k)
                veng(i) = gm(i,k,jc)*Vzz(i,k)
#endif
                boxvol  = fx*dxu(i)*dzt(k)
                engint(k,5) = engint(k,5) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,5) = zuneng(i,5) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,5) = zvneng(i,5) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,8,1,n) = termbm(k,8,1,n) + ueng(i)*boxvol
                  termbm(k,8,2,n) = termbm(k,8,2,n) + veng(i)*boxvol
                endif
2380          continue
c
c             compute ke change due to metric term globally
c             also integrate metric term, coriolis term, source term
c             & velocity over the regional volumes
c
              do 2390 i=istart,iend
                ueng(i)     = gm(i,k,jc)*Umet(i,k)
                veng(i)     = gm(i,k,jc)*Vmet(i,k)
                boxvol      = fx*dxu(i)*dzt(k)
                engint(k,4) = engint(k,4) + (uclin(i,k,jpt1)*ueng(i)
     $                        + vclin(i,k,jpt1)*veng(i))*boxvol
                zuneng(i,4) = zuneng(i,4) + ueng(i)*dzt(k)*hr(i,j)
                zvneng(i,4) = zvneng(i,4) + veng(i)*dzt(k)*hr(i,j)
                n           = n0 + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  termbm(k,9,1,n) = termbm(k,9,1,n) + ueng(i)*boxvol
                  termbm(k,9,2,n) = termbm(k,9,2,n) + veng(i)*boxvol
                  termbm(k,10,1,n) = termbm(k,10,1,n) + fV(i,k)*boxvol
                  termbm(k,10,2,n) = termbm(k,10,2,n) - fU(i,k)*boxvol
                  termbm(k,11,1,n) = termbm(k,11,1,n) + sourcu(i,k)
     $                                                *boxvol*gm(i,k,jc)
                  termbm(k,11,2,n) = termbm(k,11,2,n) - sourcv(i,k)
     $                                                *boxvol*gm(i,k,jc)
                  avwv(1,n)        = avwv(1,n) + u(i,k,jc,nc)*boxvol
                  avwv(2,n)        = avwv(2,n) + v(i,k,jc,nc)*boxvol
                  avwv(3,n)        = avwv(3,n) + p5*(w(i,k)+w(i,k-1))
     $                                *boxvol*gm(i,k,jc)
                endif
2390          continue
2400        continue
c
c           compute ke change due to wind stress globally
c
            do 2500 i=istart,iend
              ueng(i)     = gm(i,1,jc)*smf(i,1)
              veng(i)     = gm(i,1,jc)*smf(i,2)
              engint(1,7) = engint(1,7) + (uclin(i,1,jpt1)*ueng(i)
     $                      + vclin(i,1,jpt1)*veng(i))*fx*dxu(i)
              zuneng(i,7) = zuneng(i,7) + ueng(i)*hr(i,j)
              zvneng(i,7) = zvneng(i,7) + veng(i)*hr(i,j)
              dtaux       = fx*dxu(i)*ueng(i)
              dtauy       = fx*dxu(i)*veng(i)
              smflx(1,0) = smflx(1,0) + dtaux
              smflx(2,0) = smflx(2,0) + dtauy
              n         = mskhr(i,j)
              if (n .gt. 0) then
                smflx(1,n) = smflx(1,n) + dtaux
                smflx(2,n) = smflx(2,n) + dtauy
              endif
2500        continue
c
c           compute ke change due to bottom drag globally
c
            do 2600 i=istart,iend
              kz = kmu(i,j)
              if (kz .ne. 0) then
                ueng(i)   = -gm(i,1,jc)*bmf(i,1)
                veng(i)   = -gm(i,1,jc)*bmf(i,2)
                engint(kz,8) = engint(kz,8)+(uclin(i,kz,jpt1)*ueng(i)
     $                     + vclin(i,kz,jpt1)*veng(i))*fx*dxu(i)
                zuneng(i,8)  = zuneng(i,8) + ueng(i)*hr(i,j)
                zvneng(i,8)  = zvneng(i,8) + veng(i)*hr(i,j)
              endif
2600        continue
          endif
c
c-----------------------------------------------------------------------
c         compute new velocities (with incorrect vertical means)
c         also, add in remainder of coriolis term if treated implicitly
c-----------------------------------------------------------------------
c
          if (acor .eq. c0) then
            do 2800 k=1,km
              do 2790 i=istart,iend
                ua(i,k) = u(i,k,jc,nm) + c2dtuv*ua(i,k)
                va(i,k) = v(i,k,jc,nm) + c2dtuv*va(i,k)
2790          continue
2800        continue
          else
            fx    = c2dtuv*acor*c2*omega*sine(j)
            detmr = c1/(c1 + fx*fx)
            do 2900 k=1,km
              do 2890 i=istart,iend
                t1 = (ua(i,k)+fx*va(i,k))*detmr
                t2 = (va(i,k)-fx*ua(i,k))*detmr
                ua(i,k) = u(i,k,jc,nm) + c2dtuv*t1
                va(i,k) = v(i,k,jc,nm) + c2dtuv*t2
2890          continue
2900        continue
          endif
c
c-----------------------------------------------------------------------
c         subtract incorrect vertical means to get internal modes.
c         these means are due to ignoring horizontal gradients of the
c         surface pressure
c-----------------------------------------------------------------------
c
          do 3000 i=istart,iend
            ubar(i) = c0
            vbar(i) = c0
3000      continue
          do 3100 k=1,km
            do 3090 i=istart,iend
              ubar(i) = ubar(i) + ua(i,k)*dzt(k)
              vbar(i) = vbar(i) + va(i,k)*dzt(k)
3090        continue
3100      continue
          do 3200 i=istart,iend
            ubar(i) = ubar(i)*hr(i,j)
            vbar(i) = vbar(i)*hr(i,j)
3200      continue
          do 3300 k=1,km
            do 3290 i=istart,iend
              ua(i,k) = ua(i,k) - gm(i,k,jc)*ubar(i)
              va(i,k) = va(i,k) - gm(i,k,jc)*vbar(i)
3290        continue
3300      continue
c
c-----------------------------------------------------------------------
c         compute total change in k.e. for the internal mode on diag-
c         nostic timesteps. also construct d/dt of momentum and add the
c         implicit part of the coriolis term into the momentum term
c         balance. the external mode parts of both will be added after
c         the calculation for the external mode
c-----------------------------------------------------------------------
c
          if (diagts .and. eots) then
            do 3400 k=1,km
              fx = csu(j)*dyu(j)*dzt(k)/c2dtuv
#ifdef symmetry
              if (j .eq. jmtm1) fx = fx*p5
#endif
              do 3380 i=istart,iend
                tempik(i,k,1) =
     $              (uclin(i,k,jpt1)*(ua(i,k) - u(i,k,jc,nm))
     $             + vclin(i,k,jpt1)*(va(i,k) - v(i,k,jc,nm)))*fx*dxu(i)
3380          continue
              do 3390 i=istart,iend
                engint(k,1) = engint(k,1) + tempik(i,k,1)
3390          continue
3400        continue
c
            fx   = csu(j)*dyu(j)/c2dtuv
            do 3430 k=1,km
              do 3420 i=istart,iend
                n = nhreg*(mskvr(k)-1) + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  boxfac      = gm(i,k,jc)*fx*dxu(i)*dzt(k)
                  termbm(k,1,1,n) = termbm(k,1,1,n) +
     $                            (ua(i,k) - u(i,k,jc,nm))*boxfac
                  termbm(k,1,2,n) = termbm(k,1,2,n) +
     $                            (va(i,k) - v(i,k,jc,nm))*boxfac
                endif
3420          continue
3430        continue
c
            if (acor .ne. c0) then
              fx   = acor*csu(j)*dyu(j)*fcor(j)
              do 3450 k=1,km
                do 3440 i=istart,iend
                  n = nhreg*(mskvr(k)-1) + mskhr(i,j)
                  if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                    boxfac      = gm(i,k,jc)*fx*dxu(i)*dzt(k)
                    termbm(k,10,1,n) = termbm(k,10,1,n) +
     $                               (va(i,k) - v(i,k,jc,nm))*boxfac
                    termbm(k,10,2,n) = termbm(k,10,2,n) -
     $                               (ua(i,k) - u(i,k,jc,nm))*boxfac
                  endif
3440            continue
3450          continue
            endif
          endif
#ifdef skipland
        endif
3500  continue
#endif
c
c-----------------------------------------------------------------------
c     do analysis of external mode forcing on diagnostic timesteps
c-----------------------------------------------------------------------
c
      if (diagts .and. eots) then
#ifdef cyclic
        do 3600 ll=2,8
          zuneng(1,ll) = zuneng(imum1,ll)
          zvneng(1,ll) = zvneng(imum1,ll)
3600    continue
#endif
c
        do 3700 ll=2,8
#ifdef symmetry
          if (j .eq. jmtm1) then
            do 3670 i=1,imt
              zvneng(i,ll) = c0
3670        continue
          endif
#endif
          do 3680 i=2,imtm1
            tempik(i,1,1) =
     $            -p(i,j,1)*(((zvneng(i,ll)-zvneng(i-1,ll))*dyu(j)
     $            +(zvseng(i,ll)-zvseng(i-1,ll))*dyu(j-1))
     $            *dxt2r(i)*dytr(j)
     $           -((zuneng(i,ll)*dxu(i)+zuneng(i-1,ll)*dxu(i-1))*csu(j)
     $           -(zuseng(i,ll)*dxu(i)+zuseng(i-1,ll)*dxu(i-1))*csu(j-1)
     $            )*dyt2r(j)*dxtr(i))*dxt(i)*dyt(j)
3680      continue
          do 3690 i=2,imtm1
            engext(ll) = engext(ll) + tempik(i,1,1)
3690      continue
3700    continue
      endif
c
#if defined fourfil || defined firfil
c-----------------------------------------------------------------------
c     fourier filter u and v at high latitudes
c-----------------------------------------------------------------------
c
#ifdef cyclic
      do 3780 k=1,km
        ua(1,k)   = ua(imum1,k)
        va(1,k)   = va(imum1,k)
        ua(imu,k) = ua(2,k)
        va(imu,k) = va(2,k)
3780  continue
#endif
c
      call filuv (j)
c
#endif
c-----------------------------------------------------------------------
c     transfer quantities computed to the north of the present row
c     to be defined at the south in the computation of the next row
c-----------------------------------------------------------------------
c
      if (j .ne. jerow) then
        fx = csu(j)*dyu(j)*csur(j+1)*dyur(j+1)
        do 3800 k=1,km
          do 3790 i=1,imt
            fvsu(i,k) = fvn(i,k)*fx
3790      continue
3800    continue
      endif
c
      if (diagts .and. eots) then
        do 3900 ll=2,8
          do 3890 i=1,imt
            zuseng(i,ll) = zuneng(i,ll)
            zvseng(i,ll) = zvneng(i,ll)
3890      continue
3900    continue
      endif
c
#ifdef cyclic
c
c     set cyclic boundary conditions
c
      zu(imu,j) = zu(2,j)
      zv(imu,j) = zv(2,j)
      zu(1,j)   = zu(imum1,j)
      zv(1,j)   = zv(imum1,j)
c
      do 4000 k=1,km
        ua(1,k)   = ua(imum1,k)
        va(1,k)   = va(imum1,k)
        ua(imu,k) = ua(2,k)
        va(imu,k) = va(2,k)
4000  continue
c
#endif
#ifdef symmetry
c
c     set meridional component of internal mode to zero on symmetry row
c
      if (j .eq. jmtm1) then
        do 4100 k=1,km
          do 4090 i=1,imt
            va(i,k) = c0
4090      continue
4100    continue
      endif
c
#endif
c
      return
      end
