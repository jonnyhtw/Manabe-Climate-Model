      subroutine diag2
      implicit real*8 (a-h,o-z)
cfpp$ noconcur r
c
c=======================================================================
c
c     calculate & print diagnostics outside parallel region
c
c=======================================================================
c
#include "param.h"
#ifdef testcfl
# include "ccfl.h"
# ifdef cflj
# include "ccflj.h"
# endif
#endif
#include "cdiag.h"
#include "chmix.h"
#include "coord.h"
#include "cprnts.h"
#include "cregin.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "cvmix.h"
#include "docnam.h"
#include "emode.h"
#include "grdvar.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
c
#if defined congrad5pt || defined congrad9pt
c
c     spx   = the x change (longitudinal) in surface pressure
c             spx(i,j) is defined midway between xu(i) & xu(i+1)
c     spy   = the y change (latitudinal) in surface pressure
c             spy(i,j) is defined midway between yu(j) & yu(j+1)
c     sperr = the closed line integral around a "t" grid box
c             it is a measure of the error in reconstructing the
c             surface pressure.
c
c             note: the actual surface pressure is not reconstructed
c                   here. it can be constructed by assuming an
c                   arbitrary value at some "u" point on the grid and
c                   and adding "spx" & "spy" along paths to build
c                   surface pressure values relative to the arbitrary
c                   starting value.
c                   under the congrad5pt option, the error is due to
c                   "crit" & truncation due to the 5 pt operator. under
c                   the congrad9pt, the error is due only to "crit"
c                   and so may be reduced by reducing "crit" & setting
c                   "mxscan" so that convergence is reached.
c
c
      dimension scrx(imt), scry(imt)
      dimension spx(imt,jmt), spy(imt,jmt), sperr(imt,jmt)
      equivalence (spx,res), (spy,fdir), (sperr,gdir)
#endif
c
c-----------------------------------------------------------------------
c     complete and print the integrals, transports & stream function
c-----------------------------------------------------------------------
c
c     add in the external mode part of d/dt for the momentum balance
c     also add in the external mode part of the implicit coriolis term
c     and compute the surface pressure gradients
c
#if defined congrad5pt || defined congrad9pt
c
c       initialize surface pressure gradient fields
c
        do 100 j=1,jmt
          do 90 i=1,imt
            spx(i,j) = c0
            spy(i,j) = c0
90        continue
100     continue
#endif
c
#ifdef zonalavg
# if defined cpldleg3 || defined cpldleg2
        if (prdgts) then
# endif
        call zonal
# if defined cpldleg3 || defined cpldleg2
        endif
# endif
#endif
c
        fxa = c1
c
c       on mixing time steps "ptd" has been multiplied by a factor of 2
c       and the time step has to be adjusted also.
c
        if (mod(itt,nmix) .eq. 1) then
          fxa = p5
        endif
        r2dtuv = c1/c2dtuv
        do 200 j=1,jmt-1
          fddt  = csu(j)*dyu(j)/c2dtuv
          fspr  = csu(j)*dyu(j)
          atosp = acor*c2*omega*sine(j)
          f1    = atosp*csu(j)*dyu(j)
          f2    = atosp*c2dtuv
          f3    = c2dtuv/c2dtsf
          do 190 i=2,imtm1
            kz = kmu(i,j)
            if (kz .ne. 0) then
              diag1 = p(i+1,j+1,1)-p(i  ,j,1)
              diag0 = p(i  ,j+1,1)-p(i+1,j,1)
              uext  = -(diag1+diag0)*dyu2r(j)*hr(i,j)
              vext  =  (diag1-diag0)*dxu2r(i)*hr(i,j)*csur(j)
              diag3 = fxa*(ptd(i+1,j+1)-ptd(i  ,j))
              diag4 = fxa*(ptd(i  ,j+1)-ptd(i+1,j))
              dubdt = (diag3+diag4)*dyu2r(j)*hr(i,j)
              dvbdt = (diag3-diag4)*dxu2r(i)*hr(i,j)*csur(j)
              surfpx   = r2dtuv*( dubdt + f3*zu(i,j) + f2*dvbdt)
              surfpy   = r2dtuv*(-dvbdt + f3*zv(i,j) + f2*dubdt)
#if defined congrad5pt || defined congrad9pt
              spx(i,j)   = surfpx
              spy(i,j)   = surfpy
#endif
              do 180 k=1,kz
                n = nhreg*(mskvr(k)-1) + mskhr(i,j)
                if (n .gt. 0 .and. mskhr(i,j) .gt. 0) then
                  boxfac = fddt*dxu(i)*dzt(k)
                  boxspr = fspr*dxu(i)*dzt(k)
                  boxacr = f1*dxu(i)*dzt(k)
                  termbm(k,1,1,n)  = termbm(k,1,1,n)  + uext*boxfac
                  termbm(k,1,2,n)  = termbm(k,1,2,n)  + vext*boxfac
                  termbm(k,10,1,n) = termbm(k,10,1,n) + vext*boxacr
                  termbm(k,10,2,n) = termbm(k,10,2,n) - uext*boxacr
                  termbm(k,12,1,n) = termbm(k,12,1,n) - surfpx*boxspr
                  termbm(k,12,2,n) = termbm(k,12,2,n) - surfpy*boxspr
                endif
180           continue
            endif
190       continue
200     continue
#if defined congrad5pt || defined congrad9pt
c
c       construct spatial changes in surface pressure
c
        do 300 j=1,jmt-1
          do 280 i=2,imtm1
            scrx(i)  = p5*(spx(i,j) + spx(i+1,j))*csu(j)*dxt(i+1)
            scry(i)  = p5*(spy(i,j) + spy(i,j+1))*dyt(j+1)
280       continue
#ifdef cyclic
          scrx(1)   = scrx(imtm1)
          scrx(imt) = scrx(2)
          scry(1)   = scry(imtm1)
          scry(imt) = scry(2)
#else
          scrx(1)   = c0
          scrx(imt) = c0
          scry(1)   = c0
          scry(imt) = c0
#endif
          do 290 i=1,imt
            spx(i,j)   = scrx(i)
            sperr(i,j) = scry(i)
290       continue
300     continue
        do 320 j=1,jmt-1
          do 310 i=1,imt
            spy(i,j) = sperr(i,j)
310       continue
320     continue
        do 330 i=1,imt
            spy(i,jmt) = spy(i,jmt-1)
            spx(i,jmt) = spx(i,jmt-1)
330     continue
c
c       do closed line integrals around boxes. since the closed line
c       integral of a gradient should be zero, this is a measure
c       of the error involved in reconstructing the surface pressure
c       from its gradients.
c
        do 350 j=1,jmt-1
          do 340 i=2,imtm1
            sperr(i,j) = spy(i+1,j) - spy(i,j)
     $                - (spx(i,j+1) - spx(i,j))
340       continue
350     continue
c
c-----------------------------------------------------------------------
c       print a sample of the surface pressure components
c-----------------------------------------------------------------------
c
# if defined cpldleg3 || defined cpldleg2
        if (prdgts) then
# endif
c       show only a subset to keep print file from getting too large
c
        is = indp (140.0, xt, imt)
        ie = indp (160.0, xt, imt)
        js = indp (-85.0, yt, jmt)
        je = indp (-30.0, yt, jmt)
        write (stdout,'(/a50)')
     $     ' Sample longitudinal change in surface pressure '
        scl = 1000.0
        write (stdout,8000) ' spx '
     $,    itt, xt(is), xt(ie), yt(js), yt(je), scl
        call matrix (spx, imt, is, ie, -js, -je, scl)
        write (stdout,'(a50)')
     $     ' Sample latitudinal change in surface pressure '
        scl = 1000.0
        write (stdout,8000) ' spy'
     $,    itt, xt(is), xt(ie), yt(js), yt(je), scl
        call matrix (spy, imt, is, ie, -js, -je, scl)
        write (stdout,'(a60)')
     $     ' Sample closed line integral of surface pressure gradient'
ccccc   scl = 0.01
        scl = 0.0001
        write (stdout,8000) ' surf press error'
     $,    itt, xt(is), xt(ie), yt(js), yt(je), scl
        call matrix (sperr, imt, is, ie, -js, -je, scl)
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
#ifdef testcfl
c
c       monitor velocities by showing which come closest to the
c       CFL condition
c
        write (stdout,'(/60x,a12/)') ' CFL summary'

        write (stdout,'(a21,g10.3,a5,f7.2,a10,/a25,i4,a1,i4,a1,i3,a1
     $,                 3x,a18,f7.2,a1,f7.2,a1,f7.0,a1)')
     $   ' maximum U velocity (',cflum,') is ',cflup,' % the CFL'
     $,  ' at location: (i,j,k) = (',icflu,',',jcflu,',',kcflu,')'
     $,  ' (lon,lat,dpt) = (',xu(icflu),',',yu(jcflu),',',zt(kcflu),')'
   
        write (stdout,'(a21,g10.3,a5,f7.2,a10,/a25,i4,a1,i4,a1,i3,a1
     $,                 3x,a18,f7.2,a1,f7.2,a1,f7.0,a1)')
     $   ' maximum V velocity (',cflvm,') is ',cflvp,' % the CFL'
     $,  ' at location: (i,j,k) = (',icflv,',',jcflv,',',kcflv,')'
     $,  ' (lon,lat,dpt) = (',xu(icflv),',',yu(jcflv),',',zt(kcflv),')'
        write (stdout,'(a21,g10.3,a5,f7.2,a10,/a25,i4,a1,i4,a1,i3,a1
     $,                 3x,a18,f7.2,a1,f7.2,a1,f7.0,a1)')
     $   ' maximum W velocity (',cflwm,') is ',cflwp,' % the CFL'
     $,  ' at location: (i,j,k) = (',icflw,',',jcflw,',',kcflw,')'
     $,  ' (lon,lat,dpt) = (',xu(icflw),',',yu(jcflw),',',zw(kcflw),')'
        write (stdout,'(/58x,a16/)') ' end CFL summary'
#endif
c
c
c       integrate previously computed integrals vertically
c
        do 1300 k=km,1,-1
          buoy(0) = buoy(0) + buoy(k)
1300    continue
        do 1400 ll=1,8
          do 1390 k=km,1,-1
            engint(0,ll) = engint(0,ll) + engint(k,ll)
1390      continue
1400    continue
        do 1500 n=0,numreg
          if (n .gt. 0) then
            nv = (n-1)/nhreg + 1
            ks = llvreg(nv,1)
            ke = llvreg(nv,2)
            do 1420 ll=1,13
              do 1410 k=ke,ks,-1
                termbm(0,ll,1,n) = termbm(0,ll,1,n) + termbm(k,ll,1,n)
                termbm(0,ll,2,n) = termbm(0,ll,2,n) + termbm(k,ll,2,n)
1410          continue
1420        continue
          else
            ks = 1
            ke = km
          endif
          do 1480 m=1,nt
            do 1440 k=ke,ks,-1
              tvar(0,m,n) = tvar(0,m,n) + tvar(k,m,n)
1440        continue
#ifdef kdtermb
            do 1470 ll=1,9
#else
            do 1470 ll=1,8
#endif
              do 1460 k=ke,ks,-1
                termbt(0,ll,m,n) = termbt(0,ll,m,n) + termbt(k,ll,m,n)
1460          continue
1470        continue
1480      continue
1500    continue
c
c       normalize integrals by appropriate volume
c
        buoy(0) = buoy(0)/volume
        do 1600 ll=1,8
          engint(0,ll) = engint(0,ll)/volume
          engext(ll)   = engext(ll)/volume
1600    continue
        do 1700 n=0,numreg
          do 1690 m=1,nt
            tvar(0,m,n) = tvar(0,m,n)/volt(n)
            atwv(m,n)   = atwv(m,n)/volt(n)
            if (n .le. nhreg)  then
              stflx(m,n) = stflx(m,n)/volt(n)
              asst(m,n)  = asst(m,n)/areat(n)
            endif
#ifdef kdtermb
            do 1680 ll=1,9
#else
            do 1680 ll=1,8
#endif
              termbt(0,ll,m,n) = termbt(0,ll,m,n)/volt(n)
1680        continue
1690      continue
1700    continue
        do 1800 n=0,numreg
          if (n .le. nhreg) then
            smflx(1,n) = smflx(1,n)/volu(n)
            smflx(2,n) = smflx(2,n)/volu(n)
          endif
          if (n .gt. 0) then
            avwv(1,n) = avwv(1,n)/volu(n)
            avwv(2,n) = avwv(2,n)/volu(n)
            avwv(3,n) = avwv(3,n)/volu(n)
            do 1790 ll=1,13
              termbm(0,ll,1,n) = termbm(0,ll,1,n)/volu(n)
              termbm(0,ll,2,n) = termbm(0,ll,2,n)/volu(n)
1790        continue
          endif
1800    continue
c
c       2nd, compute residual terms
c
        plicin = engint(0,1) - engint(0,2) - engint(0,3) - engint(0,4)
     $         - engint(0,5) - engint(0,6)
        plicex = engext(1) - engext(2) - engext(3) - engext(4)
     $         - engext(5) - engext(6)
c
c       3rd, print the integrals
c
        write (stdout,9100)
     $                   'Globally averaged (entire ocean) Energetics'
     $,                     ' for ts =', itt, stamp, volu(0), areau(0)
        write (stdout,9101) ' time rate of change ',engint(0,1)
     $,                      engint(0,1), engext(1), engext(1)
        write (stdout,9101) ' horizontal advection', engint(0,2)
     $,                      engint(0,2), engext(2), engext(2)
        write (stdout,9101) ' vertical advection  ',engint(0,3)
     $,                      engint(0,3), engext(3), engext(3)
        write (stdout,9101) ' horizontal friction ',engint(0,4)
     $,                      engint(0,4), engext(4), engext(4)
        write (stdout,9101) ' vertical friction   ',engint(0,5)
     $,                      engint(0,5), engext(5), engext(5)
        write (stdout,9101) ' pressure forces     ',engint(0,6)
     $,                      engint(0,6), engext(6), engext(6)
        write (stdout,9101) ' ficticious sources  ',plicin
     $,                       plicin, plicex, plicex
        write (stdout,9101) ' work by wind        ',engint(0,7)
     $,                      engint(0,7), engext(7), engext(7)
        write (stdout,9101) ' bottom drag         ',engint(0,8)
     $,                      engint(0,8), engext(8), engext(8)
        buoerr = buoy(0) - engint(0,6) - engext(6)
        enleak = engint(0,2) + engint(0,3) + engext(2) + engext(3)
        write (stdout,9110) buoy(0), buoy(0), buoerr, buoerr, enleak
     $,                     enleak
c
        if (iotrmb .ne. stdout) then
          write (stdout, '(/a35,i4,a9,i10,a33/)')
     $    ' ==> Term balances written to unit ', iotrmb,' for ts ='
     $,   itt, stamp
        endif
        n = 0
        taux = smflx(1,n)*volu(n)/areau(n)
        tauy = smflx(2,n)*volu(n)/areau(n)
        write (iotrmb,10106)
     $                   'Globally averaged (entire ocean) quantities'
     $,                     ' for ts =', itt, stamp, volt(n), areat(n)
     $,                     volu(0), areau(n)
        write (iotrmb,10104)
        write (iotrmb,10098) n, ' smf(1)   = ', taux,' dynes/cm**2   '
        write (iotrmb,10098) n, ' smf(2)   = ', tauy,' dynes/cm**2   '
        do 1810 m=1,nt
          write (iotrmb,10098) n, ustf(m,1), stflx(m,n)*volt(n)/areat(n)
     $,                        ustf(m,2)
1810    continue
        write (iotrmb,10098) n, ' tot heat = ', (atwv(1,n)*volt(n))
     $,                                        ' deg C * cm**3 '
        write (iotrmb,10098) n, ' sst      = ',asst(1,n)
     $,                                       ' deg C         '
c
        do 1900 n=0,numreg
#ifdef kdtermb
c
c         supress instantaneous regional term balance writes when
c         the option to produce time averages has been selected
c
          if (n .GT. 0) goto 1875
c
#endif
          if (n .eq. 1) then
            write (iotrmb,10050)
     $  'Regional averaged Momentum & Tracer Term Balances for ts =   '
     $,     itt, stamp
          endif
          iv = 0
          if (n .gt. 0) then
            iv = (n-1)/nhreg + 1
            ih = n - (iv-1)*nhreg
            write (iotrmb,10100) 'Momentum terms averaged over region #'
     $,          n, ': ', hregnm(ih), vregnm(iv), volu(n), areau(n)
            write (iotrmb,10104)
            write (iotrmb,10101) n,' dU/dt   = ', termbm(0,1,1,n)
     $,                            ' dV/dt   = ', termbm(0,1,2,n)
            write (iotrmb,10101) n,' -Px     = ', termbm(0,2,1,n)
     $,                            ' -Py     = ', termbm(0,2,2,n)
            write (iotrmb,10101) n,' -surf Px= ', termbm(0,12,1,n)
     $,                            ' -surf Py= ', termbm(0,12,2,n)
            write (iotrmb,10101) n,' -UUx    = ', termbm(0,3,1,n)
     $,                            ' -UVx    = ', termbm(0,3,2,n)
            write (iotrmb,10101) n,' -VUy    = ', termbm(0,4,1,n)
     $,                            ' -VVy    = ', termbm(0,4,2,n)
            write (iotrmb,10101) n,' -WUz    = ', termbm(0,5,1,n)
     $,                            ' -WVz    = ', termbm(0,5,2,n)
            write (iotrmb,10101) n,'  UVtan  = ', termbm(0,13,1,n)
     $,                            ' -UUtan  = ', termbm(0,13,2,n)
#ifdef biharmonic
            write (iotrmb,10101) n,'  Uxxxx  = ', termbm(0,6,1,n)
     $,                            '  Vxxxx  = ', termbm(0,6,2,n)
            write (iotrmb,10101) n,'  Uyyyy  = ', termbm(0,7,1,n)
     $,                            '  Vyyyy  = ', termbm(0,7,2,n)
#else
            write (iotrmb,10101) n,'  Uxx    = ', termbm(0,6,1,n)
     $,                            '  Vxx    = ', termbm(0,6,2,n)
            write (iotrmb,10101) n,'  Uyy    = ', termbm(0,7,1,n)
     $,                            '  Vyy    = ', termbm(0,7,2,n)
#endif
            write (iotrmb,10101) n,'  Uzz    = ', termbm(0,8,1,n)
     $,                            '  Vzz    = ', termbm(0,8,2,n)
            write (iotrmb,10101) n,'  metric = ', termbm(0,9,1,n)
     $,                            '  metric = ', termbm(0,9,2,n)
            write (iotrmb,10101) n,'  fV     = ', termbm(0,10,1,n)
     $,                            ' -fU     = ', termbm(0,10,2,n)
            write (iotrmb,10101) n,'  source = ', termbm(0,11,1,n)
     $,                            '  source = ', termbm(0,11,2,n)
            erru = c0
            errv = c0
            do 1870 lll=2,13
              erru = erru + termbm(0,lll,1,n)
              errv = errv + termbm(0,lll,2,n)
1870        continue
            write (iotrmb,10101) n,'  error  = ', termbm(0,1,1,n) - erru
     $,                            '  error  = ', termbm(0,1,2,n) - errv
            write (iotrmb,10101) n,'  ubar   = ', avwv(1,n)
     $,                            '  vbar   = ', avwv(2,n)
     $,                            '  wbar   = ', avwv(3,n)
          endif
          if (iv .eq. 1) then
            write (iotrmb,10101) n,'  surf Uz= ', smflx(1,n)
     $,                            '  surf Vz= ', smflx(2,n)
          endif
c
#ifdef kdtermb
1875      continue
c
#endif
          if (n .eq. 0) then
            write (iotrmb,10051)
     $   'Global averaged (entire ocean) Tracer Term Balances for ts = '
     $,     itt, stamp
#ifndef kdtermb
          else
            write (iotrmb,10100) 'Tracer terms averaged over region   #'
     $,          n, ': ', hregnm(ih), vregnm(iv), volt(n), areat(n)
#endif
          endif
c
          do 1880 m=1,nt
            dchg            = termbt(0,2,m,n) + termbt(0,3,m,n) +
     $                        termbt(0,4,m,n) + termbt(0,5,m,n) +
     $                        termbt(0,6,m,n) + termbt(0,7,m,n) +
     $                        termbt(0,8,m,n)
#ifdef kdtermb
     $                      + termbt(0,9,m,n)
            tfilt(m)        = termbt(0,1,m,n) - dchg
#else
            tconv(m)        = termbt(0,1,m,n) - dchg
#endif
1880      continue
#ifdef kdtermb
c
c         supress instantaneous regional term balance writes when
c         the option to produce time averages has been selected
c
            if (n .GT. 0) goto 1895
c
#endif
          maxm = (nt-1)/7 + 1
          do 1890 mloop=1,maxm
            ms = (mloop-1)*7 + 1
            me = min(ms + 6,nt)
            write (iotrmb,10103) (trname(m),m=ms,me)
            write (iotrmb,10102) n,' dT/dt   = ', (termbt(0,1,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -UTx    = ', (termbt(0,2,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -VTy    = ', (termbt(0,3,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -WTz    = ', (termbt(0,4,m,n)
     $,                                           m=ms,me)
#if defined biharmonic && !defined isopycmix
            write (iotrmb,10102) n,'  Txxxx  = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tyyyy  = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
#else
# ifdef isopycmix
            write (iotrmb,10102) n,'  Tisox  = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tisoy  = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
# else
            write (iotrmb,10102) n,'  Txx    = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tyy    = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
# endif
#endif
            write (iotrmb,10102) n,'  Tzz    = ', (termbt(0,7,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  source = ', (termbt(0,8,m,n)
     $,                                           m=ms,me)
#ifdef kdtermb
            write (iotrmb,10102) n,'  convect= ', (termbt(0,9,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  filter+= ', (tfilt(m)
     $,                                           m=ms,me)
#else
            write (iotrmb,10102) n,'  convect= ', (tconv(m)
     $,                                           m=ms,me)
#endif
            write (iotrmb,10102) n,'  chg var= ', (tvar(0,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  tbar   = ', (atwv(m,n)
     $,                                           m=ms,me)
            if (iv .eq. 1) then
              write (iotrmb,10102) n,'  surflx = ', (stflx(m,n)
     $,                                           m=ms,me)
              taux = smflx(1,n)*volu(n)/areau(n)
              tauy = smflx(2,n)*volu(n)/areau(n)
              write (iotrmb,10105) ' Regionally averaged quantities:'
              write (iotrmb,'(1x)')
              write (iotrmb,10098) n
     $,          ' smf(1)   = ', taux,' dynes/cm**2   '
              write (iotrmb,10098) n
     $,          ' smf(2)   = ', tauy,' dynes/cm**2   '
              do 1820 m=1,nt
                write (iotrmb,10098) n
     $,          ustf(m,1), stflx(m,n)*volt(n)/areat(n), ustf(m,2)
1820          continue
              if (ms .eq. 1) then
                write (iotrmb,10098) n
     $,          ' tot heat = ', (atwv(1,n)*volt(n)),' deg C * cm**3 '
                write (iotrmb,10098) n
     $,          ' sst      = ', asst(1,n),' deg C         '
              endif
            endif
1890      continue
#ifdef kdtermb
1895      continue
#endif
1900    continue
c
#ifdef kdtermb
c
c  write term balance header file
c
      write(iointb) nt, km, nhreg, nvreg
c
c  write out info on time level and region sizes
c
      write(iointb) itt, stamp, expnam, trname
     $,            hregnm, dzt, llvreg
     $,            volt, volu, areau, areat
c
c  write out region momentum term balances
c
      write(iointb) (((termbm(0,itype,m,n),itype=1,13),
     $                         m=1,2),n=1,numreg)
c
c  write out region tracer term balances
c
      write(iointb) (((termbt(0,itype,m,n),itype=1,9),
     $                         m=1,nt),n=0,numreg)
      write(iointb) atwv
      write(iointb) asst
      write(iointb) stflx
      write(iointb) ustf
c
      write(stdout,1432) stamp, iointb
1432  format(' Instantaneous Term Balances Written On ',a32,
     $        ' To Internal Unit ',i3,/)
c
#endif
c-----------------------------------------------------------------------
c       print the northward transport of heat and salt
c-----------------------------------------------------------------------
c
#ifndef kdtermb
        write (stdout,8195)
#endif
c
c       convert heat transport to petawatts,
c       salt trnspt to 10**10 cm**3/sec
c
#ifdef testnorth
        do 2010 nsectr=0,2
          do 2000 l=1,jmt
            do 1990 ll=1,8
              ttn(ll,l,1,nsectr)=ttn(ll,l,1,nsectr)*4.186e-15
              ttn(ll,l,2,nsectr)=ttn(ll,l,2,nsectr)*1.e-10
# ifdef agetracers
              ttn(ll,l,ntage1,nsectr)=ttn(ll,l,ntage1,nsectr)*1.e-16
# endif
1990        continue
2000      continue
2010    continue
#else
        do 2000 l=1,jmt
          do 1990 ll=1,8
            ttn(ll,l,1)=ttn(ll,l,1)*4.186e-15
            ttn(ll,l,2)=ttn(ll,l,2)*1.e-10
# ifdef agetracers
            ttn(ll,l,ntage1)=ttn(ll,l,ntage1)*1.e-16
# endif
1990      continue
2000    continue
#endif
#ifdef kdtermb
      write(iointb) ttn
# ifdef testnorth
        write(iointb) tmt
# endif
#else
# ifdef testnorth
        do 2110 nsectr=0,2
          write(stdout,2117) nornam(nsectr)
          do 2100 j=2,jmtm2
            l=jmt-j
            write (stdout,8196) l, (ttn(i,l,1,nsectr),i=1,8),
     $                             (ttn(i,l,2,nsectr),i=1,8)
2100      continue
2110    continue
2117    format(//' Northward Transports ==== ',a16)
# else
        do 2100 j=2,jmtm2
          l=jmt-j
          write (stdout,8196) l, (ttn(i,l,1),i=1,8), (ttn(i,l,2),i=1,8)
2100    continue
# endif
#endif
#if defined cpldleg3 || defined cpldleg2
        if (prdgts) then
#endif
        scl = 1.e12
        write (stdout,8194) scl
        js = indp (slatxy, yt, jmt)
        je = indp (elatxy, yt, jmt)
# ifdef testnorth
        write(stdout,2118) nornam(0)
        call matrix (tmt(1,1,0), jmt, js, je, 1, km, scl)
        write(stdout,2118) nornam(1)
        call matrix (tmt(1,1,1), jmt, js, je, 1, km, scl)
        write(stdout,2118) nornam(2)
        call matrix (tmt(1,1,2), jmt, js, je, 1, km, scl)
2118    format(//' Meridional Mass Transport ==== ',a16)
# else
        call matrix (tmt, jmt, js, je, 1, km, scl)
# endif
#if defined cpldleg3 || defined cpldleg2
        endif
#endif
c
c-----------------------------------------------------------------------
c       print the stream function
c-----------------------------------------------------------------------
c
#if defined cpldleg3 || defined cpldleg2
        if (prdgts) then
#endif
        scl=1.e12
        is = indp (slonxy, xt, imt)
        ie = indp (elonxy, xt, imt)
        js = indp (slatxy, yt, jmt)
        je = indp (elatxy, yt, jmt)
        write (stdout,8000) ' stream function (sverdrups)'
     $, itt, xt(is), xt(ie), yt(js), yt(je), scl
        call matrix (p(1,1,1), imt, is, ie, -js, -je, scl)
#ifdef neptune
        scl=1.e12
        write (stdout,8000) ' neptune effect stream functn (sverdrups)'
     $, itt, xt(is), xt(ie), yt(js), yt(je), scl
        call matrix (pnept(1,1), imt, is, ie, -js, -je, scl)
#endif
#if defined cpldleg3 || defined cpldleg2
        endif
#endif
c
      return
8000  format(1x,a28,1x,'ts=',i7
     $,', lon:',f6.2,' ==> ',f6.2,',   lat:',f6.2,' ==> ',f6.2
     $,', scaling=',1pg10.3)
9100  format(/,1x,13('**********')
     $/1x,a43,a9,i10,a33/1x,'ocean volume =',e16.9,' cm**3'
     $, ', ocean surface area =',e16.9,' cm**2'//' work by:',14x
     $,     'internal mode                         external mode'/)
9101  format(a21,2(1pe15.6, ' (',z16,' hex)'))
9110  format(/' work by buoyancy forces =',1pe14.6, ' (',z16,' hex)'/
     $,      ' energy conversion error =',1pe14.6, ' (',z16,' hex)'/
     $,      ' nonlinear exchange error=',1pe14.6, ' (',z16,' hex)'/)
10050 format(/,1x,13('**********'),/1x,a61,i10,a33/)
10051 format(/,1x,a61,i10,a33/)
10098 format(1x,'(',i4,')',a12,(1pe16.8,a15))
10100 format(/1x,a37,i5,a2,a40,a20/1x,'ocean volume =',e16.9,' cm**3'
     $, ', horizontal ocean area =',e16.9,' cm**2'/)
10101 format(1x,'(',i4,')',a11,1pe15.7, 2x, a11,1pe15.7, 2x, a11
     $,      1pe15.7)
10102 format(1x,'(',i4,')',a11,7(1pe15.7))
10103 format(' Region',11x,8a15)
10104 format(' Region')
10105 format(/8x,a32)
10106 format(/,1x,13('**********')
     $/1x,a43,a9,i10,a33/1x,'ocean "t" volume =',e16.9,' cm**3 '
     $, ', ocean "t" surface area =',e16.9,' cm**2 '
     $/ 1x,'ocean "u" volume =',e16.9,' cm**3, '
     $, ' ocean "u" surface area =',e16.9,' cm**2'/)
8195  format(/,6x,'northward transport of heat (x10**15 watts)'
     $       ,22x,'northward transport of salt (x10**10 cm**3/sec)',/,
     $       3x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     $       'diffus   total'))
8196  format(i4,8f8.3,1x,8f8.3)
8194  format(/,' meridional mass transport: scaling=',g10.3)
      end
