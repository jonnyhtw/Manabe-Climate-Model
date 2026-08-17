      subroutine diag (j)
      implicit real*8 (a-h,o-z)
cfpp$ noconcur r
c
c=======================================================================
c
c     initialize, calculate & printout of diagnostics
c
c=======================================================================
c
#include "param.h"
#ifdef biharmonic
# include "cbihar.h"
#endif
#ifdef testcfl
# include "ccfl.h"
#endif
#include "cdiag.h"
#include "chmix.h"
#ifdef isopycmix
# include "cisop.h"
#endif
#include "coord.h"
#ifdef nlhmix
# include "cnlmix.h"
#endif
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
#include "slabs.h"
#include "switch.h"
c
      if (j .eq. 1) then
c
        write (stdout,9300) ' MOM"s diagnostic calculations for ts ='
     $,  itt, stamp
c
c-----------------------------------------------------------------------
c       initialize various quantities when j = 1
c-----------------------------------------------------------------------
c
#ifdef testcfl
        cflup = c0
        cflvp = c0
        cflwp = c0
#endif
        do 100 n=1,nt
          ustf(n,1) = ' unknown  '
          ustf(n,2) = ' unknown units '
          if (n .eq. 1) then
            ustf(n,1) = ' stf(1)   = '
            ustf(n,2) = ' cal/cm**2/sec '
          endif
          if (n .eq. 2) then
           ustf(n,1) = ' stf(2)   = '
           ustf(n,2) = ' gm/cm**2/sec  '
          endif
100     continue
#if defined cpldleg3 || defined cpldleg2
        ustf(3,1) = ' tracer 1='
        ustf(3,2) = '               '
#endif
        do 200 k=0,km
          buoy(k) = c0
200     continue
        do 300 ll=1,8
          do 290 k=0,km
            engint(k,ll) = c0
290       continue
300     continue
        do 400 ll=1,8
          engext(ll) = c0
          do 390 i=1,imt
            zuseng(i,ll) = c0
            zvseng(i,ll) = c0
390       continue
400     continue
c
        do 500 n=0,numreg
          if (n .gt. 0) then
            avwv(1,n) = c0
            avwv(2,n) = c0
            avwv(3,n) = c0
            do 450 ll=1,13
              do 440 k=0,km
                termbm(k,ll,1,n) = c0
                termbm(k,ll,2,n) = c0
440           continue
450         continue
          endif
          do 490 m=1,nt
            atwv(m,n) = c0
            do 460 k=0,km
              tvar(k,m,n) = c0
460         continue
#ifdef kdtermb
            do 480 ll=1,9
#else
            do 480 ll=1,8
#endif
              do 470 k=0,km
                termbt(k,ll,m,n) = c0
470           continue
480         continue
490       continue
500     continue
        do 520 n=0,nhreg
          smflx(1,n) = c0
          smflx(2,n) = c0
          do 510 m=1,nt
            stflx(m,n) = c0
            asst(m,n)  = c0
510       continue
520     continue
c
#ifdef testnorth
        do 610 nsectr=0,2
          do 600 jj=1,jmt
            do 580 m=1,ntmin2
              do 570 ll=1,8
              ttn(ll,jj,m,nsectr) = c0
570           continue
580         continue
            do 590 k=1,km
CCCCTEMP          tmt(jj,k)=c0
              tmt(jj,k,nsectr)=c0
590         continue
600       continue
610     continue
#else
        do 600 jj=1,jmt
          do 580 m=1,ntmin2
            do 570 ll=1,8
              ttn(ll,jj,m) = c0
570         continue
580       continue
          do 590 k=1,km
            tmt(jj,k)=c0
590       continue
600     continue
#endif
c
#ifdef zonalavg
c
c       Set initial zonal means to zero
c
        call zonalj(j)
#endif
c
      elseif (j .lt. jmt) then
c
c-----------------------------------------------------------------------
c       print the solution at specified latitudes
c-----------------------------------------------------------------------
c
        do 700 jlat=1,nlatpr
          jj = indp (prlat(jlat), yt, jmt)
          if (jj .eq. j .and. prlat(jlat) .le. yt(jmt)) then
            if (jlat .eq. 1) write(stdout,9000)
            is = indp (prslon(jlat), xt, imt)
            ie = indp (prelon(jlat), xt, imt)
            ks = indp (prsdpt(jlat), zt, km)
            ke = indp (predpt(jlat), zt, km)
            fx = 1.0e-2
            do 680 m=1,nt
              scl = c1
              if (m .eq. 2) scl=1.e-3
              write (stdout,9100) trname(m), itt
     $,        j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (t(1,1,jc,nc,m), imt, is, ie, ks, ke, scl)
680         continue
#ifdef cyclic
c
c           set cyclic boundary condition on w before printing
c
            do 690 k=0,km
              w(1  ,k) = w(imtm1,k)
              w(imt,k) = w(2    ,k)
690         continue
#endif
            scl = 1.e-3
            write (stdout,9100)  'w velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (w, imt, is, ie, ks, ke, scl)
c
            scl = c1
            write (stdout,9100) 'u velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (u(1,1,jc,nc), imt, is, ie, ks, ke, scl)
c
            scl = c1
            write (stdout,9100) 'v velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (v(1,1,jc,nc), imt, is, ie, ks, ke, scl)
          endif
700     continue
#ifdef zonalavg
c
c-----------------------------------------------------------------------
c        compute zonal means for each tracer and for u, v, w
# ifdef seaice
c        and sea ice
# endif
c
         call zonalj(j)
c
c-----------------------------------------------------------------------
c
#endif
c
c-----------------------------------------------------------------------
c       compute the northward transport of each tracer quantity
c       as well as the zonally integrated meridional mass transport
c-----------------------------------------------------------------------
c
c    the following cfpp$ directive fixes a bug in fpp on the cray
c
cfpp$ novector r
        if(j .lt. jmtm1) then
#ifdef testnorth
c
c        the 1313 do loop allows northward tracer and meridional mass
c        transport statistics to be computed over Global, Atlantic and
c        Indo-Pacific regions according to mask "normsk"
c
         do 1313 nsectr=0,2
#endif
          do 800 k=1,km
            vbr(k)=c0
            do 790 m=1,nt
#ifdef testnorth
              tbrs(k,m)=tbrn(k,m,nsectr)
              tbrn(k,m,nsectr)=c0
# else
              tbrs(k,m)=tbrn(k,m)
              tbrn(k,m)=c0
# endif
790         continue
800       continue
          if(j.gt.2) go to 900
          do 890 m=1,nt
            do 880 k=1,km
              tbrs(k,m)=c0
880         continue
890       continue
          do 870 k=1,km
            totdx=c0
            do 860 i=2,imtm1
#ifdef testnorth
              if (nsectr .GT. 0 .AND. normsk(i,j) .NE. nsectr) goto 860
#endif
              totdx=totdx+dxt(i)*fm(i,k,jc)
              do 850 m=1,nt
                tbrs(k,m)=tbrs(k,m)+t(i,k,jc,nc,m)*fm(i,k,jc)*dxt(i)
850           continue
860         continue
            if(totdx.ne.c0) then
              do 840 m=1,nt
                tbrs(k,m)=tbrs(k,m)/totdx
840           continue
            endif
870       continue
900       continue
          ccttj = asubh*dyur(j)
          do 1000 k=1,km
            totdx=c0
            do 990 i=2,imtm1
#ifdef testnorth
              if (nsectr .GT. 0 .AND. normsk(i,j) .NE. nsectr) goto 990
#endif
              totdx=totdx+dxt(i)*fm(i,k,jp)
              vbr(k)=vbr(k)+v(i,k,jc,nc)*dxu(i)*csu(j)
              do 980 m=1,nt
#ifdef testnorth
                tbrn(k,m,nsectr)=tbrn(k,m,nsectr)+
     $              t(i,k,jp,nc,m)*fm(i,k,jp)*dxt(i)
#else
                tbrn(k,m)=tbrn(k,m)+t(i,k,jp,nc,m)*fm(i,k,jp)*dxt(i)
#endif
980           continue
990         continue
            if(totdx.ne.c0) then
              do 970 m=1,nt
#ifdef testnorth
                tbrn(k,m,nsectr)=tbrn(k,m,nsectr)/totdx
#else
                tbrn(k,m)=tbrn(k,m)/totdx
#endif
970           continue
            endif
#ifdef testnorth
            if(k.eq.1) tmt(j,1,nsectr)=vbr(1)*dzt(1)
            if(k.gt.1) tmt(j,k,nsectr)=tmt(j,k-1,nsectr)+vbr(k)*dzt(k)
            do 960 m=1,nt
              ttn(1,j,m,nsectr)=ttn(1,j,m,nsectr)+
     $                   vbr(k)*(tbrn(k,m,nsectr)+tbrs(k,m))*p5*dzt(k)
              do 950 i=2,imtm1
              if (nsectr .GT. 0 .AND. normsk(i,j) .NE. nsectr) goto 950
                ttn(6,j,m,nsectr)=ttn(6,j,m,nsectr)+
     $             (v(i,k,jc,nc)*dxu(i)+v(i-1,k,jc,nc)*dxu(i-1))*
     $             (t(i,k,jc,nc,m)+t(i,k,jp,nc,m))*csu(j)*p25*dzt(k)
# ifdef isopycmix
                ttn(7,j,m,nsectr)=ttn(7,j,m,nsectr)-
     $              esav(i,k,m)*fm(i,k,jc)*fm(i,k,jp)*dxt(i)*dzt(k)
# else
#  ifdef biharmonic
                ttn(7,j,m,nsectr)=ttn(7,j,m,nsectr)-
     $             ccttj*fm(i,k,jc)*fm(i,k,jp)*
     $             (del2(i,k,m,jpt3)-del2(i,k,m,jpt2))*dxt(i)*csu(j)
     $             *dzt(k)
#  endif
#  ifdef nlhmix
                ttn(7,j,m,nsectr)=ttn(7,j,m,nsectr)-
     $             (cctj(i,k)/cct(j))*dyur(j)*fm(i,k,jc)*fm(i,k,jp)*
     $             (t(i,k,jp,nc,m)-t(i,k,jc,nc,m))*dxt(i)*csu(j)*dzt(k)
#  endif
#  ifdef consthmix
                ttn(7,j,m,nsectr)=ttn(7,j,m,nsectr)-
     $             ccttj*fm(i,k,jc)*fm(i,k,jp)*
     $             (t(i,k,jp,nc,m)-t(i,k,jc,nc,m))*dxt(i)*csu(j)*dzt(k)
     $             *hdcsck(k)
#  endif
# endif
950           continue
960         continue
#else
            if(k.eq.1) tmt(j,1)=vbr(1)*dzt(1)
            if(k.gt.1) tmt(j,k)=tmt(j,k-1)+vbr(k)*dzt(k)
            do 960 m=1,nt
              ttn(1,j,m)=ttn(1,j,m)+vbr(k)*(tbrn(k,m)+tbrs(k,m))*
     $                   p5*dzt(k)
              do 950 i=2,imtm1
                ttn(6,j,m)=ttn(6,j,m)+
     $                  (v(i,k,jc,nc)*dxu(i)+v(i-1,k,jc,nc)*dxu(i-1))*
     $             (t(i,k,jc,nc,m)+t(i,k,jp,nc,m))*csu(j)*p25*dzt(k)
# ifdef isopycmix
                ttn(7,j,m)=ttn(7,j,m)-
     $              esav(i,k,m)*fm(i,k,jc)*fm(i,k,jp)*dxt(i)*dzt(k)
# else
# ifdef biharmonic
                ttn(7,j,m)=ttn(7,j,m)-ccttj*fm(i,k,jc)*fm(i,k,jp)*
     $             (del2(i,k,m,jpt3)-del2(i,k,m,jpt2))*dxt(i)*csu(j)
     $             *dzt(k)
# endif
# ifdef nlhmix
                ttn(7,j,m)=ttn(7,j,m)-(cctj(i,k)/cct(j))*dyur(j)
     $                     *fm(i,k,jc)*fm(i,k,jp)*
     $             (t(i,k,jp,nc,m)-t(i,k,jc,nc,m))*dxt(i)*csu(j)*dzt(k)
# endif
# ifdef consthmix
                ttn(7,j,m)=ttn(7,j,m)-ccttj*fm(i,k,jc)*fm(i,k,jp)*
     $             (t(i,k,jp,nc,m)-t(i,k,jc,nc,m))*dxt(i)*csu(j)*dzt(k)
     $             *hdcsck(k)
# endif
#endif
950           continue
960         continue
#endif
1000      continue
          do 1100 m=1,nt
            do 1090 i=2,imtm1
#ifdef testnorth
              if (nsectr .GT. 0 .AND. normsk(i,j) .NE. nsectr) goto 1090
#endif
              totdz=c0
              vbrz=c0
              tbrz=c0
              do 1080 k=1,km
                if (kmt(i,j) .ge. k .and. kmt(i,j+1) .ge. k) then
                  vbrz=vbrz+(v(i,k,jc,nc)*dxu(i)+v(i-1,k,jc,nc)*
     $                       dxu(i-1))*dzt(k)
                  tbrz=tbrz+(t(i,k,jc,nc,m)+t(i,k,jp,nc,m))*dzt(k)
                  totdz=totdz+dzt(k)
                endif
1080          continue
              if(totdz.eq.c0) go to 1070
                tbrz=tbrz/totdz
#ifdef testnorth
                ttn(3,j,m,nsectr)=ttn(3,j,m,nsectr)+vbrz*tbrz*csu(j)*p25
                ttn(5,j,m,nsectr)=ttn(5,j,m,nsectr)-
     $                    (smf(i,1)*dxu(i)+smf(i-1,1)*
     $                    dxu(i-1))*(t(i,1,jc,nc,m)+t(i,1,jp,nc,m)-
     $                    tbrz)*csu(j)/
     $                    (8.0*omega*sine(j))
#else
                ttn(3,j,m)=ttn(3,j,m)+vbrz*tbrz*csu(j)*p25
                ttn(5,j,m)=ttn(5,j,m)-(smf(i,1)*dxu(i)+smf(i-1,1)*
     $                    dxu(i-1))*(t(i,1,jc,nc,m)+t(i,1,jp,nc,m)-
     $                    tbrz)*csu(j)/
     $                    (8.0*omega*sine(j))
#endif
1070          continue
1090        continue
1100      continue
#ifdef testnorth
          do 1200 m=1,nt
            ttn(2,j,m,nsectr)=ttn(6,j,m,nsectr)-ttn(1,j,m,nsectr)
            ttn(4,j,m,nsectr)=ttn(6,j,m,nsectr)-ttn(3,j,m,nsectr)-
     $                 ttn(5,j,m,nsectr)
            ttn(8,j,m,nsectr)=ttn(6,j,m,nsectr)+ttn(7,j,m,nsectr)
1200      continue
1313     continue
#else
          do 1200 m=1,nt
            ttn(2,j,m)=ttn(6,j,m)-ttn(1,j,m)
            ttn(4,j,m)=ttn(6,j,m)-ttn(3,j,m)-ttn(5,j,m)
            ttn(8,j,m)=ttn(6,j,m)+ttn(7,j,m)
1200      continue
#endif
        endif
      endif
c
      return
9000  format(/' Zonal section printouts at specified latitudes:'/)
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     $,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     $,', scaling=',1pg10.3)
9300  format(/25x,a39,i10,a33/)
      end
