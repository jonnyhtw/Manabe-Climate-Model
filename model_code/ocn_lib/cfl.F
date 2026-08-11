      subroutine cfl (j)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
#ifdef testcfl
#include "param.h"
#include "accel.h"
#include "ccfl.h"
#ifdef cflj
#include "ccflj.h"
#endif
#include "coord.h"
#include "ctmngr.h"
#include "docnam.h"
#include "grdvar.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
c
      cl = cflcrt * p5
      do 100 k=1,km
        dtmax = max(dtuv, dtts*dtxcel(k))
        f1    = dtmax*dyur(j)
        f2    = dtmax*dzwr(k)
#ifdef cflj
        cfluj(j,k) = c0
        cflvj(j,k) = c0
        cflwj(j,k) = c0
        imaxu(j,k) = 0
        imaxv(j,k) = 0
        imaxw(j,k) = 0
#endif
        do 90 i=2,imtm1
          cflu = abs(dtmax*dxur(i)*u(i,k,jc,nc))
          cflv = abs(f1*v(i,k,jc,nc))
          cflw = abs(f2*w(i,k))
#ifdef cflj
          if (cflu.gt.cfluj(j,k)) then
            cfluj(j,k) = cflu
            imaxu(j,k) = i
          endif
          if (cflv.gt.cflvj(j,k)) then
            cflvj(j,k) = cflv
            imaxv(j,k) = i
          endif
          if (cflw.gt.cflwj(j,k)) then
            cflwj(j,k) = cflw
            imaxw(j,k) = i
          endif
#endif
          if (cflu .ge. cl .or. cflv .ge. cl .or. cflw .ge. cl) then
            write (stdout,'(/,a33,i4,a1,i3,a1,i3,a1,/)')
     $       ' ==> cfl exceeded at coordinate (i,j,k) = (',i
     $,      ',',j,',',k,')'
            umax  = p5*dxu(i)/dtmax
            pcflu = abs(100.0*u(i,k,jc,nc)/umax)
            vmax  = p5*dyu(j)/dtmax
            pcflv = abs(100.0*v(i,k,jc,nc)/vmax)
            wmax  = p5*dzw(k)/dtmax
            pcflw = abs(100.0*w(i,k)/wmax)
            write (stdout,'(a15,f8.2,a21,g15.8,a1)')
     $       ' u reached ', pcflu,' % of the CFL limit (',umax,')'
            write (stdout,'(a15,f8.2,a21,g15.8,a1)')
     $       ' v reached ', pcflv,' % of the CFL limit (',vmax,')'
            write (stdout,'(a15,f8.2,a21,g15.8,a1)')
     $       ' w reached ', pcflw,' % of the CFL limit (',wmax,')'
            is = max(1,i-3)
            ie = min(imt,i+3)
            ks = max(1,k-3)
            ke = min(km,k+3)
            scl = c0
            fx  = 1.0e-2
            write (stdout,9100)'u velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (u(1,1,jc,nc), imt, is, ie, ks, ke, scl)
c
            write (stdout,9100) 'v velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (v(1,1,jc,nc), imt, is, ie, ks, ke, scl)
c
            write (stdout,9100)  'w velocity', itt
     $,      j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
            call matrix (w, imt, is, ie, ks, ke, scl)
c
            do 80 m=1,nt
              write (stdout,9100) trname(m), itt
     $,        j, yt(j), xt(is), xt(ie), fx*zt(ks), fx*zt(ke), scl
              call matrix (t(1,1,jc,nc,m), imt, is, ie, ks, ke, scl)
80          continue
#ifdef cflj
            write(9009) j
9009        format(/,1x,'at j=',i4,'  cfluj(j,k)  imaxu(j,k)'
     $,                            '  cflvj(j,k)  imaxv(j,k)'
     $,                            '  cflwj(j,k)  imaxw(j,k)')
            do 120 kw=1,km
              write(stdout,9010) cfluj(j,kw),imaxu(j,kw)
     $,                          cflvj(j,kw),imaxv(j,kw)
     $,                          cflwj(j,kw),imaxw(j,kw)
9010          format(10x,6e12.5)
120         continue
#endif
            stop 'CFL'
          endif
90      continue
100   continue
      if (diagts .and. eots) then
        do 200 k=1,km
          dtmax = max(dtuv, dtts*dtxcel(k))
          vmax  = p5*dyu(j)/dtmax
          wmax  = p5*dzw(k)/dtmax
          do 190 i=2,imtm1
            umax  = p5*dxu(i)/dtmax
            if (abs(100.0*u(i,k,jc,nc)/umax) .gt. cflup) then
              cflup = abs(100.0*u(i,k,jc,nc)/umax)
              cflum  = u(i,k,jc,nc)
              icflu = i
              jcflu = j
              kcflu = k
            endif
            if (abs(100.0*v(i,k,jc,nc)/vmax) .gt. cflvp) then
              cflvp = abs(100.0*v(i,k,jc,nc)/vmax)
              cflvm  = v(i,k,jc,nc)
              icflv = i
              jcflv = j
              kcflv = k
            endif
            if (abs(100.0*w(i,k)/wmax) .gt. cflwp) then
              cflwp = abs(100.0*w(i,k)/wmax)
              cflwm  = w(i,k)
              icflw = i
              jcflw = j
              kcflw = k
            endif
190       continue
200     continue
      endif
9100  format(1x,a12,1x,'ts=',i10,1x,',j=',i3,', lat=',f6.2
     $,', lon:',f6.2,' ==> ',f6.2,', depth(m):',f6.1,' ==> ',f6.1
     $,', scaling=',1pg10.3)
#endif
      return
      end

