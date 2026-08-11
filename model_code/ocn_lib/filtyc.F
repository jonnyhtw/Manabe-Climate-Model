      subroutine filtyc (j)
      implicit real*8 (a-h,o-z)
c
#ifdef seaice
# if defined fourfil || defined firfil
#  ifdef multitasking
cfpp$ noconcur r
#  endif
c
c=======================================================================
c     subroutine filtyc sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of  sea ice at the specifiied high latitude row "j".
c=======================================================================
c
#include "param.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#include "seayce.h"
#include "slabs.h"
#include "temp.h"
c
#  ifdef fourfil
c
c---------------------------------------------------------------------
c     fourier filter sea ice at high latitudes
c---------------------------------------------------------------------
c
      if((j.gt.jft1.and.j.lt.jft2) .or. j.lt.jfrst) goto 101
      jj=j-jfrst+1
      if (j.ge.jft2) jj=jj-jskpt+1
c
      iredo=0
      do 100 l=1,lsegf
        if(istf(jj,l,1).eq.0) go to 90
        is=istf(jj,l,1)
        ie=ietf(jj,l,1)
        im=ie-is+1
#    ifdef cyclic
        if(im.ne.imtm2.or.kmt(1,j).lt.1) then
          m=1
          n=nint(im*cst(j)*cstr(jft0))
        else
          m=3
          n=nint(im*cst(j)*cstr(jft0)*0.5)
        endif
#    else
        m=1
        n=nint(im*cst(j)*cstr(jft0))
#    endif
        ism1=is-1
        iea=ie
        if(ie.ge.imt) iea=imtm1
        do 40 i=is,iea
          tempik(i-ism1,1,1)=yca(i,j)
40      continue
        if(ie.ge.imt) then
          ieb=ie-imtm2
          ii=imtm1-is
          do 50 i=2,ieb
            tempik(i+ii,1,1)=yca(i,j)
50        continue
        endif
c
        call filtr (tempik(1,1,1),im,m,n,iredo)
c
        do 60 i=is,iea
          yca(i,j)=tempik(i-ism1,1,1)
60      continue
        if(ie.ge.imt) then
          do 70 i=2,ieb
            yca(i,j)=tempik(i+ii,1,1)
70        continue
        endif
90      continue
100   continue
101   continue
#  endif
c
#  ifdef firfil
c
c-----------------------------------------------------------------------
c     filter tracers at high latitudes with symmetric finite impulse
c     response filter "fir"
c-----------------------------------------------------------------------
c
      if ((j .le. jft1 .or. j .ge. jft2) .and. j .ge. jfrst) then
        jj = j - jfrst + 1
        if (j .ge. jft2) jj = jj - jskpt + 1
        numfil = numflt(jj)
        call filfir (yca(1,j), fm(1,1,jc), tempik(1,1,1), 1, 0,
     $                 numfil)
      endif
#  endif
c
# endif
#endif
      return
      end
