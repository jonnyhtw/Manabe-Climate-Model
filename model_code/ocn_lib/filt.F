      subroutine filt (j)
      implicit real*8 (a-h,o-z)
#if defined fourfil || defined firfil
#  ifdef multitasking
cfpp$ noconcur r
#  endif
c
c=======================================================================
c     subroutine filt sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of tracers at the specifiied high latitude row "j".
c=======================================================================
c
#include "param.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#include "slabs.h"
#include "temp.h"
c
#  ifdef fourfil
c
c---------------------------------------------------------------------
c     fourier filter tracers at high latitudes
c---------------------------------------------------------------------
c
      if((j.gt.jft1.and.j.lt.jft2) .or. j.lt.jfrst) goto 101
      jj=j-jfrst+1
      if (j.ge.jft2) jj=jj-jskpt+1
c
c    if previous strips were of same length, do not recompute
c    fourier coeffs
c
      isave=0
      ieave=0
      do 100 l=1,lsegf
        do 90 k=1,km
          if(istf(jj,l,k).eq.0) go to 90
          is=istf(jj,l,k)
          ie=ietf(jj,l,k)
          iredo=0
          if(is.ne.isave .or. ie.ne.ieave) then
            iredo=-1
            isave=is
            ieave=ie
            im=ie-is+1
#    ifdef cyclic
            if(im.ne.imtm2.or.kmt(1,j).lt.k) then
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
          endif
          do 80 mm=1,nt
            idx=iredo+mm
            ism1=is-1
            iea=ie
            if(ie.ge.imt) iea=imtm1
            do 40 i=is,iea
              tempik(i-ism1,k,1)=ta(i,k,mm)
40          continue
            if(ie.ge.imt) then
              ieb=ie-imtm2
              ii=imtm1-is
              do 50 i=2,ieb
                tempik(i+ii,k,1)=ta(i,k,mm)
50            continue
            endif
c
            call filtr (tempik(1,k,1),im,m,n,idx)
c
            do 60 i=is,iea
              ta(i,k,mm)=tempik(i-ism1,k,1)
60          continue
            if(ie.ge.imt) then
              do 70 i=2,ieb
                ta(i,k,mm)=tempik(i+ii,k,1)
70            continue
            endif
80       continue
90      continue
100   continue
101   continue
#  endif
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
        do 200 mm=1,nt
          call filfir (ta(1,1,mm), fm(1,1,jc), tempik(1,1,1), km, 0,
     $                 numfil)
200     continue
      endif
#  endif
c
#endif
      return
      end
