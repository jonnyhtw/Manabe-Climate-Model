      subroutine filz (j)
      implicit real*8 (a-h,o-z)
#if defined fourfil || defined firfil
c
c=======================================================================
c
c     filter vorticity at high latitudes
c
c=======================================================================
c
#include "param.h"
#include "cdiag.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "levind.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#ifdef multitasking
      dimension tempik(imt,km,1)
#else
#include "temp.h"
#endif
c
#ifdef firfil
      dimension ubar(imt), vbar(imt)
c
#endif
c=======================================================================
c
#  ifdef fourfil
        if ((j .le. jfu1 .or. j .ge. jfu2) .and. j .ge. jfrst) then
c
c---------------------------------------------------------------------
c         fourier filter ztd at high latitudes
c---------------------------------------------------------------------
c
          if (j .eq. jfu2) go to 101
c
            jj = j - jfrst + 1
            if (j .ge. jfu2) jj = jj - jskpu + 1
c
            do 100 l=1,lsegf
              is=iszf(jj,l)
              if(is.eq.0) go to 100
              ie=iezf(jj,l)
              do 80 ii=is,ie
                i=mod(ii-2,imtm2)+2
                tempik(ii+1-is,1,1)=ztd(i,j)
80            continue
              im=ie-is+1
#    ifdef cyclic
              if(im.ne.imtm2) then
                 m=1
                 n=nint(im*cst(j)*csur(jfu0))
              else
                 m=3
                 n=nint(im*cst(j)*csur(jfu0)*p5)
              endif
#    else
              m=1
              n=nint(im*cst(j)*csur(jfu0))
#    endif
c
              call filtr(tempik(1,1,1),im,m,n,0)
c
              do 90 ii=is,ie
                i=mod(ii-2,imtm2)+2
                ztd(i,j)=tempik(ii+1-is,1,1)
90            continue
100         continue
101       continue
        endif
#  endif
#  ifdef firfil
c---------------------------------------------------------------------
c       filter vorticity at high latitudes with a finite impulse
c       response filter "fir"
c---------------------------------------------------------------------
c
        if ((j .le. jfu1 .or. j .ge. jfu2) .and. j .ge. jfrst) then
c
          if (j .ne. jfu2) then
            jj = j - jfrst + 1
            if (j .ge. jfu2) jj = jj - jskpu + 1
            numfil = numflu(jj)
c
c           construct appropriate land/sea mask
c
            do 200 i=2,imtm1
              ubar(i) = c0
200         continue
c
            do 300 l=1,lseg
              isv = isz(j,l)
              if (isv .ne. 0) then
                iev = iez(j,l)
                do 290 i=isv,iev
                  ubar(i) = c1
290             continue
              endif
300         continue
c
#    ifdef cyclic
            ubar(1)   = ubar(imtm1)
            ubar(imt) = ubar(2)
#    else
            ubar(1)   = c0
            ubar(imt) = c0
#    endif
#    ifdef islands
c
c           ensure that ztd is zero next to islands
c
            do 400 i=1,imt
              tempik(i,1,1) = ztd(i,j)*ubar(i)
400         continue
c
            call filfir (tempik(1,1,1), ubar, vbar, 1, 0, numfil)
c
            do 500 i=1,imt
              if (ubar(i) .gt. c0) ztd(i,j) = tempik(i,1,1)
500         continue
#    else
            call filfir (ztd(1,j), ubar, vbar, 1, 0, numfil)
c
#    endif
          endif
        endif
#  endif
c
#endif
      return
      end
