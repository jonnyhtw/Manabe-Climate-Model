      subroutine filuv (j)
      implicit real*8 (a-h,o-z)
#if defined fourfil || defined firfil
c
c=====================================================================
c     filuv sets up input needed for fourier filtering
c     (when the "fourfil" ifdef is defined) -or- symmetric finite
c     impulse response filtering (when the "firfil" ifdef is defined)
c     of baroclinic velocities at the specifiied high latitude row "j".
c=====================================================================
c
#  ifdef multitasking
cfpp$ noconcur r
#  endif
#include "param.h"
#include "cpolar.h"
#include "emode.h"
#include "grdvar.h"
#include "index.h"
#include "scalar.h"
#include "slabs.h"
#include "switch.h"
#include "temp.h"
c
#ifdef firfil
      dimension  ubar(imt), vbar(imt)
#endif
c
c=====================================================================
c
#  ifdef fourfil
c
c---------------------------------------------------------------------
c     fourier filter u and v at high latitudes
c---------------------------------------------------------------------
c
      if((j.gt.jfu1 .and. j.lt.jfu2) .or. j.lt.jfrst) goto 701
      jj = j - jfrst + 1
      if (j.ge.jfu2) jj = jj - jskpu + 1
      fx = -c1
      if (phi(j).gt.0.) fx = c1
      isave = 0
      ieave = 0
c
      do 200 l=1,lsegf
        do 190 k=1,km
          if (isuf(jj,l,k) .eq. 0) goto 190
          is = isuf(jj,l,k)
          ie = ieuf(jj,l,k)
          iredo = 1
          if (is.ne.isave .or. ie.ne.ieave) then
            iredo = 0
            im = ie - is + 1
            isave = is
            ieave = ie
#    ifdef cyclic
            if (im .ne. imtm2) then
              m = 2
              n = nint(im*csu(j)*csur(jfu0))
            else
              m = 3
              n = nint(im*csu(j)*csur(jfu0)*p5)
            endif
#    else
           m = 2
           n = nint(im*csu(j)*csur(jfu0))
#    endif
          endif
          ism1 = is - 1
          iea = ie
          if(ie .ge. imu) iea = imum1
          do 150 i=is,iea
            tempik(i-ism1,k,1) = -fx*ua(i,k)*spsin(i) - va(i,k)*spcos(i)
            tempik(i-ism1,k,2) =  fx*ua(i,k)*spcos(i) - va(i,k)*spsin(i)
150       continue
          if (ie .ge. imu) then
            ieb = ie - imum2
            ii  = imum1 - is
            do 160 i=2,ieb
              tempik(i+ii,k,1) = -fx*ua(i,k)*spsin(i) - va(i,k)*spcos(i)
              tempik(i+ii,k,2) =  fx*ua(i,k)*spcos(i) - va(i,k)*spsin(i)
160        continue
          endif
c
          call filtr (tempik(1,k,1), im, m, n, iredo)
c
          call filtr (tempik(1,k,2), im, m, n, 1)
c
          do 170 i=is,iea
            ua(i,k) = fx*(-tempik(i-ism1,k,1)*spsin(i)
     $                   + tempik(i-ism1,k,2)*spcos(i))
            va(i,k) = -tempik(i-ism1,k,1)*spcos(i)
     $                   - tempik(i-ism1,k,2)*spsin(i)
170       continue
          if (ie .ge. imt) then
            do 180 i=2,ieb
              ua(i,k) = fx*(-tempik(i+ii,k,1)*spsin(i)
     $                     + tempik(i+ii,k,2)*spcos(i))
              va(i,k) = -tempik(i+ii,k,1)*spcos(i)
     $                     - tempik(i+ii,k,2)*spsin(i)
180         continue
          endif
190     continue
200   continue
c
      do 300 i=1,imt
        tempik(i,1,1) = c0
        tempik(i,1,2) = c0
300   continue
c
      do 400 k=1,km
        do 390 i=1,imt
          tempik(i,1,1) = tempik(i,1,1) + ua(i,k)*dzt(k)
          tempik(i,1,2) = tempik(i,1,2) + va(i,k)*dzt(k)
390     continue
400   continue
c
      do 500 i=1,imt
        tempik(i,1,1) = tempik(i,1,1)*hr(i,j)
        tempik(i,1,2) = tempik(i,1,2)*hr(i,j)
500   continue
c
      do 600 k=1,km
        do 590 i=1,imt
          ua(i,k) = ua(i,k) - tempik(i,1,1)
          va(i,k) = va(i,k) - tempik(i,1,2)
590     continue
600   continue
c
      do 700 k=1,km
        do 690 i=1,imt
          ua(i,k) = ua(i,k)*gm(i,k,jc)
          va(i,k) = va(i,k)*gm(i,k,jc)
690     continue
700   continue
c
701   continue
#  endif
#  ifdef firfil
c-----------------------------------------------------------------------
c     filter u and v at high latitudes with a finite impulse
c     response filter "fir"
c-----------------------------------------------------------------------
c
      if ((j .le. jfu1 .or. j .ge. jfu2) .and. j .ge. jfrst) then
        jj = j - jfrst + 1
        if (j .ge. jfu2) jj = jj - jskpu + 1
        numfil = numflu(jj)
c
        fx = -c1
        if (phi(j) .gt. c0) fx = c1
c
        do 720 k=1,km
          do 710 i=1,imt
            tempik(i,k,1) = -fx*ua(i,k)*spsin(i) - va(i,k)*spcos(i)
            tempik(i,k,2) =  fx*ua(i,k)*spcos(i) - va(i,k)*spsin(i)
710       continue
720     continue
c
        call filfir (tempik(1,1,1), gm(1,1,jc), ua, km, 1, numfil)
        call filfir (tempik(1,1,2), gm(1,1,jc), va, km, 1, numfil)
c
        do 740 k=1,km
          do 730 i=1,imt
            ua(i,k) = fx*(-tempik(i,k,1)*spsin(i)
     $                   + tempik(i,k,2)*spcos(i))
            va(i,k) = -tempik(i,k,1)*spcos(i)
     $                   - tempik(i,k,2)*spsin(i)
730       continue
740     continue
c
c       extract fictitious vertical mean
c
        do 800 i=2,imtm1
          ubar(i) = c0
          vbar(i) = c0
800     continue
c
        do 900 k=1,km
          do 890 i=2,imtm1
            ubar(i) = ubar(i) + ua(i,k)*dzt(k)
            vbar(i) = vbar(i) + va(i,k)*dzt(k)
890       continue
900     continue
c
        do 1000 i=2,imtm1
          ubar(i) = ubar(i)*hr(i,j)
          vbar(i) = vbar(i)*hr(i,j)
1000    continue
c
        do 1100 k=1,km
          do 1090 i=2,imtm1
            ua(i,k) = ua(i,k) - gm(i,k,jc)*ubar(i)
            va(i,k) = va(i,k) - gm(i,k,jc)*vbar(i)
1090      continue
1100    continue
c
      endif
#  endif
c
#endif
      return
      end
