      subroutine delsqr (j, jjm, jjc, jjp)
      implicit real*8 (a-h,o-z)
#ifdef biharmonic
# ifdef multitasking
cfpp$ noconcur r
# endif
c
c=======================================================================
c
c     compute del**2 for prognostic variables on row "j"
c     pointer (index) "jjm" points to the "j-1" row in memory
c     pointer (index) "jjc" points to the "j"   row in memory
c     pointer (index) "jjp" points to the "j+1" row in memory
c
c=======================================================================
c
#include "param.h"
#include "cbihar.h"
#include "chmix.h"
#include "grdvar.h"
#include "slabs.h"
c
c-----------------------------------------------------------------------
c     del**2 of tracers
c     pointer (index) "jpt3" points to the position within "del2"
c                      corresponding to "jjc"
c-----------------------------------------------------------------------
c
      do 100 m=1,nt
        do 70 k=1,km
          do 60 i=2,imt
            aux2(i,k) = dxu2r(i-1)*(t(i,k,jjc,nm,m)-t(i-1,k,jjc,nm,m))
60        continue
70      continue
c
        do 90 k=1,km
          do 80 i=2,imtm1
            del2(i,k,m,jpt3) = bbt(j)*dxt4r(i)*
     $           (fm(i+1,k,jjc)*aux2(i+1,k) - fm(i-1,k,jjc)*aux2(i,k))
     $         +  cct(j)*fm(i,k,jjp)*(t(i,k,jjp,nm,m)-t(i,k,jjc,nm,m))
     $         +  ddt(j)*fm(i,k,jjm)*(t(i,k,jjm,nm,m)-t(i,k,jjc,nm,m))
80        continue
#ifdef cyclic
          del2(1,k,m,jpt3)   = del2(imtm1,k,m,jpt3)
          del2(imt,k,m,jpt3) = del2(2,k,m,jpt3)
#else
          del2(1,k,m,jpt3)   = c0
          del2(imt,k,m,jpt3) = c0
#endif
90      continue
100   continue
c
c-----------------------------------------------------------------------
c     del**2 of momentum
c     pointer (index) "jpt3" points to the position within "del2"
c                      corresponding to "jjc"
c-----------------------------------------------------------------------
c
      do 200 k=1,km
        do 190 i=2,imt
          aux3(i,k) = dxt4r(i)*(u(i,k,jjc,nm)-u(i-1,k,jjc,nm))
          aux4(i,k) = dxt4r(i)*(v(i,k,jjc,nm)-v(i-1,k,jjc,nm))
190     continue
200   continue
c
c     zonal component
c
      m = nt + 1
      do 300 k=1,km
        do 290 i=2,imtm1
          del2(i,k,m,jpt3) = bbu(j)*(dxu2r(i)*(aux3(i+1,k) - aux3(i,k)))
     $                     + ccu(j)*(u(i,k,jjp,nm) - u(i,k,jjc,nm))
     $                     + ddu(j)*(u(i,k,jjm,nm) - u(i,k,jjc,nm))
     $                     + ggu(j)*u(i,k,jjc,nm) - hhu(j)*
     $                      dxu2r(i)*(v(i+1,k,jjc,nm) - v(i-1,k,jjc,nm))
290     continue
#ifdef cyclic
        del2(1,k,m,jpt3)   = del2(imtm1,k,m,jpt3)
        del2(imt,k,m,jpt3) = del2(2,k,m,jpt3)
#else
        del2(1,k,m,jpt3)   = c0
        del2(imt,k,m,jpt3) = c0
#endif
300   continue
c
c     meridional component
c
      m = nt + 2
      do 400 k=1,km
        do 390 i=2,imtm1
          del2(i,k,m,jpt3) = bbu(j)*(dxu2r(i)*(aux4(i+1,k) - aux4(i,k)))
     $                     + ccu(j)*(v(i,k,jjp,nm) - v(i,k,jjc,nm))
     $                     + ddu(j)*(v(i,k,jjm,nm) - v(i,k,jjc,nm))
     $                     + ggu(j)*v(i,k,jjc,nm) + hhu(j)*
     $                      dxu2r(i)*(u(i+1,k,jjc,nm) - u(i-1,k,jjc,nm))
390     continue
#ifdef cyclic
        del2(1,k,m,jpt3)   = del2(imtm1,k,m,jpt3)
        del2(imt,k,m,jpt3) = del2(2,k,m,jpt3)
#else
        del2(1,k,m,jpt3)   = c0
        del2(imt,k,m,jpt3) = c0
#endif
400   continue
#ifdef tcvmix
c
c-----------------------------------------------------------------------
c     del**2 of turbulent kenetic energy
c     pointer (index) "jpt3" points to the position within "del2"
c                      corresponding to "jjc"
c-----------------------------------------------------------------------
c
c       q2
c
        m = nt + 3
        do 500 k=1,km
          do 490 i=2,imt
            aux2(i,k) = dxu2r(i-1)*(q2(i,k,jjc,nm)-q2(i-1,k,jjc,nm))
490       continue
500     continue
c
        do 600 k=1,km
          do 590 i=2,imtm1
            del2(i,k,m,jpt3) = bbt(j)*dxt4r(i)*
     $           (fm(i+1,k,jjc)*aux2(i+1,k) - fm(i-1,k,jjc)*aux2(i,k))
     $         + cct(j)*fm(i,k,jjp)*(q2(i,k,jjp,nm)-q2(i,k,jjc,nm))
     $         + ddt(j)*fm(i,k,jjm)*(q2(i,k,jjm,nm)-q2(i,k,jjc,nm))
590       continue
# ifdef cyclic
          del2(1,k,m,jpt3)   = del2(imtm1,k,m,jpt3)
          del2(imt,k,m,jpt3) = del2(2,k,m,jpt3)
# else
          del2(1,k,m,jpt3)   = c0
          del2(imt,k,m,jpt3) = c0
# endif
600     continue
# ifdef leq
c
c       q2l
c
        m = nt + 4
        do 700 k=1,km
          do 690 i=2,imt
            aux2(i,k) = dxu2r(i-1)*(q2l(i,k,jjc,nm)-q2l(i-1,k,jjc,nm))
690       continue
700     continue
c
        do 800 k=1,km
          do 790 i=2,imtm1
            del2(i,k,m,jpt3) = bbt(j)*dxt4r(i)*
     $           (fm(i+1,k,jjc)*aux2(i+1,k) - fm(i-1,k,jjc)*aux2(i,k))
     $         + cct(j)*fm(i,k,jjp)*(q2l(i,k,jjp,nm)-q2l(i,k,jjc,nm))
     $         + ddt(j)*fm(i,k,jjm)*(q2l(i,k,jjm,nm)-q2l(i,k,jjc,nm))
790       continue
#  ifdef cyclic
          del2(1,k,m,jpt3)   = del2(imtm1,k,m,jpt3)
          del2(imt,k,m,jpt3) = del2(2,k,m,jpt3)
#  else
          del2(1,k,m,jpt3)   = c0
          del2(imt,k,m,jpt3) = c0
#  endif
800     continue
# endif
#endif
#endif
      return
      end
