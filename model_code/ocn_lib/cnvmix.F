      subroutine cnvmix (j)
      implicit real*8 (a-h,o-z)
#ifdef constvmix
#ifdef multitasking
cfpp$ noconcur r
#endif
c=====================================================================
c
c     set up coefficients for constant vertical diffusion
c
c     (vertical diffusion of tracers and momentum is constant in time
c     but is allowed to vary in grid space in a predetermined way
c     through the use of scaling factors)
# ifdef blvmix
c
c     vertical tracer diffusion coefficients are set as a function of
c     depth, following the method of Bryan and Lewis (JGR, 1979)
# endif
c
c=====================================================================
c
#include "param.h"
#include "cvmix.h"
#ifdef implicitvmix
# include "levind.h"
# include "index.h"
# include "slabs.h"
# include "temp.h"
#endif
c
      do 100 k=1,km
        do 90 i=1,imt
          vdc(i,k) = fkph*vdcsci(i)*vdcscj(j)*vdcsck(k)
          vvc(i,k) = fkpm*vvcsci(i)*vvcscj(j)*vvcsck(k)
90      continue
100   continue
c
#ifdef implicitvmix
c
c         ks=1: compare lev. 1 to 2; 3 to 4; etc.
c         ks=2: compare lev. 2 to 3; 4 to 5; etc.
c
          do 200 ks=1,2
c
c           find density for entire slab for stability determination
c
            call statec (t(1,1,jc,nm,1), t(1,1,jc,nm,2), tempik(1,1,1)
     $,                 ks)
# ifdef skipland
            do 190 lstrip=1,lseg
              istart = ist(j,lstrip)
              iend   = iet(j,lstrip)
              if (istart .ne. 0) then
# endif
c
c               set "heavy water" in land to insure stability there
c
                dense = 1.e15
                do 150 i=istart,iend
                  k = kmt(i,j) + 1
                  if (k .le. km) then
                    tempik(i,k,1) = dense
                  endif
150            continue
c
c               if unstable, set mixing to the limits
c
                do 180 m=1,nt
                  do 170 k=ks,kmm1,2
                    do 160 i=istart,iend
                      if (tempik(i,k,1) .gt. tempik(i,k+1,1)) then
                        vdc(i,k) = vdclim
                        vvc(i,k) = vvclim
                      endif
160                continue
170              continue
180            continue
# ifdef skipland
              endif
190        continue
# endif
200       continue
#endif
#endif
      return
      end
