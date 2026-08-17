      subroutine setlsm (lsts, lsuv)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c     using the horizontal region masks "mskhr", this routine creates
c     landsea masks of zero for land and 1 for sea for both the
c      ts & uv grids.
c=======================================================================
c
#include "param.h"
#include "cregin.h"
c
c     dimension land/masks for ts & uv grids
c
      dimension lsts(imt,jmt), lsuv(imt,jmt)
c=======================================================================
c
c     compute "lsts" from horizontal region masks "mskhr"
c
      do 100 j=1,jmt
        do 90 i=1,imt
          if (mskhr(i,j) .ne. 0) then
            lsts(i,j) = 1
          else
            lsts(i,j) = 0
          endif
90      continue
100   continue
c
c     compute the uv grids "lsuv" from the ts grids "lsts"
c     (avoid using 1st and last columns, instead set
c     cyclic conditions)
c
      do 200 jrow=1,jmt-2
        do 190 i=2,imt-2
          lsuv(i,jrow) = lsts(i,jrow)*lsts(i+1,jrow)
     $                  *lsts(i,jrow+1)*lsts(i+1,jrow+1)
190     continue
        lsuv(imt-1,jrow) = lsts(imt-1,jrow)*lsts(imtm1,jrow+1)
     $                    *lsts(2,jrow)*lsts(2,jrow+1)
#ifdef cyclic
        lsuv(imt  ,jrow) = lsuv(2    ,jrow)
        lsuv(1    ,jrow) = lsuv(imt-1,jrow)
#endif
200   continue
c
      do 300 i=2,imt-2
        lsuv(i,jmt-1) = lsts(i,jmtm2)*lsts(i+1,jmtm2)
300   continue
#ifdef cyclic
      lsuv(imt  ,jmt-1) = lsuv(2    ,jmt-1)
      lsuv(1    ,jmt-1) = lsuv(imt-1,jmt-1)
#endif
c
      do 400 jrow=1,jmt
        jj = jmt + 1 - jrow
        write(stdout,9401) jj, (lsts(ii,jj),ii=1,imt)
400   continue
c
      do 500 jrow=1,jmt
        jj = jmt + 1 - jrow
        write(stdout,9401) jj, (lsuv(ii,jj),ii=1,imt)
500   continue
c
      return
9401  format(1x,i3,1x,100(i1))
      end


