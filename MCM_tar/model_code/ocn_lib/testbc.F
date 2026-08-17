      subroutine testbc (bc, ls, flag)
      implicit real*8 (a-h,o-z)
#ifdef bctest
c
c=======================================================================
c
c     compute zonal averages of surface boundary conditions 
c     for a quick test
c
c=======================================================================
c
#include "param.h"
c
c   variables:
c     ls   =  land/masks used in testing of surface boundary conditions
c             that were read in in "getvc"
c     bc   =  12 monthly boundary condition fields
c     flag =  flag for land values
c
      dimension bc(imt,jmt,12)
      dimension ls(imt,jmt)
c
c=======================================================================
c
      do 100 mnth = 1,12
        write(stdout,9101) mnth
        do 90 jrow=2,jmtm1
          sum = c0
          count = c0
          do 80 i=2,imtm1
            if (ls(i,jrow) .ne. 0) then
              sum = sum + bc(i,jrow,mnth)
              count = count + c1
            endif
80        continue
          if (count .gt. c0) sum = sum/count
          write(stdout,9102) jrow, sum, count
90      continue
100   continue
c
#endif
      return
9101  format(/' Check of surface boundary conditions for month ',i2)
9102  format(' jrow=',i3,5x,'average= ',e12.6,5x,'over ',f6.0,' pts')
      end
