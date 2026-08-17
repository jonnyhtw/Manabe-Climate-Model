      subroutine intrpo
      implicit real*8 (a-h,o-z)
c                                                                    
#include "param.h"
c
      parameter ( ix=imtm2/2, il=jmtm2)
      parameter ( nbc=3 )
c                                                                   
#include "ascoma.h"
c
#include "ctask.h"
c
#include "ctmngr.h"
c
#include "grdvar.h"
c
#include "slabs.h"
c
      common /ocom/ bc(imt,jmt,nbc)
c
      dimension bci(ix,jmt), toa(ix)                             
c                                                                  
c.......................................................................
c                                                                      
c  East-west interpolation from ocean grid to atmos grid.  
c   Assumes atmos gridbox is centered at Greenwich meridian,          
c       and ocean gridbox interface is at Greenwich meridian.         
c                                                                   
c                     i=    1       2       3  .....
c           Atmos grid: |   x   |   x   |   x   |......
c           Ocean grid: | o | o | o | o | o | o |.......
c                     i=  2 | 3   4   5   6   7 .....  
c                           |
c                           |<-- 0 degrees longitude
      do 10 m=1,nbc
        do 20 jj=1,jmt                                                   
          do 30 i=1,ix                                                     
            bci(i,jj) = (bc(2*i+1,jj,m) + bc(2*i,jj,m)) * p5                
  30      continue                                                        
          if (m .eq. 2) then
             do 40 i=1,ix
               toa(i) = (ta(2*i+1,1,1) + ta(2*i,1,1)) * p5
               if (toa(i) .gt. -2.0) bci(i,jj) = c0
  40         continue
          endif
  20    continue
c                                                                      
        nmonth = 1
c  Reverse j for atmos grid:                                         
        if (m .lt. nbc) then
          do 50 jj=1,il                                                  
            do 50 i=1,ix                                              
              bca(i,jj,m,nmonth) = bci(i,jmt-jj)                             
  50      continue
        endif
        if (m .eq. nbc) then
          do 60 jj=1,il                                                  
            do 60 i=1,ix                                              
              smask(i,jj) = bci(i,jmt-jj)                             
  60      continue
        endif
  10  continue
c                                                                     
      return                                                         
      end                                                           
