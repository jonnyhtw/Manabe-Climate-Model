      subroutine bcatm(j)          
      implicit real*8 (a-h,o-z)
C                                 
#include "param.h"
c
      parameter (iy=jmtm2/2, ix=imtm2/2)                  
      parameter (nhem=2, il=nhem*iy)                  
c                                                                      
#include "ascoma.h"
c
#include "ctask.h"
c
#include "ctmngr.h"
c
#include "grdvar.h"
c
#ifdef seaice
#include "seayce.h"
#endif
c
#include "slabs.h"
c
#include "switch.h"
c
      common /ocom/ bc(imt,jmt,3)  
      common /ocom/ tbci, ybci, sump
#ifdef cpldleg2
     $             ,mts
#endif
c                                            
c.......................................................................
c                                            
#ifdef cpldleg2
      if ( first .and. j.eq.jsrow ) then
#endif
#ifdef cpldleg3
      if ( j .eq. jsrow+1 ) then
#endif
        do 10 m=1,3
          do 10 jj=1,jmt
            do 10 i=1,imt
              bc(i,jj,m) = c0  
10      continue
        tbci = c0
        ybci = c0
        sump = c0
#ifdef cpldleg2
        mts  = 0
#endif
      endif
c
c   save unfiltered sst and ice for atmos boundary data,
c    (and fm array to check the land-sea mask in subroutine atmos):
      do 20 i=1,imt
#ifdef cpldleg2
        bc(i,j,1) = bc(i,j,1) + ta(i,1,1) * fm(i,1,jc)
# ifdef seaice
        bc(i,j,2) = bc(i,j,2) + yca(i,j) * fm(i,1,jc)
# endif
#endif
#ifdef cpldleg3
        bc(i,j,1) = ta(i,1,1) * fm(i,1,jc)
# ifdef seaice
        bc(i,j,2) = yca(i,j)  * fm(i,1,jc)
# endif
#endif
        bc(i,j,3) = fm(i,1,jc)
  20  continue
#ifdef cpldleg2
      if ( j .eq. jerow ) mts = mts + 1
#endif
c
c   compute area integral for consistency check:
      do 30 i=2,imt-1
        fctr = fm(i,1,jc) * cst(j) * dxt(i) * dyt(j)
        tbci = tbci + ta(i,1,1) * fctr
#ifdef seaice
        ybci = ybci + yca(i,j) * fctr
#endif
        sump = sump + fctr
  30  continue
c
#ifdef cpldleg2
      if ( endmon .and. j.eq.jerow ) then
        do 40 jj=1,jmt
          do 40 i=1,imt
            bc(i,jj,1) = bc(i,jj,1) / float(mts)
            bc(i,jj,2) = bc(i,jj,2) / float(mts)
40      continue
#endif
#ifdef cpldleg3
      if ( j .eq. jerow ) then
#endif
        tbci = tbci / sump
        ybci = ybci / sump
c
        do 50 m=1,3                                          
          do 50 i=1,imt                                        
            bc(i, 1 ,m) = bc(i,  2  ,m)                             
            bc(i,jmt,m) = bc(i,jmt-1,m)                              
  50    continue                                                    
c                                                        
c   Interpolate from ocean to atmos grid
        call intrpo
c
        write(6,51) itt, month, tbci, ybci
#ifdef cpldleg2
     $             ,mts
#endif
  51    format(1x,'Boundary data for atmos computed at ts=',i10,
#ifdef cpldleg3
     $            '  month=',i3,'    ta=',e16.8,'  yce=',e16.8)
#endif
#ifdef cpldleg2
     $            '  month=',i3,'    ta=',e16.8,'  yce=',e16.8,
     $            '  mts=',i3)
c
        do 60 m=1,3
          do 60 jj=1,jmt
            do 60 i=1,imt
              bc(i,jj,m) = c0
60      continue
        tbci = c0
        ybci = c0
        sump = c0
        mts = 0
#endif
      endif
c                                                         
      return
      end
