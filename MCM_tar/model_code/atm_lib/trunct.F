      subroutine trunct (sp,kn)
#ifdef triang
c
c  subroutine trunc truncates rhomboidal spectral fields to triangular
c    entry points: initrc - initializes /truncm/
c
#include "parm.h"
c
#include "parmjx.h"
      parameter ( mx2=2*mx, jmtr2=jx*(jx+1)-2 )
c
      common /truncm/ zerof(mx2,jx)
c
      dimension sp(mx2,jx,kx),cmp(jmtr2,kx)
c
      do 4 k=1,kn
       do 3 n=1,jx
        do 3 m=1,mx2
         sp(m,n,k)=zerof(m,n)*sp(m,n,k)
  3    continue
  4   continue
      return
c
      entry initrc
c
      do 20 n=1,jx
       do 20 m=1,mx2
        zerof(m,n) =1.0
   20 continue
      do 5 n=2,jx
       inz=mx2-2*n+3
       do 2 m=inz,mx2
  2     zerof(m,n)=0.0
  5   continue
#endif
c
      return
      end
