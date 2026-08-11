      subroutine nor1st (nunit)
      implicit real*8 (a-h,o-z)
#ifdef testnorth
c
c=======================================================================
c
c     Subroutine nor1st is i/o routine for the masks used in northward 
c     stuff
c
c=======================================================================
c
#include "param.h"
#include "cdiag.h"
#include "levind.h"
c
c  nunit = unit to be written to or read from
c
c
      character *1 chrow(imt), chnum(0:9)
      data (chnum(n),n=0,9) /'0','1','2','3','4','5','6','7','8','9'/
c
      dimension ncol(imt)
c
c-----------------------------------------------------------------------
c       read in (write out) horizontal masks
c       set linel to length of desired formatted printout line
c-----------------------------------------------------------------------
c
        linemx = 100
        linel  = 105
        line   = linel - 5
        if (line .gt. linemx) line = linemx
        nwr = (imt/line) + 1
c
        do 100 i=1,imt
          ncol(i) = mod(i,10)
100     continue
c
        kntnor = 0
        read (nunit,9099)
        read (nunit,9099)
        do 1000 n=1,nwr
          ia = 1 + (line*(n-1))
          ib = ia + line - 1
          if (ib .gt. imt) ib = imt
          write(stdout,9021) (ncol(i),i=ia,ib)
          read (nunit,9099)
          do 990 jj=1,jmt
            jjj = jmt - jj + 1
            read (nunit,9022) jr , (normsk(i,jjj),i=ia,ib)
            if (jr .ne. jjj) then
              write (stdout,999) nunit, jjj, jr
              write (stderr,999) nunit, jjj, jr
              stop 'normsk  '
            endif
            write(stdout,9022)jjj,(normsk(i,jjj),i=ia,ib)
            do 980 i=ia,ib
              if (i.NE.1 .AND. i.NE.imt) then
                if (kmt(i,jjj).GT.0 .AND. normsk(i,jjj).EQ.0) then
                  write(stdout,9116) i,jjj,kmt(i,jjj),normsk(i,jjj)
                elseif (kmt(i,jjj).EQ.0 .AND. normsk(i,jjj).NE.0) then
                  write(stdout,9116) i,jjj,kmt(i,jjj),normsk(i,jjj)
                else
                  kntnor = kntnor + 1
                endif
              endif
980         continue
990       continue
1000    continue
c
        write(stdout,9123) kntnor
c
#endif
      return
c
999   format(/' error => bad j-row when reading normsk from unit ',
     $       i3,/'   expected',i4,'    read in',i4)
9021  format('  i=>',100(i1))
9022  format(1x,i3,1x,100(i1))
9116  format('CHECK normsk at i,j=',2i4,'  kmt=',i3,'  normsk=',i3)
9123    format(' A total of ',i6,' normsk points checked out OK')     
9099  format(1x)
      end
