      subroutine findex(kxx, jjmax, kmax, jf1, jf2, imax, isf, ief)
      implicit real*8 (a-h,o-z)
#ifdef fourfil
c
c=======================================================================
c
c  findex finds and prints starting and ending indices
c         for filtering, where:
c             kxx  = field of maximum levels for the quantity
c                    being filtered
c             jjmax = number of rows to be filtered
c             kmax  = maximum number of levels to be filtered
c             jf1   = last row in the south to be filtered
c             jf2   = first row in the north to be filtered
c             imax  = maximum i index to be filtered
c             isf   = returned values of starting indices
c             ief   = returned values of ending indices
c
c=======================================================================
c
#include "param.h"
#include "pfil.h"
#include "switch.h"
c
      dimension kxx(imt,jmt)
     $        ,isf(jjmax,lsegf,kmax),ief(jjmax,lsegf,kmax)
     $        ,iis(lsegf+1),iie(lsegf+1)
c
c=======================================================================
c     begin executable code to find start and end indices
c=======================================================================
c
      jj = 0
      do 100 j = jfrst,jmtm1
        if (j.gt.jf1 .and. j.lt.jf2) goto 100
        jj = jj+1
        do 90 k = 1,kmax
           do 60 l = 1,lsegf+1
              iis(l) = 0
              iie(l) = 0
60         continue
           l = 1
           if (kxx(2,j).ge.k) then
              iis(1) = 2
           endif
           do 70 i = 2,imax-1
              if (kxx(i-1,j).lt.k .and. kxx(i,j).ge.k) then
                 iis(l) = i
              endif
              if (kxx(i,j).ge.k .and. kxx(i+1,j).lt.k) then
                 if (i.ne.iis(l) .or. (i.eq.2 .and. kxx(1,j).ge.k))then
                    iie(l) = i
                    l = l+1
                 else
                    iis(l) = 0
                 endif
              endif
70         continue
           if (kxx(imax-1,j).ge.k .and. kxx(imax,j).ge.k) then
              iie(l) = imax-1
              l = l+1
           endif
           lm = l-1
#ifdef cyclic
           if (lm.gt.1) then
              if (iis(1).eq.2 .and. iie(lm).eq.imax-1
     $                        .and. kxx(1,j).ge.k) then
                 iis(1) = iis(lm)
                 iie(1) = iie(1) + imax-2
                 iis(lm) = 0
                 iie(lm) = 0
                 lm = lm-1
              endif
           endif
#endif
           if (lm .gt. lsegf) then
              write (stdout,99) j, k
              write (stderr,99) j, k
              stop ' findex '
           endif
           do 80 l = 1,lsegf
              isf(jj,l,k) = iis(l)
              ief(jj,l,k) = iie(l)
80         continue
90      continue
100   continue
c
c     print them
c
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
      llast=lsegf
      if (llast .gt. 11) llast=11
      jj=jj+1
      do 200 j=jmtm1,jfrst,-1
        if (j.gt.jf1 .and. j.lt.jf2) go to 200
        jj = jj-1
        if (kmax .gt. 1) then
          write (stdout,9001) j
          do 190 k=1,kmax
            write (stdout,9002) k,(isf(jj,l,k),ief(jj,l,k),l=1,llast)
190       continue
        else
          write (stdout,9002) j,(isf(jj,l,1),ief(jj,l,1),l=1,llast)
        endif
200   continue
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
c
99    format (/' error => lsegf not large enough.  j=',i4,'    k=',i3)
9001  format (/' indices for row ',i3,':')
9002  format (1x,i9,3x,11(i5,i4))
#endif
      return
      end
