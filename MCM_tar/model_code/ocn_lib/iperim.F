      subroutine iperim (kmz)
      implicit real*8 (a-h,o-z)
#ifdef islands
c
c=======================================================================
c
c     find island perimeter coordinates & analyze for errors in land
c     mass distribution
c
c     inputs:
c
c     kmz   = scratch area
c     kmt   = number of levels on the "t" grid
c     kmu   = number of levels on the "u" grid
c
c     output:
c
c     all island variables in "index.h"
c
c=======================================================================
c
      character*1 dir, blank
      logical prblm, fatal
#include "param.h"
#include "coord.h"
#include "index.h"
#include "levind.h"
#if defined cpldleg3 || defined cpldleg2
# include "switch.h"
#endif
      dimension kmz(imt,jmt)
      dimension ncol(imt)
c
      data blank /' '/
c
c-----------------------------------------------------------------------
c     statement functions
c-----------------------------------------------------------------------
c
      iwest(i) = mod(i-2+imt,imt) + 1
      ieast(i) = mod(i,imt) + 1
c
c-----------------------------------------------------------------------
c     save the original "kmt" field (to be restored later)
c     and initialize "chisl" field
c-----------------------------------------------------------------------
c
      do 100 j=1,jmt
        do 90 i=1,imt
          kmz(i,j) = kmt(i,j)
90      continue
100   continue
c
c-----------------------------------------------------------------------
c     search for & knockout isolated bays ... "t" grid boxes at the
c     surface which can't be influenced by advection
c-----------------------------------------------------------------------
c
      do 400 j=2,jmtm1
        do 390 i=2,imtm1
          if (kmt(i,j) .ne. 0) then
            if (kmu(i,j) .eq. 0 .and. kmu(i-1,j) .eq. 0 .and.
     $        kmu(i,j-1) .eq. 0 .and. kmu(i-1,j-1) .eq. 0) then
              kmt(i,j) = 0
            endif
          endif
390     continue
400   continue
c
c-----------------------------------------------------------------------
c     find the coordinates of the perimeter of each island
c-----------------------------------------------------------------------
c
      fatal = .false.
      do 600 n=1,nisle
c
c-----------------------------------------------------------------------
c       convert island latitude & longitude "seed" values to nearest
c       model points incices.
c-----------------------------------------------------------------------
c
        i = indp (alonis(n), xt, imt)
        j = indp (alatis(n), yt, jmt)
        if (kmt(i,j) .ne. 0) then
          write (stdout,'(a17,i2,a18,5x,a5,f8.2,a6,f8.2)')
     $    ' => island seed #',n,' is an ocean point'
     $,   ' lon=',alonis(n),', lat=',alatis(n)
          stop '=>iperim'
        endif
c
c-----------------------------------------------------------------------
c       determine a first perimeter point
c-----------------------------------------------------------------------
c
        if (j .le. jmt/2) then
          jend = jmt
          jd   = 1
        else
          jend = 1
          jd   = -1
        endif
        do 500 jj=j,jend,jd
          if (kmt(i,jj) .ne. 0) then
            j = jj
            go to 510
          endif
500     continue
510     continue
        i0  = i
        j0  = j
        dir = blank
c
c-----------------------------------------------------------------------
c       scan for the rest of the perimeter
c-----------------------------------------------------------------------
c
        if (n .eq. 1) then
          iofs(n)  = 0
        else
          iofs(n)  = nippts(n-1) + iofs(n-1)
        endif
        lenp       = iofs(n) + 1
        iperm(lenp) = i
        jperm(lenp) = j
c
        do 580 m=1,maxipp
          iw  = iwest(i)
          ie  = ieast(i)
          js  = max(1,j-1)
          jn  = min(jmt,j+1)
          if (iw .eq. 1) then
            iw  = imt-1
          endif
          if (ie .eq. imt) then
            ie  = 2
          endif
c
c-----------------------------------------------------------------------
c         analyze the geometry to verify that the island and any other
c         land mass are separated by at least 2 ocean "t" grid boxes
c-----------------------------------------------------------------------
c
          prblm = .false.
c
          if ((kmt(ie,jn) .eq. 0 .and. kmt(iw,js) .eq. 0) .or.
     $        (kmt(ie,js) .eq. 0 .and. kmt(iw,jn) .eq. 0)) then
            if (kmt(iw,j) .ne. 0 .and. kmt(ie,j) .ne. 0 .and.
     $          kmt(i,jn) .ne. 0 .and. kmt(i,js) .ne. 0 ) prblm = .true.
            if (kmt(iw,j) .ne. 0 .and. kmu(iw,j) .eq. 0 .and.
     $          kmu(iw,js) .eq. 0) prblm = .true.
            if (kmt(ie,j) .ne. 0 .and. kmu(i,j) .eq. 0 .and.
     $          kmu(i,js) .eq. 0) prblm = .true.
            if (kmt(i,jn) .ne. 0 .and. kmu(i,j) .eq. 0 .and.
     $          kmu(iw,j) .eq. 0) prblm = .true.
            if (kmt(i,js) .ne. 0 .and. kmu(i,js) .eq. 0 .and.
     $          kmu(iw,js) .eq. 0) prblm = .true.
          endif
          if (prblm) then
            write (stdout,9200) n, i, j, xt(i), yt(j)
            fatal = .true.
          endif
c
c-----------------------------------------------------------------------
c         look west (only if not comming from west)
c-----------------------------------------------------------------------
c
          k2 = kmt(iw,jn)  * kmt(iw,js)
          k3 = kmt(i,jn)   * kmt(i,js)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(iw,j) .ne. 0 .and. k .eq. 0 .and. dir .ne. 'e') then
            if (lenp .lt. maxipp) then
              if (iw .ne. i0 .or. j .ne. j0) then
                dir          = 'w'
                lenp         = lenp + 1
                i            = iw
                iperm(lenp)  = iw
                jperm(lenp)  = j
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>west'
            endif
          endif
c
c-----------------------------------------------------------------------
c         look north (only if not comming from north)
c-----------------------------------------------------------------------
c
          k2 = kmt(ie,jn)  * kmt(iw,jn)
          k3 = kmt(ie,j)   * kmt(iw,j)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(i,jn) .ne. 0 .and. k .eq. 0 .and. jn .ne. j
     $        .and. dir .ne. 's') then
            if (lenp .lt. maxipp) then
              if (i .ne. i0 .or. jn .ne. j0) then
                dir          = 'n'
                lenp         = lenp + 1
                j            = jn
                iperm(lenp)  = i
                jperm(lenp)  = jn
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>north'
            endif
          endif
c
c-----------------------------------------------------------------------
c         look east (only if not comming from east)
c-----------------------------------------------------------------------
c
          k2 = kmt(ie,jn)  * kmt(ie,js)
          k3 = kmt(i,jn)   * kmt(i,js)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(ie,j) .ne. 0 .and. k .eq. 0 .and. dir .ne. 'w') then
            if (lenp .lt. maxipp) then
              if (ie .ne. i0 .or. j .ne. j0) then
                dir          = 'e'
                lenp         = lenp + 1
                i            = ie
                iperm(lenp)  = ie
                jperm(lenp)  = j
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>east'
            endif
          endif
c
c-----------------------------------------------------------------------
c         look south (only if not comming from south)
c-----------------------------------------------------------------------
c
          k2 = kmt(ie,js)  * kmt(iw,js)
          k3 = kmt(ie,j)   * kmt(iw,j)
          if (k2 .eq. 0 .or. k3 .eq. 0) then
            k = 0
          else
            k = 1
          endif
          if (kmt(i,js) .ne. 0 .and. k .eq. 0 .and. js .ne. j
     $        .and. dir .ne. 'n') then
            if (lenp .lt. maxipp) then
              if (i .ne. i0 .or. js .ne. j0) then
                dir          = 's'
                lenp         = lenp + 1
                j            = js
                iperm(lenp)  = i
                jperm(lenp)  = js
                go to 580
              else
                go to 590
              endif
            else
              write (stdout,'(a23)') '=> maxipp limit reached'
              stop '=>south'
            endif
          endif
580     continue
590     continue
        nippts(n) = lenp - iofs(n)
# if defined cpldleg3 || defined cpldleg2
        if (pintro) then
# endif
        write (stdout,9000) n, nippts(n)
        ls = iofs(n) + 1
        le = lenp
        write (stdout,9100) (iperm(l),jperm(l), l=ls,le)
# if defined cpldleg3 || defined cpldleg2
        endif
# endif
600   continue
c
c-----------------------------------------------------------------------
c     restore the original "kmt" field
c-----------------------------------------------------------------------
c
      do 700 j=1,jmt
        do 690 i=1,imt
          kmt(i,j) = kmz(i,j)
690     continue
700   continue
c
c-----------------------------------------------------------------------
c     if a fatal error was encountered print out error message and stop
c     otherwise load perimeter field into "kmz" and print it out
c-----------------------------------------------------------------------
c
      if (fatal) then
c
        write (stdout,*)
     $   ' Fatal error(s) in "iperim" land mass analysis'
        write (stdout,*)
     $    ' The "kmt" field must be fixed up before continuing ...'
        stop '=>iperim'
c
      else
c
c-----------------------------------------------------------------------
c     initialize scratch "kmz" field
c-----------------------------------------------------------------------
c
        do 800 j=1,jmt
          do 790 i=1,imt
            kmz(i,j) = 0
790       continue
800     continue
c
c-----------------------------------------------------------------------
c       identify perimeter points with island number
c-----------------------------------------------------------------------
c
        do 900 ni=1,nisle
          ls = iofs(ni) + 1
          le = ls + nippts(ni) - 1
          do 890 m=ls,le
            kmz(iperm(m),jperm(m)) = ni
890       continue
900     continue
c
c-----------------------------------------------------------------------
c       print out island perimeter field
c-----------------------------------------------------------------------
c
        linemx = 100
        linel  = 105
        line   = linel - 5
        if (line .gt. linemx) line = linemx
        nwr = (imt/line) + 1
c
        do 950 i=1,imt
          ncol(i) = mod(i,10)
950     continue
c
# if defined cpldleg3 || defined cpldleg2
       if (pintro) then
# endif
       write (stdout,9300)
        do 1000 n=1,nwr
          ia = 1 + (line*(n-1))
          ib = ia + line - 1
          if (ib .gt. imt) ib = imt
          write(stdout,9400) (ncol(i),i=ia,ib)
          do 990 jj=1,jmt
            jjj = jmt - jj + 1
            write(stdout,9500)jjj,(kmz(i,jjj),i=ia,ib)
990       continue
1000    continue
# if defined cpldleg3 || defined cpldleg2
       endif
# endif
c
      endif
c
9000  format (/' Island # ',i2,' perimeter length is ',i4, ' points.'
     $, ' (i,j) coordinates are:'/)
9100  format (10(1x,'(',i4,',',i4,') '))
9200  format (///1x, ' ==> Error on island # ',i2,' at (i,j) = (',i4,','
     $,  i4,') ... (lon,lat) = (',f7.2,',',f7.2,')'/
     $, 1x,'    There must be at least two ocean points separating '
     $,    'an island and another land mass.'/)
9300  format (/,' Island perimeter locations '/)
9400  format('  i=>',100(i1))
9500  format(1x,i3,1x,100(i1))
c
#endif
      return
c
      end
