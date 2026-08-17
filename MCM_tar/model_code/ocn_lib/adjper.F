      subroutine adjper(per)
      implicit real*8 (a-h,o-z)
c
c   This routine is valid for the r30 ocean model only!
#ifdef ocnr30
c
c At some coastal ocean gridpoints, the surface boundary data p-e+r
c  is very large because of river runoff. Then spread the large
c  p-e+r over the neighboring pairs of gridpoints.
c  Annual mean (p-e+r) greater than 5 cm/day is considered to be
c  very large.
c
#include "param.h"
c
#ifdef amazon
# include "camazon.h"
#endif
c
      dimension per(imt,jmt)
c
c <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c   Description of the argument list for subr distr:
ccccc   call distr(per, ib,jb, i1,j1, i2,j2, iflag)
ccccc        ib,jb is i,j coords of gridpoints with large p-e+r
ccccc        i1,j1 is gridpoint of southwest corner of spread region
ccccc        i2,j2 is gridpoint of northeast corner of spread region
ccccc        iflag=(0,1)  (don't, do) print per before & after adjusting
c
c  case 1: Africa's Congo River at i=8&9, j=39
      call distr(per,  8,39,  6,37,  9,40, 1)
c
c  case 2: Africa's Niger River at i=6&7, j=43
      call distr(per,  6,43,  4,41,  7,44, 1)
c
c  case 3: Amazon River at i=168&169, j=42
      call distr(per, 168,42, iamzns,jamzns, iamzne,jamzne, 1)
c
c  case 4: Nile River at i=18&19, j=56
      call distr(per, 18,56, 16,56, 21,56, 0)
c
c  case 5: China's Yangtze River at i=68&69, j=56
      call distr(per, 68,56, 68,55, 69,57, 0)
c
      return
      end
c
      subroutine distr(per,ib,jb,i1,j1,i2,j2,iflag)
      implicit real*8 (a-h,o-z)

c
#include "param.h"
c
#include "grdvar.h"
c
#include "levind.h"
c
      dimension rkmt(imt,jmt)
      dimension per(imt,jmt)
      dimension wtj(jmt)
c
ccccc*****
      dimension perw(imt)
ccccc*****
c
c <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c    compute area weights, assumes dxt is constant
c                              and dzt is constant in the region.
      do 10 j=1,jmt
        wtj(j) = cst(j) * dyt(j)
10    continue
c
      write(stdout,9991) ib, jb
9991  format(1x,'Adjusting per boundary data at i,j=',2i3)
c
ccccc*****
      if (iflag .eq. 1) then
c   write the data before adjusting per
c   units of per are cm/sec; convert to cm/day for convenience
c    (sign of per is changed to -(p-e+r) in subr bcocn.F)
      do 109 j=j2,j1,-1
        do 108 i=i1,i2
          perw(i) = per(i,j) * 86400.
108     continue
        write(stdout,9901) i1, j, (perw(iw),iw=i1,i2)
9901    format(1x,'  For (i,j) =',i3,i3,6f12.5)
109   continue
      endif
ccccc*****
c
#ifdef testaper
c  compute area integral to check code
      perint = c0
      wtsum = c0
      do 101 jj=j1,j2
        do 101 ii=i1,i2
          if (kmt(ii,jj) .ne. 0) then
            perint = perint + per(ii,jj)*86400.*wtj(jj)
            wtsum = wtsum + wtj(jj)
          endif
101   continue
      perint = perint/wtsum
#endif
c
      pts = c0
      spread = c0
      wtsum = c0
      do 110 j=j1,j2
        do 110 i=i1,i2
          rkmt(i,j) = c0
          if (kmt(i,j) .ne. 0) then
            spread = spread + per(i,j) * wtj(j)
            pts = pts + c1
            rkmt(i,j) = c1
            wtsum = wtsum + wtj(j)
          endif
110   continue
c
      if (pts .gt. c0) then
        spread = spread / wtsum
      else
        write(stdout,9921) ib, jb
9921    format(/,' *** IN SUBR ADJPER, PTS = 0 for I,J=',2i5,' ***',/)
      endif
      do 120 j=j1,j2
        do 120 i=i1,i2
          per(i,j) = spread * rkmt(i,j)
120   continue
c
ccccc*****
      if (iflag .eq. 1) then
c   write the data after adjusting per
      perwi = per(ib,jb) * 86400.
      write(stdout,9911) perwi
9911  format(1x,' After adjusting, the values are:',f9.5)
ccc   write(stdout,9911)
ccc9911  format(1x,'  After adjusting, the values are:')
ccc   do 112 j=j2,j1,-1
ccc     do 111 i=i1,i2
ccc       perw(i) = per(i,j) * 86400.
111     continue
ccc     write(stdout,9901) i1, j, (perw(iw),iw=i1,i2)
112   continue
      endif
ccccc*****
c
#ifdef testaper
c  compute area integral to check code
      perout = c0
      wtsum = c0
      do 191 jj=j1,j2
        do 191 ii=i1,i2
          if (kmt(ii,jj) .ne. 0) then
            perout = perout + per(ii,jj)*86400.*wtj(jj)
            wtsum = wtsum + wtj(jj)
          endif
191   continue
      perout = perout/wtsum
      write(stdout,9931) perint, perout, pts, (wtj(jw),jw=j1,j2)
9931  format(1x,' Area integral before & after spreading:',2e14.6,
     $          '  pts=',f3.0,'  wtj=',4e12.5)
#endif
      return
c
#endif
      end
