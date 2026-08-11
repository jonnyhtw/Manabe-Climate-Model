      subroutine topog
      implicit real*8 (a-h,o-z)
c
c-----------------------------------------------------------------------
c     construct an "idealized" world ... piece by piece
c     note: this topography will map into arbitrary model resolution
c-----------------------------------------------------------------------
c
#include "param.h"
#include "coord.h"
#include "grdvar.h"
#include "levind.h"
#if defined ocnr15 || defined ocnr30
# include "iounit.h"
c
      dimension rkmt (imt)
#endif
c
#ifdef writekmt
      dimension ind(imt)
c
#endif
c     define the model topography on the "t" grid using "kmt"
c     set maximum level everywhere except on boundary
c
      do 100 j=2,jmtm1
        do 90 i=2,imtm1
          kmt(i,j) = km
90      continue
100   continue
c
      do 200 j=1,jmt
        kmt(imt,j) = 0
        kmt(1,j) = 0
200   continue
      do 300 i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1) = 0
300   continue
#if defined ocnr15 
c
c      set t,s levels to values for ocnr15 runs
c
      open(iokmt,file='cpldfkmt_r15')
#endif
#if defined ocnr30 
# if defined km18   
c
c      set t,s levels to values for 18-levels
c
      open(iokmt,file='cpldfkmtk18_r30')
# endif
#endif
c
#if defined ocnr15 || defined ocnr30
      do 315 j=1,jmt
        read(iokmt,15) jj, (rkmt(i),i=1,imt)
        if (jj .ne. j) then
          write (stdout,99) j, jj
          write (stderr,99) j, jj
          stop ' topog  '
        endif
        do 314 i=1,imt
          kmt(i,j) = rkmt(i)
 314    continue
 315  continue
      close (iokmt)
# ifdef testnorth
c  ensure that northernmost and southermost rows are land. 
      do 339 i=1,imt
        kmt(i,jmt) = 0
        kmt(i,1) = 0
339   continue
# endif
c
#else
c
c     antarctica
c
      call setkmp (-90.0, 0.0, c360, -80.0, 0.0, c360, 0)
      call setkmp (-80.0, c360-25.0, c360, -70.0, c360, c360, 0)
      call setkmp (-80.0, 0.0, c360, -70.0, 0.0, 170.0, 0)
      call setkmp (-80.0, c360-135.0, c360-60.0, -68.0, c360-75.0,
     $              c360-60.0, 0)
      call setkmp (-70.0, 0.0, 155.0, -67.0, 50.0, 145.0, 0)
c
c     australia
c
      call setkmp (-35.0, 116.0, 120.0, -31.0, 114.0, 130.0, 0)
      call setkmp (-38.0, 140.0, 151.0, -31.0, 130.0, 151.0, 0)
      call setkmp (-31.0, 115.0, 153.0, -20.0, 113.0, 149.0, 0)
      call setkmp (-20.0, 113.0, 149.0, -11.0, 131.0, 143.0, 0)
c
c     south america
c
      call setkmp (-50.0, c360-74.0, c360-68.0, -40.0, c360-73.0,
     $              c360-62.0, 0)
      call setkmp (-40.0, c360-73.0, c360-62.0, -20.0, c360-70.0,
     $              c360-40.0, 0)
      call setkmp (-20.0, c360-70.0, c360-40.0, -16.0, c360-81.0,
     $              c360-35.0, 0)
      call setkmp (-16.0, c360-81.0, c360-35.0, 0.0, c360-80.0,
     $             c360-50.0, 0)
      call setkmp (0.0, c360-80.0, c360-50.0, 11.0, c360-75.0,
     $             c360-60.0, 0)
c
c     central america
c
      call setkmp (6.0, c360-78.0, c360-75.0, 20.0, c360-105.0,
     $             c360-97.0, 0)
      call setkmp (20.0, c360-105.0, c360-97.0, 30.0, c360-115.0,
     $             c360-94.0, 0)
c
c     north america
c
      call setkmp (25.0, c360-82.0, c360-80.0, 30.0, c360-85.0,
     $             c360-81.0, 0)
      call setkmp (30.0, c360-115.0, c360-80.0, 40.0, c360-124.0,
     $             c360-74.0, 0)
      call setkmp (40.0, c360-124.0, c360-74.0, 50.0, c360-124.0,
     $             c360-57.0, 0)
      call setkmp (50.0, c360-124.0, c360-57.0, 60.0, c360-140.0,
     $             c360-64.0, 0)
      call setkmp (60.0, c360-165.0, c360-64.0, 65.0, c360-140.0,
     $             c360-64.0, 0)
      call setkmp (65.0, c360-140.0, c360-64.0, 70.0, c360-162.0,
     $             c360-72.0, 0)
      call setkmp (70.0, c360-162.0, c360-140.0, 72.0, c360-157.0,
     $             c360-157.0, 0)
      call setkmp (70.0, c360-130.0, c360-70.0, 75.0, c360-120.0,
     $             c360-80.0, 0)
c
c     greenland
c
      call setkmp (60.0, c360-45.0, c360-45.0, 75.0, c360-58.0,
     $             c360-19.0, 0)
c
c     africa
c
      call setkmp (-35.0, 19.0, 28.0, 6.0, 8.0, 50.0, 0)
      call setkmp (6.0, 0.0, 50.0, 18.0, 0.0, 56.0, 0)
      call setkmp (18.0, 0.0, 56.0, 26.0, 0.0, 59.0, 0)
      call setkmp (6.0, c360-10.0, c360, 18.0, c360-18.0, c360, 0)
      call setkmp (18.0, c360-18.0, c360, 26.0, c360-15.0, c360, 0)
c
c     northern africa &  europe & asia
c
      call setkmp (26.0, c360-15.0, c360, 40.0, c360-7.0, c360, 0)
      call setkmp (40.0, c360-7.0, c360, 50.0, c360, c360, 0)
c
      call setkmp (26.0, 0.0, 126.0, 40.0, 0.0, 122.0, 0)
      call setkmp (40.0, 0.0, 130.0, 50.0, 0.0, 140.0, 0)
      call setkmp (50.0, 0.0, 140.0, 60.0, 8.0, 140.0, 0)
      call setkmp (60.0, 8.0, 163.0, 65.0, 13.0, 180.0, 0)
      call setkmp (65.0, 13.0, 188.0, 70.0, 20.0, 180.0, 0)
      call setkmp (70.0, 70.0, 180.0, 75.0, 90.0, 100.0, 0)
c
c     add an "idealized" ridge in the atlantic & pacific
c
      level = indp (3500.e2, zw, km)
c
      call setkmp (-20.0, c360-20.0, c360-10.0, 30.0, c360-45.0,
     $              c360-35.0, level)
      call setkmp (30.0, c360-45.0, c360-35.0, 60.0, c360-20.0,
     $              c360-30.0, level)
c
      call setkmp (-60.0,c360-100.0, c360-130.0, 40.0, c360-160.0,
     $              180.0, level)
      level = indp (2000.e2, zw, km)
      call setkmp (-50.0, c360-120.0, c360-120.0, 30.0, 190.0, 190.0,
     $             level)
#endif
#ifdef cyclic
c
c     set cyclic boundary conditions
c
      do 400 j=1,jmt
        kmt(1,j)   = kmt(imtm1,j)
        kmt(imt,j) = kmt(2,j)
400   continue
#endif
c
#ifdef writekmt
      do 320 i=1,imt
        ind(i) = i
320   continue
      is = 1
350   continue
      ie = ie + 40
      if (ie .gt. imt) ie = imt
      write(6,9011) (ind(i),i=is,ie)
9011  format(/,3x,40i3)
      do 330 j=1,jmt
        jj=jmt+1-j
        write(stdout,9030) jj,(kmt(i,jj),i=is,ie)
9030    format(1x,i2,40i3)
330     continue
      is = is + 40
      if (is .lt. imt) go to 350
#endif
c
      return
#if defined ocnr15 || defined ocnr30
  15  format(i3,/,(10f5.1))
  99  format(/' error => problem when reading in kmt'
     $       /' expected j= ',i3,'  got ',i5)
#endif
c
      end
