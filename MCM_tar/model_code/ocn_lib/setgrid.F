      subroutine cmesh (dx, apb, amb, wid, n, ndir)
      implicit real*8 (a-h,o-z)
c
c
c=======================================================================
c
c                   c-o-s-i-n-e   m-e-s-h
c
c     calculate variable width spacings "dx" in domain of width "wid"
c     which is imbedded between fine & coarse resolution domains on
c     either side. the variation of "dx" between fine & coarse
c     resolutions follows a cosine.
c
c     inputs:
c
c     apb=a+b = coarse resolution
c     amb=a-b = fine resolution
c     wid = width between coarse & fine resolution
c     ndir = (1,-1) = spacing (inc,dec) with increasing index
c
c     outputs:
c
c     dx = variable width spacings
c     n  = number of variable spaces calculated to span wid
c        = 2 * wid / ( coarse + fine resoultion)
c
c=======================================================================
c
#include "param.h"
#include "scalar.h"
#include "switch.h"
c
      dimension dx(1000)
c
      a     = p5*(apb+amb)
      b     = p5*(apb-amb)
      c     = (wid/a) + 0.001
      n     = c
      rn    = c1/float(n)
c
      do 100 ii=1,n
        if (ndir .eq. 1)  then
          i = n+1-ii
        elseif (ndir .eq. -1) then
          i = ii
        else
          write (stdout, 9009) ndir
          write (stderr, 9009) ndir
          stop ' cmesh1'
        endif
        dx(i) = a + b*cos(pi*rn*(ii-p5))
100   continue
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
      write (stdout,9000) apb, amb, wid, c
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
c
      frac = c-int(c)
      if (frac .gt. 0.01*a) then
        write (stdout, 9010)
        write (stderr, 9010)
        stop ' cmesh  '
      endif
      return
9000  format (/1x,'max res=',e11.4,'  min res=',e11.4,'  width=',e11.4,
     $'  number of grid boxes=',e14.7)
9009  format (/' error => "ndir" =',i8,' but should be 1 or -1 '/)
9010  format (/' error => non integral number of boxes in "wid"'/)
      end
      subroutine grids
      implicit real*8 (a-h,o-z)
c
c
c=======================================================================
c
c     calculate all grid related items
c
c=======================================================================
c
c
#include "param.h"
#ifndef implicitvmix
#include "accel.h"
#endif
#include "coord.h"
#include "grdvar.h"
#include "scalar.h"
#include "switch.h"
#ifdef worldocn
# ifdef ocnr15
c
c  e-w grid orientation
c
c                    -3.75        3.75
c  atmosphere   -7.5    |    0.0    |    7.5
c             ===48=====|=====1=====|=====2=
c                /|\    |    /|\    |    /|\
c               / | \   |   / | \   |   / |
c              /  |  \  |  /  |  \  |  /  |
c             _|__|__|__|__|__|__|__|__|__|_
c       ocean  t  v  t  v  t  v  t  v  t  v
c             96 96  1  1  2  2  3  3  4  4
c                   97 97 98 98
c
c-----------------------------------------------------------------------
# endif
# ifdef ocnr30
c
c  e-w grid orientation
c
c                       -1.875           1.875
c  atmosphere    -3.75     |      0.0      |      3.75
c             ====96=======|=======1=======|=======2=
c                / | \     |     / | \     |     / | \
c               /  |  \    |    /  |  \    |    /  |
c              /   |   \   |   /   |   \   |   /   |
c             _| __| __| __| __| __| __| __| __| __|_
c       ocean  t   v   t   v   t   v   t   v   t   v
c            192 192   1   1   2   2   3   3   4   4
c                    193 193 194 194
c
c-----------------------------------------------------------------------
# endif
#endif
c
           xwid(1) = 360.0 + (360.0/imtm2)*c2
           xmin(1) = 360.0/imtm2
           xmax(1) = xmin(1)
           idir(1) = 1
           stlon = c0 - (360.0/imtm2)*c2
c
c
c-----------------------------------------------------------------------
c     set width & height of grid boxes (in degrees)
c-----------------------------------------------------------------------
c
      call mesh (dxt, dyt)
c
c-----------------------------------------------------------------------
c     set latitudes & longitudes of "t" & "u,v" grid points (in degrees)
c-----------------------------------------------------------------------
c
      xt(1) = stlon + p5*dxt(1)
      xu(1) = stlon + dxt(1)
      do 100 i=2,imt
        xt(i) = xt(i-1) +  p5*(dxt(i) + dxt(i-1))
        xu(i) = xu(i-1) +  dxt(i)
100   continue
c
#ifdef worldocn
c     set "t" & "u,v" latitudes for ocnr15 or ocnr30 world ocean cases
c
      jmtr2 = jmt/2
      yu(jmtr2) = c0
      do 150 j=1,jmtr2-1
        yu(jmtr2+j) = yu(jmtr2+j-1) +  dyt(jmtr2+j)
        yu(jmtr2-j) = -yu(jmtr2+j)
150   continue
      yu(jmt) = yu(jmt-1) +  dyt(jmt)
c
      yt(1) = yu(1) - p5*dyt(1)
      do 200 j=2,jmt
        yt(j) = yt(j-1) +  p5*(dyt(j) + dyt(j-1))
200   continue
#else
      yt(1) = stlat + p5*dyt(1)
      yu(1) = stlat + dyt(1)
      do 200 j=2,jmt
        yt(j) = yt(j-1) +  p5*(dyt(j) + dyt(j-1))
        yu(j) = yu(j-1) +  dyt(j)
200   continue
#endif
c
c     convert to centimeters
c
      do 300 j=2,jmtm1
        dyt(j) = dyt(j)*radius/radian
300   continue
      dyt(1)   = dyt(2)
      dyt(jmt) = dyt(jmtm1)
c
      do 400 i=2,imtm1
        dxt(i) = dxt(i)*radius/radian
400   continue
      dxt(1)   = dxt(2)
      dxt(imt) = dxt(imtm1)
#ifdef cyclic
c
c     set cyclic conditions on dxt
c
      dxt(1)=dxt(imtm1)
      dxt(imt)=dxt(2)
#endif
c
c-----------------------------------------------------------------------
c     compute auxiliary arrays based upon the spacing specified above
c-----------------------------------------------------------------------
c
      do 500 k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
500   continue
c
      dzw(0) = p5*dzt(1)
      zw(1)  = dzt(1)
      do 600 k=2,km
        dzw(k-1) = p5*(dzt(k-1)+dzt(k))
        zw(k)    = zw(k-1)+dzt(k)
600   continue
c
      dzw(km)   = p5*dzt(km)
      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)
      zt(1)     = dzw(0)
c
      do 700 k=1,km
        dzwr(k-1)  = c1/dzw(k-1)
        dzw2r(k-1) = p5/dzw(k-1)
        zt(k+1)  = zt(k) + dzw(k)
700   continue
c
      do 800 k=1,km
#if defined implicitvmix || defined isopycmix
        dztur(k) = c1/(dzw(k-1)*dzt(k))
        dztlr(k) = c1/(dzw(k)*dzt(k))
#endif
#ifdef tcvmix
        dzwur(k) = c1/(dzt(k)*dzw(k))
        if (k .lt. km) dzwlr(k) = c1/(dzt(k+1)*dzw(k))
#endif
        dztr(k)  = c1/dzt(k)
800   continue
c
      phi(1)   = yu(1)/radian
      phit(1)  = phi(1) - p5*dyt(1)/radius
      sumdy    = phi(1)
      dyu(jmt) = dyt(jmt)
c
      do 900 j=1,jmt
        if (j .ne. jmt) dyu(j) = p5*(dyt(j)+dyt(j+1))
        dytr(j)  = c1/dyt(j)
        dyt2r(j) = p5/dyt(j)
        dyt4r(j) = p25/dyt(j)
        dyur(j)  = c1/dyu(j)
        dyu2r(j) = p5/dyu(j)
        dyu4r(j) = p25/dyu(j)
        if (j .ne. jmt) sumdy = sumdy + dyt(j+1)/radius
        if (j .ne. jmt) phi(j+1) = sumdy
        if (j .ne. 1) phit(j) = p5*(phi(j-1)+phi(j))
        cst(j)  = cos(phit(j))
        csu(j)  = cos(phi (j))
        sine(j) = sin(phi(j))
        cstr(j) = c1/cst(j)
        csur(j) = c1/csu(j)
        tng(j)  = sine(j)/csu(j)
#ifndef nohilats
        tanra(j)= tng(j)/radius
#endif
900   continue
c
      dxu(imt)=dxt(imt)
#ifdef cyclic
      dxu(imt)=p5*(dxt(2)+dxt(3))
#endif
c
      do 1000 i=1,imt
        if (i .ne. imt) dxu(i) = p5*(dxt(i)+dxt(i+1))
        dxtr(i)  = c1/dxt(i)
        dxt2r(i) = p5/dxt(i)
        dxt4r(i) = p25/dxt(i)
        dxur(i)  = c1/dxu(i)
        dxu2r(i) = p5/dxu(i)
        dxu4r(i) = p25/dxu(i)
1000  continue
c
#ifndef implicitvmix
c     for convection code, compute values needed to include effects
c     of tracer timestep acceleration on effective layer thicknesses
c
      do 1100 k=1,km
        dztxcl(k) = dzt(k)/dtxcel(k)
1100  continue
c
      do 1200 k=1,kmm1
        dzwxcl(k) = c1/(dztxcl(k)+dztxcl(k+1))
1200  continue
      dzwxcl(km) = c0
c
#endif
c-----------------------------------------------------------------------
c     print grid related arrays
c-----------------------------------------------------------------------
c
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
      write (stdout,9101)
      write (stdout,9001) dzt
      write (stdout,9102)
      write (stdout,9001) dzw
      write (stdout,9103)
      write (stdout,9001) zw
      write (stdout,9104)
      write (stdout,9001) zt
      write (stdout,9105)
      write (stdout,9001) yt
      write (stdout,9106)
      write (stdout,9001) yu
      write (stdout,9107)
      write (stdout,9001) xt
      write (stdout,9108)
      write (stdout,9001) xu
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
c
       return
c
9001  format (1x,10g13.5)
9101  format (/,' "t" grid box thickness (cm): "dzt(k) k=1..km"')
9102  format (/,' "w" grid box thickness (cm): ',
     $        '"dzw(k) k=0..km"')
9103  format (/,' depth to "t" box bottom (cm): "zw(k) k=1..km"')
9104  format (/,' depth to "t" & "u,v" grid points (cm): "zt(k) ',
     $        'k=1..kmp1"')
9105  format (/,' latitude of t,s points (deg): "yt(j) j=1,jmt"')
9106  format (/,' latitude of u,v points (deg): "yu(j) j=1,jmt"')
9107  format (/,' longitude of t,s points (deg): "xt(i) i=1,imt"')
9108  format (/,' longitude of u,v points (deg): "xu(i) i=1,imt"')
c
       end
      function indp(c, ca, ia)
c
c
c=======================================================================
c
c     indp = index of nearest data point
c
c     inputs:
c
c     c      = arbitrary data point...same units as ca
c     ca     = array of data points  (must be monotonically increasing)
c     ia     = dimension of ca
c
c     output:
c
c     indp =  index of nearest data point to "c"
c             if "c" is outside the domain of "ca" then indp is set to
c             1 or ia depending on  whether c(1) or c(ia) is closest
c
c             note: if "ca" is dimensioned ca(0:ia) in the calling
c                   program, then the returned index should be reduced
c                   by one to account for the zero base.
c
c     example:
c
c     let model depths be defined by the following:
c     parameter (km=5)
c     dimension z(km)
c     data z /5.0, 10.0, 50.0, 100.0, 250.0/
c
c     k1 = indp (12.5, z, km)
c     k2 = indp (0.0, z, km)
c
c     k1 would be set to 2, & k2 would be set to 1 so that
c     z(k1) would be the nearest data point to 12.5 and z(k2) would
c     be the nearest data point to 0.0
c
c=======================================================================
c
c
#ifdef multitasking
cfpp$ noconcur r
#endif
      dimension ca(ia)
c
      if (c .lt. ca(1) .or. c .gt. ca(ia)) then
        if (c .lt. ca(1))  indp = 1
        if (c .gt. ca(ia)) indp = ia
        return
      else
        do 10 i=2,ia
          if (c .le. ca(i)) then
            indp = i
            if (ca(i)-c .gt. c-ca(i-1)) indp = i-1
            go to 11
          endif
10      continue
11      continue
      endif
      return
      end
      subroutine mesh (dx, dy)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     construct grid spacing "dx" & "dy" by building up "nxres"
c     and "nyres" sub domains
c
c     inputs:
c
c     see "coord.h"
c
c
c     outputs:
c
c     dx = array of longitudinal widths of "t" boxes
c     dy = array of latitudinal heights of "t" boxes
c     see "coord.h" for others
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "switch.h"
      dimension dx(imt), dy(jmt)
c
c-----------------------------------------------------------------------
c     construct longitudinal domains
c-----------------------------------------------------------------------
c
      m  = 1
      xs = stlon
      do 30 l=1,nxres
#if defined cpldleg3 || defined cpldleg2
        if (pintro) then
#endif
        write (stdout,'(/,1x,a22,i2)') 'longitude section # ', l
#if defined cpldleg3 || defined cpldleg2
        endif
#endif
        call cmesh (dx(m), xmax(l), xmin(l), xwid(l), n, idir(l))
        sum = c0
        do 20 i=m,m+n-1
          sum = sum + dx(i)
20      continue
        xe = xs + sum
#if defined cpldleg3 || defined cpldleg2
        if (pintro) then
#endif
        write (stdout,'(/,1x,a21,g14.7,a3,g14.7)')
     $       'longitudinal domain: ', xs, '==>', xe
        write (stdout,'(/,1x,a5)')  ' dxt='
        write (stdout,'(1x,10g12.6)') (dx(m+i-1),i=1,n)
#if defined cpldleg3 || defined cpldleg2
        endif
#endif
        xs = xe
        m = m + n
30    continue
c
c-----------------------------------------------------------------------
c     construct latitudinal domains
c-----------------------------------------------------------------------
c
      m  = 1
      ys = stlat
c
      do 100 l=1,nyres
#if defined cpldleg3 || defined cpldleg2
        if (pintro) then
#endif
        write (stdout,'(/,1x,a21,i2)') 'latitude section # ', l
#if defined cpldleg3 || defined cpldleg2
        endif
#endif
c
#ifdef worldocn
c       get gaussian latitudes to match atmospheric model grid (ocean 
c       model's u,v pts match gaussian half lats, in order to have t,s
c       box north-south dimensions match those of spectral model, 
c       thereby conserving heat and water fluxes across air-sea 
c       interface)
c
        n = jmt
# if defined cpldleg3 || defined cpldleg2
        if (pintro) then
# endif
# ifdef ocnr15
        write (stdout,'(/,1x,a28)') 'using r15 gaussian latitudes'
# endif
# ifdef ocnr30
        write (stdout,'(/,1x,a28)') 'using r30 gaussian latitudes'
# endif
# if defined cpldleg3 || defined cpldleg2
        endif
# endif
c
        call speclt ( dy )
c
        stlat = -90.0 - p5*dy(1)
        ys = stlat
#else
        call cmesh (dy(m), ymax(l), ymin(l), ywid(l), n, jdir(l))
#endif
c
        sum = c0
        do 90 j=m,m+n-1
          sum = sum + dy(j)
90      continue
        ye = ys + sum
# if defined cpldleg3 || defined cpldleg2
        if (pintro) then
# endif
        write (stdout,'(/,1x,a20,g14.7,a3,g14.7)')
     $         'latitudinal domain: ', ys, '==>', ye
        write (stdout,'(/,1x,a5)')  ' dyt='
        write (stdout,'(1x,10g12.6)') (dy(m+j-1),j=1,n)
# if defined cpldleg3 || defined cpldleg2
        endif
# endif
        ys = ye
        m = m + n
100   continue
c
      return
      end
