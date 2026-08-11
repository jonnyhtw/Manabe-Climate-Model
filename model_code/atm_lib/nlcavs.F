      subroutine nlcavs
#if defined mixed_layer && defined ice_dynamics
c
c  subroutine to calculate and apply non-linear drag cavitating
c  fluid velocity correction in spherical C-grid
c
c  G.M. Flato - Institute of Ocean Sciences, Sidney, BC
c
c*******************************************************
c
c  9/8/93
c
c  modified to allow for irregular latitude spacing
c
c  C. F. Ip - GFDL, Princeton, NJ
c
c  9/17/93
c
c  modified to use with global climate model
c
c  C. F. Ip - GFDL, Princeton, NJ
c*******************************************************
c
#include "parm.h"
c
#include "variceco.h"
c
      iterat=1
      maxit=1500
      toler=0.0005
c
  100 continue
c     print *, 'iterat = ', iterat
c
c  update velocity arrays
c
      do j=1,np
      do i=1,nl
         uice(i,j,2)=uice(i,j,1)*uvmx(i,j)
         vice(i,j,2)=vice(i,j,1)*uvmy(i,j)
      enddo
      enddo
c
      do j=1,np-1
      kkm=1
c
c iterate more for cells near the poles
c
      if (ihem .eq. 1) then
         if(j.eq.np-1) kkm=6
         if(j.eq.np-2) kkm=3
      end if
c
      if (ihem .eq. 2) then
         if (j .eq. 3) kkm = 6
         if (j .eq. 4) kkm = 3
      end if
c
      do kkt=1,kkm
      do i=1,nl
c     do k=2,nl+1
c        i = k
c     if (i .eq. nl+1) then
c        i = 1
c     end if
      i1=i+1
      if(i1.gt.nl) i1=i1-nl
c
c  do nothing if not an open water grid cell
c
      if(heffm(i+1,j+1).eq.0.) go to 101
c
c  do not correct outflow grid cells
c
      if(out(i+1,j+1).eq.1.) go to 101
c
c do not correct for isolated open water grid cell
c
      uvms = uvmx(i,j) + uvmx(i1,j) + uvmy(i,j) + uvmy(i,j+1)
      if (uvms .eq. 0.) go to 101
c
c set latitudes and differences because of irregular spacings (cfi)
c
      phi1 = arlat(j)
      phi2 = arlat(j+1)
      phi12 = arlatc(j+1)
      dphi = abs(phi1-phi2)
      dphic1 = abs(arlatc(j)-arlatc(j+1))
      dphic2 = abs(arlatc(j+1)-arlatc(j+2))
c     print *, i, j, i1
c     print *, phi1, phi2, phi12
c     print *, dphi, dphic1, dphic2
c     print *, i, j, Axi(i,j), Ayi(i,j), uice(i,j,1), vice(i,j,1)
c     print *, i1, j, Axi(i1,j), Ayi(i1,j)
c     print *, dlam
c     print *, '---------------------------------'
c
c.o   phi1=phi0+float(j)*dphid
c.o   phi2=phi0+float(j+1)*dphid
c.o   phi12=(phi1+phi2)/2.
c
c compute corresponding cosines
c
      cosp1 = cos(phi1)
      cosp2 = cos(phi2)
      cosp12 = cos(phi12)
c
c.o   cosp1=cos(phi1*.0174533)
c.o   cosp2=cos(phi2*.0174533)
c.o   cosp12=cos(phi12*.0174533)
c
c  calculate divergence (before)
c
c     print *, '1 ', uice(2,2,1), vice(2,2,1)
      D=(1./(R*cosp12))*((1./dphi)*(vice(i,j+1,1)*cosp2
     <-vice(i,j,1)*cosp1)+(1./dlam)*(uice(i1,j,1)-uice(i,j,1)))
c
c  calculate pressure increment
c
c     print *, '2 ', uice(2,2,1), vice(2,2,1)
c     pq = (((uvmx(i,j)/Axi(i,j)+uvmx(i1,j)/Axi(i1,j))
c    <     /(cosp12*dlam**2))+(cosp1*uvmy(i,j)/Ayi(i,j)/dphic1
c    <     +cosp2*uvmy(i,j+1)/Ayi(i,j+1)/dphic2)/dphi)
c     if (pq .eq. 0.) then
c        print *, 'pq 0', i, j, i1, j+1
c        stop
c     end if
c     print *, '3 ', uice(2,2,1), vice(2,2,1), D
c
c adjust because of irregular latitudes (cfi)
c
      delp1=(-D*R*R*cosp12)/(((uvmx(i,j)/Axi(i,j)+uvmx(i1,j)/Axi(i1,j))
     <     /(cosp12*dlam**2))+(cosp1*uvmy(i,j)/Ayi(i,j)/dphic1
     <     +cosp2*uvmy(i,j+1)/Ayi(i,j+1)/dphic2)/dphi)
c.o  <     /(cosp12*dlam**2))+(cosp1*uvmy(i,j)/Ayi(i,j)
c.o  <     +cosp2*uvmy(i,j+1)/Ayi(i,j+1))/(dphi**2))
c     print *, '4 ', uice(2,2,1), vice(2,2,1), delp
      if(D.le.0.) delp=delp1
c
c  now test pressure against strength
c
c     print *, '5 ', uice(2,2,1), vice(2,2,1), delp
      pdiff=xmpress(i+1,j+1)-(presur(i+1,j+1)+delp)
      if(pdiff.lt.0.) delp=xmpress(i+1,j+1)-presur(i+1,j+1)
c
      if(D.gt.0..and.presur(i+1,j+1).le.0.) then
         delp=0.
      endif
      if(D.gt.0..and.presur(i+1,j+1).gt.0.) then
         delp=(max(delp1,-presur(i+1,j+1)))
      endif
c     print *, '6 ', i, j, are(i+1,j+1,1), hef(i+1,j+1,1)
c     print *, D, xmpress(i+1,j+1), presur(i+1,j+1), pdiff, delp
c     print *, '-------------------------------------'
c
c  calculate and apply velocity corrections
c
      du1=(-1./(Axi(i,j)*R*cosp12))*(delp/dlam)
      du2=(1./(Axi(i1,j)*R*cosp12))*(delp/dlam)
      dv1=(-1./(Ayi(i,j)*R))*(delp/dphi)
      dv2=(1./(Ayi(i,j+1)*R))*(delp/dphi)
      uice(i,j,1)=(uice(i,j,1)+du1)*uvmx(i,j)
      uice(i1,j,1)=(uice(i1,j,1)+du2)*uvmx(i1,j)
      vice(i,j,1)=(vice(i,j,1)+dv1)*uvmy(i,j)
      vice(i,j+1,1)=(vice(i,j+1,1)+dv2)*uvmy(i,j+1)
      presur(i+1,j+1)=presur(i+1,j+1)+delp
c
  101 continue
      enddo
      enddo
      enddo
c
c  calculate max error
c
      xmaxe=0.
      do j=1,np
      do i=1,nl
         uerr=abs(uice(i,j,1)-uice(i,j,2))
         verr=abs(vice(i,j,1)-vice(i,j,2))
         if(uerr.ge.xmaxe) then
            xmaxe=uerr
            iu=i
            ju=j
            xmaxu=uerr
         endif
         if(verr.ge.xmaxe) then
            xmaxe=verr
            iv=i
            jv=j
            xmaxv=verr
         endif
      enddo
      enddo
c
c  compare with tolerance
c
      if(xmaxe.gt.toler) then
         iterat=iterat+1
         if(iterat.lt.maxit) then
            go to 100
         else
            print*,'    **** max. iterations exceeded in NLCAVS'
            print*,' iu  ju xmaxu ',iu,ju,xmaxu
            print*,' iv  jv xmaxv ',iv,jv,xmaxv
         endif
      endif
c
      print*,'# of iteration in NLCAVS ',iterat
      print*,'max error ',xmaxe
c
#endif
      return
      end
