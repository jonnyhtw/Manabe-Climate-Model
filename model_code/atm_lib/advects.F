      subroutine advects
#if defined mixed_layer && defined ice_dynamics
c
c  subroutine to do advection in spherical C-grid
c  using upstream differencing
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
c  9/17/83
c
c  modified to allow use of global climate model
c
c  C. F. Ip - GFDL, Princeton, NJ
c*******************************************************
c
#include "parm.h"
c
c.ocinclude "iohist.h"
c
#include "variceco.h"
c
      do j=1,np
      do i=1,nl
         hef(i+1,j+1,2)=hef(i+1,j+1,1)
         are(i+1,j+1,2)=are(i+1,j+1,1)
      enddo
      enddo
c
c
      do j=1,np-1
      do i=1,nl
      i0=i
      if(i0.eq.1) i0=nl+1
      i1=i+1
      if(i1.gt.nl) i1=i1-nl
      i2=i+2
      if(i2.gt.nl+1) i2=i2-nl
c
c use latitude arrays and grid length (in radians) (cfi)
c
      phi1 = arlat(j)
      phi2 = arlat(j+1)
      dphi = abs(phi1-phi2)
c
c.o   phi1=phi0+float(j)*dphid
c.o   phi2=phi0+float(j+1)*dphid
c
c use latitude array for center of grid point (cfi)
c
      phi12 = arlatc(j+1)
c
c.o   phi12=(phi1+phi2)/2.
c
c compute corresponding cosines (cfi)
c
      cosp1 = cos(phi1)
      cosp2 = cos(phi2)
      cosp12 = cos(phi12)
c
c.o   cosp1=cos(phi1*.0174533)
c.o   cosp2=cos(phi2*.0174533)
c.o   cosp12=cos(phi12*.0174533)
c
c  if not an open water grid cell, do nothing
c
      if(heffm(i+1,j+1).gt.0.1) then
c
         um=uice(i,j,1)
         vm=vice(i,j,1)
         ump=uice(i1,j,1)
         vmp=vice(i,j+1,1)
c
         if(um.ge.0.) qinx=um*hef(i0,j+1,2)*heffm(i0,j+1)
         if(um.lt.0.) qinx=um*hef(i+1,j+1,2)*heffm(i+1,j+1)
         if(vm.ge.0.) qiny=vm*hef(i+1,j,2)*heffm(i+1,j)
         if(vm.lt.0.) qiny=vm*hef(i+1,j+1,2)*heffm(i+1,j+1)
c
         if(ump.ge.0.) qinxp=-ump*hef(i+1,j+1,2)*heffm(i+1,j+1)
         if(ump.lt.0.) qinxp=-ump*hef(i2,j+1,2)*heffm(i2,j+1)
         if(vmp.ge.0.) qinyp=-vmp*hef(i+1,j+1,2)*heffm(i+1,j+1)
         if(vmp.lt.0.) qinyp=-vmp*hef(i+1,j+2,2)*heffm(i+1,j+2)
c
         qijh=hef(i+1,j+1,2)*heffm(i+1,j+1)
c
c  now do advection of open water fraction
c
         if(um.ge.0.) qax=um*are(i0,j+1,2)*heffm(i0,j+1)
         if(um.lt.0.) qax=um*are(i+1,j+1,2)*heffm(i+1,j+1)
         if(vm.ge.0.) qay=vm*are(i+1,j,2)*heffm(i+1,j)
         if(vm.lt.0.) qay=vm*are(i+1,j+1,2)*heffm(i+1,j+1)
c
         if(ump.ge.0.) qaxp=-ump*are(i+1,j+1,2)*heffm(i+1,j+1)
         if(ump.lt.0.) qaxp=-ump*are(i2,j+1,2)*heffm(i2,j+1)
         if(vmp.ge.0.) qayp=-vmp*are(i+1,j+1,2)*heffm(i+1,j+1)
         if(vmp.lt.0.) qayp=-vmp*are(i+1,j+2,2)*heffm(i+1,j+2)
c
         qija=are(i+1,j+1,2)*heffm(i+1,j+1)
c
c  update thickness array
c
      if(i.eq.1.and.j.eq.8) then
         f1=-qiny*cosp1*R*dlam*1.D-6
      endif
      if(i.eq.36.and.j.eq.8) then
         f2=-qiny*cosp1*R*dlam*1.D-6
      endif
         hef(i+1,j+1,1)=qijh+(deltat/(R*cosp12))
     <      *((qinx+qinxp)/dlam+(qiny*cosp1+qinyp
     <      *cosp2)/dphi)
     <      *heffm(i+1,j+1)
c
c  update area array
c
         are(i+1,j+1,1)=qija+(deltat/(R*cosp12))
     <      *((qax+qaxp)/dlam+(qay*cosp1+qayp*cosp2)/dphi)
     <      *heffm(i+1,j+1)
c
      endif
c
      enddo
      enddo
c
c.o   write(intgl,888)iday,f1+f2
  888 format(i4,f16.8)
c
#endif
      RETURN
      END
