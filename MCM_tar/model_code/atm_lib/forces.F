      subroutine forces
#if defined mixed_layer && defined ice_dynamics
c
c
c  subroutine to calculate free drift velocity
c  in B-grid and obtain forcing fields in C-grid. B-grid
c  velocity field is used as an initial guess to speed 
c  iterative solution in FREES.
c
c  G.M. Flato - Institute of Ocean Sciences, Sidney, BC
c
c***********************************************************
c
c  7/23/93
c
c  modified to distinguish between air drag and ocean drag
c  turning angle
c
c  C. F. Ip - GFDL, Princeton, NJ
c
c  9/17/93
c
c  modified to allow use with global climate model
c
c  C. F. Ip - GFDL, Princeton, NJ
c************************************************************
c
#include "parm.h"
c
#include "variceco.h"
c
      diff = (theta-phi)*piice/180.
      sintmp = sin(diff)
      costmp = cos(diff)
c
c.o   sinphi=0.4226D+00
c.o   cosphi=0.9063D+00
c.o   fcor=1.46D-04
c.o   rhoice=0.92*1.D+03
c
c  calculate linear free drift
c
c
c  cfi
c
c  the solution can be computed analytically
c  using complex notation approach
c
c  linear means linear water drag (still quadratic wind drag)
c
c  cfi
c
      do j=1,np
      do i=1,nl
c     do k = 2, nl+1
c        i = k 
c        if (i .eq. nl+1) then
c           i = 1
c        end if
         i0=i
         if(i0.eq.1) i0=nl+1
         heffa=(hef(i0,j,1)+hef(i+1,j,1)+hef(i0,j+1,1)
     <        +hef(i+1,j+1,1))/4.
	 delta = atan((sintmp*dwatn(i,j)+rhoice*heffa*fcorp(j)
     >     *cosphi)/(costmp*dwatn(i,j)+rhoice*heffa*fcorp(j)*sinphi))
c.o      delta=57.29575*atan((rhoice*heffa*fcor*cosphi)/
c.o  <         (dwatn(i,j)+rhoice*heffa*fcor*sinphi))
c
c use surface wind stress as input (not geostrophic wind)
c
         Ui = (dairn(i,j)*cos(delta))/(dwatn(i,j)*costmp
c        Ug=sqrt(gairx(i,j)**2+gairy(i,j)**2)
c        Ui = (dairn(i,j)*Ug*cos(delta))/(dwatn(i,j)*costmp
     >     +rhoice*heffa*fcorp(j)*sinphi)
c.o      Ui=(dairn(i,j)*Ug*cos(delta*.0174533))/
c.o  <      (dwatn(i,j)+rhoice*heffa*fcor*sinphi)
c
c use surface wind stress as input (not geostrophic wind)
c
         gamma=atan2(tauyin(i,j),tauxin(i,j))
c        gamma=atan2(gairy(i,j),gairx(i,j))
c.o      gamma=57.29575*atan2(gairy(i,j),gairx(i,j))
         uice(i,j,1) = (Ui*cos(gamma-delta)+gwatx(i,j))
         vice(i,j,1) = (Ui*sin(gamma-delta)+gwaty(i,j))
c.o      uice(i,j,1)=(Ui*cos((gamma-delta)*.0174533)+gwatx(i,j))
c.o      vice(i,j,1)=(Ui*sin((gamma-delta)*.0174533)+gwaty(i,j))
      enddo
      enddo
c
c  average forces to C-Grid
c
         do j=1,np-1
         do i=1,nl
c     do k = 2, nl+1
c        i = k 
c        if (i .eq. nl+1) then
c           i = 1
c        end if
            i0=i
            if(i0.eq.1) i0=nl+1
            i1=i+1
            if(i1.gt.nl) i1=i1-nl
            i2=i+2
            if(i2.gt.nl+1) i2=i2-nl
            Aij=dwatn(i,j)*costh
c.o         Aij=dwatn(i,j)*cosphi
            Bij=((hef(i+1,j+1,1)+hef(i0,j+1,1)+hef(i0,j,1)
     <      +hef(i+1,j,1))/4.*rhoice*fcorp(j))+(dwatn(i,j)*sinth)
c.o         Bij=((hef(i+1,j+1,1)+hef(i0,j+1,1)+hef(i0,j,1)
c.o  <      +hef(i+1,j,1))/4.*rhoice*fcor)+(dwatn(i,j)*sinphi)
            Ai1j=dwatn(i1,j)*costh
c.o         Ai1j=dwatn(i1,j)*cosphi
            Bi1j=((hef(i2,j+1,1)+hef(i+1,j+1,1)+hef(i+1,j,1)
     <      +hef(i2,j,1))/4.*rhoice*fcorp(j))+(dwatn(i1,j)*sinth)
c.o         Bi1j=((hef(i2,j+1,1)+hef(i+1,j+1,1)+hef(i+1,j,1)
c.o  <      +hef(i2,j,1))/4.*rhoice*fcor)+(dwatn(i1,j)*sinphi)
            Aij1=dwatn(i,j+1)*costh
c.o         Aij1=dwatn(i,j+1)*cosphi
            Bij1=((hef(i+1,j+2,1)+hef(i0,j+2,1)+hef(i0,j+1,1)
     <      +hef(i+1,j+1,1))/4.*rhoice*fcorp(j+1))+(dwatn(i,j+1)*sinth)
c.o         Bij1=((hef(i+1,j+2,1)+hef(i0,j+2,1)+hef(i0,j+1,1)
c.o  <      +hef(i+1,j+1,1))/4.*rhoice*fcor)+(dwatn(i,j+1)*sinphi)
c
            XX(i,j)=((Aij*uice(i,j,1)-Bij*vice(i,j,1))+(Aij1
     <              *uice(i,j+1,1)-Bij1*vice(i,j+1,1)))/2.
            YY(i,j)=((Aij*vice(i,j,1)+Bij*uice(i,j,1))+(Ai1j
     <              *vice(i1,j,1)+Bi1j*uice(i1,j,1)))/2.
         enddo
         enddo
c
#endif
      return
      end
