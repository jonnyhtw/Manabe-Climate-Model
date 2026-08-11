      subroutine frees
#if defined mixed_layer && defined ice_dynamics
c
c  subroutine to calculate free drift velocity field
c  in C-grid by under-relaxation
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
c  9/16/93
c
c  modified to allow use of global climate model
c
c  C. F. Ip - GFDL, Princeton, NJ
c*******************************************************
c
#include "parm.h"
c
#include "variceco.h"
c
#include "ctlblk.h"
c
      nmax=500
      tol=0.00001
      w=0.4
c
c  relaxation scheme to find velocity field given pressure field
c
      do kk=1,nmax
c
      do j=1,np
      do i=1,nl
         uice(i,j,2)=uice(i,j,1)*uvmx(i,j)
         vice(i,j,2)=vice(i,j,1)*uvmy(i,j)
      enddo
      enddo
c
      do j=1,np
      do i=1,nl
c     do k=2,nl+1
c      	 i = k
c        if (i .eq. nl+1) then
c           i = 1
c        end if
         jm1=j-1
         if(j.eq.1) jm1=1
         im1=i-1
         if(i.eq.1) im1=nl
         i0=i
         if(i0.eq.1) i0=nl+1
         i1=i+1
         if(i1.gt.nl) i1=i1-nl
c
c calculate latitude difference between centers of grid (cfi) 
c and cosine of center of grid
c
         dphi = abs(arlatc(j+1)-arlatc(j))
         cosp = cos(arlatc(j+1))
c
c.o      cosp=cos((phi0+(float(j)+0.5)*dphid)*.0174533)
         dpdx=(presur(i+1,j+1)-presur(i0,j+1))/(R*dlam*cosp)
         dpdy=(presur(i+1,j+1)-presur(i+1,j))/(R*dphi)

         if(j+1.le.np) then
         Bvbar=0.25*Bxi(i,j)*(vice(im1,j,1)+vice(i,j,1)
     <         +vice(im1,j+1,1)+vice(i,j+1,1))
         else
         Bvbar=0.5*Bxi(i,j)*(vice(im1,j,1)+vice(i,j,1))
         endif

         Bubar=0.25*Byi(i,j)*(uice(i1,j,1)+uice(i,j,1)
     <         +uice(i,jm1,1)+uice(i1,jm1,1))
         R1=dpdx-XX(i,j)-Bvbar
         R2=dpdy-YY(i,j)+Bubar
         uice(i,j,1)=(w*(-1./Axi(i,j))*R1
     <               +(1.-w)*uice(i,j,2))*uvmx(i,j)
         vice(i,j,1)=(w*(-1./Ayi(i,j))*R2
     <               +(1.-w)*vice(i,j,2))*uvmy(i,j)
      enddo
      enddo
c
      errm=0.
      do j=1,np
      do i=1,nl
         erru=abs(uice(i,j,1)-uice(i,j,2))
         errv=abs(vice(i,j,1)-vice(i,j,2))
         if(erru.gt.errm)errm=erru
         if(errv.gt.errm)errm=errv
      enddo
      enddo
      if(errm.le.tol) go to 999
c
      enddo
      if (kk.ge.nmax) then
        print*,'max number of iterations exceeded in frees'
        print *, errm
      end if
c
  999 print*,'# of iterations in frees =',kk
#endif
      return
      end
