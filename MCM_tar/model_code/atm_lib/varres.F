c
c subroutine to reset everything for each hemisphere
c
      subroutine varres
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
#include "variceco.h"
c
      do j = 1, iy
         arlat(j) = 0.
         arlatc(j) = 0.
         fcorp(j) = 0.
         fcorc(j) = 0.
      end do
      arlatc(iy+1) = 0.
      fcorc(iy+1) = 0.
c
      do j = 1, iy
      do i = 1, ix
         uvmx(i,j) = 0.
	 uvmy(i,j) = 0.
	 uvmb(i,j) = 0.
	 xx(i,j) = 0.
	 yy(i,j) = 0.
	 axi(i,j) = 0.
	 ayi(i,j) = 0.
	 bxi(i,j) = 0.
	 byi(i,j) = 0.
	 dwatn(i,j) = 0.
	 dwatx(i,j) = 0.
	 dwaty(i,j) = 0.
	 dairn(i,j) = 0.
	 dairx(i,j) = 0.
	 dairy(i,j) = 0.
	 gairx(i,j) = 0.
	 gairy(i,j) = 0.
	 gwatx(i,j) = 0.
	 gwaty(i,j) = 0.
	 tauxin(i,j) = 0.
	 tauyin(i,j) = 0.
      end do
      end do
c
      do j = 1, iy+1
      do i = 1, ix+1
         out(i,j) = 0.
	 heffm(i,j) = 0.
	 tauxo(i,j) = 0.
	 tauyo(i,j) = 0.
	 xmpress(i,j) = 0.
         presur(i,j) = 0.
         fxi(i,j) = 0.
         fyi(i,j) = 0.
      end do
      end do
c
      do j = 1, iy
      do i = 1, ix
      do l = 1, 3
         uice(i,j,l) = 0.
         vice(i,j,l) = 0.
	 hef(i,j,l) = 0.
	 are(i,j,l) = 0.
      end do
      end do
      end do
c
      do j = 1, iy
      do l = 1, 3
         hef(ix+1,j,l) = 0.
	 are(ix+1,j,l) = 0.
      end do
      end do
c
      do i = 1, ix
      do l = 1, 3
         hef(i,iy+1,l) = 0.
	 are(i,iy+1,l) = 0.
      end do
      end do
c
#endif
      return
      end
