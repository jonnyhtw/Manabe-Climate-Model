      subroutine invtri (z, topbc, botbc, vmc, tdt, kmz, mask, is, ie
     $,                  j)
      implicit real*8 (a-h,o-z)
#if defined implicitvmix || defined isopycmix
c
c=======================================================================
c     solve the vertical diffusion equation implicitly using the
c     method of inverting a tridiagonal matrix as described in
c     richtmyer and morton.
c     this routine assums that the variables are defined at grid points
c     and the top and bottom b.c. are flux conditions.
c
c     inputs:
c     z         = right hand side terms
c     topbc     = top boundary condition
c     botbc     = bottom boundary condition
c     vmc       = vertical mixing coeff
c     tdt       = 2 * timestep
c     kmz       = level indicator
c     mask      = land - sea mask
c     is        = index of starting longitude
c     ie        = index of ending longitude
c     j         = latitude row
c
c     outputs:
c     z         = returned solution
c
c=======================================================================
c
#include "param.h"
#include "cvmix.h"
#include "grdvar.h"
      dimension z(imt,km)
      dimension topbc(imt), botbc(imt)
      dimension vmc(imt,km)
      dimension kmz(imt,jmt), tdt(km)
      real mask(imt,km)
c
      dimension a(imt,km), b(imt,km), c(imt,km), d(imt,km)
      dimension e(imt,0:km), f(imt,0:km), g(imt)
c
      do 100 k=2,km
        do 90 i=is,ie
          a(i,k)   = vmc(i,k-1)*dztur(k)*tdt(k)*aidif
          c(i,k)   = vmc(i,k)*dztlr(k)*tdt(k)*aidif
          b(i,k)   = c1 + a(i,k) + c(i,k)
          d(i,k)   = z(i,k)
          e(i,k-1) = c0
          f(i,k-1) = c0
90      continue
100   continue
c
c     b. c. at top
c
      k = 1
      do 300 i=is,ie
        a(i,k)   = dztr(k)*tdt(k)*aidif
        c(i,k)   = vmc(i,k)*dztlr(k)*tdt(k)*aidif
        b(i,k)   = c1 + c(i,k)
        d(i,k)   = z(i,k)
        e(i,k-1) = c0
        f(i,k-1) = c0
300   continue
c
c     b. c. at bottom
c
      do 400 i=is,ie
        kz = kmz(i,j)
        if (kz .ne. 0) then
          b(i,kz) = c1 + a(i,kz)
          c(i,kz) = dztr(kz)*tdt(kz)*aidif
          e(i,kz) = c0
          f(i,kz) = -botbc(i)
        endif
400   continue
c
c-----------------------------------------------------------------------
c     now invert
c-----------------------------------------------------------------------
c
      do 500 k=km,1,-1
        do 490 i=is,ie
          if (k .le. kmz(i,j)) then
            g(i)     = c1/(b(i,k)-c(i,k)*e(i,k))
            e(i,k-1) = a(i,k)*g(i)
            f(i,k-1) = (d(i,k)+c(i,k)*f(i,k))*g(i)
          endif
490     continue
500   continue
c
c     b.c. at surface
c
      do 600 i=is,ie
        z(i,1) = (e(i,0)*topbc(i) + f(i,0))*mask(i,1)
600   continue
c
      do 700 k=2,km
        do 690 i=is,ie
          z(i,k) = (e(i,k-1)*z(i,k-1) + f(i,k-1))*mask(i,k)
690     continue
700   continue
#endif
      return
      end
