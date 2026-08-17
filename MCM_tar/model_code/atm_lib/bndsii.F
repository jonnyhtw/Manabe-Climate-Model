      subroutine bndsii
#if defined mixed_layer && defined ice_dynamics
c
c  subroutine to generate boundary masks for spherical
c  grid cavitating fluid model
c
c  G.M. Flato - Institute of Ocean Sciences, Sidney, BC
c
c  modified to be used in global climate model
c  all velocity masks are computed from thickness mask
c
c  Chi Ip - GFDL, Princeton, NJ
c*******************************************************
c
#include "parm.h"
c
#include "variceco.h"
c
      do j=1,np
      do i=1,nl
         uvmx(i,j)=0.
         uvmy(i,j)=0.
      enddo
      enddo
c
c  calculate velocity mask in C-Grid
c
      do j=1,np
      do i=1,nl
         i0=i
         if(i0.eq.1) i0=nl+1
         htx=heffm(i0,j+1)+heffm(i+1,j+1)
         hty=heffm(i+1,j)+heffm(i+1,j+1)
         if(htx.gt.1.5) uvmx(i,j)=1.0
         if(hty.gt.1.5) uvmy(i,j)=1.0
      enddo
      enddo
c
c  calculate velocity mask in B-Grid
c
      do j=1,np
      do i=1,nl
         if(j.ne.1.and.j.ne.np) then
         i1=i-1
         if(i1.lt.1)i1=nl
            tj=uvmx(i,j)+uvmy(i,j)+uvmx(i,j-1)+uvmy(i1,j)
            if(tj.lt.3.5) then
               uvmb(i,j)=0.
            else
               uvmb(i,j)=1.
            endif
         else
            uvmb(i,j)=0.
         endif
      enddo
      enddo
c
#endif
      return
      end
