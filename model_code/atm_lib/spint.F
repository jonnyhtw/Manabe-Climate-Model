      subroutine spint(spv,ps,glmn,scale)
c
c  this subroutine computes in spectral space the covariance of spectral
c   variables spv and ps times a scale factor (scale) and returns result
c   to variable glmn.  generally this is a global mean pressure-weighted
c   integral of spv.  (see notes of nov. 24, 1981 and sept. 13, 1983.)
c   factors of 1/2 and 2 result from spherical harmonic normalization
c   conventions.  (see notes of january 27, 1982.)
c
c  note: variables spv and ps must be even functions in a hemispheric
c        model or else this routine will return zero or a meaningless
c        result.
c
#include "parm.h"
#include "parmjx.h"
#include "parmkp.h"
c
      complex spv(mx,jx,kx),ps(mx,jx)
c
#include "sgwgt1.h"
c
      glmn=0.0
      do 2 k=1,kx
      do 2 n=1,jtrun1,jin
      glmn=glmn + real (ps(1,n)*conjg(spv(1,n,k)))*dq(k)
  2   continue
      do 3 k=1,kx
      do 3 n=1,jtrun1,jin
      do 3 m=2,mx
      glmn=glmn+2.*real(ps(m,n)*conjg(spv(m,n,k)))*dq(k)
  3   continue
      glmn=.5*scale*glmn
      return
      end
