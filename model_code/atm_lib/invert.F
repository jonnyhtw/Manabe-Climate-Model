      subroutine invert(exmat,n,nmt)
c
c   called from: kernal(impvar)
c   subroutines called: matinv
c
c   this subroutine performs an inversion of nmt n x n matrices.
c   the input array exmat is stored in a smaller working array a.
c   the inverted matrix is transferred from the same working array
c   back to exmat.
c
#include "parm.h"
c
      parameter (m=kx)
c   m is the maximum value of n from the argument list.
c   the number of model levels kx must not exceed the parameter m.
c
      dimension a(m,m), exmat(nmt,n,n)
c
c
c  <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      if ( n .le. kx ) go to 5
      write(6,6)
      stop
c
    5 continue
      do 20 lp=1,nmt
      do 10 j=1,n
      do 10 i=1,n
   10 a(i,j)=exmat(lp,i,j)
c
      call matinv(a,n,det)
c
      do 30 j=1,n
      do 30 i=1,n
   30 exmat(lp,i,j)=a(i,j)
   20 continue
c
    6 format(10x,'stop in subroutine invert. the work array a is too',
     *       ' small.')
c
      return
      end
