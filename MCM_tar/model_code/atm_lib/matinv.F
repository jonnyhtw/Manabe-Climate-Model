      subroutine  matinv(a, n, det)
c
c   called from: invert
c   subroutines called: none
c
c   this subroutine inverts n x n nonsingular matrix and finds
c   its determinant.
c
c   arguments -
c        a   : input matrix of dimension n x n to be inverted
c        n   : dimension of matrix a (must be equal to kx)
c        det : determinant of matrix a
c
c   this subroutine uses the method of elementary
c   tranformations to find the inverse of matrix a. the input
c   matrix a is destroyed in computation and the inverse matrix
c   (a-inverse) is stored in a. for numerical accuracy, elementary
c   transformations are performed in conjunction with the
c   'pivotal' element method.
c
c   the maximum size, n, of the matrix 'a' to be inverted is
c   limited by the value of the parameter nmax. however, this can be
c   changed easily by changing the parameter statement for nmax.
c   if the input matrix is singular (determinant less
c   than 1.0e-40), an error message is printed out and the
c   program is terminated.
c
#include "parm.h"
c
      parameter ( nmax = kx, nmax2 = 2 * nmax )
c
      dimension a(n,n)
      dimension ac(nmax2,nmax), dd(nmax2), h(nmax2), temp(nmax2,nmax)
c
c  <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      if ( n .gt. nmax) then
        write(6,3457)
        stop
      endif
c
      np1 = n + 1
      n2  = 2.0 * n
c
c   initialize arrays
c
      det = 1.0
      do 1 j=1,n
        do 1 i=1,n
          ac(i,j) = a(i,j)
    1 continue
c
      do 3 j=1,n
        do 2 i=np1,n2
          ac(i,j) = 0.0
    2   continue
        ac(j+n,j) = 1.0
    3 continue
c
      do 10 k=1,n
c
c   find the largest element in the row
c
        do 32 j=k,n
          h(j) = ac(k,j)
   32   continue
        m   = n - k + 1
        vmx = 0.0
        l   = k
        do 31 lh=1,m
          index = k + lh - 1
          test  = abs(h(index))
          if (test .ge. vmx) then
            vmx = test
            l   = index
          endif
   31   continue
c
c   interchange columns if the largest element is not the
c   diagonal element.
c
        if ((k-l) .lt. 0) then
          do 33 i=k,n2
            dd(i)   = ac(i,k)
            ac(i,k) = ac(i,l)
            ac(i,l) = dd(i)
   33     continue
          det = -det
        endif
c
c   divide the column by the largest element
c
        det = det * ac(k,k)
#ifndef a_32_bit_machine
        if (abs(det) .lt. 1.0e-40) then
#endif
#ifdef a_32_bit_machine
c  For 32 bit machines, make the number larger so that there 
c  is not a compiler error message.
        if (abs(det) .lt. 1.0e-30) then
#endif
          write(6,3456) det
          stop
        endif
c
        d = 1.0 / ac(k,k)
        do 60 i=k,n2
          h(i) = ac(i,k) * d
   60   continue
c
        do 61 j=1,n
          do 61 i=k,n2
            temp(i,j) = h(i) * ac(k,j)
   61   continue
c
        do 62 j=1,n
          do 62 i=k,n2
            ac(i,j) = ac(i,j) - temp(i,j)
   62   continue
c
        do 63 i=k,n2
          ac(i,k) = h(i)
   63   continue
   10 continue
c
      do 12 j=1,n
        do 11 i=1,n
          a(i,j) = ac(i+n,j)
   11   continue
   12 continue
c
c
 3456 format(/,10x,' determinant of a =',e15.4,/,10x,'singular matrix',
#ifndef a_32_bit_machine
     *       ' a - determinant is less than 1.0e-40., stop')
#endif
#ifdef a_32_bit_machine
     *       ' a - determinant is less than 1.0e-30., stop')
#endif
 3457 format(/,10x,' dimension of matrix a is greater than nmax. it ',
     *       'can not be inverted. change the parameter statement in',
     *       ' subr. matinv.')
c
      return
      end
