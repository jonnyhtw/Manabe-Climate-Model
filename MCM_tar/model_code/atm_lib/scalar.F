      subroutine scalar(ftvar,s,kmx)
c
c  this subroutine computes the sum and difference of the partial dot
c    products of the legendre polynomials "colyl" and a spherical
c    harmonic variable (svble) by the equations below, where the first
c    partial product involves the odd n=1 to jx indices and the
c    second, the even indices:
c
c  ftvar(m,k,1 or 2) = summation of colyl(m,n)*svble(m,n,k)
c    for n=1,jx,2  plus (or minus) the summation of
c    colyl(mn,n)*svble(m,n,k) for n=2,jx,2; k=1,kmx.
c
c  in other words, the sum is stored at the third index = 1
c    and the difference at index = 2.
c
c  see notes of january 27, 1982 and march 3, 1984.
c
#include "parm.h"
#include "parmjx.h"
#include "clnpol.h"
c
      dimension ftvar(mx2,kmx,nhem), polvar(mx4),c(mx2,jx)
      dimension s(mx2,jx,kmx)
c
      equivalence(colyl,c)
c
c  below are statement functions for the various resolutions of the
c   model.  statements functions are used because they are very fast
c   on the cray ymp.
c
#ifdef r15
      f15(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)
#endif
#ifdef r21
      f21(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)
#endif
#ifdef r30
      f30(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)+s(i,23,k)*c(i,23)+
     4         s(i,25,k)*c(i,25)+s(i,27,k)*c(i,27)+s(i,29,k)*c(i,29)+
     5         s(i,31,k)*c(i,31)
#endif
#ifdef r63
      f63(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)+s(i,23,k)*c(i,23)+
     4         s(i,25,k)*c(i,25)+s(i,27,k)*c(i,27)+s(i,29,k)*c(i,29)+
     5         s(i,31,k)*c(i,31)+s(i,33,k)*c(i,33)+s(i,35,k)*c(i,35)+
     6         s(i,37,k)*c(i,37)+s(i,39,k)*c(i,39)+s(i,41,k)*c(i,41)+
     7         s(i,43,k)*c(i,43)+s(i,45,k)*c(i,45)+s(i,47,k)*c(i,47)+
     8         s(i,49,k)*c(i,49)+s(i,51,k)*c(i,51)+s(i,53,k)*c(i,53)+
     9         s(i,55,k)*c(i,55)+s(i,57,k)*c(i,57)+s(i,59,k)*c(i,59)+
     1         s(i,61,k)*c(i,61)+s(i,63,k)*c(i,63)
#endif
#ifdef t21
      f21(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)
#endif
#ifdef t42
      f42(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)+s(i,23,k)*c(i,23)+
     4         s(i,25,k)*c(i,25)+s(i,27,k)*c(i,27)+s(i,29,k)*c(i,29)+
     5         s(i,31,k)*c(i,31)+s(i,33,k)*c(i,33)+s(i,35,k)*c(i,35)+
     6         s(i,37,k)*c(i,37)+s(i,39,k)*c(i,39)+s(i,41,k)*c(i,41)+
     7         s(i,43,k)*c(i,43)
#endif
#ifdef t84
      f84(i,k)=s(i, 1,k)*c(i, 1)+s(i, 3,k)*c(i, 3)+s(i, 5,k)*c(i, 5)+
     1         s(i, 7,k)*c(i, 7)+s(i, 9,k)*c(i, 9)+s(i,11,k)*c(i,11)+
     2         s(i,13,k)*c(i,13)+s(i,15,k)*c(i,15)+s(i,17,k)*c(i,17)+
     3         s(i,19,k)*c(i,19)+s(i,21,k)*c(i,21)+s(i,23,k)*c(i,23)+
     4         s(i,25,k)*c(i,25)+s(i,27,k)*c(i,27)+s(i,29,k)*c(i,29)+
     5         s(i,31,k)*c(i,31)+s(i,33,k)*c(i,33)+s(i,35,k)*c(i,35)+
     6         s(i,37,k)*c(i,37)+s(i,39,k)*c(i,39)+s(i,41,k)*c(i,41)+
     7         s(i,43,k)*c(i,43)+s(i,45,k)*c(i,45)+s(i,47,k)*c(i,47)+
     8         s(i,49,k)*c(i,49)+s(i,51,k)*c(i,51)+s(i,53,k)*c(i,53)+
     9         s(i,55,k)*c(i,55)+s(i,57,k)*c(i,57)+s(i,59,k)*c(i,59)+
     1         s(i,61,k)*c(i,61)+s(i,63,k)*c(i,63)+s(i,65,k)*c(i,65)+
     2         s(i,67,k)*c(i,67)+s(i,69,k)*c(i,69)+s(i,71,k)*c(i,71)+
     3         s(i,73,k)*c(i,73)+s(i,75,k)*c(i,75)+s(i,77,k)*c(i,77)+
     4         s(i,79,k)*c(i,79)+s(i,81,k)*c(i,81)+s(i,83,k)*c(i,83)+
     5         s(i,85,k)*c(i,85)
#endif
c
c  compute products of cloned polynomials "colyl" and the variable
c    "svble" at level "k" for each m=1,mx2 and n=1,jx combination.
c    using the statement functions above.  also note that two
c    rows are computed to get longer vectors (the mx4).  finially
c    if jx is odd, the last jx spectral variable is summed.
c
      do 12 k=1,kmx
#ifdef r15
         do 20 mn=1,mx4
          polvar(mn)=f15(mn,k)
20       continue
         do 21 mn=1,mx2
          polvar(mn)=polvar(mn)+s(mn,17,k)*c(mn,17)
21       continue
#endif
c
#ifdef r21
         do 22 mn=1,mx4
          polvar(mn)=f21(mn,k)
 22      continue
         do 23 mn=1,mx2
          polvar(mn)=polvar(mn)+s(mn,23,k)*c(mn,23)
 23      continue
#endif
c
#ifdef r30
         do 24 mn=1,mx4
          polvar(mn)=f30(mn,k)
 24      continue
#endif
#ifdef r63
         do 28 mn=1,mx4
          polvar(mn)=f63(mn,k)
 28      continue
         do 29 mn=1,mx2
          polvar(mn)=polvar(mn)+s(mn,65,k)*c(mn,65)
 29      continue
#endif
#ifdef t21
         do 22 mn=1,mx4
          polvar(mn)=f21(mn,k)
 22      continue
         do 23 mn=1,mx2
          polvar(mn)=polvar(mn)+s(mn,23,k)*c(mn,23)
 23      continue
#endif
c
#ifdef t42
         do 24 mn=1,mx4
          polvar(mn)=f42(mn,k)
 24      continue
#endif
c
#ifdef t84
         do 24 mn=1,mx4
          polvar(mn)=f84(mn,k)
 24      continue
#endif
c
c  compute sum and difference fourier components corresponding to
c    northern and southern hemispheres respectively from sums of the
c    odd and even n=1 to jx indices and store in ftvar.
       do 16 m=1,mx2
        ftvar(m,k,1)=polvar(m)+polvar(m+mx2)
  16   continue
        if(nhem.eq.2) then
         do 18 m=1,mx2
          ftvar(m,k,2)=polvar(m)-polvar(m+mx2)
 18      continue
        endif
 12   continue
c
      return
c
      entry prescl
c     ===== ======
      return
      end
