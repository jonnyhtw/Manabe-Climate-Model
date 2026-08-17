      common / psbcom2 /   tdty(mx2,jx,kx,maxcpu)
#ifdef rad
      common / psbcom2 / rdteny(mx2,jx,kx,maxcpu)
#endif
      common / psbcom2 /    psy(mx2,jx, 1,maxcpu)
#ifndef grdmx
      common / psbcom2 / qmxdty(mx2,jx,kx,maxcpu)
#endif
c
c  Definitions:
c qmxdty, tdty, rdteny, are tendencies from slabb which are
c stored as a function of iy so that they can be summed in
c such an order to make the answers (checksums) reproduce
c when the code is multitasked.
