c
      common / psacom2 / vordty(mx2,jx,kx,maxcpu)
      common / psacom2 / divdty(mx2,jx,kx,maxcpu)
      common / psacom2 /    eky(mx2,jx,kx,maxcpu)
      common / psacom2 /   tdty(mx2,jx,kx,maxcpu)
#ifndef grdmx
      common / psacom2 / qmxdty(mx2,jx,kx,maxcpu)
#endif
      common / psacom2 /    psy(mx2,jx, 1,maxcpu)
      common / psacom2 /  psdty(mx2,jx, 1,maxcpu)
#ifdef dirnal
      common / psacom2 /  qdoty(mx2,jx,kx,maxcpu)
      common / psacom2 /  somgy(mx2,jx,kx,maxcpu)
#endif
c
c  Definitions:
c vordty, tdty, and friends are tendencies from slaba which are
c stored as a function of iy so that they can be summed in
c such an order to make the answers (checksums) reproduce
c when the code is multitasked.

