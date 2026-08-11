      common / trans  / transm(ix,0:kx,0:kx), transz(ix,0:kx)
      common / trans  / absrbx(ix,0:kx,2)

#ifndef unitsk
      common / trnstr / glbtrs(ix,0:kx,0:kx,il)
      common / trnstr / glbtrz(ix,0:kx,     il)
      common / trnstr / glbabs(ix,0:kx,2   ,il)
#endif

#ifdef dirnal
      common / trans  / rheat(kp,ix), sing(ix), song(ix)
#endif
