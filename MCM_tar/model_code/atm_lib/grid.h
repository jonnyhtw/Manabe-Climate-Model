c
c  *** warning  ****
c  the order of the variables in the common block / grid /
c   is linked to the variable order in the common blocks
c   / ift / and / fftstr /.  change the order at your own
c   risk.  adding variables to the end of the common blocks
c   should be safe.
c
      common / grid  / gvor   (ix,kx,nhem), gdiv   (ix,kx,nhem)
      common / grid  / gucos  (ix,kx,nhem), gvcos  (ix,kx,nhem)
      common / grid  / cgt    (ix,kx,nhem), gps    (ix,   nhem)
      common / grid  / gpsdx  (ix,   nhem), gpsdy  (ix,   nhem)
      common / grid  / gzs    (ix,   nhem)
      common / grid  / gqmix  (ix,kx,nhem)
#ifdef rad
      common / grid  / grad   (ix,kx,nhem)
#endif
#ifdef dirnal
      common / grid  / omg    (ix,kx,nhem)
#endif
c
      common / grid  / agrid  (ix,kx,nhem)
      common / grid  / gvbrps (ix,   nhem)
      common / grid  / bgrid  (ix,kx,nhem)
      common / grid  / gteq   (ix,kx,nhem)
      common / grid  / gtu    (ix,kx,nhem), gtv    (ix,kx,nhem)
      common / grid  / gekcos (ix,kx,nhem)
      common / grid  / gpseq  (ix,   nhem)
      common / grid  / gqmxeq (ix,kx,nhem)
      common / grid  / gqmixu (ix,kx,nhem), gqmixv (ix,kx,nhem)
