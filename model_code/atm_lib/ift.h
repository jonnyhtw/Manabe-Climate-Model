c
c  *** warning  ****
c  the order of the variables in the common block / ift /
c   is linked to the variable order in the common block
c   / grid /.  change the order at your own risk.  adding
c   variables to the end of the common blocks
c   should be safe.
c
      common / ift    /  aift(mx,kx,nhem)
      common / ift    / spift(mx,   nhem)
      common / ift    /  bift(mx,kx,nhem)
      common / ift    / tift (mx,kx,nhem)
      common / ift    / tuift(mx,kx,nhem), tvift(mx,kx,nhem)
      common / ift    / ekift(mx,kx,nhem)
      common / ift    / psift(mx,   nhem)
      common / ift    / qeift(mx,kx,nhem)
      common / ift    / quift(mx,kx,nhem), qvift(mx,kx,nhem)
      common / ift    / zaift(mx,kx,nhem), zbift(mx,kx,nhem)
#ifdef dirnal
      common / ift    / qdift(mx,kx,nhem), ogift(mx,kx,nhem)
#endif
