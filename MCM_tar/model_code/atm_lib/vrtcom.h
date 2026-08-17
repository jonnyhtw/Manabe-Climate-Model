c
      common / vrtcom / ks, kv, ksm, kvp, kvm, kxmv
      common / vrtcom / kll
c  
c  All these variables are used in vertical diffusion (VERTDF).
c
c  Definitions:
c     ks   = number of levels to be diffused
c     ksm  = ks - 1
c     kv   = kx - ks + 1, index of top level to be diffused
c     kvp  = kv + 1
c     kvm  = kv - 1
c     kxmv = kx - kv
c     kll  = the top of the boundary layer for the gravity wave comp.
