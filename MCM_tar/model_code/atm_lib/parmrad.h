c
c   kx=number of atmospheric model layers
c   nb=number of spectral bands
c   ng=number of absorbing gases
c   ncc=number of clouds in the simulation plus one dummy cloud
c
      parameter (k2=kx+2)
      parameter (nb=19, ng=3, ngp=ng +1)
      parameter (kpx2=kp*2, kpx2m=kpx2-1, kpx2p=kpx2+1)
#ifdef fixcld
      parameter (ncc=4, lx=ncc+1, nc=ncc+1)
#endif
#ifdef prdcld
      parameter (lx=k2, nc=k2)
#endif
