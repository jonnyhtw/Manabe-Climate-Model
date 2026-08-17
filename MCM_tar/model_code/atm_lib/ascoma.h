c   Atmosphere surface boundary condition variables:
c
c   bca(i,j,1,m)   ocean surface temperature   degrees celsius  
c   bca(i,j,2,m)   sea ice                     cm
c   smask(i,j)     land-sea mask               1=sea; 0=land
c
#if defined cpldleg1 || defined cpldleg2
# if defined seanal || defined seasonal
      common / ascoma / bca(ix,il,2,12), smask(ix,il)
# else
      common / ascoma / bca(ix,il,2,1), smask(ix,il)
# endif
#endif
#ifdef cpldleg3
      common / ascoma / bca(ix,il,2,1), smask(ix,il)
#endif
