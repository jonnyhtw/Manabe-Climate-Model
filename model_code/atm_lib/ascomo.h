c   Ocean surface boundary condition variables:
c
c   bco(i,j,1,m)   heat flux                         cal/cm**2/sec
c   bco(i,j,2,m)   negative freshwater flux (e-p-r)  cm/sec
c
c   if defined icecpld then (cfi)
c   bco(i,j,3,m)   salt flux from ice code
c
c   bco(i,j,3,m)   snow flux (snow + glacier runoff)   cm/sec
c   bco(i,j,4,m)   wind stress, x-component          dynes/cm**2
c   bco(i,j,5,m)   wind stress, y-component          dynes/cm**2
#ifdef windbc
c   bco(i,j,6,m)   wind speed                        cm/sec      
c   bco(i,j,7,m)   (wind speed)**3  see vertdf       (cm/sec)**3   
c
c   nbo = 7 if wind speed is included
      parameter (nbo = 7)
#else
      parameter (nbo = 5)
#endif
c
#ifdef cpldleg1
ccc   common / ascomo / is not needed for leg1 of coupled model
#endif
c
#ifdef cpldleg2
# if defined seanal || defined seasonal
      common / ascomo / bco(imt,jmt,nbo,2)
# else
      common / ascomo / bco(imt,jmt,nbo,1)
# endif
#endif
#if defined cpldleg3 || defined icecpld
      common / ascomo / bco(imt,jmt,nbo,1)
#endif
#if defined cpldleg2 || defined cpldleg3
      common / ascomo / tbco(nbo), sumt, sumu
#endif
