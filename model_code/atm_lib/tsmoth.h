#ifdef consrv_ps
      double precision alpha, om2al
c
#endif
      common / tsmoth / alpha     , om2al
c
c  alpha is the Robert filter constant, set in atmos or namelist
c  om2al is (1-2*alpha); these are used in linear
c

