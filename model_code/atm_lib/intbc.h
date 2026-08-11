#if defined cpldleg3 || defined cpldleg1 || defined ice_dynamics
      common /intbc/ oyv(jmt), wto(jmt), wta(il), ay(il)
c
c     oyv is latitudes on the ocean u,v grid in radians
c     wto is area weights for the latitudes on the ocean t,s grid
c     wta is area weights for the latitudes on the atmos grid
c     ay  is latitudes on the atmos grid in radians
c
#endif
