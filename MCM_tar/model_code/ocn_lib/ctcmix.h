c====================== include file "ctcmix.h"=========================
c
c    variables for turbulence closure scheme
c
c    inputs:
c           initialized in include "dtcmix.h" & namelist
c
c    derived quantities:
c    tke     = turbulent kinetic energy tendency
c    advq2   = advection of tke
c    hdq2    = horizontal diffusion of tke
c    sprod   = shear production
c    bprod   = buoyancy production
c    dissp   = dissipation of tke
c
c
#ifndef multitasking
# define task
#endif
c
cSGI
cc    task common   /ctcmix/
           common   /ctcmix/
     $               bprod(imt,km), sprod(imt,km), dissp(imt,km)
     $,             tke(imt,km), bprol(imt,km)
     $,             advq2(imt,km), fvsq(imt,km)
     $,             hdq2(imt,km)
     $,             el(imt,km,nslabs), elmax(imt,km)
     $,             eeq(imt,km), ffq(imt,km)
     $,             gmp1(imt,km), gh(imt,km)
     $,             sm(imt,km), sh(imt,km)
#ifdef leq
     $,             tkel(imt,km)
     $,             advq2l(imt,km)
     $,             hdq2l(imt,km)
#endif
      common /ctcmx1/
     $              a1, b1, a2, b2, cc1, e1, e2, e3
     $,             vdqlim
     $,             alpha, vk, sq
     $,             rhoh2o
     $,             small, big, deps
     $,             p032, p1, p154, p4, p48
     $,             c15, c100
     $,             jpr1,ipr1,jpr2,ipr2
c
