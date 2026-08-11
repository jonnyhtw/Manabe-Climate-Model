c====================== include file "chmix.h" =========================
c
c     horizontal mixing coefficients
c
c     ah  = lateral eddy diffusivity (cm**2/sec)
#if defined consthmix || defined isopycmix
c     hdcsck = multiplicative scale factors applied to "ah" for
c              "consthmix" & "biharmonic" cases and to "ahisop" for 
c              "isopycmix" cases.
c              these scale factors allow for variations in the z
c              direction only (horizontal diffusivity for tracers is
c              thus allowed to vary in k-level grid space, in a
c              predetermined, non-time dependent fashion)
#endif
c     am  = lateral eddy viscosity (cm**2/sec)
#ifdef consthmix
c     hvcsck = multiplicative scale factors applied to "am" for
c              "consthmix" cases, allowing horizontal viscosity to vary
c              in k-level grid space, in a predetermined, non-time
c              dependent fashion
#endif
c     aq  = lateral turbulent kinetic energy diffusivity (cm**2/sec)
c
      common /chmix1/ am, ah, asubh, asubm
#ifdef tcvmix
     $,               aq
#endif
#if defined consthmix || defined isopycmix
     $,               hdcsck(km), hvcsck(km)
#endif
     $,               bbu(jmt), ccu(jmt), ddu(jmt), ggu(jmt), hhu(jmt)
     $,               bbt(jmt), cct(jmt), ddt(jmt)
c
