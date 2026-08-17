c====================== include file "cvmix.h" =========================
c
c     variables used for vertical diffusion
c
c     inputs: (set through namelist)
c
c     fkph   = base vertical diffusion coefficient (cm**2/sec)
#ifdef constvmix
c     vdcsci, vdcscj, vdcsck = multiplicative scale factors applied
c              to "fkph" for "constvmix" cases.  these scale factors
c              allow for variations in the x, y & z directions,
c              respectively (vertical diffusivity for tracers is thus
c              allowed to vary in grid space, in a predetermined,
c              non-time dependent fashion)
#endif
c     fkpm   = base vertical viscosity coefficient (cm**2/sec)
#ifdef constvmix
c     vvcsci, vvcscj, vvcsck = multiplicative scale factors applied
c              to "fkpm" for "constvmix" cases, allowing vertical
c              viscosity to vary in grid space, in a predetermined,
c              non-time dependent fashion
#endif
c     bvdc   = background vertical diffusion constant
c     bvvc   = background vertical viscosity constant
c     vvclim = vertical viscosity coefficient limit
c     vdclim = vertical diffusion coefficient limit
c
c     aidif  = coefficient for implicit time differencing for
c              vertical diffusion. aidif=1 gives the fully implicit
c              case. aidif=0 gives the fully explicit case
c              note: not used unless "implicitvmix" or "isopycmix"
c                    is enabled
#ifdef constvmix
# ifdef blvmix
c            
c     to set "vdc" as a function of depth using the same functional
c     relationship as in Bryan & Lewis (JGR, pp 2503-2517, 1979) the
c     following four variables are set through namelist.
c     (the variations are set through "vdcsck(k)" in program ocean)
c
c     afkph  = vertical diffusion coefficient at depth "zfkph" and the
c              midpoint of the "vdc" range (cm**2/sec)
c     dfkph  = range of the vertical diffusion coefficient between the
c              surface and great depths (cm**2/sec) i.e., "vdc" values
c              at the surface will approach "afkph"-("dfkph"/2) and
c              at great depths "vdc" approaches "afkph"+("dfkph"/2)
c     sfkph  = scaling for the sharpness of the crossover in "vdcsck"
c              near "zfkph. (cm**-1)
c              the smaller the value of "sfkph", the sharper
c              the vertical gradient of "vdcsck" near  "zfkph".
c              (one can think of 1/sfkph as the vertical extent of the
c              transition zone between thermocline and deep water
c              vertical mixing zones)
c     zfkph  = depth at which transition between surface and deep fkph
c              takes place and the depth where "vdc"="afkph" (cm)
c
c         Note: to duplicate the Bryan & Lewis values, blkdta sets
c               afkph=0.8, dfkph=1.05, sfkph=4.5e-5 & zfkph=2500.0e2
# endif
#endif
c
c
c     derived quantities:
c
c     vvc  = vertical viscosity coeff
c     vdc  = vertical diffusion coeff 
c            (vvc & vdc are calculated at the zw(k) depths)
c
#ifndef tcvmix
# ifndef multitasking
#  define task
# endif
cSGI
cc    task common /cvmix/ vvc(imt,km), vdc(imt,km)
           common /cvmix/ vvc(imt,km), vdc(imt,km)
#endif
c
      common /cvmix1/ fkph, fkpm, bvdc, bvvc, vvclim, vdclim
#ifdef constvmix
     $,               vdcsci(imt), vdcscj(jmt), vdcsck(km)
     $,               vvcsci(imt), vvcscj(jmt), vvcsck(km)
# ifdef blvmix
     $,               afkph, dfkph, sfkph, zfkph
# endif
#endif
     $,               aidif
c
