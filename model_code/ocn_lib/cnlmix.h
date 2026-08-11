c====================== include file "cnlmix.h" ========================
c
c     arrays used for non-linear horizontal viscosity
c
c     dten     = tension rate of strain
c     dtens    = tension rate of strain at southern row
c     dshear   = shear rate of strain
c     dshers   = shear rate of strain at southern row
c     dpure    = pure rate of deformation
c     dpures   = pure rate of deformation at southern row
c     dpuret   = pure rate of deformation
c     dkmc     = non-linear equivalent of am
c     effk     = effective  wavenumber in x direction
c     effky    = effective  wavenumber in y direction
c     bbtj     = nonlinear diffusion coefficient for tracers
c     cctj     = nonlinear diffusion coefficient for tracers
c     ddtj     = nonlinear diffusion coefficient for tracers
c
c
c
      common /cnlmix/ dten(imt,km),   dtens(imt,km) 
     $,               dshear(imt,km), dshers(imt,km)
     $,               dpure(imt,km),  dpures(imt,km), dpuret(imt,km)
     $,               dkmc(imt,km)
     $,               dkhcx(imt,km), dkhcy(imt,km)
     $,               bbtj(imt,km), cctj(imt,km), ddtj(imt,km)
     $,               effk(imt), effky(jmt)
c
