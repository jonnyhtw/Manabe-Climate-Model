c
      integer adsk, bdsk
c
      common / iohist / itap, irs, irsf
      common / iohist / adsk, bdsk
      common / iohist / kdsk
      common / iohist / iunt, ibd
      common / iohist / iunitr,iunitw
      common / iohist / ialb, idsk
      common / iohist / izsvar
      common / iohist / irunof
      common / iohist / intgl
      common / iohist / ictocn
      common / iohist / iycthk, iycrst
      common / iohist / isstyc, iflux, iflxan, iflxac
c
c  itap   =    input; spectral atm restart tape
c  irs    =   output; spectral atm history tape
c  irsf   =   output; spectral atm restart tape
c  adsk   =   output; /aux/ data from PASSA, SLABA
c  bdsk   =   output; /aux/ data from PASSB, SLABB
c  kdsk   = internal; used to store solar radiation tendencies (diurnal)
c  iunt   =    input; mixed and fixed ocean restart tape
c  ibd    =    input; mixed and fixed ocean restart tape
c  iunitr =    input; fixed SST and YCE data
c  iunitw = internal; fixed SST and YCE data fourier comp.
c  ialb   =    input; surface bare soil albedo field
c  idsk   = internal; used to store longwave trans. func. from rad code
c  izsvar =    input; subgrid topog variance, used for grav wave drag
c  irunof =    input; indexes for moving runoff from land to ocean(cpld)
c  intgl  =in/output; unit # for time series file of model integrals
c  ictocn = internal; ocn control block used in cpld models
c  isstyc =in/output; instan. sst and yce fields from ocean (cpld)
c  iflux  = internal; fluxes from atm to ocn on ocn grid (cpld)
c  iflxan = internal; fluxes from atm to ocn on ocn grid (cpld) ann model only
c  iflxac =   output; fluxes from atm to ocn on ocn grid (cpld) leg 1 only
c  iycthk =in/output; ice thickness for cav. fluid sea ice model
c  iycrst =    input; data needed for restart of cav. fluid sea ice model
c
c
