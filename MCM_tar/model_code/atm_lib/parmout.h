c set conditional moisture parameters
#ifdef moist
      parameter( ivf=4, icv=5 )
#endif
#ifdef dry
      parameter( ivf=3, icv=2 )
#endif
c
c set radiation parameter
#ifdef rad
# ifndef cloud_forcing
      parameter ( irm=6 )
# endif
# ifdef cloud_forcing
      parameter ( irm=12 )
# endif
#endif
#ifndef rad
      parameter ( irm=0 )
#endif
c
c irm   = the # of radiation variables output
c ivf   = # of vertically diffused variables output
c icv   = # of convection variables output
c 

