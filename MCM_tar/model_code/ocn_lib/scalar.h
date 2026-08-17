c====================== include file "scalar.h" ========================
c
c     various scalar quantities:
c
c     dtts   = time step for density & tracers (in seconds)
c     dtuv   = time step for baroclinic velocity (in seconds)
c     dtsf   = time step for barotropic velocity (in seconds)
c     c2dtts = 2*dtts   
c     c2dtuv = 2*dtuv   
c     c2dtsf = 2*dtsf 
c     area   = surface area of ocean (cm**2)  
c     volume = volume of ocean (cm**3)
c     omega  = earth's rotation rate (radians/sec)
c     radius = earth's radius (cm)
c     grav   = earth's gravitational acceleration (cm/sec**2)  
c     cdbot  = bottom drag coefficient
c     ncon   = number of  passes through convective code in tracer
c
cSGI scalar chg to scalro to avoid conflict with atm scalar
      common /scalro/ dtts, dtuv, dtsf, c2dtts, c2dtuv, c2dtsf
     $,               area, volume, omega, radius, grav, cdbot
      common /scalri/ ncon
c
c     bring in non dimensional constants
c
#include "ndcon.h"
c
