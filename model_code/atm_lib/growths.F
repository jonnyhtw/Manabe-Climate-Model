c
      subroutine growths
#if defined mixed_layer && defined ice_dynamics
c
c calculate thermodynamic growth of two-level thickness model
c based on the Hibler [1979] ice model
c
#include "parm.h"
c
#include "param.h"
c
c location of heat flux adjustment
c
#include "bcs.h"
c
c location of idto which is length (sec) of atmospheric leg
c
#include "ctlblk.h"
c
c input surface temperature to be used in land mask calculation
c
#include "grdfld.h"
c
c input oceanic heat flux and snowfall
c
#include "qa.h"
c 
c location of liquid equivalent ice and ocean temperature
c
#include "surfo.h"
c
#include "varice.h"
c
#include "variceco.h"
c
#ifdef icecpld
c
c location of first layer depth of ocean model
c
#include "mixdepth.h"
c
#endif
c
c hint is initial ice thickness entering thermodynamic calculation
c
      dimension hint(ix,il)
c
c armin is minimum compactness (computational purpose)
c one is maximum compactness
c
      data armin /0.2/, one /1./
c
c rhoic is density of sea ice (in g/cm^3)
c alhfice is latent heat of fusion of sea ice (in cal/g)
c
      data rhoic /1./, alhfice /80./
c     data rhoic /0.92/, alhfice /79.67/
c
c alhfrz is volumetric heat of fusion of water (in cal/cm^3)
c
      data alhfrz /80./
c
c cminsec is conversion factor from minute to sec
c
      data cminsec /60./
c
c co is ocean volumetric heat capacity (cal/cm^3/K) 
c olandtm is temperature of ocean over land points
c
      data co /1./, olandtm /290./
c
#ifndef icecpld
c
c depml is mix-layer depth (cm)
c
      data depml /5000./
#endif
c
c time step length (in sec) in thermodynamic calculation of ice model
c
c     data dtice /86400./
c
c cyear is computational days in a year
c
      data cyear /365./ 
c
c h0 is demarcation thickness (in cm) between thick and thin ice
c
      h0north = 50.
      h0south = 50.
c
c set up mixed-layer depth (in cm)
c
#ifdef icecpld
c
c use first layer depth of ocean model
c
      depml = dzti(1)
c
#endif
      do 9 j = 1, il
      do 9 i = 1, ix
c
c set up mask and reset oceanic heat flux
c
         if (ts(i,j) .gt. 0.0) then
            heffmg(i,j) = 0.
         else
            heffmg(i,j) = 1.
         end if
c
c get current value of heat flux adjustment
#ifdef q_flux
         qbot(i,j) = bcqflx(i,j)
#else
         qbot(i,j) = 0.
#endif
c
c temporary store liquid equivalent ice 
c
         yceol(i,j) = yce(i,j)
c
c record initial ice thickness used for calculation of heat flux to ocean
c
         hint(i,j) = heff(i,j,1)
c
   9  continue
      do 11 j = 1, il
  11  continue
c
c store heat fluxes (in cal/cm^2/min) from atmosphere over ice and water
c
      do 12 j = 1, il
      do 12 i = 1, ix
         qficec(i,j) = qfice(i,j)
         qfowc(i,j) = qfow(i,j)
  12  continue
c
c calculate total heat flux (in cal/cm^2/min) from atmosphere
c
      do 14 j = 1, il
      do 14 i = 1, ix
       	 htemp = heff(i,j,2)*100000.
       	 ar(i,j) = amin1(area(i,j,2),htemp)
         qflc(i,j) = qfice(i,j)*ar(i,j) + qfow(i,j)*(1.-ar(i,j))
c        print *, '----------------------------'
c        print *, i, j, qflc(i,j), qao(i,j), qfice(i,j), qfow(i,j)
c        print *, i, j, heff(i,j,2), ar(i,j), area(i,j,2)
c        print *, '----------------------------'
  14  continue
c     stop
c
c checking growth rate
c
c     do 27 j = 1, il
c     do 27 i = 1, ix
c        htemp = heff(i,j,2)*10000.
c        ar(i,j) = amin1(area(i,j,2),htemp)
c        sice(i,j) = qfice(i,j)*ar(i,j) + qfow(i,j)*(1.-ar(i,j))
c        write (6,26) sice(i,j), qao(i,j), heff(i,j,1), ar(i,j)
c  26    format (4f16.5)
c  27 continue
c     stop
c
      do 19 j = 1, il
      do 19 i = 1, ix
c
c put minimum on ice compactness used in compactness growth 
c
      area(i,j,2) = amax1(area(i,j,2),armin)
c
c change compactness to zero if ice thickness is zero
c
      htemp = heff(i,j,2)*10000.
      ar(i,j) = amin1(area(i,j,2),htemp)
c
c account for oceanic heat flux from bottom (in cal/min/cm^2)
c
      qfowi(i,j) = -qbot(i,j) - qfow(i,j)
      qficei(i,j) = -qbot(i,j) - qfice(i,j)
c
c account of snowfall minus evaporation (cm/timestep)
c
      qfowi(i,j) = qfowi(i,j) + sno(i,j)*alhfrz*cminsec
      qficei(i,j) = qficei(i,j) + sno(i,j)*alhfrz*cminsec
c
c convert to ice thickness (in cm) per timestep (in min)
c
      fow(i,j) = qfowi(i,j)/rhoic/alhfice*cliqsol
      fice(i,j) = qficei(i,j)/rhoic/alhfice*cliqsol
  19  continue
c
      do 29 j = 1, il
      do 29 i = 1, ix
c
c ice thickness total growth rate
c
      sice(i,j) = fice(i,j)*ar(i,j) + fow(i,j)*(1.-ar(i,j))
  29  continue
c
c separate ice growth rates (in cm/timestep)
c
      do 36 j = 1, il
      do 36 i = 1, ix
         siceic(i,j) = fice(i,j)*ar(i,j)*dticet/cminsec
         siceow(i,j) = fow(i,j)*(1.-ar(i,j))*dticet/cminsec
  36  continue
c
      do 39 j = 1, il
      do 39 i = 1, ix
c
c distinguish between northern and southern hemisphere
c
      if (j .le. iy) then
         h0 = h0north
      else
         h0 = h0south
      end if
c
c ice compactness total growth rate
c
      sarea(i,j) = 0.
      if (fow(i,j) .gt. 0.) then
         sarea(i,j) = sarea(i,j) + fow(i,j)/h0*(1.-area(i,j,2))
      end if
      if (sice(i,j) .lt. 0.) then
         htemp = heff(i,j,1) + 0.001
         sarea(i,j) = sarea(i,j) + area(i,j,2)/2./htemp*sice(i,j)
      end if
  39  continue
c
      do 49 j = 1, il
      do 49 i = 1, ix
c
c convert to ice thickness (dtice in sec)
c
      area(i,j,1) = area(i,j,1) + sarea(i,j)*dticet/cminsec
      heff(i,j,1) = heff(i,j,1) + sice(i,j)*dticet/cminsec
c
c ice (in m) or area growth in one ice time step
c
      sareac(i,j) = sarea(i,j)*dticet/cminsec
      sicec(i,j) = sice(i,j)*dticet/cminsec
c
c acccount for heat (in negative ice) in the mixed layer to melt sea ice
c
      hml =  co*depml*(toint(i,j)-tfro)/rhoic/alhfice*cliqsol
      heff(i,j,1) = heff(i,j,1) - hml
c
c get available heat (in negative ice) in the mixed layer for history
c
      hmlc(i,j) = hml
c
c adjust for additional heat (in negative ice) pass on to mixed layer
c
      if (heff(i,j,1) .lt. 0.) then
         heatml(i,j) = -heff(i,j,1)
         heff(i,j,1) = 0.
      else
         heatml(i,j) = 0.
      end if
  49  continue
c
c calculate heat and water fluxes entering the ocean 
c
      do 59 j = 1, il
      do 59 i = 1, ix
c
c heat flux to ocean in ice unit
c
         hocean = sice(i,j)*dticet/cminsec - (heff(i,j,1)-hint(i,j))
c        hocean = hmlc(i,j) - heatml(i,j)
c
c convert heat flux to ocean to units of cal/cm^2/min
c
         qoceani(i,j) = -hocean*rhoic*alhfice/cliqsol/dticet*cminsec
c
c        qtol = 1.e-10
c        qoabs = abs(qoceani(i,j))
c        if (qoabs .lt. qtol) then
c           qoceani(i,j) = 0.
c        end if
c
c water budget includes snow and sublimation
c
         wocean = heff(i,j,1) - hint(i,j) - sno(i,j)*cliqsol*dticet
c
c convert water flux to ocean to units of g/cm^2/min
c
         poceani(i,j) = -wocean*rhoic/cliqsol/dticet*cminsec
c
  59  continue
c
c     if (ktau .eq. 921540) then
c        stop
c     end if
c
      do 61 j = 1, il
      do 61 i = 1, ix
c
c calculate ocean mixed layer temperature within a time step
c
      toint(i,j) = toint(i,j) + qoceani(i,j)*dticet/cminsec/co/depml
      toint(i,j) = toint(i,j)*heffmg(i,j)+olandtm*(1.-heffmg(i,j))
c     print *, i, j, toint(i,j)
c     to(i,j) = tox
c     to(i,j) = tfro + heatml(i,j)/co/depml*rhoic*alhfice/cliqsol
c     to(i,j) = to(i,j)*heffmg(i,j) + olandtm*(1.-heffmg(i,j))
c
c set compactness to zero when there is no ice
c
      areatemp = heff(i,j,1)/0.0001
      area(i,j,1) = amin1(area(i,j,1),areatemp)
c
c convert solid sea ice to liquid equivalent ice
c
      yce(i,j) = heff(i,j,1)/cliqsol*heffmg(i,j)
c     yce(i,j) = heff(i,j,1)*heffmg(i,j)
c
c put upper and lower limit on compactness
c
      area(i,j,1) = amin1(one,area(i,j,1))
      area(i,j,1) = amax1(armin,area(i,j,1))
  61  continue
c
      do 69 j = 1, il
      do 69 i = 1, ix
c
c prepare ice thickness and compactness for next time step
c
      heff(i,j,1) = heff(i,j,1)*heffmg(i,j)
      area(i,j,1) = area(i,j,1)*heffmg(i,j)
      heff(i,j,2) = heff(i,j,1)
      area(i,j,2) = area(i,j,1)
c
c set up ice thickness and compactness use in surface heat budget
c
      area(i,j,2) = amax1(armin,area(i,j,2))
c
c for no lead calculations
c
c     if (heff(i,j,2) .gt. 0.) then
c        area(i,j,2) = 1.
c        area(i,j,1) = 1.
c     end if
      hice(i,j) = heff(i,j,2)/area(i,j,2)
c
c for no lead calculatins
c
c     hice(i,j) = heff(i,j,2)
c     hice(i,j) = heff(i,j,2)*cliqsol/area(i,j,2)
c     if (area(i,j,2) .eq. armin) then
c        hice(i,j) = 0.
c     end if
  69  continue
c
#endif
      return
      end
