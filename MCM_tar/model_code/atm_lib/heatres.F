c
      subroutine heatres
#if defined restore && defined ice_dynamics
c
c        called from : ocean
c       entry points : none
c subroutines called : none
c
c subroutine to restore heat content to observation values (flux adjustment)
c requires results from subroutine bciol
c
#include "parm.h"
c
c location of observation sst and sea ice
c
#include "bcs.h"
c
c location of time for restart
c
#include "ctlblk.h"
c
c input surface temperature to be used in land mask calculation
c
#include "grdfld.h"
c
c freezing point of sea water
c
#include "param.h"
c
c location of model sst and sea ice
c
#include "surfo.h"
c
#include "varice.h"
c
#include "variceco.h"
c
#ifdef icecpld
c
c depth of first few levels of ocean model
c
#include "mixdepth.h"
#endif
c
      dimension sealnd(ix,il)
c
c co is ocean volumetric heat capacity (cal/cm^3/K)
c olandtm is temperature of ocean over land points
c
      data co /1./, olandtm /290.0/
c
#ifndef icecpld
c
c depml is mix-layer depth (cm)
c
      data depml /5000./
#endif
c
c rhoice is density of sea ice (in g/cm^3)
c alhfice is latent heat of fusion of sea ice (in cal/g)
c
      data rhoic /1./, alhfice /80./
c
c time step length (in sec) of ice model
c
      data dtice /86400./
c
c cminsec is conversion factor from minute to sec
c
      data cminsec /60./
c
#ifdef icecpld
c
c use first level of ocean model as mixed-layer depth
c
      depml = dzti(1)
#endif
c
c rtime is restoring time (in sec); 20 days < rtime < 60 days
c
      rtime = 30.*86400.
c
c set rtime based on Billy's initial experiments
c
      rtime = 10.*86400.
c
c     go to 40
c
c sealnd is landsea mask, 0 for land; 1 for sea
c
      do 10 j = 1, il
      do 10 i = 1, ix
         if (ts(i,j) .gt. 0.0) then
            sealnd(i,j) = 0.
         else
            sealnd(i,j) = 1.
         end if
  10  continue
c
      do 20 j = 1, il
      do 20 i = 1, ix
c
c reset restoring heat flux term
c
         qres(i,j) = 0.
c
         if (sealnd(i,j) .eq. 1.) then
c
c hicm is heat content calculated from the ice dynamics model
c
            hicm = co*depml*(to(i,j)-tfro)
     >             - rhoic*alhfice*heff(i,j,2)/cliqsol
c
c hstar is heat content calculated from observed sst and sea ice
c
            hstar = co*depml*(bcto(i,j)-tfro) 
     >              - rhoic*alhfice*bcyce(i,j)*bcconc(i,j)
c
c qres is restoring heat flux (in cal/cm^2/sec)
c
            qres(i,j) = (hstar-hicm)/rtime
c
c restore heat content 
c
            hnew = hicm + dtice*qres(i,j)
c
c           print *, 'before'
c           print *, i, j, to(i,j), hice(i,j), area(i,j,2)
c
c restoring sst and ice thickness from new heat content
c
            if (hnew .gt. 0.) then
               to(i,j) = hnew/co/depml + tfro
               hice(i,j) = 0.
               area(i,j,2) = 0.
            else
               to(i,j) = tfro
c
c restore ice compactness
c
               ares = (bcconc(i,j)-area(i,j,2))/rtime
c
c areaco is constraint placed on ice compactness for next time step
c
               areaco = area(i,j,2) + dtice*ares
               area(i,j,2) = areaco
c
               if (bcconc(i,j) .eq. 0.) then
c                 areaco = area(i,j,2)
c                 areaco = 0.2
c                 area(i,j,2) = areaco
                  print *, 'bcconc is zero'
                  print *, i,j,area(i,j,2),heff(i,j,2),hice(i,j)
                  print *, ts(i,j),to(i,j),bcto(i,j),bcyce(i,j)
           	  print *, '--------------------------'
c              else
c                 areaco = bcconc(i,j)
c                 area(i,j,2) = bcconc(i,j)
               end if
c
               hice(i,j) = -hnew/rhoic/alhfice/areaco*cliqsol
            end if
c
c           print *, 'after'
c           print *, i, j, to(i,j), hice(i,j), area(i,j,2)
c           print *, i, j, bcto(i,j), bcyce(i,j), bcconc(i,j)
c           print *, '-----------------'
         end if
c
c fixed surface conditions
c
c        hice(i,j) = bcyce(i,j)*cliqsol*sealnd(i,j)
c        area(i,j,2) = bcconc(i,j)*sealnd(i,j)
c        to(i,j) = bcto(i,j)
  20  continue
c
      do 30 j = 1, il
      do 30 i = 1, ix
c
c put land mask on to and hice (set ocean temp to olandtm for land points)
c
         to(i,j) = to(i,j)*sealnd(i,j) + olandtm*(1.-sealnd(i,j))
         hice(i,j) = hice(i,j)*sealnd(i,j)
         area(i,j,2) = area(i,j,2)*sealnd(i,j)
c
c convert restoring heat flux to cal/cm^2/min
c
	 qres(i,j) = qres(i,j)*cminsec*sealnd(i,j)
c
c reset other ice variables after restoration of heat content
c
         heff(i,j,2) = hice(i,j)*area(i,j,2)
         yce(i,j) = heff(i,j,2)/cliqsol
         heff(i,j,1) = heff(i,j,2)
         area(i,j,1) = area(i,j,2)
         ar(i,j) = area(i,j,2)
c        print *, i, j, to(i,j), ts(i,j), sealnd(i,j)
c        print *, hice(i,j), heff(i,j,2), yce(i,j), area(i,j,2)
c        print *, heff(i,j,1), area(i,j,1), ar(i,j), qres(i,j)
c        print *, '------------------------------------------'
c
  30  continue
  40  continue
c
c write out ice restart file with after restoring heat content
c
      call rstrti
c
#endif
      return
      end
