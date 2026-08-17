c
c subroutine to prepare inputs from ice model to ocean model
c
      subroutine bcice
#if defined icecpld
c
#include "parm.h"
c
      parameter (imt=2*ix+2, jmt=il+2)
c
c location of boundary condition quantities to ocean model
c
#include "ascomo.h"
c
c location of surface temperature
c
#include "grdfld.h"
c
c location of quantities used in interpolation
c
#include "intbc.h"
c
c location of precipitation
c
#include "qa.h"
c
#include "varice.h"
c
      dimension zin(ix,il)
      dimension bci(imt,jmt,nbo)
      dimension tauxoc(ix,il), tauyoc(ix,il)
      dimension seamk(imt,jmt), seamko(imt,jmt)
c
c conversion factor from minute to second
c
      data cms /60./
c
c average salinity in ice
c
      data avesali /0.035/
c
c set up land mask
c
      do 10 j = 1, il
      do 10 i = 1, ix
         if (ts(i,j) .gt. 0.) then
            zin(i,j) = 1.0
         else
            zin(i,j) = 0.0
         end if
  10  continue
c
      do 20 j = 1, il
      do 20 i = 1, ix
         tauxoc(i,j) = tauxic(i,j) + fxice(i,j)
         tauyoc(i,j) = tauyic(i,j) + fyice(i,j)
	 tauxic(i,j) = tauxoc(i,j)
	 tauyic(i,j) = tauyoc(i,j)
  20  continue
c
      lgrid = 0
      call intrpa(qocean,bci(1,1,1),lgrid,zin,seamk)
      call intrpa(wocean,bci(1,1,3),lgrid,zin,seamk)
      lgrid = 1
      call intrpa(tauxoc,bci(1,1,4),lgrid,zin,seamko)
      call intrpa(tauyoc,bci(1,1,5),lgrid,zin,seamko)
c
c     do 30 j = 1, jmt 
c     do 30 i = 1, imt
c
c convert heat flux from cal/cm^2/min to cal/cm^2/sec
c
         bci(i,j,1) = bci(i,j,1)/cms
c
c convert water flux to salt flux g/cm^2/sec
c (positive means decrease in salinity)
c
         bci(i,j,3) = -bci(i,j,3)/cms*avesali
c
  30  continue
c
c correct for water flux to ocean and add in ice melt contribution
c
c     do 40 j = 1, jmt
c     do 40 i = 1, imt
         bci(i,j,2) = -bci(i,j,2) + bci(i,j,3)/avesali
  40  continue
c
c put fluxes to boundary variables
c
c     do 50 n = 1, nbo
      do 50 j = 1, jmt
      do 50 i = 1, imt
c        bco(i,j,n,1) = bci(i,j,n)
         bco(i,j,1,1) = bci(i,j,1)
         bco(i,j,2,1) = bci(i,j,2)
         bco(i,j,3,1) = bci(i,j,3)
         bco(i,j,4,1) = bci(i,j,4)
         bco(i,j,5,1) = bci(i,j,5)
  50  continue
#endif
      return
      end
