c
c subroutine to generate ice code dimensions, constants and masks
c
      subroutine bndsi
#if defined mixed_layer && defined ice_dynamics
c
#include "parm.h"
c
      parameter (jmt=il+2)
c
c obtain latitudes from climate model
c
#include "astro.h"
c
c obtain temperature used as mask
c
#include "grdfld.h"
c
c obtain latitudes for grid center and grid corners directly
c
#include "intbc.h"
c
c put in constant 
c
#include "variceco.h"
c
      piice = acos(-1.)
c
c rhoice is density of sea ice (kg/m^3) in dynamic calculation
c
      rhoice = 920.
c
c r is radius (in m) of earth
c
      r = 6378000.
c
c phi is wind turning angle
c
      phi = 25.*piice/180.
      if (ihem .eq. 2) then
         phi = -phi
      end if
c
c no wind turning angle if using wind stress vectors
c
      phi = 0.
c
      sinphi = sin(phi)
      cosphi = cos(phi)
c
c theta is water turning angle
c
      theta = 25.*piice/180.
      if (ihem .eq. 2) then
         theta = -theta
      end if
      sinth = sin(theta)
      costh = cos(theta)
c
c pstari is ice strength (empirical) of solid sea ice at 1 meter
c
      pstari = 275000.
c
c set latitude points and Coriolis parameters
c
      alimit = alat(jlimit) - 1.
c     alimit = 20.
      j = 1
  12  continue
      if (alat(j) .lt. alimit) then
         go to 18
      end if
      j = j + 1
      go to 12
  18  continue
c
c compute the number
c
c     np = j + 3
      np = j + 4
c
c testing the output
c
c     np = 16
c
      if (np .gt. iy+1) then
         print *, 'ice field extend beyond one hemisphere'
         np = iy
      end if
c     print *, alimit, np, jlimit, iy
#ifndef icecpld
c
c declare the values of arlat and arlatc (in radians)
c arlat is corner point, arlatc is center point
c
      if (ihem .eq. 1) then
         do 27 j = 1, np
            arlatc(np+1-j) = alat(j)*piice/180.
  27     continue
         arlatc(np+1) = 0.5*(piice/2.+arlatc(np))
      end if
c
      if (ihem .eq. 2) then
         do 29 j = 2, np+1
            arlatc(j) = alat(il-j+2)*piice/180.
c           arlatc(np+1-j) = -arlatc(np+1-j)
  29     continue
         arlatc(1) = 0.5*(-piice/2.+arlatc(2))
c        arlatc(np+1) = -arlatc(np+1)
      end if
c
      do 39 j = 1, np
         arlat(j) = 0.5*(arlatc(j)+arlatc(j+1))
  39  continue
c     do j = 1, np
c        print *, j, arlat(j)*180./piice, arlatc(j)*180./piice
c     end do
#endif
#ifdef icecpld
c
c use latitude values from coupled atmosphere ocean model
c
      if (ihem .eq. 1) then
         do 26 j = 1, np
            arlatc(np+1-j) = ay(il+1-j)
            arlat(np+1-j) = oyv(il+1+1-j)
  26     continue
         arlatc(np+1) = 2*(piice/2.-arlatc(np)) + arlatc(np)
      end if
c
      if (ihem .eq. 2) then
         do 28 j = 2, np+1
            arlatc(j) = ay(j-1)
            arlat(j-1) = oyv(j-1)
  28     continue
         arlatc(1) = 2.*(-piice/2.-arlatc(2)) + arlatc(2)
         arlat(1) = 0.5*(arlatc(1)+arlatc(2))
      end if
c     print *, '----------------'
c     do j = 1, np
c        print *, j, arlat(j)*180./piice, arlatc(j)*180./piice
c     end do
c
#endif
c
c compute varying Coriolis parameters 
c fcorp is corner point, fcorc is center point
c
      omega = 2.*piice/86400.
      do 41 j = 1, np
         fcorp(j) = 2.*omega*sin(arlat(j))
	 fcorc(j) = 2.*omega*sin(arlatc(j))
  41  continue
      fcorc(np+1) = 2.*omega*sin(arlatc(np+1))
c
c set longitude points 
c
      nl = ix
      dlamd = 360./float(ix)
      dlam = dlamd*piice/180.
c
c refresh outflow cell locations
c
      do 44 j = 1, np
      do 44 i = 1, nl
         out(i,j) = 0.
  44  continue
c
c set up ice thickness mask for northern or southern hemisphere
c
      do 45 j = 1, np-1
      do 45 i = 1, nl
         inew = i 
c        if (i .eq. nl-1) then
c           inew = 1
c        end if
c        if (i .eq. nl) then
c           inew = 1
c        end if
c
c northern hemisphere
c
         k = j
c
c southern hemisphere
c
         if (ihem .eq. 2) then
c           k = il + 1 - j
            k = il - np + j + 1
         end if
c
         if (ts(inew,k) .ge. 0.) then
            heffm(i,np+1-j) = 0.
         else
            heffm(i,np+1-j) = 1.
         end if
  45  continue
c
      do 46 i = 1, nl+1
         heffm(i,1) = 0.
         heffm(i,np+1) = 0.
c
c out are cells of outflow where ice strength is zero
c
	 if (ihem .eq. 1) then 
            out(i,2) = 1.
         end if
         if (ihem .eq. 2) then
            out(i,np) = 1.
c           out(i,2) = 1.
         end if
  46  continue
c
c cyclic condition
c
      do 49 j = 1, np+1
         heffm(nl+1,j) = heffm(1,j)
  49  continue
c
c     do j = 1, np+1
c        k = np + 1 - j + 1
c        write(6,59) (heffm(i,k),i=1,nl+1)
c 59     format (48f2.0)
c     end do
c
#endif
      return
      end
