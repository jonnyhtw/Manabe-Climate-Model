c
      subroutine sphercav
#ifdef ice_dynamics
c********************************************************************
c
c          Spherical C-grid Cavitating Fluid Sea Ice Model
c               (Version with island at North Pole)
c
c  Ref: Flato, G.M. and W.D. Hibler III, "Modeling Pack Ice as a 
c       Cavitating Fluid", J. Phys. Oceanogr., 22(6), 626-651, 1992.
c
c  written by: G.M. Flato, Institute of Ocean Sciences, Sidney, BC
c             (formerly at Dartmouth College, Hanover, NH)
c             E-Mail: flato@ios.bc.ca
c
c********************************************************************
c
c  Original code: 29/Mar/90
c
c  revisions:18/Dec/91 (jim waugh) replaced non-standard functions
c                      cosd, atan2d, atand, and sind with standard
c                      functions cos, atan2, atan, and sin.
c                      also, from within subroutine bnds remove the
c                      DATA statements that were initializing COMMON
c                      variables and place them in a BLOCK DATA
c                      subprogram where the fortran77 standard specifies
c                      that they should go.
c            20/July/93 (greg flato) modified to use non-linear drag 
c		       parameterization. Misc. code clean-up.
c
c            7/22/93   (chi ip) distinguish between wind drag
c                      and water drag turning angles
c
c            9/8/93    (chi ip) modified to allow for irregular
c                      latitude spacing
c
c            9/16/93   (chi ip) modified to a subroutine to allow
c                      for use with global climate model
c
c********************************************************************
c
#include "parm.h"
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
#include "iohist.h"
c
c location of liquid equivalent ice for initialization
c
#include "surfo.h"
c
#include "varice.h"
c
#include "variceco.h"
c
      common/six/tmix(ix+1,il+1),tice(ix+1,il+1),yneg(ix+1,il+1),
     <fw(ix+1,il+1)
c
      dimension fh(ix+1,il+1),heffa(ix+1,il)
      dimension sealn(ix,il)
c
c record intermediate ice thickness and compactness
c
      dimension heffint(ix,il), areaint(ix,il)
c
c growth of thickness and compactness within a time step
c
      dimension gicead(ix,il), gicether(ix,il) 
      dimension garead(ix,il), garether(ix,il)
c
      character*8 monthsi(12)
      dimension mo(12)
c
      data monthsi / 'jan.    ', 'feb.    ', 'mar.    ',
     >              'apr.    ', 'may     ', 'june    ',
     >              'july    ', 'aug.    ', 'sep.    ',
     >              'oct.    ', 'nov.    ', 'dec.    ' /
      DATA MO/31,28,31,30,31,30,31,31,30,31,30,31/
c
c cmcm is conversion from meter to centimeter
c cscgs is conversion of stress to cgs units
c
      data cmcm /100./, cscgs /10./
c
c velcap is ice velocity limit (m/s)
c
      data velcap /0.00001/
#ifdef mixed_layer
c
c execute the main body only for interactive ocean
c otherwise, read the ice thickness and compactness from
c common block / bcs /
c
c get days used in monthly average
c
      monthc = 0
      do i = 1, 12
         if (month .eq. monthsi(i)) then
            monthc = i 
         end if
      end do
c
c iday calculated in atmospheric leg is ahead of ice day by 1
c
      if (iday .eq. 1) then
         monthc = monthc - 1
      end if
c
c put cyclic condition on months of the year
c
      if (monthc .eq. 0) then
         monthc = 12
      end if
c
c daysm is the number of days in this month
c
      daysm = float(mo(monthc))
c
      idyc = iday - 1
      if (iday .eq. 1) then
         idyc = mo(monthc)
      end if
c
      print *, '---------------------------------------'
      print *, 'entering dynamic thermodynamic ice code.'
      print *, 'model date   ', monthsi(monthc), idyc, ', ', iyear
c
c  specify parameters
c
      fcor=1.46D-04
c
c ntime is number of time step used within a day
c
      ntime = 2
c
c deltat, dticet are dynamic, thermodynamic time step length (in sec)
c
      deltat = 86400./float(ntime)
      dticet = 86400./float(ntime)
c
c prepare for time stepping
c
      do j = 1, il
      do i = 1, ix
c
c reset cummulative average variables
c
	 giceadv(i,j) = 0.
         gicetherm(i,j) = 0.
         gareadv(i,j) = 0.
         garetherm(i,j) = 0.
         pocean(i,j) = 0.
         qocean(i,j) = 0.
c
c toint is the ocean temperature used in the ice model
c
         toint(i,j) = to(i,j)
      end do
      end do
c
c begin time stepping within ice model 
c
      nti = 1
  19  continue
c
      print *, '---------------------'
      print *, 'time step no ', nti
c
c indicator for calculation of hemispheres
c inhem of 0 is two hemispheres, inhem of 1 is north, ihem of 2 is south
c
      inhem = 0
      ihem = inhem
      if (inhem .eq. 0) then
         ihem = 1
      end if
c
c get ice compactness and ice thickness (in cm) fields
c
      do j = 1, il
      do i = 1, ix
c
c h,a(i,j,3) keeps the ice thickness (in cm), compactness before ice model
c
         heff(i,j,3) = heff(i,j,1)
         area(i,j,3) = area(i,j,1)
c
c h,aint are intermediate ice thickness and compactness
c
         heffint(i,j) = heff(i,j,1)
         areaint(i,j) = area(i,j,1)
c
      end do
      end do
c     do j = 1, il
c     do i = 1, ix
c        tauxic(i,j) = float(i)
c        tauyic(i,j) = float(j)
c     end do
c     end do
c
#ifdef icecpld
         call windcurr
#endif
c
  39  continue
c
      print *, '------------------'
      print *, 'hemisphere no is ', ihem
c
c  set up parameters and boundary masks
c
      call varres
      call bnds
      call bndsi
      call bndsii
c
c put in thickness and compactness for dynamic section
c
      call dynthick
c
c change ice thickness from centimeter to meter
c
      do j = 1, np+1
      do i = 1, nl+1
         hef(i,j,1) = hef(i,j,1)/cmcm
      end do
      end do
c
c put in ice velocities for dynamic section
c
      call dynveloc
c
c take care of outflow region
c
      do j = 1, np+1
      do i = 1, nl+1
         if (out(i,j) .eq. 1.) then
            hef(i,j,1) = 0.
            hef(i,j,2) = 0.
            are(i,j,1) = 0.
            are(i,j,2) = 0.
         end if
      end do
      end do
c
      do j=1,np+1
      do i=1,nl+1
         yneg(i,j)=0.
         tice(i,j)=273.0D+00
      enddo
      enddo
c
#ifndef icecpld
c
c  read in annual average current field
c
      call currs
c
#endif
c
c***********************************************************
c
c  begin ice velocity 
c
c***********************************************************
c
c  calculate strength and initialize pressure field to 0 
c
      do j=1,np+1
      do i=1,nl+1
         xmpress(i,j)=pstari*hef(i,j,1)*exp(-20.
     <                *(1.-are(i,j,1)))
c
c  cfi
c
c  prevent using out of bound element in array
c
c        presur(i,j) = 0.
c.o      presur(i+1,j+1)=0.
      enddo
      enddo
c
#ifndef icecpld
c
c  read wind field or wind stress field
c
      call winds
c
#endif
c
c change wind stress from cgs to mks units (N/m^2)
c
      do j = 1, np
      do i = 1, nl
         tauxin(i,j) = tauxin(i,j)/cscgs
         tauyin(i,j) = tauyin(i,j)/cscgs
      end do
      end do
c
      nt = 0
      tautol = 0.
      if (ihem .eq. 4) then
      do j = np, 1, -1
      do i = 1, nl
         tauspeed = sqrt(tauxin(i,j)**2+tauyin(i,j)**2)
         if (tauspeed .ne. 0.) then
            tautol = tautol + tauspeed
            nt = nt + 1
         end if
         print *, '------------------------'
         print *, i+1,j+1,tauxo(i+1,j+1),tauyo(i+1,j+1)
         print *, i,j,tauxin(i,j),tauyin(i,j)
      end do
      end do
      stop
      end if
c     if (nt .ne. 0) then
c        print *, 'average stress', nt, tautol/nt
c     end if
c
c  put minimum value on suface wind stress components to avoid
c  problems with inverse trig function in FORCES
c
      do j=1,np
      do i=1,nl
         if(abs(tauxin(i,j)).lt.0.00005) tauxin(i,j)=0.00005
         if(abs(tauyin(i,j)).lt.0.00005) tauyin(i,j)=0.00005
      enddo
      enddo
c
c  put minimum value on wind velocity components to avoid
c  problems with inverse trig function in FORCES
c
      do j=1,np
      do i=1,nl
         if(abs(gairx(i,j)).lt.0.00005) gairx(i,j)=0.00005
         if(abs(gairy(i,j)).lt.0.00005) gairy(i,j)=0.00005
      enddo
      enddo
c
c...Modified Euler time step for non-linear drag terms
c
      do kme=1,2
c
c  average velocities here for second iteration
c
      if(kme.ge.2) then
         do j = 1, np
         do i = 1, nl
            uice(i,j,1) = 0.5*(uice(i,j,1)+uice(i,j,3))
            vice(i,j,1) = 0.5*(vice(i,j,1)+vice(i,j,3))
         end do
         end do
      endif      
c 
c  update velocity array
c
      do j = 1, np
      do i = 1, nl     
         uice(i,j,3) = uice(i,j,1)
         vice(i,j,3) = vice(i,j,1)    
      end do
      end do
c
c  calculate drag coefficients at B-grid locations where forcing is
c  specified (averaging to C-grid is done in FORCES)
c
      do j = 1, np
      do i = 1, nl
         i0 = i
         if(i0.eq.1) i0=nl+1
         im1=i-1
         if(im1.eq.0) im1=nl
         jm1=j-1
         if(jm1.eq.0) jm1=1
         uavb=0.5*(uice(i,j,1)+uice(i,jm1,1))
         vavb=0.5*(vice(i,j,1)+vice(im1,j,1))
c
c use surface wind stress vectors as air drag
c
         dairn(i,j)=sqrt(tauxin(i,j)**2+tauyin(i,j)**2)
c        dairn(i,j)=1.56D-03*sqrt(gairx(i,j)**2+gairy(i,j)**2)
         dwatn(i,j)=5.5D+00*sqrt((uavb-gwatx(i,j))**2
     <      +(vavb-gwaty(i,j))**2)
c
c put minimum on drag coefficient to avoid spuriously large velocities
c when near zero drag occurs.
c
c        if(dwatn(i,j).lt.0.25) dwatn(i,j)=0.25
c
c*** Option: linear drag coefficients ***
c  Ref: Flato, G.M. and W.D. Hibler III, "On a simple sea-ice dynamics
c   model for climate studies", Annals of Glaciology, 14, 72-77, 1990.
c  NOTE: for linear drag calculation, need only one loop over this time-
c   stepping scheme; i.e. replace statement ~35 lines up by "do kme=1,1"
c
c         dairn(i,j)=0.01256
cajb The following line was "un-commented" during August of year 11
cajb (during the first austral winter of the restoring run).
          dwatn(i,j)=0.6524
c         
c.o      dwatx=dwatn(i,j)
c.o      dwaty=dwatn(i,j)
c
c use linear water drag if near zero drag occurs
c
         if(dwatn(i,j).lt.0.25) dwatn(i,j)=0.6524
c
         Axi(i,j)=dwatn(i,j)*costh
         Bxi(i,j)=((hef(i+1,j+1,1)+hef(i0,j+1,1))/2.*rhoice*fcorc(j+1))
     <      +(dwatn(i,j)*sinth)
         Ayi(i,j)=dwatn(i,j)*costh
         Byi(i,j)=((hef(i+1,j+1,1)+hef(i+1,j,1))/2.*rhoice*fcorp(j))
     <      +(dwatn(i,j)*sinth)
      end do
      end do
c     print *, 'ihem = ', ihem
c     do j = 1, np
c     do i = 1, nl
c        print *, i,j,hef(i,j,1),are(i,j,1),xmpress(i,j),presur(i,j)
c        print *, axi(i,j),ayi(i,j),uice(i,j,1),vice(i,j,1)
c        print *, '------------------------------------'
c     end do
c     end do
c
c  calculate forcing fields
c
      if(kme.eq.1) call forces      
c
c  calculate free drift velocity field
c
      call frees
c     do j = 1, np
c     do i = 1, nl
c     write (iycthk, *) i,j,uice(i,j,1),vice(i,j,1),hef(i,j,1)
c     end do
c     end do
c     stop
c
c  calculate corrected velocity field
c
      call nlcavs
c
c...end of Modified Euler time step
c
      end do
c
c testing for zero velocity ice movement
c
c     print *, '-----------------------------------'
c     print *, 'ihem = ', ihem
c     do j = 1, np
c     do i = 1, nl
c        uice(i,j,1) = 0.09
c        vice(i,j,1) = 0.
c        print *, i,j, uice(i,j,1), vice(i,j,1)
c     end do
c     end do
c     print *, '-----------------------------------'
c
c  do advection using upstream differencing
c
      call advects
c
c zero outflow regions
c
      do j = 1, np+1
      do i = 1, nl+1
         if (out(i,j) .eq. 1.) then
            hef(i,j,1) = 0.
            are(i,j,1) = 0.
	    presur(i,j) = 0.
         end if
      end do
      end do
c
      do j = 1, np+1
         hef(1,j,1) = hef(nl+1,j,1)
         are(1,j,1) = are(nl+1,j,1)
         presur(1,j) = presur(nl+1,j)
      end do
c
c change ice thickness to units of centimeters
c
      do j = 1, np+1
      do i = 1, nl+1
         hef(i,j,1) = hef(i,j,1)*cmcm
      end do
      end do
c
c zero out small ice velocities
c
      do j = 1, np
      do i = 1, nl
         uspe = abs(uice(i,j,1))
         vspe = abs(vice(i,j,1))
         if (uspe .lt. velcap) then
            uice(i,j,1) = 0.
         end if
         if (vspe .lt. velcap) then
            vice(i,j,1) = 0.
         end if
      end do
      end do
c
c  interpolate velocity field back to grid center for plotting 
c
      do j = 1, np
      do i = 1, nl
         i0 = i - 1
         if (i0 .lt. 1) i0 = nl
         jm1 = j - 1
         if (jm1 .eq. 0) jm1=1
         uice(i,j,2) = heffm(i,j)*0.5*(uice(i0,jm1,1)
     <              +uice(i,jm1,1))
         vice(i,j,2) = heffm(i,j)*0.5*(vice(i0,j,1)
     <              +vice(i0,jm1,1))
      end do
      end do
c
c calculate internal ice force
c
      call fint
c
c convert interal ice force to cgs units
c
      do j = 1, np
      do i = 1, nl
         fxi(i,j) = fxi(i,j)*cscgs
         fyi(i,j) = fyi(i,j)*cscgs
      end do
      end do
c
c prepare ice thickness for thermodynamic growth
c
      call thermthick
c
c put ice velocities for each hemisphere in global ice variables
c
      call thermveloc
c
c     print *, '-----------------------------------'
c     print *, 'ihem = ', ihem
c     do j = 1, il
c     do i = 1, ix
c        print *, i,j,heff(i,j,1),presu(i,j)
c     end do
c     end do
c     print *, '-----------------------------------'
c
      jlimit = 0
      jinit = 1
      je = iy
      js = 1
      if (ihem.eq.2) then
      jinit = il 
      je = iy + 1
      js = -1
      end if
      do j = jinit, je, js
      do i = 1, ix
      jl = j
      if (ihem .eq. 2) then
      jl = il - j + 1
      end if
      if (heff(i,j,1) .gt. 0.) then
      jlimit = jl
      end if
      end do
      end do
      print *, 'number of latitude grid points is ', np
c
c check if more hemisphere calculations are needed
c
      if (inhem .lt. 1) then
         if (ihem .lt. 2) then
            ihem = ihem + 1
            go to 39
         end if
      end if
c
c****************************************************
c
c end of dynamic section
c
c****************************************************
c
c advected thickness and compactness in ice growth
c
c.o   do j=1,np
c.o   do i=1,nl
c.o      heff(i,j,2)=heff(i,j,1)
c.o      area(i,j,2)=area(i,j,1)
c.o   enddo
c.o   enddo
c
c  do thermodynamic growth/decay
c
  49  continue
c     do j = 1, il
c     do i = 1, ix
c     print *, '------------------------------------------------'
c     print *, i, j, uic(i,j,1), vic(i,j,1), presu(i,j)
c     print *, taux(i,j), tauy(i,j), windvx(i,j), windvy(i,j)
c     print *, '------------------------------------------------'
c     end do
c     end do
c     stop
c
c calculate ice thickness and compactness from advection
c
      do j = 1, il
      do i = 1, ix
         gicead(i,j) = heff(i,j,1) - heffint(i,j)
         garead(i,j) = area(i,j,1) - areaint(i,j)
         heffint(i,j) = heff(i,j,1)
         areaint(i,j) = area(i,j,1)
      end do
      end do
c
c
      call growths
c
c set up for next time step in thermodynamic calculation
c
      do j = 1, il
      do i = 1, ix
         heff(i,j,2) = heff(i,j,1)
         area(i,j,2) = area(i,j,1)
c
c set compactness to zero when there is no ice 
c
         areatemp = heff(i,j,1)/0.00001
c
c area(i,j,2) is used in surface temperature weighting 
c
	 area(i,j,2) = amin1(area(i,j,2),areatemp)
      end do
      end do
c
c get ice thickness and compactness from thermodynamics and ice model
c
      do j = 1, il
      do i = 1, ix
         gicether(i,j) = heff(i,j,1) - heffint(i,j)
         garether(i,j) = area(i,j,1) - areaint(i,j)
         gicet(i,j) = heff(i,j,1) - heff(i,j,3)
         garet(i,j) = area(i,j,1) - area(i,j,3)
      end do
      end do
c
c convert gice, gare quantities to monthly values to be average
c
      do j = 1, il
      do i = 1, ix
         giceadv(i,j) = gicead(i,j)*daysm + giceadv(i,j)
         gareadv(i,j) = garead(i,j)*daysm + giceadv(i,j)
         gicetherm(i,j) = gicether(i,j)*daysm + gicetherm(i,j)
         garetherm(i,j) = garether(i,j)*daysm + gicetherm(i,j)
      end do
      end do
c
c prepare ice compactness field for ice history 
c
      do j = 1, il
      do i = 1, ix
         ar(i,j) = area(i,j,2)
      end do
      end do
c
c prepare flux inputs to ocean (with one-day time step)
c
      do j = 1, il
      do i = 1, ix
         qocean(i,j) = qoceani(i,j)/float(ntime) + qocean(i,j)
         pocean(i,j) = poceani(i,j)/float(ntime) + pocean(i,j)
      end do
      end do
c
c prepare internal ice stress to ocean (with one-day time step)
c
      do j = 1, il
      do i = 1, ix
         fxice(i,j) = fxicei(i,j)/float(ntime) + fxice(i,j)
         fyice(i,j) = fyicei(i,j)/float(ntime) + fyice(i,j)
      end do
      end do
c
c add up ice velocities at ice locations
c
      do j = 1, il
      do i = 1, ix
         if (heff(i,j,1) .gt. 0.) then
            uicea(i,j) = uicea(i,j) + uicec(i,j)*100.
            vicea(i,j) = vicea(i,j) + vicec(i,j)*100.
            uicec(i,j) = uicec(i,j)*100.
            vicec(i,j) = vicec(i,j)*100.
         end if
	 if (heff(i,j,1) .le. 0.) then
            uicec(i,j) = 0.
       	    vicec(i,j) = 0.
         end if
      end do
      end do
c
c  print out velocity and thickness fields at end of month
c
c.o   kmonth = 0
c.o   do kc = 1,12
c.o      kmonth=kmonth+mo(kc)
c.o      if(iday.eq.kmonth) then
c
      if (iday .eq. 1) then
         kc = monthc
         rmo = mo(kc)*1.0
c        write (iycthk) ((uicea(i,j)/rmo,vicea(i,j)/rmo,i=1,ix), j=1,il)
         do j = 1, il
         do i = 1, ix
	    uicea(i,j) = uicea(i,j)/rmo
	    vicea(i,j) = vicea(i,j)/rmo
         end do
         end do
         write (iycthk) uicea
         write (iycthk) vicea
c
         do j = 1, il
         do i = 1, ix
               uicea(i,j) = 0.
               vicea(i,j) = 0.
         end do
         end do
       end if
c.o         write(10,100) ((10.*heff(i,np+1-j,1),i=2,nl+1),j=1,np)
c.o         write(10,101)
c.o         write(10,100) ((presur(i,np+1-j)/1000.,i=2,nl+1),j=1,np)
c.o         do j=1,np
c.o         do i=1,nl+1
c.o            i1=i+1
c.o            if(i.eq.nl+1) i1=2
c.o            heffa(i,j)=heff(i1,j+1,1)
c.o            if(heffm(i1,j+1).lt.0.1) heffa(i,j)=-1.
c.o         enddo
c.o         enddo
c.o         write(11,102) ((heffa(i,j),i=1,nl+1),j=1,np)
c
c.o      endif
c.o   enddo
c
 1000 continue
 1001 continue
c
      if (nti .lt. ntime) then
         nti = nti + 1
         go to 19
      end if
c
      print *, '---------------------------------------'
c
#endif
c
#ifdef fixed_sst
c
c sealn is landsea mask, 0 for land; 1 for sea
c
          do j = 1, il
	  do i = 1, ix
            if (ts(i,j) .gt. 0.0) then
               sealn(i,j) = 0.
            else
               sealn(i,j) = 1.
            end if
          end do
	  end do
	do j = 1, il
          do i = 1,ix
            hice(i,j) = bcyce(i,j)*cliqsol*sealn(i,j)
            area(i,j,2) = bcconc(i,j)*sealn(i,j)
            heff(i,j,2) = hice(i,j)*bcconc(i,j)
c
            yceol(i,j) = yce(i,j)
            yce(i,j) = bcyce(i,j)*bcconc(i,j)*sealn(i,j)
c
          enddo
        enddo
c
c preparing variables for ocean and ice history
c
        do j = 1, il
        do i = 1, ix
           heff(i,j,1) = heff(i,j,2)
           area(i,j,1) = area(i,j,2)
           ar(i,j) = bcconc(i,j)*sealn(i,j)
        end do
        end do
#endif
c
#ifndef restore
c
c write ice restart if necessary
c
      call rstrti
#endif
c
#ifdef icecpld
c
c prepare linking to ocean model
c
c     call bcice
#endif
c
c change format statement suitable for cray YMP at GFDL
c
  100 format(2f3.0)
c 100 format(36f3.0)
  101 format(///)
  102 format(2g10.4)
c 102 format(37g10.4)
  103 format(2g16.5)
c
#endif
      return
      end
