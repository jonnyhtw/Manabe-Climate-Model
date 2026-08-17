      subroutine indati
c
c         called from: atmos
c        entry points: rstrti
c  subroutines called: none
c
c  this routine is very similar to indato.  it reads the initial
c   ice restart tape.  the entry point rstrto, writes the ice
c   model data to a history file.
c
#include "parm.h"
c
      parameter (iunt = 1)
c
c location of time for restart
c
#include "ctlblk.h"
c
c contains surface temperature
c
#include "grdfld.h"
c
c contains ocean temperature and liquid equivalent ice
c
#include "surfo.h"
c
#include "varice.h"
c
#include "variceco.h"
c
      dimension pland(ix,il)
c
c set ice-liquid equivalent multipler
c
      cliqsol = 1.06
#ifndef first_spinup
c
c skip control block record
      read (iunt)
c
      call rdice(iunt)
c
      print *, 'ice restart data read in ', month, ' ', iday, ' ', iyear
#endif
c
#ifdef first_spinup
c
c ice variables not presented 
c
      print *, '------------------------------------------------'
      print *, 'initialize ice variables using atmospheric data!'
c
c temperature on top of ice is treated as surface temperature
c
      do j = 1, il
      do i = 1, ix
         pland(i,j) = 0.
         tsice(i,j) = -ts(i,j)
c
c distinquish between land and ocean point
c
         if (ts(i,j) .gt. 0.) then
            pland(i,j) = 1.
         end if
      end do
      end do
c
#ifdef init_ice
c
c hinit is initial ice thickness
c
      hinit = 300.
c
c specify an uniform ice sheet
c
      inorth = int(il/8)
      isouth = il - inorth + 1
c
c northern hemisphere
c
      do j = 1, inorth
      do i = 1, ix
         if (pland(i,j) .eq. 0.) then
            yce(i,j) = hinit
         end if
      end do
      end do
c
c southern hemisphere
c
      do j = isouth, il
      do i = 1, ix
         if (pland(i,j) .eq. 0.) then
            yce(i,j) = hinit
         end if
      end do
      end do
#endif
c
c set up ice thickness and compactness
c
      do j = 1, il
      do i = 1, ix
         heff(i,j,1) = yce(i,j)*cliqsol*(1.-pland(i,j))
         if (yce(i,j).gt.0..and.pland(i,j).ne.1..and..not.swamp) then
            area(i,j,1) = 1.
         end if
c     print *, i,j,yce(i,j),heff(i,j,1),area(i,j,1),ts(i,j),pland(i,j)
      end do
      end do
c
#endif
c
c fix all variables using restart data
c
      do j = 1, il
      do i = 1, ix
         heff(i,j,2) = heff(i,j,1)
         area(i,j,2) = area(i,j,1)
c
c set up ice thickness for surface heat budget
c
	 if (area(i,j,2) .ne. 0.) then
         hice(i,j) = heff(i,j,2)/area(i,j,2)
         end if
c
c set compactness to zero when there is no ice (for temperature weighting)
c
         areatemp = heff(i,j,1)/0.0001
         area(i,j,2) = amin1(area(i,j,2),areatemp)
c
c input secondary ice velocities
c
         uic(i,j,2) = uic(i,j,1)
         vic(i,j,2) = vic(i,j,1)
      end do
      end do
c
      rewind (iunt)
      return
c
      entry rstrti
c     ===== ======
c
c  called from sphercav
c
c
      rewind (iunt)
c
      call wrtctl(iunt,0)
c
      call wrtice(iunt)
c
      print *, 'ice restart file for ', month, ' ', iday, ' ', iyear
c
      close (iunt)
c
      return
      end
