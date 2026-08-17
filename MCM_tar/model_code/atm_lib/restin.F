c
c subroutine to read in restart for ice variables
c
      subroutine restin
#ifdef ice_dynamics
c
#include "parm.h"
c
#include "iohist.h"
c
#include "varice.h"
c
#include "variceco.h"
c
c read in previous end of the month ice variables
c
c     open (unit=iycrst, file='ice.restart', form='unformatted')
c     open (unit=iycthk, file='ice.thick.extent', form='formatted')
      open (unit=iycthk, file='ice.aver', form='unformatted')
      read (iycrst,end=17) ((heff(i,j,1),area(i,j,1),presu(i,j),
     >           i=1,ix), j=1,il)
      read (iycrst,end=17) ((uic(i,j,1),vic(i,j,1),i=1,ix), j=1,il+2)
c
      print *, 'finishing reading ice restart'
 17   continue
      print *, 'no data in ice restart file!'
 18   continue
c
#endif
      return
      end
