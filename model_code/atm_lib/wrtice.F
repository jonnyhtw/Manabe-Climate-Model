      subroutine wrtice(iunit)
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
#include "surfo.h"
c
c This subroutine writes sea ice restart data.
c
c
c write in ice restart file
c
      write (iunit) ((area(i,j,1),i=1,ix),j=1,il)
      write (iunit) ((heff(i,j,1),i=1,ix),j=1,il)
      write (iunit) tsice
c
c restart variables for ice dynamics calculations
c
#ifdef mixed_layer
      write (iunit) ((uic(i,j,1),i=1,ix),j=1,il)
      write (iunit) ((vic(i,j,1),i=1,ix),j=1,il)
      write (iunit) presu
#endif
c
      return
      end
