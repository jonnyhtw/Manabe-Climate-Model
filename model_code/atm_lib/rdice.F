      subroutine rdice(iunit)
c
#include "parm.h"
c
#include "varice.h"
c
#include "variceco.h"
c
#include "surfo.h"
c
c This subroutine reads sea ice restart data.
c
#ifndef first_spinup
c
c read in ice restart file
c
      read (iunit) ((area(i,j,1),i=1,ix),j=1,il)
      read (iunit) ((heff(i,j,1),i=1,ix),j=1,il)
      read (iunit) tsice
c
c restart variables for ice dynamics calculations
c
#ifdef mixed_layer
      read (iunit) ((uic(i,j,1),i=1,ix),j=1,il)
      read (iunit) ((vic(i,j,1),i=1,ix),j=1,il)
      read (iunit) presu
#endif
c
#endif
c
      return
      end
