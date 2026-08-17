      subroutine rdocn(iunit)
c
#include "parm.h"
c
#include "anlocn.h"
c
#include "ctlblk.h"
c
#include "grdfld.h"
c
#include "mrstr.h"
c
#include "qa.h"
c
#include "surfo.h"
c
#ifdef ice_dynamics
c
c input ice velocities in ocean history (cfi)
c
#include "varice.h"
#endif
c
#include "zmanl.h"
c
c This subroutine reads ocean restart data.
c
      read (iunit) to
      read (iunit) yce
      read (iunit) qao
      read (iunit) sno
      read (iunit) ycemlt
c     read (iunit) zmto,zmyce,zmqo,zmsno
c
#ifdef ice_dynamics
c
c add in ice velocities variables (cfi)
c
#ifndef first_spinup
         read (iunit) tauxic
         read (iunit) tauyic
         read (iunit) ar
c
c additional ice quatities from ice dynamics model (cfi)
c
#ifdef mixed_layer
         read (iunit) qficec
         read (iunit) qfowc
         read (iunit) uicec
         read (iunit) vicec
         read (iunit) giceadv
         read (iunit) gicetherm
         read (iunit) gareadv
         read (iunit) garetherm
c
c        read (iunit) siceic
c        read (iunit) siceow
c        read (iunit) hmlc
c        read (iunit) sicec
c        read (iunit) sareac
c        read (iunit) qflc
c        read (iunit) gicet
c        read (iunit) garet
c        read (iunit) windvx
c        read (iunit) windvy
c
c restoring heat flux (cfi)
c
#ifdef restore
	 read (iunit) qres
#endif
c
#endif
c
#endif
c
#endif
c
#ifdef q_flux
#ifndef first_spinup
      read (iunit) qbot
#endif
#endif
c
      read (iunit) zmto,zmyce,zmqo,zmsno
c
      return
      end
