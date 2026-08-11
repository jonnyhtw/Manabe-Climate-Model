      subroutine wrtocn(iunit)
c
c This subroutine reads ocean restart data.
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
      write(iunit) to
      write(iunit) yce
      write(iunit) qao
      write(iunit) sno
      write(iunit) ycemlt
c
#ifdef ice_dynamics
c
c put ice variables in output (cfi)
c
      write (iunit) tauxic
      write (iunit) tauyic
      write (iunit) ar
c     write (iunit) fxice
c     write (iunit) fyice
c
c additional ice quantites output from ice dynamics model (cfi)
c
#ifdef mixed_layer
         write (iunit) qficec
         write (iunit) qfowc
         write (iunit) fxice
         write (iunit) fyice
c        write (iunit) uicec
c        write (iunit) vicec
         write (iunit) giceadv
         write (iunit) gicetherm
         write (iunit) gareadv
         write (iunit) garetherm
c        write (iunit) siceic
c        write (iunit) siceow
c        write (iunit) hmlc
c        write (iunit) sicec
c        write (iunit) sareac
c        write (iunit) qflc
c        write (iunit) gicet
c        write (iunit) garet
c        write (iunit) windvx
c        write (iunit) windvy
c
c restoring heat flux (cfi)
c
#ifdef restore
	 write (iunit) qres
#endif
c
#endif
c
#endif
c
#ifdef q_flux
      write(iunit) qbot
#endif
c
      write(iunit) zmto,zmyce,zmqo,zmsno
c
      return
      end
