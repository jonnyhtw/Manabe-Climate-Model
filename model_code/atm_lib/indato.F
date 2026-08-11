      subroutine indato
c
c         called from: airsea
c        entry points: rstrto
c  subroutines called: none
c
c  this routine is very similar to indata.  it reads the initial
c   ocean restart tape.  the entry point rstrto, writes the mixed
c   layer ocean data to a history file.
c
#include "parm.h"
c
c
#include "anlocn.h"
c
#include "ctlblk.h"
c
#include "grdfld.h"
c
#include "iohist.h"
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
      namelist / ocean_in /  irto
c
c
      data irto    /  0 /
c
  201 format (/' ocean restart file read in indato at ktau = ',i8,
     *        '  the timestep specified is',i10)
c
      read (5,ocean_in)
      write(6,ocean_in)
c
 100  continue
c
cccc  read ocean data from 55
c
      call rdctl(iunt)
c
      call rdocn(iunt)
c
      write(6,201) ktau,irto
c
      if(irto .ne. ktau) go to 100
      rewind iunt
c
c
c  save ocean variables in control block to use in indata to
c  update them in the ctlblk
      ktauo2  = ktauo
      ifcst2  = ifcsto
      idto2   = idto
      kpero2  = kpero
      mono2   = mono
      oyear2  = oyear
      apo2    = apo
      ito2    = ito
      its2    = its
c
#ifdef mixed_layer
      mlayer = .true.
#endif
c
#ifdef fixed_sst
      mlayer = .false.
#endif
c
c change qflx when qflux model is specified (cfi)
c
#ifdef q_flux
      qflx = .true.
#endif
c
      return
c
      entry rstrto
c     ===== ======
c
c  called from ocean
c
c
      iunit = ibd
c
 510  continue
      call wrtctl(iunit, 0)
c
      call wrtocn(iunit)
c
      if (ktau  .lt. ktaust) return
      if (iunit .eq.   iunt) return
c
      iunit = iunt
      go to 510
c
      end
