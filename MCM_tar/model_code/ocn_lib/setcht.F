      subroutine setcht (vname, units, hgrid, vgrid, klev, nv)
      implicit real*8 (a-h,o-z)
#ifdef selfdoc
c
#include "param.h"
#include "histjr.h"
c
      character * 16 vname, units
      character * 2  hgrid, vgrid
c
c---------------------------------------------------------------------
c
      nvtot = iblk
      nvchk = nvchk + 1
c
      if (nvchk .GT. nvtot) then
        write(stdout) ' SELFDOC PROBLEM: nvchk, nvtot =',nvchk,nvtot
        write(stderr) ' SELFDOC PROBLEM: nvchk, nvtot =',nvchk,nvtot
        stop 'setcht1'
      endif
c
      if (nvchk .NE. nv) then
        write(stdout) ' SELFDOC PROBLEM: nvchk, nv =',nvchk,nv
        write(stdout) ' vname, units, klev =',vname,units,klev
        write(stderr) ' SELFDOC PROBLEM: nvchk, nv =',nvchk,nv
        write(stderr) ' vname, units, klev =',vname,units,klev
        stop 'setcht2'
      endif
c
      varnam(nv) = vname
      varunt(nv) = units
      varhgr(nv) = hgrid
      varvgr(nv) = vgrid
      klevar(nv) = klev
c
ccccc      write(6,66) nv, vname, units, hgrid, vgrid, klev
c
#endif
      return
ccccc66    format(1x,i4,2x,a16,2x,a16,2x,a2,2x,a2,2x,i3)
      end
