      subroutine aquit(nmbr)
c
c subroutine to stop program abnormally through a division by zero.
c
CCCCinclude "parm.h"
c
c
#include "ctlblk.h"
c
c
c
      data a / 0.0 /
c
      write(6,1) nmbr
      iunit  = 6
      ifrmat = 1
ccc   call wrtctl(iunit, ifrmat)
c
#ifdef __sgi
ccc   call trace_back_stack_and_print()
ccc   call abort()
#endif

      b=1.0/a
c
      write(6,3) a,b
c
      stop
10    continue
      return
c
    1 format('1          error stop at statement number ',i9)
    3 format('  in quit.F (aquit), a=',o23,' b=',o23)
c
      end
