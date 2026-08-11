      program readr
c
c  This program reads the control block and outputs the
c  ktau, month, and year into 3 files.
c
       character*3 tempmon
c
#include "ctlblk.h"
c

       read (5,410) iunit
       write(6,411) iunit
  410  format(i4)
  411  format(' fortran unit to be read is ',i10)
c
c  call rdctl to read control block
c
       call rdctl(iunit)
c
       write(6,400) ktau, month, iyear
  400  format(' ktau = ',i10,'  month =',a8,'  iyear =',i10)
c
c  write ktau, month, iyear to units to communicate with shell.
c
       open(21,form='formatted',file='irt.data')
       open(22,form='formatted',file='month.data')
       open(23,form='formatted',file='iyear.data')
c
       write(21,401) ktau
  401  format(i10)
       tempmon = month(1:3)
       write(22,402) tempmon
  402  format(a3)
       write(23,403) iyear
  403  format(i10)
c
       stop
       end
