      subroutine indata
c
c  called from: atmos
c  entry point: rstart
c  subroutines called: swell, march (from entry rstart)
c
c  indata stores the model data at the beginning of the run by either
c  reading a restart tape or generating initial data.
c
c
#include "parm.h"
#include "parmjx.h"
#include "parmkp.h"
#include "parmout.h"
#include "parmdsk.h"
c
      complex vor  , div  , t  , psln
      complex vorm , divm , tm , psm
      complex vordt, divdt, tdt, psdt
      complex qmix , qmixm, qmxdt
      complex ucos , vcos , psdx, psdy, zs, ek
#ifdef rad
      complex radten
#endif
c
#ifdef seanal
      logical save
#endif
c
#ifdef force_month_end
      character*8 msave
c
#endif
c
c
#include "anal1.h"
c
#include "anal2.h"
c
#ifdef covcld
#include "cloudcx.h"
#endif
c
#include "cloudp.h"
c
#include "ctlblk.h"
c
#include "grdfld.h"
c
#include "gauss.h"
c
#include "gridw.h"
c
#include "iohist.h"
c
#include "lats.h"
c
#include "mrstr.h"
c
#include "nxstps.h"
c
#include "param.h"
c
#include "qa.h"
c
#include "radten.h"
c
#include "sfharm.h"
c
#include "surfo.h"
c
#ifdef print_zstar
      dimension izstar(ix,il),la(ix),lb(ix)
#endif
      dimension kocen(iclt)
c
#include "water.h"
c
#include "wimplt.h"
c
#ifdef force_month_end
# ifdef seanal
      data months / 0 /
# endif
#endif
#ifdef seanal
      data save /.true./
#endif
      data in/'n'/,is/'s'/
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c  zero vertical velocities to insure that boundary values vanish
      do 35 k=1,kp
       do 35 i=1,ix
#ifdef impcit
        gwexpt(i,k)=0.0
#endif
        gwdot(i,k)=0.0
 35   continue
c
   36 continue
      call rdctl(itap)
      bndry  = .false.
      hybrid = .false.
      diurnl = .false.
#ifdef annlmn
      seasnl = .false.
#endif
      swamp  = .false.
c
#ifdef dirnal
      diurnl = .true.
#endif
#ifdef seanal
      seasnl = .true.
#endif
c
      call rdatm(itap)
c
      write(6,9005) tmean
      write(6,9006) tmnshr
c
      write(6,61) ktau,irt
c
c continue reading restart data if necessary until desired ktau is found
c
      if(ktau.ne.irt) go to 36
c
c  update ocean variables in control block with those from indato
      ktauo  = ktauo2
      ifcsto = ifcst2
      idto   = idto2
      kpero  = kpero2
      mono   = mono2
      oyear  = oyear2
      apo    = apo2
      ito    = ito2
      its    = its2
c
      bndry  = .false.
      hybrid = .false.
      diurnl = .false.
#ifdef annlmn
      seasnl = .false.
#endif
      swamp  = .false.
c
#ifdef dirnal
      diurnl = .true.
#endif
#ifdef seanal
      seasnl = .true.
#endif
#ifdef mixed_layer
      mlayer = .true.
#endif
#ifdef fixed_sst
      mlayer = .false.
#endif
      call aspace (intgl)
c
#ifdef print_zstar
c
c  for land-sea map at beginning of run envoke this ifdef
      write (6,121) karray(1)
      do 123 i=1,ix
       lones(i) = i/10
       ltens(i) = mod(i,10)
 123  continue
c     izstar=0
      do 83 j=1,il
       do 83 i=1,ix
        if (ts(i,j) .lt. 0.0) then
          izstar(i,j) = 1
        else
          izstar(i,j) = 0
        endif
   83 continue
      ixmax = min0 ( ix , 120 )
      do 119 j=1,il
       jl   = min0(j,il+1-j)
       lat  = in
       if(j.gt.iy) lat = is
       write (6,120) j, ylat(j), lat, (izstar(i,j),i=1,ixmax)
 119  continue
      write (6,122) lones
      write (6,125) ltens
#endif
c
c
      return
c
  120 format(1x,i3,f6.2,a1,1x,120i1)
  121 format(//11x,' spectral model land-sea points (truncation =',i3,
     *  ')'////)
  122 format(/12x,120i1)
  125 format(12x,120i1)
c
c
      entry rstart
c     ===== ======
c
c  write restart data for tape.
c  rstart is called from atmos.
c
      lrs   = .false.
#ifdef seanal
      lrad  = .false.
      save  = lrad
# ifdef force_month_end
      msave = month
# endif
c
c  compute current day, month, year of the model for the control block
c
      call march
c
      lrad = save
# ifdef force_month_end
      if ( month  .ne. msave )  months = months + 1
      if ( nmonth .ne. 0     )  then
        if(months .eq. nmonth)  ktaust = ktau
      endif
# endif
#endif
c
c
      call wrtctl(irs, 0)
c
      call wrtatm(irs)
c
c
      write(6,60) irs,ktau,day
c
c
      if(ktau .ne. ktaust) then
c
        return
c
      else
c
c rewind history unit at end of run
        rewind irs
c
c  write restart data on irsf 
        call wrtctl(irsf, 0)
c
        call wrtatm(irsf)
c
        write(6,61) irsf,ktau
c
c
      endif
c
      return
c
   60 format(/,' history file written on unit ',i2,' at tau = ',i8,
     *          '  day = ',f11.3/)
   61 format(/' restart file read at tau = ',i10,
     *       '     the timestep specified  = ',i10/)
   70 format(' restart file written on unit ',i2,' at tau = ',i8,
     *          '  day = ',f11.3/)
9005  format(/,' tmean  =',14f8.2)
9006  format(  ' tmnshr =',15f8.2)
c
      end
