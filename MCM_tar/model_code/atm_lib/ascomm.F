      subroutine ascomm
c
#include "parm.h"
#include "parmjx.h"
#include "anal2.h"
#include "ctlblk.h"
#include "gauss.h"
#include "grdfld.h"
#include "qa.h"
#include "surfc.h"
#ifdef ice_dynamics
c
c variables in ice used in climate model (cfi)
c
#include "varice.h"
#endif
c
      dimension sea(ix)
c
c
      if (kper.eq.0) return
      xi    = 1.0/float(ix)
      ak    = 1.0/float(kper)
      aks   = 1.0/float(kper*idt)
      kpert = kper
      kper  = 0
      do 26 j=1,il
      do 26 i=1,ix
       qao(i,j)=qao(i,j)*ak
       sno(i,j)=sno(i,j)*aks
#ifdef ice_dynamics
c
c correctly calculate net heat flux over ice and open water (cfi)
c
       qfice(i,j) = qfice(i,j)*ak
       qfow(i,j) = qfow(i,j)*ak
c
c calculate average wind stress vectors for ice model (cfi)
c
       tauxic(i,j) = tauxic(i,j)*ak
       tauyic(i,j) = tauyic(i,j)*ak
       windvx(i,j) = windvx(i,j)*ak
       windvy(i,j) = windvy(i,j)*ak
#endif
 26   continue
c
      write (6,101) day,mon,kpert
  101 format('0ocean b.c. data computed in ascomm at day=',
     *       f9.3,'  mon=',i3,'  kper=',i3)
c
c  compute global mean heat flux and snow flux into the ocean
c   and store into anal2 for the anals code.
c
      tqao = 0.0
      tsno = 0.0
      twt  = 0.0
c
      do 10 j=1,il
       tqoj = 0.0
       tsnj = 0.0
       seaj = 0.0
       ls = min0(j,2*iy+1-j)
       wt = 0.5*weight(ls)*cosqr(ls)*float(3-nhem)
       do 11 i=1,ix
        if (ts(i,j).gt.0.0) then
          sea(i) = 0.0
        else
          sea(i) = 1.0
        endif
        tqoj = tqoj + qao(i,j)*sea(i)
        tsnj = tsnj + sno(i,j)*sea(i)
        seaj = seaj +          sea(i)
   11  continue
       tqao = tqao + tqoj*xi*wt
       tsno = tsno + tsnj*xi*wt
       twt  = twt  + seaj*xi*wt
   10 continue
      tqao = tqao / twt
      tsno = tsno / twt
      samqao = samqao + tqao
      samsno = samsno + tsno
c
      return
      end
