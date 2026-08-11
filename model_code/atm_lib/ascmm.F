      subroutine ascmm
c
#include "parm.h"
#include "parmjx.h"
#include "anlocn.h"
#include "ctlblk.h"
#include "gauss.h"
#include "grdfld.h"
#include "iohist.h"
#include "qa.h"
#include "surfc.h"
c
      dimension zinc(il),rainwt(il)
c
#include "zmanl.h"
c
c
      write (6,250)
c
      xi     = 1.0/float(ix)
      wt     = 0.0
      somqao = 0.0
      somsno = 0.0
      somto  = 0.0
      somyce = 0.0
      ycnh   = 0.0
      tonh   = 0.0
      qonh   = 0.0
      snnh   = 0.0
      wtnh   = 0.0
      wtsh   = 0.0
      qosh   = 0.0
      snsh   = 0.0
      tosh   = 0.0
      ycsh   = 0.0
c
c     get zonal means
c
      do 210 j=1,il
       zinc(j)=0.0
       do 220 i=1,ix
        if (ts(i,j) .gt. 0.0) then
          zin(i)=1.0
        else
          zin(i)=0.0
        endif
        zinc(j)=zinc(j)+1.0-zin(i)
220    continue
       if(zinc(j).ne.0.) then
         zmqo (j)=zmqo (j)/(zinc(j)*float(kpero))
         zmsno(j)=zmsno(j)/(zinc(j)*float(kpero))
         zmto (j)=zmto (j)/(zinc(j)*float(kpero))
         zmyce(j)=zmyce(j)/(zinc(j)*float(kpero))
       else
         zmqo (j)=0.0
         zmsno(j)=0.0
         zmto (j)=0.0
         zmyce(j)=0.0
       endif
c
c  print zonal means here
c
       write (6,260) j,zinc(j),zmqo(j),zmsno(j),zmyce(j),zmto(j)
c
  210 continue
c
c     compute global averages
c
      do 200 ll=1,iy
      do 200 jj=1,jlp,il
       l1=max0(ll,jj-ll)
       rainwt(ll) = 0.5*weight(ll)*cosqr(ll)*float(3-nhem)
       somqao     = somqao+zmqo (l1)*rainwt(ll)*zinc(l1)*xi
       somsno     = somsno+zmsno(l1)*rainwt(ll)*zinc(l1)*xi
       somyce     = somyce+zmyce(l1)*rainwt(ll)*zinc(l1)*xi
       somto      = somto +zmto (l1)*rainwt(ll)*zinc(l1)*xi
       wt         = wt  +            rainwt(ll)*zinc(l1)*xi
200   continue
c
c     compute hemispheric averages
      do 470 ll=1,iy
       qonh = qonh+zmqo (ll)*rainwt(ll)*zinc(ll)*xi
       tonh = tonh+zmto (ll)*rainwt(ll)*zinc(ll)*xi
       ycnh = ycnh+zmyce(ll)*rainwt(ll)*zinc(ll)*xi
       snnh = snnh+zmsno(ll)*rainwt(ll)*zinc(ll)*xi
       wtnh = wtnh+           rainwt(ll)*zinc(ll)*xi
       l1=ilp-ll
       qosh = qosh+zmqo (l1)*rainwt(ll)*zinc(l1)*xi
       tosh = tosh+zmto (l1)*rainwt(ll)*zinc(l1)*xi
       ycsh = ycsh+zmyce(l1)*rainwt(ll)*zinc(l1)*xi
       snsh = snsh+zmsno(l1)*rainwt(ll)*zinc(l1)*xi
       wtsh = wtsh+           rainwt(ll)*zinc(l1)*xi
 470  continue
c
      do 280 j=1,il
       zmqo (j)=0.0
       zmto (j)=0.0
       zmyce(j)=0.0
       zmsno(j)=0.0
 280  continue
      somqao = somqao/(wtnh+wtsh)
      somto  = somto /(wtnh+wtsh)
      somyce = somyce/(wtnh+wtsh)
      somsno = somsno/(wtnh+wtsh)
      qonh = qonh/wtnh
      qosh = qosh/wtsh
      tonh = tonh/wtnh
      tosh = tosh/wtsh
      ycnh = ycnh/wtnh
      ycsh = ycsh/wtsh
      snnh = snnh/wtnh
      snsh = snsh/wtsh
      wtg  = wtnh+wtsh
c
      write(6,485) qosh,tosh,ycsh,snsh
      write(6,495) qonh,tonh,ycnh,snnh
      write(6,475) wtnh,wtsh,wtg
      write(6,465) somto , somyce
      write(6,466) somqao, somsno
c
c     write (intgl,100) tosh, ycsh, qosh, snsh
c     write (intgl,100) tonh, ycnh, qonh, snnh
c     write (intgl,100) somto,somyce,somqao,somsno
c
      return
c
c100  format (4e16.8)
 250  format(6x,'j',5x,'zinc',5x,'zmqo',16x,'zmsno',15x,'zmyce',
     *       15x,'zmto')
260   format(5x,i2,5x,f4.1,5e20.8)
485   format(10x,'southern hemi values for qo,to,yce, and sno resp',
     *                        4(1x,e16.8))
495   format(10x,'northern hemi values for qo,to,yce, and sno resp',
     *                        4(1x,e16.8))
475   format(1x,' the sum of northern hemi weights=',e16.8,
     *          '  the south hemi weights=',e16.8,
     *          ' and their sums=',e16.8)
465   format(11x,'global mean to=',f11.6,'global mean yce=',e16.8)
466   format(11x,'global mean qo=',f11.6,'global mean sno=',e16.8)
c
c
      end
