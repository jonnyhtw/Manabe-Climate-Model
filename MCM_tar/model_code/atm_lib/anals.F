       subroutine anals
c
c this subroutine computes running annual means for :
c     sea sfc temp       = samto,
c     sea yce            = samyce,
c     ocn heat flux      = samqao,
c     short wave top     = samsw,
c     long wave top      = samlw,
c     snow flux into ocn = samsno,
c     mean atm temp      = samta
#if defined cpldleg3 || defined cpldleg1
c     pme flux into ocn  = sampme
#endif
c  from the atmosphere side of the airsea model (sam prefix).
c
#include "parm.h"
#include "anal1.h"
c
#include "anal2.h"
c
#include "ctlblk.h"
c
#include "iohist.h"
c
c
c if kanal = 0 return
c this means problems
      if (kanal .eq. 0) return
c
c compute monthly means
      rmt    = 1.0/float(kanal)
      samta  = samta *rmt
      samsw  = samsw *rmt
      samlw  = samlw *rmt
      samto  = samto *rmt
      samqao = samqao*rmt
      samyce = samyce*rmt
      samsno = samsno*rmt
#if defined cpldleg3 || defined cpldleg1
      sampme = sampme*rmt
#endif
c
      write(6,405) kanal,samta ,samsw ,samlw,samto,
     *                   samyce,samsno,samqao
#if defined cpldleg3 || defined cpldleg1
     $                  ,sampme
#endif
c
      write (intgl,415) day, mon, iyear
c
      write (intgl,416) samta , samsw , samlw, samto,
     *                  samyce, samqao, samsno
#if defined cpldleg3 || defined cpldleg1
     $                 ,sampme
#endif
c
c store global monthly averages into proper month
      toan (mon) = samto
      ycean(mon) = samyce
      qoan (mon) = samqao
      snoan(mon) = samsno
      swan (mon) = samsw
      rlwan(mon) = samlw
      taan (mon) = samta
#if defined cpldleg3 || defined cpldleg1
      pmean(mon) = sampme
#endif
c
c  compute annual means
      ttoa = 0.0
      tyce = 0.0
      tqo  = 0.0
      tsno = 0.0
      tsw  = 0.0
      tlw  = 0.0
      tta  = 0.0
#if defined cpldleg3 || defined cpldleg1
      tpme = 0.0
#endif
      rmt=1.0/float(mnth)
      do 10 i=1,mnth
       ttoa = ttoa + toan (i)*rmt
       tyce = tyce + ycean(i)*rmt
       tsw  = tsw  + swan (i)*rmt
       tlw  = tlw  + rlwan(i)*rmt
       tqo  = tqo  + qoan (i)*rmt
       tta  = tta  + taan (i)*rmt
       tsno = tsno + snoan(i)*rmt
#if defined cpldleg3 || defined cpldleg1
       tpme = tpme + pmean(i)*rmt
#endif
10    continue
c
      write(6,410) ttoa,tqo,tyce,tsno,tsw,tlw,tta
#if defined cpldleg3 || defined cpldleg1
     $            ,tpme
#endif
c
      write (intgl,416)  tta, tsw, tlw, ttoa, tyce, tqo, tsno
#if defined cpldleg3 || defined cpldleg1
     $                  ,tpme
#endif
c
c zero summed vbls
      kanal  = 0
      samta  = 0.0
      samsw  = 0.0
      samlw  = 0.0
      samyce = 0.0
      samto  = 0.0
      samqao = 0.0
      samsno = 0.0
#if defined cpldleg3 || defined cpldleg1
      sampme = 0.0
#endif
c
c  ingrement atm month counter
      mon = mon+1
      if (mon .eq. 13) mon = 1
c
      return
c
  405 format(1x,' kanal=',i3,' averaged means: tbar=',f8.4,' swtop=',
     * f8.6,' lwtop=',f8.6,' tocn=',f8.4,' tyce=',f8.4,/,30x,
#if defined cpldleg3 || defined cpldleg1
     * ' tsno=',e16.8,' tqo=',e16.8,' tpme=',e16.8)
#else
     * ' tsno=',e16.8,' tqo=',e16.8)
#endif
  410 format(/,1x,' the running annual means(ie, the last 12 averages) '
     * ,/15x,' toan=',f8.4,/,15x,' qoan=',e16.8,/,15x,'ycean=',
     * e16.8,/,15x,'snoan=',e16.8,/,15x,' swan=',e16.8,/,15x,' lwan=',
#if defined cpldleg3 || defined cpldleg1
     * e16.8,/,14x,'tbaran=',f8.4,/,14x,' pmean=',e16.8)
#else
     * e16.8,/,14x,'tbaran=',f8.4)
#endif
  415 format (f10.1,2i10)
#if defined cpldleg3 || defined cpldleg1
  416 format (10x,8e15.7)
#else
  416 format (10x,7e16.8)
#endif
c
      end
