      subroutine bcocn
#if defined cpldleg3 || defined cpldleg1
c
#include "parm.h"
c
#include "parmjx.h"
c
      parameter ( imt=2*ix+2, jmt=il+2 )
c
#include "anal2.h"
c
#include "ascomo.h"
c
#include "ctlblk.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "iohist.h"
c
#include "intbc.h"
c
#include "qa.h"
c
#ifdef ice_dynamics
c
c variables in ice used in climate model (cfi)
c
#include "varice.h"
#endif
c
      dimension zqao(il), zper(il), zsno(il), ztx(il), zty(il)
#ifdef windbc
      dimension zwind(il), zwin3(il)
#endif
      dimension zinc(il), wt(il)
      dimension bci(imt,jmt,nbo)
      dimension seamk(imt,jmt), seamko(imt,jmt)
      dimension zin(ix,il)
#ifdef cpldleg1
      dimension bcout(imt,jmt)
#endif
#ifdef annlmn
      dimension bcold(imt,jmt,nbo)
#endif
c
      character*16 jform
      data jform /'     j      zin '/
c                  123456789a123456
#if !defined cpldleg1
      character*16 bform
      data bform /'                '/
#endif
c
c  <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
      if (kper .eq. 0) return
#ifdef cpldleg1
      if (kpermo .eq. 0) return
#endif
c
      aix = 1.0 / float(ix)
      ak  = 1.0 / float(kper)
      aks = 1.0 / float(kper*idt)
c
#ifdef cpldleg1
c   accumulate data for monthly averages in arrays qaot, etc.
c    use arrays qao, etc to store daily means for the history tape
      do 50 j=1,il
        do 50 i=1,ix
          qaot (i,j) = qaot (i,j) + qao (i,j)      
          pert (i,j) = pert (i,j) + per (i,j)      
          snot (i,j) = snot (i,j) + sno (i,j)      
          tauxt(i,j) = tauxt(i,j) + taux(i,j)      
          tauyt(i,j) = tauyt(i,j) + tauy(i,j)      
# ifdef windbc
          windst(i,j) = windst(i,j) + winds(i,j)      
          wind3t(i,j) = wind3t(i,j) + wind3(i,j)      
# endif
50    continue
# ifdef bco48
c   write accumulated monthly boundary data for ocean to unit iflxac
c    when restarting the model at less than 1 month intervals
      rewind iflxac
      write(iflxac) qaot, pert, snot, tauxt, tauyt
# ifdef windbc
     $         ,windst, wind3t
# endif
# endif
c
      write(6,9050) ktau, kper, kpermo
9050  format(' accumulating ocean surface boundary data in bcocn',
     $       ' at ktau=',i10,'  kper=',i7,'  kpermo=',i7/)
#else
      kpert = kper
#endif
      kper  = 0
c
#ifdef cpldleg1
c   compute daily means of ocean boundary data before writing history,
c    then compute zonal and global means for samqao, etc. and return.
#endif
      do 100 j=1,il
        do 100 i=1,ix
          qao (i,j) = qao  (i,j) * ak
          per (i,j) = per  (i,j) * aks
          sno (i,j) = sno  (i,j) * aks
          taux(i,j) = taux (i,j) * ak
          tauy(i,j) = tauy (i,j) * ak
#ifdef windbc
          winds(i,j)= winds(i,j) * ak
          wind3(i,j)= wind3(i,j) * ak
#endif
#ifdef ice_dynamics
c
c correctly calculate net heat flux over ice and open water (cfi)
c
          tauxic(i,j) = tauxic(i,j)*ak
          tauyic(i,j) = tauyic(i,j)*ak
          windvx(i,j) = windvx(i,j)*ak
          windvy(i,j) = windvy(i,j)*ak
#endif
  100 continue
c
#ifdef cpldleg1
      if (monend .eq. 1) then
c   compute monthly means of ocean boundary data at the end of the month
        akm  = 1.0 / float(kpermo)
        aksm = 1.0 / float(kpermo*idt)
        do 110 j=1,il
          do 110 i=1,ix
            qaot (I,J) = qaot (I,J) * akm
            pert (I,J) = pert (I,J) * aksm
            snot (I,J) = snot (I,J) * aksm
            tauxt(I,J) = tauxt(I,J) * akm
            tauyt(I,J) = tauyt(I,J) * akm
# ifdef windbc
            windst(i,j)= windst(i,j) * akm
            wind3t(i,j)= wind3t(i,j) * akm
# endif
  110   continue
c
        kpert = kpermo
        kpermo = 0
      endif
#endif
c
#ifdef annlmn
c   for annual mean model, compute the constant for time filter
                         tmcnst = 14.0 * 86400.0                         
      if (day .gt. 28.0) tmcnst = 28.0 * 86400.0                         
c
c   no time filter first 14 days for annual mean model
                         wtnew = float(kpert * idt) / tmcnst             
      if (day .lt. 14.0) wtnew = 1.0                                     
                         wtold = 1.0 - wtnew                             
      write(6,9015) tmcnst, wtnew, wtold                                  
c
c   for annual mean model, read previous (not interpolated) boundary
c    data and apply time filter using wtold and wtnew
c
      rewind iflxan                                                          
      read (iflxan) bcold                                                    
      do 220 j=1,jmt                                                  
        do 220 i=1,imt                                                
          qao (i,j) = (qao (i,j) * wtnew) + (bcold(i,j,1) * wtold)        
          per (i,j) = (per (i,j) * wtnew) + (bcold(i,j,2) * wtold)        
          sno (i,j) = (sno (i,j) * wtnew) + (bcold(i,j,3) * wtold)        
          taux(i,j) = (taux(i,j) * wtnew) + (bcold(i,j,4) * wtold)        
          tauy(i,j) = (tauy(i,j) * wtnew) + (bcold(i,j,5) * wtold)        
# ifdef windbc
          winds(i,j)= (winds(i,j)* wtnew) + (bcold(i,j,6) * wtold)        
          wind3(i,j)= (wind3(i,j)* wtnew) + (bcold(i,j,7) * wtold)        
# endif
220   continue                                                           
      rewind iflxan                                                          
      write (iflxan) qao, per, sno, taux, tauy                               
# ifdef windbc
     $          ,winds, wind3
# endif
#endif
c
c   compute zonal means   
c
      do 120 j=1,il
        zinc(j) = 0.0
        zqao(j) = 0.0
        zper(j) = 0.0
        zsno(j) = 0.0
        ztx (j) = 0.0
        zty (j) = 0.0
#ifdef windbc
        zwind(j)= 0.0
        zwin3(j)= 0.0
#endif
        do 140 i=1,ix
         if (ts(i,j) .gt. 0.0) then
           zin(i,j) = 1.0
         else
           zin(i,j) = 0.0
         endif
140     continue
c
c   note that seamk is stored south-to-north for ocean grid
        do 150 i=1,ix
          seamk(2*i+1,jmt-j) = zin(i,j)
          seamk(2*i  ,jmt-j) = zin(i,j)
          zinc(j) = zinc(j) + (1.0 - zin(i,j))
150     continue
        seamk(  1,jmt-j) = seamk(imt-1,jmt-j)
        seamk(imt,jmt-j) = seamk(    2,jmt-j)
c
#ifdef cpldleg1
        if (monend .eq. 0) then
#endif
        do 160 i=1,ix
          zqao(j) = zqao(j) + (qao (i,j) * (1.0 - zin(i,j)))
          zper(j) = zper(j) + (per (i,j) * (1.0 - zin(i,j)))
          zsno(j) = zsno(j) + (sno (i,j) * (1.0 - zin(i,j)))
          ztx (j) = ztx (j) + (taux(i,j) * (1.0 - zin(i,j)))
          zty (j) = zty (j) + (tauy(i,j) * (1.0 - zin(i,j)))
# ifdef windbc
          zwind(j)= zwind(j)+ (winds(i,j)* (1.0 - zin(i,j)))
          zwin3(j)= zwin3(j)+ (wind3(i,j)* (1.0 - zin(i,j)))
# endif
160     continue
#ifdef cpldleg1
        endif
#endif
c
#ifdef cpldleg1
        if (monend .eq. 1) then
          do 170 i=1,ix
           zqao(j) = zqao(j) + (qaot (i,j) * (1.0 - zin(i,j)))
           zper(j) = zper(j) + (pert (i,j) * (1.0 - zin(i,j)))
           zsno(j) = zsno(j) + (snot (i,j) * (1.0 - zin(i,j)))
           ztx (j) = ztx (j) + (tauxt(i,j) * (1.0 - zin(i,j)))
           zty (j) = zty (j) + (tauyt(i,j) * (1.0 - zin(i,j)))
# ifdef windbc
           zwind(j)= zwind(j)+ (windst(i,j)* (1.0 - zin(i,j)))
           zwin3(j)= zwin3(j)+ (wind3t(i,j)* (1.0 - zin(i,j)))
# endif
170       continue
        endif
#endif
c
        if (zinc(j) .ne. 0.0) then
          zqao(j) = zqao(j) / zinc(j)
          zper(j) = zper(j) / zinc(j)
          zsno(j) = zsno(j) / zinc(j)
          ztx (j) = ztx (j) / zinc(j)
          zty (j) = zty (j) / zinc(j)
#ifdef windbc
          zwind(j)= zwind(j)/ zinc(j)
          zwin3(j)= zwin3(j)/ zinc(j)
#endif
        endif
c
120   continue
c
      do 180 I=1,imt
        seamk(i,  1) = seamk(i,    2)
        seamk(i,jmt) = seamk(i,jmt-1)
180   continue
c
c   compute global means        
c
      tqao = 0.0
      tper = 0.0
      tsno = 0.0
      ttx  = 0.0
      tty  = 0.0
#ifdef windbc
      twin = 0.0
      twin3= 0.0
#endif
      wtts = 0.0
      do 200 ll=1,iy
        do 200 jj=1,jlp,il
          l1 = max0(ll,jj-ll)
          wt(ll) = 0.5 * weight(ll) * cosqr(ll) * float(3-nhem)
          tqao  = tqao + (zqao(l1) * wt(ll) * zinc(l1) * aix)
          tper  = tper + (zper(l1) * wt(ll) * zinc(l1) * aix)
          tsno  = tsno + (zsno(l1) * wt(ll) * zinc(l1) * aix)
          ttx   = ttx  + (ztx (l1) * wt(ll) * zinc(l1) * aix)
          tty   = tty  + (zty (l1) * wt(ll) * zinc(l1) * aix)
#ifdef windbc
          twin  = twin + (zwind(l1)* wt(ll) * zinc(l1) * aix)
          twin3 = twin3+ (zwin3(l1)* wt(ll) * zinc(l1) * aix)
#endif
          wtts  = wtts + (           wt(ll) * zinc(l1) * aix)
200   continue
      tqao = tqao / wtts
      tper = tper / wtts
      tsno = tsno / wtts
      ttx  = ttx  / wtts
      tty  = tty  / wtts
#ifdef windbc
      twin = twin / wtts
      twin3= twin3/ wtts
#endif
c
c   store global mean heat flux and sno over ocean for anals
c
      samqao = samqao + tqao
      samsno = samsno + tsno
      sampme = sampme + tper
c
#ifdef cpldleg1
c
      if (monend .eq. 0) return
c
#endif
c
      write(6,9010) ktau, day, mon, kpert
c
c   print zonal means        
c
#if !defined cpldleg1
      ndayp = day + 0.01
      if (mod(ndayp,30) .eq. 0) then
#endif
        write(6,9020) jform
        do 210 j=1,il
          write(6,9025)
     $          j, zinc(j), zqao(j), zper(j), zsno(j), ztx(j), zty(j)
#ifdef windbc
     $           ,zwind(j), zwin3(j)
#endif
  210   continue
#if !defined cpldleg1
      else
        write(6,9020) bform
      endif
#endif
c
c   print global means       
c
      write(6,9030) tqao, tper, tsno, ttx, tty
#ifdef windbc
     $             ,twin, twin3
#endif
c
c   interpolate to ocean grid
c
c   lgrid=0 for ocean t,s grid (qao, per, sno);
c   lgrid=1 for u,v grid (taux, tauy)
      lgrid = 0
#ifdef cpldleg1
      call intrpa(qaot,  bci(1,1,1), lgrid, zin, seamk)
      call intrpa(pert , bci(1,1,2), lgrid, zin, seamk)
      call intrpa(snot , bci(1,1,3), lgrid, zin, seamk)
      lgrid = 1
      call intrpa(tauxt, bci(1,1,4), lgrid, zin, seamko)
      call intrpa(tauyt, bci(1,1,5), lgrid, zin, seamko)
# ifdef windbc
      lgrid = 0
      call intrpa(windst,bci(1,1,6), lgrid, zin, seamk)
      call intrpa(wind3t,bci(1,1,7), lgrid, zin, seamk)
# endif
#endif
#ifdef cpldleg3
      call intrpa(qao,  bci(1,1,1), lgrid, zin, seamk)
      call intrpa(per , bci(1,1,2), lgrid, zin, seamk)
      call intrpa(sno , bci(1,1,3), lgrid, zin, seamk)
      lgrid = 1
      call intrpa(taux, bci(1,1,4), lgrid, zin, seamko)
      call intrpa(tauy, bci(1,1,5), lgrid, zin, seamko)
# ifdef windbc
      lgrid = 0
      call intrpa(winds,bci(1,1,6), lgrid, zin, seamk)
      call intrpa(wind3,bci(1,1,7), lgrid, zin, seamk)
# endif
#endif
c
c   compute means on the ocean grid
c
      wtts = 0.0
      wtuv = 0.0
      tqao = 0.0
      tper = 0.0
      tsno = 0.0
      ttx  = 0.0
      tty  = 0.0
#ifdef windbc
      twin = 0.0
      twin3= 0.0
#endif
      do 260 j=2,jmt-1
        do 260 i=2,imt-1
          wtts = wtts + (             wto(j) * (1.0 - seamk(i,j)))
          tqao = tqao + (bci(i,j,1) * wto(j) * (1.0 - seamk(i,j)))
          tper = tper + (bci(i,j,2) * wto(j) * (1.0 - seamk(i,j)))
          tsno = tsno + (bci(i,j,3) * wto(j) * (1.0 - seamk(i,j)))
          ttx  = ttx  + (bci(i,j,4) * cos(oyv(j)) * seamko(i,j))
          tty  = tty  + (bci(i,j,5) * cos(oyv(j)) * seamko(i,j))
#ifdef windbc
          twin = twin + (bci(i,j,6) * wto(j) * (1.0 - seamk(i,j)))
          twin3= twin3+ (bci(i,j,7) * wto(j) * (1.0 - seamk(i,j)))
#endif
          wtuv = wtuv + (             cos(oyv(j)) * seamko(i,j))
260   continue
      tqao = tqao / wtts
      tper = tper / wtts
      tsno = tsno / wtts
      ttx  = ttx / wtuv
      tty  = tty / wtuv
#ifdef windbc
      twin = twin / wtts
      twin3= twin3/ wtts
#endif
      write(6,9035) tqao, tper, tsno, ttx, tty
#ifdef windbc
     $             ,twin, twin3
#endif
c
      do 240 j=1,jmt
        do 240 i=1,imt
c   convert interpolated qao to ocean units: cal/cm2/min to cal/cm2/sec
          bci(i,j,1) =  bci(i,j,1) / 60.
c   change sign of per (positive per must decrease salinity)
          bci(i,j,2) = -bci(i,j,2)
c   change signs of wind stress
          bci(i,j,4) = -bci(i,j,4)
          bci(i,j,5) = -bci(i,j,5)
240   continue
      tqao =  tqao / 60.
      tper = -tper
      ttx  = -ttx
      tty  = -tty
c
#ifdef resetps
c   set per=0. and sno=0. for ocean model in first ** atm. days <-----
c
      if (day .lt. 100.0) then
        do 280 j=1,jmt
          do 280 i=1,imt
            bci(i,j,2) = 0.0
            bci(i,j,3) = 0.0
280     continue
        write(6,9090) day
9090    format(/,'*****set per=0. and sno=0. for ocean at day=',F8.3/)
        tper = 0.0
        tsno = 0.0
      endif
#endif
c
c   print only the integrals having the units (or sign) changed
      write(6,9040) tqao, tper, ttx, tty
c
#ifdef cpldleg3
      do 300 n=1,nbo
        do 300 j=1,jmt
          do 300 i=1,imt
            bco(i,j,n,1) = bci(i,j,n)
300   continue
#endif
c
#ifdef cpldleg1
c   for leg1 of coupled model, write bc data on unit iflux
      rewind iflux
      do 320 n=1,nbo
        do 330 j=1,jmt
          do 330 i=1,imt
            bcout(i,j) = bci(i,j,n)
330     continue
        write(iflux) bcout
320   continue
#endif
c
#ifdef seanal
c   for seasonal model, annual mean wind stress is computed
c    in the ocean model
#endif
c
c
9010  format(/,' Ocean surface boundary data computed at ktau=',i10,
     $            '  day=',f12.3,'  month=',i3,'  kper=',i6,/)
#ifdef annlmn
9015  format(55x,'tmcnst=',f12.3,'  wtnew=',f7.4,'  wtold=',f7.4)  
#endif
9020  format(a16,5x,'qao ',11x,'per ',11x,'sno ',
#ifndef windbc
     $          11x,'taux',11x,'tauy')
#endif
#ifdef windbc
     $          11x,'taux',11x,'tauy',11x,'wind',11x,'wind**3')
#endif
9025  format(     5x,i2,5x,f4.1,7e15.7)
9030  format(' on atmos grid: ',7e15.7)
9035  format(' on ocean grid: ',7e15.7)
9040  format(' in ocean units:',2e15.7,15x,2e15.7,/)
c
#endif
      return
      end
