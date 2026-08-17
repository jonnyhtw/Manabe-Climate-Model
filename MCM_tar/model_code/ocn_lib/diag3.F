      subroutine diag3
      implicit real*8 (a-h,o-z)
cfpp$ noconcur r
c
#ifdef kdtermb
c
#include "param.h"
#include "cdiag.h"
#include "coord.h"
#include "cregin.h"
#include "ctmngr.h"
#include "docnam.h"
#include "grdvar.h"
#include "iounit.h"
# if defined cpldleg3 || defined cpldleg2
# include "switch.h"
# endif
c
c=======================================================================
c
c     time average regional term balances and write them out
c
c     instantaneous values were written to internal unit "iointb" on 
c     diagnostic timesteps by subroutine diag2.  here those values are 
c     read in, time averaged, and written out to unit "iotrm2"
c
c=======================================================================
c
c     locally dimension arrays used in time averaging term balances
c
      dimension termdt(9,nt,0:numreg),termdm(13,2,numreg)
      dimension atwvd(nt,0:numreg), asstd(nt,0:nhreg), avwadv(3,numreg)
      dimension stflxd(nt,0:nhreg)
c
c     locally dimension array used in time averaging northward tracer
c     transports
# ifdef testnorth
      dimension ttnd(8,jmt,nt,0:2)
      dimension tmtd(jmt,km,0:2)
# else
      dimension ttnd(8,jmt,nt)
# endif
c
      character*32 stampp
c
c     cnterm = count of number of instantaneous time levels being 
c              averaged
c
c=======================================================================
c
c     rewind the internal i/o unit to which instantaneous term balance
c     info was previously written and initialize variables and arrays
c     to be used in the time averaging process
c
      rewind (iointb)
c
      cnterm = c0
c
      do 100 itype=1,13
        do 90 m=1,2
          do 80 n=1,numreg
            termbm(0,itype,m,n) = c0
 80       continue
 90     continue
100   continue
c
# ifdef testnorth
       do 130 j=1,jmt
         do 120 k=1,km
           do 110 nsec=0,2
             tmt(j,k,nsec) = c0
110       continue
120     continue
130   continue
# endif
      do 200 m=1,nt
        do 160 j=1,jmt
          do 150 itype=1,8
# ifdef testnorth
            do 140 nsectr=0,2
             ttn(itype,j,m,nsectr) = c0
140         continue
# else
            ttn(itype,j,m) = c0
# endif
150       continue
160     continue
        do 190 n=0,numreg
          atwv(m,n) = c0
          asst(m,n) = c0
          stflx(m,n)= c0
          do 180 itype=1,9
            termbt(0,itype,m,n) = c0
180       continue
190     continue
200   continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
201   continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     read in term balance info one time level at a time, and
c     accumulate the values
c
      read (iointb,end=499) nnt
c
      cnterm = cnterm + 1.0
c
      read (iointb) iitt, stampp
      if (cnterm .eq. c1) stamp1 = stampp
      read (iointb) (((termdm(itype,m,n),itype=1,13),
     $                         m=1,2),n=1,numreg)
      read (iointb) (((termdt(itype,m,n),itype=1,9),
     $                         m=1,nt),n=0,numreg)
      read (iointb) ((atwvd(m,n),m=1,nt),n=0,numreg)
      read (iointb) ((asstd(m,n),m=1,nt),n=0,nhreg)
      read (iointb) ((stflxd(m,n),m=1,nt),n=0,nhreg)
      read (iointb) ((ustf(m,i),m=1,nt),i=1,2)
c
# ifdef testnorth
      read (iointb) ((((ttnd(i,j,m,nsectr),i=1,8),j=1,jmt),m=1,nt)
     $              ,nsectr=0,2)
      read (iointb) (((tmtd(j,k,nsectr),j=1,jmt),k=1,km),nsectr=0,2)
# else
      read (iointb) (((ttnd(i,j,m),i=1,8),j=1,jmt),m=1,nt)
# endif
c
      do 300 itype=1,13
        do 290 m=1,2
          do 280 n=1,numreg
            termbm(0,itype,m,n) = termbm(0,itype,m,n)
     $                          + termdm(itype,m,n)
280       continue
290     continue
300   continue
c
# ifdef testnorth
       do 330 j=1,jmt
         do 320 k=1,km
           do 310 nsec=0,2
             tmt(j,k,nsec) = tmt(j,k,nsec) + tmtd(j,k,nsec)
310       continue
320     continue
330   continue
# endif
      do 400 m=1,nt
        do 360 j=1,jmt
          do 350 itype=1,8
# ifdef testnorth
            do 340 nsectr=0,2
              ttn(itype,j,m,nsectr) = ttn(itype,j,m,nsectr) +
     $             ttnd(itype,j,m,nsectr)
340         continue
# else
            ttn(itype,j,m) = ttn(itype,j,m) + ttnd(itype,j,m)
# endif
350       continue
360     continue
        do 390 n=0,numreg
          atwv(m,n) = atwv(m,n)  + atwvd(m,n)
          asst(m,n) = asst(m,n)  + asstd(m,n)
          stflx(m,n)= stflx(m,n) + stflxd(m,n)
          do 380 itype=1,9
            termbt(0,itype,m,n) = termbt(0,itype,m,n)
     $                          + termdt(itype,m,n)
380       continue
390     continue
400   continue
c
      goto 201 
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
499   continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      rewind (iointb)
c
c     time average the term balances
c
      if (cnterm .EQ. c0) goto 799
c
      do 600 itype=1,13
        do 590 m=1,2
          do 580 n=1,numreg
            termbm(0,itype,m,n) = termbm(0,itype,m,n) / cnterm
580       continue
590     continue
600   continue
c
# ifdef testnorth
       do 630 j=1,jmt
         do 620 k=1,km
           do 610 nsec=0,2
             tmt(j,k,nsec) = tmt(j,k,nsec) / cnterm
610       continue
620     continue
630   continue
c
# endif
      do 700 m=1,nt
        do 660 j=1,jmt
          do 650 itype=1,8
# ifdef testnorth
            do 640 nsectr=0,2
              ttn(itype,j,m,nsectr) = ttn(itype,j,m,nsectr) / cnterm
640         continue
# else
            ttn(itype,j,m) = ttn(itype,j,m) / cnterm
# endif
650       continue
660     continue
        do 690 n=0,numreg
          atwv(m,n) = atwv(m,n)  / cnterm
          asst(m,n) = asst(m,n)  / cnterm
          stflx(m,n)= stflx(m,n) / cnterm
          do 680 itype=1,9
            termbt(0,itype,m,n) = termbt(0,itype,m,n) / cnterm
680       continue
690     continue
700   continue
c
799   continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     write out time averaged term balances to external disk unit
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c   set limits to prevent ieee conversion underflow & overflow
c
      itype = 13
      length = (km+1) * itype * 2 * numreg
      call ieeeok(termbm, length, .false.)
c
#ifdef kdtermb
      itype = 9
#else
      itype = 8
#endif
      length = (km+1) * itype * nt * (numreg+1)
      call ieeeok(termbt, length, .false.)
c
c
c     write term balance header file
      write(stdout,123) iotrm2
123   format(' WRITING TERM BALANCES TO UNIT iotrm2 =',i3)
c
      write(iotrm2) nt, km, jmt, nhreg, nvreg, cnterm
      write(iotrm2) yu
c
c     write out info on time level and region sizes
c
      write(iotrm2) iitt, stamp1, stampp, expnam, trname
     $,            hregnm, dzt, llvreg
     $,            volt, volu, areau, areat
c
c     write out regional momentum term balances
c
      write(iotrm2) (((termbm(0,itype,m,n),itype=1,13),
     $                         m=1,2),n=1,numreg)
c
c      write out regional tracer term balances
c
      write(iotrm2) (((termbt(0,itype,m,n),itype=1,9),
     $                         m=1,nt),n=0,numreg)
      write(iotrm2) atwv
      write(iotrm2) asst
      write(iotrm2) stflx
      write(iotrm2) ustf
c
      write(iotrm2) ttn
#ifdef testnorth
      write(iotrm2) tmt
#endif
c
      write(stdout,1533) stamp, iotrm2, cnterm, stamp1, stampp
c
# if defined cpldleg3 || defined cpldleg2
      if (prterm) then
# endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     write out formatted time averaged term balances to paper unit
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
       do 1900 n=1,numreg
          write (iotrmb,10050)
     $ 'Regional Term Balances (per sec): Time Avg. Ending on ts =  '
     $,    itt, stamp
            iv = (n-1)/nhreg + 1
            ih = n - (iv-1)*nhreg
            write (iotrmb,10100) 'Momentum terms averaged over region #'
     $,          n, ': ', hregnm(ih), vregnm(iv), volu(n), areau(n)
            write (iotrmb,10104)
            write (iotrmb,10101) n,' dU/dt   = ', termbm(0,1,1,n)
     $,                            ' dV/dt   = ', termbm(0,1,2,n)
            write (iotrmb,10101) n,' -Px     = ', termbm(0,2,1,n)
     $,                            ' -Py     = ', termbm(0,2,2,n)
            write (iotrmb,10101) n,' -surf Px= ', termbm(0,12,1,n)
     $,                            ' -surf Py= ', termbm(0,12,2,n)
            write (iotrmb,10101) n,' -UUx    = ', termbm(0,3,1,n)
     $,                            ' -UVx    = ', termbm(0,3,2,n)
            write (iotrmb,10101) n,' -VUy    = ', termbm(0,4,1,n)
     $,                            ' -VVy    = ', termbm(0,4,2,n)
            write (iotrmb,10101) n,' -WUz    = ', termbm(0,5,1,n)
     $,                            ' -WVz    = ', termbm(0,5,2,n)
            write (iotrmb,10101) n,'  UVtan  = ', termbm(0,13,1,n)
     $,                            ' -UUtan  = ', termbm(0,13,2,n)
#ifdef biharmonic
            write (iotrmb,10101) n,'  Uxxxx  = ', termbm(0,6,1,n)
     $,                            '  Vxxxx  = ', termbm(0,6,2,n)
            write (iotrmb,10101) n,'  Uyyyy  = ', termbm(0,7,1,n)
     $,                            '  Vyyyy  = ', termbm(0,7,2,n)
#else
            write (iotrmb,10101) n,'  Uxx    = ', termbm(0,6,1,n)
     $,                            '  Vxx    = ', termbm(0,6,2,n)
            write (iotrmb,10101) n,'  Uyy    = ', termbm(0,7,1,n)
     $,                            '  Vyy    = ', termbm(0,7,2,n)
#endif
            write (iotrmb,10101) n,'  Uzz    = ', termbm(0,8,1,n)
     $,                            '  Vzz    = ', termbm(0,8,2,n)
            write (iotrmb,10101) n,'  metric = ', termbm(0,9,1,n)
     $,                            '  metric = ', termbm(0,9,2,n)
            write (iotrmb,10101) n,'  fV     = ', termbm(0,10,1,n)
     $,                            ' -fU     = ', termbm(0,10,2,n)
            write (iotrmb,10101) n,'  source = ', termbm(0,11,1,n)
     $,                            '  source = ', termbm(0,11,2,n)
            erru = c0
            errv = c0
            do 1870 lll=2,13
              erru = erru + termbm(0,lll,1,n)
              errv = errv + termbm(0,lll,2,n)
1870        continue
            write (iotrmb,10101) n,'  error  = ', termbm(0,1,1,n) - erru
     $,                            '  error  = ', termbm(0,1,2,n) - errv
            write (iotrmb,10101) n,'  ubar   = ', avwv(1,n)
     $,                            '  vbar   = ', avwv(2,n)
     $,                            '  wbar   = ', avwv(3,n)
c
          if (iv .eq. 1) then
            write (iotrmb,10101) n,'  surf Uz= ', smflx(1,n)
     $,                            '  surf Vz= ', smflx(2,n)
          endif
c
          if (n .eq. 0) then
            write (iotrmb,10051)
     $   'Global averaged (entire ocean) Tracer Term Balances for ts = '
     $,     itt, stamp
          else
            write (iotrmb,10100) 'Tracer terms averaged over region   #'
     $,          n, ': ', hregnm(ih), vregnm(iv), volt(n), areat(n)
          endif
c
          maxm = (nt-1)/7 + 1
          do 1890 mloop=1,maxm
            ms = (mloop-1)*7 + 1
            me = min(ms + 6,nt)
            write (iotrmb,10103) (trname(m),m=ms,me)
            write (iotrmb,10102) n,' dT/dt   = ', (termbt(0,1,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -UTx    = ', (termbt(0,2,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -VTy    = ', (termbt(0,3,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,' -WTz    = ', (termbt(0,4,m,n)
     $,                                           m=ms,me)
#if defined biharmonic && !defined isopycmix
            write (iotrmb,10102) n,'  Txxxx  = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tyyyy  = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
#else
# ifdef isopycmix
            write (iotrmb,10102) n,'  Tisox  = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tisoy  = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
# else
            write (iotrmb,10102) n,'  Txx    = ', (termbt(0,5,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  Tyy    = ', (termbt(0,6,m,n)
     $,                                           m=ms,me)
# endif
#endif
            write (iotrmb,10102) n,'  Tzz    = ', (termbt(0,7,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  source = ', (termbt(0,8,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  convect= ', (termbt(0,9,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  filter+= ', (tfilt(m)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  chg var= ', (tvar(0,m,n)
     $,                                           m=ms,me)
            write (iotrmb,10102) n,'  tbar   = ', (atwv(m,n)
     $,                                           m=ms,me)
            if (iv .eq. 1) then
              write (iotrmb,10102) n,'  surflx = ', (stflx(m,n)
     $,                                           m=ms,me)
              taux = smflx(1,n)*volu(n)/areau(n)
              tauy = smflx(2,n)*volu(n)/areau(n)
              write (iotrmb,10105) ' Regionally averaged quantities:'
              write (iotrmb,'(1x)')
              write (iotrmb,10098) n
     $,          ' smf(1)   = ', taux,' dynes/cm**2   '
              write (iotrmb,10098) n
     $,          ' smf(1)   = ', taux,' dynes/cm**2   '
              write (iotrmb,10098) n
     $,          ' smf(2)   = ', tauy,' dynes/cm**2   '
              do 1820 m=1,nt
                write (iotrmb,10098) n
     $,          ustf(m,1), stflx(m,n)*volt(n)/areat(n), ustf(m,2)
1820          continue
              if (ms .eq. 1) then
                write (iotrmb,10098) n
     $,          ' tot heat = ', (atwv(1,n)*volt(n)),' deg C * cm**3 '
                write (iotrmb,10098) n
     $,          ' sst      = ', asst(1,n),' deg C         '
              endif
            endif
1890      continue
1900    continue
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
c
c-----------------------------------------------------------------------
c       print the northward transport of tracers
c-----------------------------------------------------------------------
c
# ifdef testnorth
        do 2150 nsectr=0,2
        write(stdout,1217) nornam(nsectr)
1217    format(//40x,' Northward Transports ==== ',a16)
#endif
        ntimes = (nt+1)/2
        do 2100 ntim=1,ntimes
          nt1 = 1 + (ntim-1)*2
          nt2 = nt1 + 1
          if (ntim.EQ.1) then
            write (stdout,8195)
          else
            write (stdout,8197) nt1, nt2
          endif
c
# ifdef testnorth
          do 2100 j=2,jmtm2
            l=jmt-j
            if (nt2 .gt. nt) then
              write (stdout,8196) l, (ttn(i,l,nt1,nsectr),i=1,8)
            else
              write (stdout,8196) l, (ttn(i,l,nt1,nsectr),i=1,8),
     $                               (ttn(i,l,nt2,nsectr),i=1,8)
            endif
2100      continue
# else
        do 2100 j=2,jmtm2
          l=jmt-j
          if (nt2 .gt. nt) then
            write (stdout,8196) l, (ttn(i,l,nt1),i=1,8)
          else
            write (stdout,8196) l, (ttn(i,l,nt1),i=1,8),
     $                             (ttn(i,l,nt2),i=1,8)
          endif
2100    continue
# endif
# ifdef testnorth
2150    continue
# endif
c
#endif
      return
1533  format(/' REGIONAL TERM BALANCES AVERAGED ON ',a32,
     $        ' AND WRITTEN TO UNIT ',i3/
     $        ' (',f6.0,' TIME LEVELS AVERAGED,',t35,
     $        ' 1st SAMPLE IN MEAN =',a32,/,t35,
     $        'LAST SAMPLE IN MEAN =',a32,' )'/)
8195  format(/,30x,'>> Time Averaged Global Northward Tracer Transports'
     $       ,' <<'/,5x,'northward transport of heat (x10**15 watts)'
     $       ,22x,'northward transport of salt (x10**10 cm**3/sec)',/,
     $       2x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     $       'diffus   total'))
8196  format(i3,8f8.3,1x,8f8.3)
8197  format(/,5x,'northward transport of tracer number ',i2,'    '
     $       ,22x,'northward transport of tracer number ',i2,/,
     $       2x,2(3x,'x mean  x eddy  z mean  z eddy   ekman tot adv  ',
     $       'diffus   total'))
10050 format(///,1x,13('**********'),/1x,a61,i10,a33/)
10051 format(///,1x,a61,i10,a33/)
10098 format(1x,'(',i4,')',a12,(1pe16.8,a15))
10100 format(/1x,a37,i5,a2,a40,a20/1x,'ocean volume =',e16.9,' cm**3'
     $, ', horizontal ocean area =',e16.9,' cm**2'/)
10101 format(1x,'(',i4,')',a11,1pe15.7, 2x, a11,1pe15.7, 2x, a11
     $,      1pe15.7)
10102 format(1x,'(',i4,')',a11,7(1pe15.7))
10103 format(' Region',11x,8a15)
10104 format(' Region')
10105 format(/8x,a32)
10106 format(///,1x,13('**********'))
      end
