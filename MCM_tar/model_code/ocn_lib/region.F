      subroutine region (nunit, wrform, wrtr, wrk)
      implicit real*8 (a-h,o-z)
#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     subroutine region is provided to allow the writing out of area
c     and volume weighted averages of tracers and surface tracer fluxes
c     under user defined horizontal regions within the model grid.
c     these quantities are computed via calls to subroutine regnj.
c     (see "cregin.h" for more details on variables, and "reg1st.F" for
c     set up and horizontal mask i/o considerations).
c
c     this version of the tracer average code is markedly different
c     from the original MOM version in some respects.  this version is
c     set up to produce output compatible with k dixon's graphical
c     analysis programs.  as a results of these changes, this routine
c     has been rewritten, and the changes do not appear as isolated
c     ifdef options.  the ifdef option "kdtravg" must be defined to
c     get these changes.                                  (kd: feb 93)
c
c     the region averages are calculated at intervals of "travg" days,
c     and are written to unit "iotavg" or a user specified unit passed
c     into this routine as "nunit".
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#ifdef cpldleg3
#include "cinlnd.h"
#endif
#include "cregin.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "docnam.h"
#include "grdvar.h"
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#include "slabs.h"
#include "switch.h"
c
c  nunit = unit to be written to or read from
c  wrform= true(false) switch for formatted(unformatted) writes
c  wrtr  = switch to write tracer information
c  wrk   = switch to write k-level averages
c  sumbt = volume weighted tracer sum for a given region
c  avgbt = volume weighted tracer average for a given region
c  sumbk = volume weighted tracer sum for a region at a k level
c  avgbk = volume weighted tracer average for a region at a k level
c  sumgt = total global volume wtd tracer sum
c  avggt = total global volume wtd tracer average
c  sumgk = total global volume wtd tracer sum at a k level
c  avggk = total global volume wtd tracer average at a k level
c  sumbf = area weighted tracer flux sum for a given region
c  avgbf = average area weighted tracer flux for a region
c  sumgf = global area weighted tracer flux sum
c  avggf = global average area weighted tracer flux
c

      logical  wrform, wrtr, wrk
c
c-----------------------------------------------------------------------
c       write out regional tracer means
c-----------------------------------------------------------------------
c
      if (wrtr) then
        if(wrform) then
         npass = nhreg/9
         if (mod(nhreg,9) .GT. 0) npass = npass + 1
ccccc     if (nhend .gt. nhreg) nhend=nhreg
          do 800 m=1,nt
           nhst = 1
           do 795 np=1,npass
            nhend= nhst + 8
            if (nhend .GT. nhreg) nhend = nhreg
            write(nunit,9004) trname(m), itt, stamp
            if(wrk) then
              if (nhst .eq. 1) then
                write(nunit,9001) (mask,mask=nhst,nhend)
              else
                write(nunit,9101) (mask,mask=nhst,nhend)
              endif
              do 790 k=1,km
                if (nhst .eq. 1) then
                  write(nunit,9002) k, avggk(k,m),
     $                           (avgbk(mask,k,m),mask=nhst,nhend)
                else
                  write(nunit,9102) k,
     $                           (avgbk(mask,k,m),mask=nhst,nhend)
                endif
790           continue
            endif
            if (nhst .eq. 1) then
              write(nunit,9003) avggt(m),(avgbt(mask,m),mask=nhst,nhend)
              write(nunit,9014) avggf(m),(avgbf(mask,m),mask=nhst,nhend)
            else
              write(nunit,9103) (avgbt(mask,m),mask=nhst,nhend)
              write(nunit,9114) (avgbf(mask,m),mask=nhst,nhend)
            endif
            nhst = nhend + 1
795        continue
800       continue
#ifdef seaice
c
          nhst = 1
          do 805 np=1,npass
           nhend= nhst + 8
           if (nhend .GT. nhreg) nhend = nhreg
           write(nunit,9004) ycname, itt, stamp
             if (nhst .eq. 1) then
               write(nunit,9001) (mask,mask=nhst,nhend)
             else
               write(nunit,9101) (mask,mask=nhst,nhend)
             endif
             k = 1
             if (nhst .eq. 1) then
               write(nunit,9002) k, avggi,
     $                        (avgbi(mask),mask=nhst,nhend)
             else
               write(nunit,9102) k,
     $                        (avgbi(mask),mask=nhst,nhend)
             endif
             write(nunit,9006)
            nhst = nhend + 1
805       continue
c
#endif
        else
c
c    write out some extra header info for unformatted files 
c
          nnt = nt
          kkm = km
          nnhreg = nhreg
c
          write(nunit) nnt, kkm, nnhreg
          write(nunit) itt, stamp, trname, hregnm, expnam,
     $      dzt, volbk, volbt, volgk, volgt, iceflg
          if(wrk) write(nunit)  avggk, avgbk
          write(nunit) avggt, avgbt
          write(nunit) avggf, avgbf
#ifdef seaice
c --------------------------------------------------------------
c   set limits to prevent ieee conversion underflow & overflow
c --------------------------------------------------------------
c
          call ieeeok(avggi,     1, .false.)
          call ieeeok(avgbi, nhreg, .false.)
c
          write(nunit) ycname, avggi, avgbi
#endif
#ifdef cpldleg3
ccccMEDtest
ccccc     write(nunit) areamed, vinlnd(1)
ccccc     write(nunit) delmed, dtrdt
ccccMEDtest
#endif
#ifdef worldocn
          write(stdout,9007) nunit, itt, stamp
#endif
        endif
      endif
c
      return
c
9001  format('    k','  All Regions ',9(1x,i6,5x))
9101  format('    k','      Regions:',9(1x,i6,5x))
9002  format(1x,i4,10(1x,e11.5))
9102  format(1x,i4,12x,9(1x,e11.5))
9003  format('  SUM',10(1x,e11.5))
9103  format('  SUM',12x,9(1x,e11.5))
9004  format(/' Volume Weighted Averages for ',a12,' on ts =',i10,a33)
#ifdef seaice
9006  format(/)
#endif
9007  format(/,1x,'Regional Tracer Averages written to unit ',i2,
     $       ' on ts =',i10,' (',a32,')'/)
9014  format(/' FLUX',10(1x,e11.5),/)
9114  format(/' FLUX',12x,9(1x,e11.5),/)
      end


      subroutine regnj (j)
      implicit real*8 (a-h,o-z)

#ifdef multitasking
cfpp$ noconcur r
#endif
c
c=======================================================================
c
c     subroutine regnj is provided to allow the computation of area
c     and volume weighted averages of tracers and surface tracer fluxes
c     under user defined horizontal regions within the model grid.
c     individual averages are calculated for each region as specified by
c     the integer region masks stored in "mskhr", and for each of the
c     "km" model levels.
c     subroutine regnj is called once per model j-row
c     (see "cregin.h" for more details on variables, and "reg1st.F" for
c     set up and horizontal mask i/o considerations).
c
c     the region averages are calculated at intervals of "travg" days,
c     and are written to unit "iotavg".
c
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "cregin.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "docnam.h"
#include "grdvar.h"
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#include "slabs.h"
c
c  sumbt = volume weighted tracer sum for a given region
c  avgbt = volume weighted tracer average for a given region
c  sumbk = volume weighted tracer sum for a region at a k level
c  avgbk = volume weighted tracer average for a region at a k level
c  sumgt = total global volume wtd tracer sum
c  avggt = total global volume wtd tracer average
c  sumgk = total global volume wtd tracer sum at a k level
c  avggk = total global volume wtd tracer average at a k level
c  sumbf = area weighted tracer flux sum for a given region
c  avgbf = average area weighted tracer flux for a region
c  sumgf = global area weighted tracer flux sum
c  avggf = global average area weighted tracer flux
c
c-----------------------------------------------------------------------
c     initialize sums for averages under horizontal regions
c     before processing the first j-row
c-----------------------------------------------------------------------
c
      if (j .eq. 2) then
        do 100 m=1,nt
          do 90 mask=1,nhreg
            sumbf(mask,m) = c0
            do 80 k=1,km
              sumbk(mask,k,m) = c0
              avgbk(mask,k,m) = c0
80          continue
90        continue
100     continue
#ifdef seaice
c
        iceflg = .TRUE.
c
        do 110 mask=1,nhreg
          sumbi(mask) = c0
110     continue
#else
        iceflg = .FALSE.
#endif
      endif
c
c-----------------------------------------------------------------------
c     accumulate sums for averages under horizontal regions
c     (compute in units of meters, rather than cm)
c     do for all valid j-rows
c-----------------------------------------------------------------------
c
      if (j .ge. 2 .and. j .le. jmtm1) then
        do 200 i=2,imtm1
          mask = mskhr(i,j)
          if (mask .ne. 0) then
            boxar = cst(j) * dxt(i) * dyt(j) * fm(i,1,jc) * 0.0001
            do 190 m=1,nt
              sumbf(mask,m) = sumbf(mask,m) + stf(i,m) * boxar
# if defined cpldleg3 && !defined cpld_restore
     $                                      + stfadj(i,m) * boxar
c    for coupled model, lump in terr & serr into regional
c    surface flux summary statistics
# endif
              do 180 k=1,km
                sumbk(mask,k,m) = sumbk(mask,k,m) +
     $               t(i,k,jc,nc,m) * boxar * dzt(k) * fm(i,k,jc) * 0.01
180           continue
190         continue
#ifdef seaice
            sumbi(mask) = sumbi(mask) + yca(i,j) * boxar * 0.01
#endif
          endif
200     continue
      endif
c
c-----------------------------------------------------------------------
c     compute averages under horizontal regions after accumulating all
c     jrows ... prepare values to be written out by subroutine region
c     when processing the last j-row
c-----------------------------------------------------------------------
c
      if (j .eq. jmtm1) then
c
c       initialize sums for horizontal regions & tracer averages
c
        do 300 m=1,nt
          sumgf(m) = c0
          avggf(m) = c0
          sumgt(m) = c0
          avggt(m) = c0
          do 280 mask=1,nhreg
            sumbt(mask,m) = c0
            avgbt(mask,m) = c0
            avgbf(mask,m) = c0
280       continue
          do 290  k=1,km
            sumgk(k,m) = c0
            avggk(k,m) = c0
290       continue
300     continue
#ifdef seaice
c
        sumgi = c0
        avggi = c0
        do 310 mask=1,nhreg
          avgbi(mask) = c0
310     continue
#endif
c
c       compute sums for tracer averages over horizontal regions
c
        do 400 m=1,nt
          do 390 mask=1,nhreg
            sumgf(m) = sumgf(m) + sumbf(mask,m)
            do 380  k=1,km
               sumbt(mask,m) = sumbt(mask,m) + sumbk(mask,k,m)
              sumgk(k,m) = sumgk(k,m) + sumbk(mask,k,m)
380         continue
            sumgt(m) = sumgt(m) + sumbt(mask,m)
390       continue
400     continue
#ifdef seaice
c
       do 410 mask=1,nhreg
         sumgi = sumgi + sumbi(mask)
410    continue
#endif
c
        do 500 k=1,km
          if (volgk(k) .gt. c0) then
            rvolgk = c1 / volgk(k)
            do 490 m=1,nt
              avggk(k,m) = sumgk(k,m) * rvolgk
              do 480 mask=1,nhreg
                if (volbk(mask,k) .gt. c0) then
                  avgbk(mask,k,m) = sumbk(mask,k,m) / volbk(mask,k)
                endif
480           continue
490         continue
          endif
500     continue
c
        rvolgt = c1 / volgt
        rareag = c1 / areag
        do 600 m=1,nt
          avggt(m) = sumgt(m) * rvolgt
          avggf(m) = sumgf(m) * rareag
          do 590 mask=1,nhreg
            if (volbt(mask) .gt. c0) then
              avgbt(mask,m) = sumbt(mask,m) / volbt(mask)
            endif
            if (areab(mask) .gt. c0) then
              avgbf(mask,m) = sumbf(mask,m) / areab(mask)
            endif
590       continue
600     continue
#ifdef seaice
c
        avggi = sumgi * rareag
ccccc   write(stdout,19) avggi, sumgi,  rareag
19      format(/,1x,'avggi, sumgi,  rareag=',3e12.5)
        do 610 mask=1,nhreg
          if (areab(mask) .gt. c0) then
            avgbi(mask) = sumbi(mask) / areab(mask)
          endif
610     continue
#endif
c
      endif
c
      return
      end
