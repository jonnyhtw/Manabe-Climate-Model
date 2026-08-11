      subroutine reg1st (nunit, wrform, wrvol, wrk, wrmask, rmask)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     Subroutine reg1st is i/o routine for the user defined horizontal
c     and vertical regional masks used in MOM.
c     It is also the i/o routine that can be used to write volume
c     information for the horizontal regions used in calculating
c     volume weighted averages of tracers and surface tracer fluxes.
c     Both formatted and unformatted i/o is supported.
c     (see "cregin.h" for more details on variables).
c
c     The "mskhr" horizontal masks are used in calculations in
c     subroutines clinic, tracer, and region.
c     The "mskvr" vertical masks are used in calculations in subroutines
c     clinic, tracer, and diag2.
c
c     when "readrmsk" is defined, "mskhr","mskvr", "hregnm" &"vregnm"
c     will be read from file "regionmasks" which is assigned to
c     "iormsk" in the same manner as this routine writes that
c     information. ("iormsk" is set in "blkdta.F")
c     alternatively, the horizontal region masks can be specified in
c     "ocean.F", and although a simple, latitude dependent default
c     setting is provided, the user should consider specifying a more
c     appropriate set of masks.
c
c=======================================================================
c
#include "param.h"
#include "cregin.h"
c
c  nunit = unit to be written to or read from
c  wrform= true(false) switch for formatted(unformatted) writes
c  wrvol = switch to write volume & area information
c  wrk   = switch to write k-level volume information
c  wrmask= switch to write horizontal region masks field
c  rmask = switch to read horizontal region masks field from specified
c          unit
c  settop, setbot = used in defining level limits for vertical regions
c
      logical  wrform, wrvol, wrk, wrmask, rmask
      logical  settop, setbot
c
#ifdef kdtravg
      character *1 chrow(imt), chnum(0:25)
      data (chnum(n),n=0,25) /'0','1','2','3','4','5','6','7','8','9',
     $                        'A','B','C','D','E','F','G','H','I','J',
     $                        'K','L','M','N','O','P'/
c
#endif
      dimension ncol(imt)
c
c-----------------------------------------------------------------------
c     write out regional volume information
c-----------------------------------------------------------------------
c
        if (wrvol) then
          if(wrform) then
            write(nunit,9000)
#ifdef kdtravg
           npass = nhreg/9
           if (mod(nhreg,9) .GT. 0) npass = npass + 1
           nhst = 1
           do 750 np=1,npass
            nhend= nhst + 8
            if (nhend .gt. nhreg) nhend=nhreg
#endif
            if(wrk) then
#ifdef kdtravg
              if (nhst .eq. 1) then
                write(nunit,9001) (mask,mask=nhst,nhend)
              else
                write(nunit,9101) (mask,mask=nhst,nhend)
              endif
              do 700 k=1,km
               if (nhst .eq. 1) then
                write(nunit,9002) k, volgk(k),
     $                           (volbk(mask,k),mask=nhst,nhend)
               else
                write(nunit,9102) k,           
     $                           (volbk(mask,k),mask=nhst,nhend)
               endif
700           continue
#else
              write(nunit,9001) (mask,mask=1,nhreg)
              do 700 k=1,km
                write(nunit,9002) k, volgk(k),
     $                           (volbk(mask,k),mask=1,nhreg)
700           continue
#endif
            endif
#ifdef kdtravg
            if (nhst .eq. 1) then
              write(nunit,9003) volgt, (volbt(mask),mask=nhst,nhend)
              write(nunit,9004) areag, (areab(mask),mask=nhst,nhend)
            else
              write(nunit,9103) (volbt(mask),mask=nhst,nhend)
              write(nunit,9104) (areab(mask),mask=nhst,nhend)
            endif
            nhst = nhend + 1
750       continue
#else 
            write(nunit,9003) volgt, (volbt(mask),mask=1,nhreg)
            write(nunit,9004) areag, (areab(mask),mask=1,nhreg)
#endif
          else
            if(wrk) write(nunit) volgk, volbk
            write(nunit) volgt, volbt
            write(nunit) areag, areab
          endif
        endif
c
c-----------------------------------------------------------------------
c       read in (write out) horizontal & vertical region masks
c       set linel to length of desired formatted printout line
c-----------------------------------------------------------------------
c
        if (wrmask .or. rmask) then
c
          if(wrform) then
            if (wrmask) write(nunit,9011)
            if (wrmask) write(nunit,9012)
     $         (' domain for hor mask #',i,'=',hregnm(i),i=1,nhreg)
            if ( rmask) then
              read (nunit,9099)
              read (nunit,9099)
            endif
            if ( rmask) read (nunit,9013) (hregnm(i),i=1,nhreg)
            linemx = 100
            linel  = 105
            line   = linel - 5
            if (line .gt. linemx) line = linemx
            nwr = (imt/line) + 1
c
            do 900 i=1,imt
              ncol(i) = mod(i,10)
900         continue
c
            do 1000 n=1,nwr
              ia = 1 + (line*(n-1))
              ib = ia + line - 1
              if (ib .gt. imt) ib = imt
              if (wrmask) write(nunit,9021) (ncol(i),i=ia,ib)
              if ( rmask) read (nunit,9099)
              do 990 jj=1,jmt
                jjj = jmt - jj + 1
#ifdef kdtravg
c    put in code to allow reading of up to 25 region masks
c
                if (wrmask) then
                  do 1680 i=ia,ib
                    chrow(i)= chnum(mskhr(i,jjj))
1680              continue
                  write (nunit,9019) jjj, (chrow(i),i=ia,ib)
                endif
                if ( rmask) then
                  read (nunit,9019) jr, (chrow(i),i=ia,ib)
                  do 1700 i=ia,ib
                    do 1690 nch=0,25
                      if (chrow(i) .EQ. chnum(nch)) mskhr(i,jjj) = nch
1690                continue
1700              continue
#else
                if (wrmask) write(nunit,9022)jjj,(mskhr(i,jjj),i=ia,ib)
                if ( rmask) then
                  read (nunit,9022) jr , (mskhr(i,jjj),i=ia,ib)
#endif
                  if (jr .ne. jjj) then
                    write (stdout,999) nunit, jjj, jr
                    write (stderr,999) nunit, jjj, jr
                    stop ' regin  '
                  endif
                endif
990           continue
1000        continue
c
            if (wrmask) write(nunit,9031)
            if (wrmask) write(nunit,9032)
     $         (' domain for ver mask #',i,'=',vregnm(i),i=1,nvreg)
            if (wrmask) write(nunit,9034) mskvr
c
            if ( rmask) then
               read (nunit,9099)
               read (nunit,9099)
            endif
            if ( rmask) read (nunit,9033) (vregnm(i),i=1,nvreg)
            if ( rmask) read (nunit,9034) mskvr
c
          else
            if (wrmask) then
              write(nunit) mskhr, mskvr
              write(nunit) hregnm, vregnm
            endif
            if (rmask) then
              read (nunit) mskhr, mskvr
              read (nunit) hregnm, vregnm
            endif
          endif
        endif
c
c       if vertical masks were read in, set level limits for defining
c       vertical regions in term balance calculations
c
        if (rmask) then
          do 1100 i=1,nvreg
            settop = .false.
            setbot = .false.
            do 1090 k=1,km
              kk = km-k+1
              if (mskvr(k)  .eq. i .and. .not. settop) then
                llvreg(i,1) = k
                settop = .true.
              endif
              if (mskvr(kk) .eq. i .and. .not. setbot) then
                llvreg(i,2) = kk
                setbot = .true.
              endif
1090        continue
1100        continue
        endif
c
      return
c
999   format(/' error => bad j-row when reading regionmasks from unit ',
     $       i3,/'   expected',i4,'    read in',i4)
9000  format(/' Horizontal regional volumes [cubic m] and areas [sq m]')
9001  format(/'   k',' All Regions  ',9(1x,i6,5x))
9002  format(1x,i3,10(1x,1pe11.5))
9003  format(' SUM',10(1x,1pe11.5))
9004  format(/'AREA',10(1x,1pe11.5))
9011  format(/' Horizontal regional mask names & domains:')
9012  format(a22,i4,a1,a40)
9013  format(27x,a40)
9021  format('  i=>',100(i1))
9022  format(1x,i3,1x,100(i1))
9031  format(/' Vertical regional mask names & domains:')
9032  format(a22,i4,a1,a20)
9033  format(27x,a20)
9034  format (1x, 42i3)
9099  format(1x)
#ifdef kdtravg
9019  format(1x,i3,1x,100(a1))
9101  format(/'   k','     Regions:',9(1x,i6,5x))
9102  format(1x,i3,12x,9(1x,1pe11.5))
9103  format(' SUM',12x,9(1x,1pe11.5))
9104  format(/'AREA',12x,9(1x,1pe11.5))
#endif
      end
