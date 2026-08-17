      subroutine getvbc (ionum, mon, bcfld)
      implicit real*8 (a-h,o-z)
#if defined worldocn && defined seasonal
c
#include "param.h"
c
# if defined r15vbc
c=======================================================================
c
c     read in vertical boundary conditions for clinic & tracer
c     for wold ocean model runs.  uses data assembled by Dixon
c     in March, 1991
c
c=======================================================================
c
c
      dimension bcfld(imt,jmt)
c
c     read in one month at a time
c
      read(ionum) rmon, rimt, rjmt, bcfld
c
      mo = rmon
      id = rimt
      jd = rjmt
c
c     check that the month and grid dimensions are as expected
c
      if (mo .ne. mon) then 
        write(stdout,9701) mo, mon 
        stop 'getvbc' 
      endif
c
      if (id .ne. imt) then
        imtwr = imt
        write(stdout,9702) id, imtwr
        stop 'getvbc '
      endif
c
      if (jd .ne. jmt) then
        jmtwr = jmt 
        write(stdout,9703) jd, jmtwr
        stop 'getvbc'
      endif
c
c
9701  format(/' ===> ERROR: Month Mismatch <==='/
     $ ' Read in ',i4,' but expected ',i3,/)
9702  format(/' ===> ERROR: I Dimension Mismatch <==='/
     $ ' Read in ',i4,' but expected ',i3,/)
9703  format(/' ===> ERROR: J Dimension Mismatch <==='/
     $ ' Read in ',i4,' but expected ',i3,/)
c
      return
# endif
c
# if defined cpldleg2
c=======================================================================
c
c    read vertical boundary conditions for leg2 of coupled model
c
c=======================================================================

c   specify the location of taux (nwsx) in the boundary data array
      parameter (nwsx=4)

#include "ascomo.h"

#include "cvbc.h"

#include "iounit.h"

#include "seasnl.h"

#include "switch.h"

      dimension bcom(imt,jmt,2)


c   word-addressable i/o is used for the boundary data.
c   subroutine odam cannot be called here, because when the diskless
c   i/o option is used, odam can only be used for the "slab" data

      if (first) then
        nblk = nwds*nbo*12
        istats=1
        call wopen (ionum, nblk/512+1, istats)
c   get data for "month behind"
        nloca = (mbhd-1)*nwds*nbo + 1
        do 2 n=1,nbo
          call getwa (ionum, bco(1,1,n,2), nloca, nwds)
          nloca = nloca + nwds
2       continue
      endif

c   get data for "month ahead"
      nloca = (mahd-1)*nwds*nbo + 1
      do 22 n=1,nbo
        call getwa (ionum, bco(1,1,n,1), nloca, nwds)
        nloca = nloca + nwds
22    continue

      return

      entry getwsa

c   read monthly wsx into bcom(1,1,1) and wsy into bcom(1,1,2)
      r12 = c1 / 12.0
      do 280 m=1,12
        nloc = (m-1)*nwds*nbo + 1
        nloctx = nloc   + (nwsx-1)*nwds
        nlocty = nloctx + nwds
        call getwa(iobc, bcom(1,1,1), nloctx, nwds)
        call getwa(iobc, bcom(1,1,2), nlocty, nwds)
        do 270 jj=1,jmt
          do 260 i=1,imt
            wsxavg(i,jj) = wsxavg(i,jj) + bcom(i,jj,1) * r12
            wsyavg(i,jj) = wsyavg(i,jj) + bcom(i,jj,2) * r12
260       continue
270     continue
280   continue

      write(stdout,9599) wsxavg(14,14), wsyavg(14,14)
9599  format(1x,'compute annual mean wind stress:',
     $          '  wsxavg, wsyavg at (14,14)=',2e13.5)

      return
# endif


# if defined cpldleg3
c=======================================================================

c    read vertical boundary conditions for coupled model

c=======================================================================
cnetcdf
      integer iret, ncid, var1_id,var2_id,start(3), count(3)

      save ncid, var1_id,var2_id

#include "netcdf.inc"

cnetcdf

#include "cvbc.h"

#include "iounit.h"

#include "seasnl.h"

#include "switch.h"

      entry getfad

c   read surface flux adjustment data

# ifdef crayio
      if (first) then
        call ostart(iobca, nwds*2*12, nwds, 1)
c   get data for "month behind"
        nlocb = (mbhd-1)*nwds*2 + 1
        call oget(iobca, nwds, nlocb, terr(1,1,2))
        nlocb = nlocb + nwds
        call oget(iobca, nwds, nlocb, serr(1,1,2))
      endif

c   get data for "month ahead"
      nloca = (mahd-1)*nwds*2 + 1
      call oget(iobca, nwds, nloca, terr(1,1,1))
      nloca = nloca + nwds
      call oget(iobca, nwds, nloca, serr(1,1,1))

# else

      if (first) then
c   get data for "month behind"
ccc        rewind iobca
ccc        do 720 monts=1,mbhd
ccc          read(iobca) ((terr(i1,j1,2),i1=1,imt),j1=1,jmt)
ccc          read(iobca) ((serr(i2,j2,2),i2=1,imt),j2=1,jmt)
ccc720     continue
cnetcdf
cread in netcdf data
        iret = nf_open('terrserr_r30.nc',nf_read, ncid)
        call check_err(iret)
        iret = nf_inq_varid(ncid, 'terr', var1_id)
        call check_err(iret)
        iret = nf_inq_varid(ncid, 'serr', var2_id)
        call check_err(iret)


        start(1) = 1
        start(2) = 1
        start(3) = mbhd
        count(1) = imt
        count(2) = jmt
        count(3) = 1
        iret = nf_get_vara_real (ncid, var1_id, start, count, 
     &                           terr(1,1,2))
        call check_err(iret)
        iret = nf_get_vara_real (ncid, var2_id, start, count, 
     &                           serr(1,1,2))
        call check_err(iret)

cnetcdf
      endif

c   get data for "month ahead"
ccc   rewind iobca
ccc   do 730 monts=1,mahd
ccc     read(iobca) ((terr(i1,j1,1),i1=1,imt),j1=1,jmt)
ccc     read(iobca) ((serr(i2,j2,1),i2=1,imt),j2=1,jmt)
ccc730   continue
        start(1) = 1
        start(2) = 1
        start(3) = mahd
        count(1) = imt
        count(2) = jmt
        count(3) = 1
        iret = nf_get_vara_real (ncid, var1_id, start, count, 
     &                           terr(1,1,1))
        call check_err(iret)
        iret = nf_get_vara_real (ncid, var2_id, start, count, 
     &                           serr(1,1,1))
        call check_err(iret)

# endif

      return
# endif

#endif
      end
