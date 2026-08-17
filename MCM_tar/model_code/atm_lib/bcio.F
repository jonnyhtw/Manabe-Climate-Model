      subroutine bcio
c
#ifdef cpldleg2
c   This param.h is from the ocean library; it contains imt & jmt:
#include "param.h"
c
      parameter ( ix=(imt-2)/2, il=jmt-2 )
c
#else
#include "parm.h"
c
      parameter ( imt=ix*2+2, jmt=il+2 )
c
#endif
#include "ascoma.h"
c
#include "ascomo.h"
c
#if !defined cpldleg2
#include "iohist.h"
#endif
c
#ifdef cpldleg2
#include "ctmngr.h"
c
#include "iounit.h"
c
#endif
#ifdef cpldleg1
#include "qa.h"
c
#endif
cnetcdf
      common / lats / ylats(il)
c lats is a atm common block
c
      include 'netcdf.inc'

      integer iret, ncid, lon_dim, lat_dim
      integer lon_id, lat_id, ivbl1_id
      integer lon_rank, lat_rank
      parameter (lon_rank=1, lat_rank=1)
      integer lon_dims(lon_rank), lat_dims(lat_rank)
      integer ivbl_rank
      parameter (ivbl_rank=2)
      integer ivbl_dims(ivbl_rank)
      integer lon_len,lat_len, ix, il
      parameter (lon_len=ix, lat_len=il)

      real*4  lon(lon_len)
      real*4  lat(lat_len)

cnetcdf
cSGI
      real*4 bcam
cSGI
      dimension bcam(ix,il), bcom(imt,jmt)
c
c   <<<<<<<<<<<<<<<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
c
c   read boundary data at the beginning of the run
c
c   note that "iobc" and "iobcw" are defined in "iounit.h" for the
c    ocean model, but they are not defined yet for the atmos model.
c    (for the coupled model, bcio.F is in the atmos library.)
#ifdef cpldleg1
c   for leg1 of coupled model: 
c    read unit iflux is not needed (but see write unit iflux below).
c    read boundary data for atmos (to and yce) from unit isstyc.
      rewind isstyc
      nm = 1
# if defined seanal || defined seasonal
      do 9 nm=1,12
# endif
        do 8 n=1,2
          read(isstyc) bcam
          do 7 j=1,il
            do 7 i=1,ix
              bca(i,j,n,nm) = bcam(i,j)
7         continue
8       continue
        iunit = isstyc
        write(6,9005) iunit, nm
9005    format(1x,' Read atmos surface boundary data from unit ',i3,
     $            ' for m=',i3)
9     continue
# ifdef bco48
c   read accumulated monthly boundary data for ocean from unit 48
c    when restarting the model at less than 1 month intervals
      rewind 48
      read(48) qaot, pert, snot, tauxt, tauyt
# ifdef windbc
     $        ,windst, wind3t
# endif
      write(6,9048) qaot(7,68)
9048  format(/,1x,' read accumulated bc data for ocean from unit 48.',
     $            '  qaot(7,68)=',e14.5,/)
# endif
#endif
c
#ifdef cpldleg2
c   for leg2 of coupled model:
c    read boundary data for ocean (qao, pme, etc.) from unit iobc
c     in subr. setvbc
c    read unit iobcw is not needed (but see write unit iobcw below).
#endif
c
#ifdef cpldleg3
c   for coupled model:
c    read unit iflux is not needed (boundary data for the ocean
c     is computed by the atmos model each day).
c    read boundary data for atmos (to and yce) from unit isstyc.
cnetcdf
ccc      rewind isstyc
ccc      do 160 n=1,2
ccc        read(isstyc) bcam
ccc        do 150 j=1,il
ccc          do 150 i=1,ix
ccc            bca(i,j,n,1) = bcam(i,j)
ccc150     continue
ccc160   continue
cread in netcdf data
      nf_read = 0
      iret = nf_open('sstyce_r30_in.nc',nf_read, ncid)
      call check_err(iret)
      iret = nf_inq_varid(ncid, "sst", var_id)
      call check_err(iret)
      iret = nf_get_var_real (ncid, var_id, bcam)
      call check_err(iret)
      do j=1,il
       do i=1,ix
        bca(i,j,1,1) = bcam(i,j)
       end do
      end do

      iret = nf_inq_varid(ncid, "yce", var_id)
      call check_err(iret)
      iret = nf_get_var_real (ncid, var_id, bcam)
      call check_err(iret)
      do j=1,il
       do i=1,ix
        bca(i,j,2,1) = bcam(i,j)
       end do
      end do
ccc   print *,' in bcio: read sst and yce '

cnetcdf
      write(6,9020)
9020  format(1x,'Read atmos surface boundary data from unit isstyc')
#endif
      return
c
c
      entry bcout
c
c   write boundary data at the end of the run
c
#ifdef cpldleg1
c   for leg1 of coupled model: 
c    write unit iflux monthly mean boundary data for ocean
c    in subr. bcocn at end of the month.
c    write unit isstyc is not needed.
#endif
c
#ifdef cpldleg2
c   for leg2 of coupled model:
c    write unit iobc is not needed.
c    compute monthly mean boundary data for atmos
c    in subr. bcatm at end of each month.
      rewind iobcw
      nmon = 1
c
      do 26 n=1,2
        do 25 j=1,il
          do 25 i=1,ix
            bcam(i,j) = bca(i,j,n,nmon)
25      continue
c
c   set limits to prevent ieee conversion underflow & overflow
        length = ix*il
        call ieeeok(bcam, length, .false.)
c
        write(iobcw) bcam
26    continue
# ifdef gtime
      write(6,9030) iobcw, month
# else
      montho = month - 1
      if (montho .eq. 0) montho = 12
      write(6,9030) iobcw, montho
# endif
9030  format(1x,'Write atmos surface boundary data to unit',i3,
     $          ' for month=',i3)
29    continue
#endif
c
#ifdef cpldleg3
c   for coupled model:
c    write unit iflux is not needed (boundary data for the ocean
c     is computed by the atmos model each day).
c    store boundary data for atmos in subr. bcatm at end of each day.

      nmon = 1

cnetcdf
ccc      rewind isstyc
c
ccc      do 46 n=1,2
ccc        do 45 j=1,il
ccc          do 45 i=1,ix
ccc            bcam(i,j) = bca(i,j,n,nmon)
ccc45      continue
cccc
ccc        write(isstyc) bcam
ccc46    continue
* enter define mode
      iret = nf_create('sstyce_r30_out.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'longitude', ix, lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'latitude', il, lat_dim)
      call check_err(iret)
* define variables
      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid, 'longitude', NF_REAL, lon_rank, lon_dims,
     &lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid, 'latitude', NF_REAL, lat_rank, lat_dims,
     &lat_id)
      call check_err(iret)
      ivbl_dims(2) = lat_dim
      ivbl_dims(1) = lon_dim
      iret = nf_def_var(ncid, 'sst', NF_REAL, ivbl_rank, ivbl_dims,
     &ivbl1_id)
      call check_err(iret)
      iret = nf_def_var(ncid, 'yce', NF_REAL, ivbl_rank, ivbl_dims,
     &ivbl2_id)
      call check_err(iret)
* assign attributes
      iret = nf_put_att_text(ncid, lon_id, 'long_name', 9,
     &'longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lon_id, 'units', 12,
     &'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'long_name', 8,
     &'latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, lat_id, 'units', 13,
     &'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ivbl1_id, 'long_name', 27,
     &'sst restart data for atm   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ivbl1_id, 'units', 8, 'K       ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ivbl2_id, 'long_name', 27,
     &'yce thickness restart atm  ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ivbl2_id, 'units', 8, 'cm      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 40,
     &'sst and yce on atm grid                 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title', 16,
     &'sst and yce data')
      call check_err(iret)
* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)

c compute alongs
      aintvl = 360.0/float(ix)
      do 24 i=1,ix
       lon(i) = (i-1)*aintvl 
ccc    if (lon(i) .lt. 0) lon(i)=lon(i)+360.
ccc    if (lon(i) .gt. 360) lon(i)=lon(i)-360.
ccc    if (i.le.10) print *,'in bcio ', i,lon(i)
   24 continue
      do 23 j=1,il
       lat(j) = -1.0*ylats(j)
   23 continue
ccc    do j=1,jmt
ccc     if (j.le.10) print *,'in bcio ', j,lat(j)
ccc    end do

      iret = nf_put_var_real(ncid, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid, lat_id, lat)
      call check_err(iret)

c put out data
      do j=1,il
       do i = 1,ix
        bcam(i,j) = bca(i,j,1,1)
       end do
      end do
      iret = nf_put_var_real(ncid, ivbl1_id, bcam)
      call check_err(iret)

      do j=1,il
       do i = 1,ix
        bcam(i,j) = bca(i,j,2,1)
       end do
      end do
      iret = nf_put_var_real(ncid, ivbl2_id, bcam)
      call check_err(iret)

      iret = nf_close(ncid)
      call check_err(iret)
ccc   print *,' wrote sstyce to sstyce_r30_out.nc'

c
#endif
      return
      end
