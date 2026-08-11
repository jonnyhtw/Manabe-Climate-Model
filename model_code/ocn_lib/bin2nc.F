      subroutine bin2nc
c
c     program fgennc

#include "netcdf.inc"

* error status return
      integer  iret
* netCDF id
ccc   integer  ncid
* dimension ids
      integer  lon_dim
      integer  lat_dim
      integer  time_dim
      integer  sdsk_dim
* dimension lengths
      integer  lon_len
      integer  lat_len
      integer  time_len

      integer imt, jmt

      parameter (lon_len = 194)
      parameter (lat_len = 82)
c
      parameter (ix = 96, il=80)
c
ccc   parameter (time_len = NF_UNLIMITED)
      parameter (time_len = 1)
c

#include "ncstuf.h"

#include "param.h"

      parameter (n11=20, n12=imt*jmt*nkflds
     $,           nsl=nslab*jmt, ntb=n11+n12+ntlev*nsl)

c
* variable ids
      integer  i, j, k, ir
      integer  lon_id
      integer  lat_id
      integer  time_id
ccc   integer  sdsk_id

* rank (number of dimensions) for each variable
      integer  lon_rank
      integer  lat_rank
      integer  time_rank
      integer  ivbl_rank
      integer  idsk_rank
      parameter (lon_rank = 1)
      parameter (lat_rank = 1)
      parameter (time_rank = 1)
      parameter (ivbl_rank = 3)
      parameter (idsk_rank = 1)
* variable shapes
      integer  lon_dims(lon_rank)
      integer  lat_dims(lat_rank)
      integer  time_dims(time_rank)
      integer  ivbl_dims(ivbl_rank)
      integer  idsk_dims(idsk_rank)
      integer  start(ivbl_rank), count(ivbl_rank)
* data variables
      real*4  lon(lon_len)
      real*4  lat(lat_len)
      real*4  time(time_len), rtime
* attribute vectors
      real*8  realval(1)
* enter define mode
      iret = nf_create('ocean_r30_out.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 's_len', ntb, sdsk_dim)
      call check_err(iret)
* define variables

      idsk_dims(1) = sdsk_dim
      
      iret = nf_def_var(ncid, 'sdisk', NF_DOUBLE, idsk_rank, idsk_dims,
     &  sdsk_id) 
      call check_err(iret)


* assign attributes
      
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 40, 
     &'ocean restart                           ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title', 16, 
     &'ocean restart   ')
      call check_err(iret)
* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)


      return
      end
ccc    
ccc   subroutine check_err(iret)
ccc   integer iret
ccc   include 'netcdf.inc'
ccc   if (iret .ne. NF_NOERR) then
ccc   print *, nf_strerror(iret)
ccc   stop
ccc   endif
ccc   end
