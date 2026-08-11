      subroutine calc_sea_level
      implicit real*8 (a-h,o-z)
#ifdef sfctopog
c
c=======================================================================
c
c     construct the diagnostic sea level (in cm) and store in sealev
c
c  originally copied from /net/mw/kd_mom/sealev.F   (17 Oct 1996)
c  debugged with K. Dixon
c=======================================================================
c
#include "param.h"
#include "coord.h"
#include "ctmngr.h"
#include "emode.h"
#include "grdvar.h"
#include "levind.h"
#include "iounit.h"
#include "scalar.h"
#include "switch.h"
c
c
c     local variables used for sea level calculation
c
c     lpn, lps,
c     lpw, lpe   = Laplacian coefficients for relaxation
c     slx        = partial zonal gradient of sea level
c     sly        = partial meridional gradient of sea level
c     sl_sor     = successive over-relaxation factor
c     tol        = error tolerance in cm
c     maxit      = maximum number of relaxation iterations
c     sl_mask    = mask (1.0,0.0) for TS points. if all of the four 
c                  surrounding UV points are land then=0.0, else=1.0
c
c     the user may alter sl_sor, tol, and maxit by modifying the
c     fortran code
c
c     the following two arrays are found in emode.h
c
c     divslh     = divergence of sea level times depth
c     sealev     = sea level anomalies [centimeters]
c
c=======================================================================
c
cnetcdf
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
      parameter (ix=(imt-2)/2, il=jmt-2)
      parameter (lon_len=imt, lat_len=jmt)

      real*4  lon(lon_len)
      real*4  lat(lat_len)
c
      common / lats / ylats(il)
c lats is a atm common block
c
cnetcdf

      real slx(imt,jmt), sly(imt,jmt)
      real lpn(imt,jmt), lps(imt,jmt)
      real lpw(imt,jmt), lpe(imt,jmt)
      real sl_mask(imt,jmt)
      real sl_sor, tol
      real lpf, err, errmax
      real ar_sum(jmt), sl_sum(jmt), totar_sum, totsl_sum
c
c=======================================================================
c
ctemp
ctemp   write (stdout, '(/, "sealev anomaly field IN" )')
ctemp   CALL MAPIT (stdout, imt, jmt, 1.0e-01, 0.0, sealev )
ctemp
c---------------------------------------------------------------------
c       initialize sea level gradient and divergence fields
c       and relaxation coefficients
c       note: sealev is initialized to zero in ocean.F the first
c             time thru ocean.  after that, the initial guess for
c             the sealev field is the result computed the last
c             timestep that calc_sea_level was called.  this reduces
c             the number of relaxation iteration needed.
c---------------------------------------------------------------------
c
        do j=1,jmt
          do i=1,imt
            slx(i,j)    = c0
            sly(i,j)    = c0
            divslh(i,j) = c0
            lpn(i,j)= c0
            lps(i,j)= c0
            lpe(i,j)= c0
            lpw(i,j)= c0
            sl_mask(i,j) = c0
          end do 
        end do 
c
c---------------------------------------------------------------------
c     compute slx and sly from zu and zv (vertcal avgerage momentum
c     terms computed in clinic)
c---------------------------------------------------------------------
c
        f3 = c1/(c2dtsf*grav)
        do j=1,jmt-1
          do i=2,imtm1
            kz = kmu(i,j)
            if (kz .ne. 0) then
              slx(i,j) = f3*zu(i,j)
              sly(i,j) = f3*zv(i,j)
            endif
            hcheck = h(i-1,j)+h(i,j)+h(i-1,j-1)+h(i,j-1)
            if (hcheck .GT. c0) sl_mask(i,j) = c1
          end do
        end do
c
# ifdef cyclic
        do j=1,jmtm1
          slx(1,j)   = slx(imt-1,j)
          slx(imt,j) = slx(2,j)
          sly(1,j)   = sly(imt-1,j)
          sly(imt,j) = sly(2,j)
        end do
# endif
c
c---------------------------------------------------------------------
c     compute five point Laplacian relaxation coefficients,
c     compute the divergence of sea_level*depth from slx and sly
c     and normalize forcing and coefficients for relaxation.
c---------------------------------------------------------------------
c
      do j=2,jmtm1
        fx=c1
        fxa=p5*cstr(j)*cstr(j)
        fxb=p5*csu(j  )*cstr(j)*dytr(j)*dyur(j  )
        fxc=p5*csu(j-1)*cstr(j)*dytr(j)*dyur(j-1)
        fxd = p5
        do i=2,imtm1
          if (sl_mask(i,j) .EQ. c1) then
            lpn(i,j)=fxb*(h(i-1,j)+h(i,j))
            lps(i,j)=fxc*(h(i-1,j-1)+h(i,j-1))
            lpe(i,j)=fxa*dxur(i)*dxtr(i)*(h(i,j)+h(i,j-1))
            lpw(i,j)=fxa*dxur(i-1)*dxtr(i)*(h(i-1,j)+h(i-1,j-1))
            lpf     =fx/(lpn(i,j)+lps(i,j)+lpe(i,j)+lpw(i,j))
c
            duh = fxd*((h(i,j)*slx(i,j) - h(i-1,j)*slx(i-1,j))*dyu(j)
     &             +(h(i,j-1)*slx(i,j-1) - h(i-1,j-1)*slx(i-1,j-1))
     &             *dyu(j-1))
            dvh = fxd*((h(i,j)*sly(i,j)*dxu(i)
     &             + h(i-1,j)*sly(i-1,j)*dxu(i-1))*csu(j)
     &             -(h(i,j-1)*sly(i,j-1)*dxu(i)
     &             + h(i-1,j-1)*sly(i-1,j-1)*dxu(i-1))*csu(j-1))
            divslh(i,j) = (duh+dvh)*cstr(j)*dxtr(i)*dytr(j)*lpf
c
            lpn(i,j)=lpn(i,j)*lpf
            lps(i,j)=lps(i,j)*lpf
            lpe(i,j)=lpe(i,j)*lpf
            lpw(i,j)=lpw(i,j)*lpf
          endif
        end do
      end do
c
# ifdef cyclic
        do j=1,jmtm1
          lpn(1,j)   = lpn(imt-1,j)
          lpn(imt,j) = lpn(2,j)
          lps(1,j)   = lps(imt-1,j)
          lps(imt,j) = lps(2,j)
          lpe(1,j)   = lpe(imt-1,j)
          lpe(imt,j) = lpe(2,j)
          lpw(1,j)   = lpw(imt-1,j)
          lpw(imt,j) = lpw(2,j)
          divslh(1,j)   = divslh(imt-1,j)
          divslh(imt,j) = divslh(2,j)
        end do
# endif
c
c---------------------------------------------------------------------
c     relax until tol exceeds errmax or maxit iterations are performed
c----------------------------------------------------------------------
c
c     the user may alter sl_sor, tol, and maxit here
c 
      sl_sor = 1.57
      tol    = 0.01
      maxit  = 800
c
      do it=1,maxit
c
        errmax = 0.0
        do j=2,jmtm1
          do i=2,imtm1
            est = lpn(i,j)*sealev(i,j+1)+lps(i,j)*sealev(i,j-1)
     &           +lpe(i,j)*sealev(i+1,j)+lpw(i,j)*sealev(i-1,j)
     &           -divslh(i,j)
            slnew = sl_sor*est+(1.0-sl_sor)*sealev(i,j)
            err = slnew - sealev(i,j)
            sealev(i,j) = slnew
            if (err .lt. 0.0) err = -err
            if (errmax .lt. err) then
              errmax = err
              ierrmax = i
              jerrmax = j
            endif
          end do
        end do
c
# ifdef cyclic
        do j=1,jmtm1
          sealev(1,j)   = sealev(imt-1,j)
          sealev(imt,j) = sealev(2,j)
        end do
c
# endif
        if (errmax .lt. tol) go to 100
c
      end do
c
100   continue
c
c
ctemp
ctemp   write (stdout, '(/, "sealev anomaly field B " )')
ctemp   CALL MAPIT (stdout, imt, jmt, 1.0e-01, 0.0, sealev )
ctmp
c---------------------------------------------------------------------
c     normalize sealev to zero global average
c----------------------------------------------------------------------
c
      totsl_sum = 0.0
      totar_sum = 0.0
c
      do j=1,jmt
        do i=1,imt
          sl_sum(j) = 0.0
          ar_sum(j) = 0.0
        end do
      end do
c
      do j=2,jmtm1
        do i=2,imtm1
          if (sl_mask(i,j) .EQ. c1) then
            ar = dxt(i)*cst(j)*dyt(j)
            sl_sum(j) = sl_sum(j) + sealev(i,j)*ar
            ar_sum(j) = ar_sum(j) + ar
          endif
        end do
      end do
c
      do j=2,jmtm1
        totsl_sum = totsl_sum + sl_sum(j)
        totar_sum = totar_sum + ar_sum(j)
      end do
      totsl_sum = totsl_sum/totar_sum
c
      do j=2,jmtm1
        do i=2,imtm1
          if (sl_mask(i,j) .EQ. c1) then
            sealev(i,j) = sealev(i,j)-totsl_sum
          endif
        end do
      end do
c
# ifdef cyclic
        do j=1,jmtm1
          sealev(1,j)   = sealev(imt-1,j)
          sealev(imt,j) = sealev(2,j)
        end do
c
# endif
      write(stdout,'(" sealev relax: it errmax i j totsl_sum =",
     $              (i8,4x,e12.6,2i8,4x,e12.6))')
     $              it, errmax, ierrmax, jerrmax, totsl_sum
c
c     write out sea level anomaly field for use as a first guess
c     when restarting the model
c
cnetcdf
ccc   rewind (iosealev)
ccc   write(iosealev) sealev
* enter define mode
      iret = nf_create('sealevel_r30_out.nc', NF_CLOBBER, ncid)
      call check_err(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'longitude', imt, lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'latitude', jmt, lat_dim)
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
      iret = nf_def_var(ncid, 'sealev', NF_DOUBLE, ivbl_rank, ivbl_dims,
     &ivbl1_id)
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
     &'sea level restart data     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, ivbl1_id, 'units', 8, 'cm      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'source', 40,
     &'sea level restart data                  ')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'title', 16,
     &'sea level data  ')
      call check_err(iret)
* leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)

c compute alongs
      aintvl = 360.0/float(imt-2)
      do 24 i=1,imt
       lon(i) = (i-1)*aintvl - aintvl/2.
ccc    if (lon(i) .lt. 0) lon(i)=lon(i)+360.
ccc    if (lon(i) .gt. 360) lon(i)=lon(i)-360.
ccc    if (i.le.10) print *,'in sealev ', i,lon(i)
   24 continue
c
      lat(1) = -90.
      lat(jmt) = 90.
      do 23 j=1,il
       lat(j+1) = ylats(j)
   23 continue
ccc    do j=1,jmt
ccc     if (j.le.10) print *,'in sealev ', j,lat(j)
ccc    end do

      iret = nf_put_var_real(ncid, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid, lat_id, lat)
      call check_err(iret)

c put out data
      iret = nf_put_var_double(ncid, ivbl1_id, sealev)
      call check_err(iret)

      iret = nf_close(ncid)
      call check_err(iret)
      print *,' wrote sealev to sealevel_r30_out.nc'
cnetcdf
c
ctemp
ctemp   write (stdout, '(/, "sealev anomaly field A " )')
ctemp   CALL MAPIT (stdout, imt, jmt, 1.0e-01, 0.0, sealev )
ctemp
ctemp   write (stdout, '(/, "sl_mask field  " )')
ctemp   CALL MAPIT (stdout, imt, jmt, 1.0, 0.0, sl_mask )
ctmp
ctemp
ctemp
c
#endif
c
      return
      end
