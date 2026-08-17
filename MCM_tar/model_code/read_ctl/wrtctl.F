      subroutine wrtctl(iunit, ifrmat)
c
c  This subroutine writes the control block
c   Definitions:
c     iunit  = the unit number for wrtting the contyrol block
c     ifrmat = formatted or unformatted data; 1 = formatted
c
#include "ctlblk.h"

cnetcdf
#include "netcdf.inc"

#include "ncstuff_atm.h"

* error status return
      integer  iret
* dimension ids
      integer  lon_dim
      integer  lat_dim
* dimension lengths
      integer  lon_len
      integer  lat_len
      integer  nwavlen
      integer  mwavlen
      integer  klevlen
      integer  klvplen
      parameter (lon_len = 96)
      parameter (lat_len = 80)
      parameter (mwavlen = 31*2)
      parameter (nwavlen = 32)
      parameter (klevlen = 18)
      parameter (klvplen = klevlen + 1)
* rank (number of dimensions) for each variable
      integer  lon_rank
      integer  lat_rank
      parameter (lon_rank = 1)
      parameter (lat_rank = 1)
      parameter (nwavrank = 1)
      parameter (mwavrank = 1)
      parameter (klevrank = 1)
      parameter (klvprank = 1)
      parameter (ir2d_rank = 2)
      parameter (ir3d_rank = 3)
      parameter (is2d_rank = 2)
      parameter (is3d_rank = 3)
* variable shapes
      integer  lon_dims(lon_rank)
      integer  lat_dims(lat_rank)
      integer  nwavdims(nwavrank)
      integer  mwavdims(mwavrank)
      integer  klevdims(klevrank)
      integer  klvpdims(klvprank)
      integer  ir2d_dims(ir2d_rank)
      integer  is2d_dims(is2d_rank)
      integer  is3d_dims(is3d_rank)
* data variables
      real*4  lon(lon_len)
      real*4  lat(lat_len)
      real*4  nwav(nwavlen)
      real*4  mwav(mwavlen)
      real*4  klev(klevlen)
* attribute vectors
      real*4  realval(1)

      integer*4 i4byte
      real*4 r4word
      real*8 r8word


cnetcdf
c
c
c
      if (ifrmat .ne. 1) then

        if (iunit .ne. 72) then
c
c  write control block to unformatted unit
      write(iunit)     expnr,           ktau,  ifcst, ifcstm,
     *                   idt,     ls,   jrow,    day, ktaust,
     *                karray,
     *                 intrs,  intpr,  intzm,  intrd,  intdk,
     *                  iday,  month,  iyear, intbal,  intsw,
     *                   nin,    lts,   igmt,    ltc,  kanal,
     *                 sprec,  sevap,  sqbar,  stbar,  sgvft,
     *                 scaht,  sgrad,    fct,    fcq,   diss,
     *                balnce,   prnt,  first,   rsrt, seasnl,
     *                  lrad, implct, resetm, resetp,  ldisk,
     *                 kdisk,    map,   sgmt,  jumps, resetd,
     *                resets, resetn,  bndry, hybrid, diurnl,
     *                 ktauo, ifcsto,   idto,  kpero,   mono,
     *                 oyear,  ioprt,   lmon,    apo, monend,
     *                   mon,   kper,     ap, kpermo,    ito,
     *                   its,    iys,    lrs,   iseg,    leg,
     *                    ih,  swamp, mlayer,   qflx, nmonth,
     *            narray
cnetcdf
        else
c  write netcdf file
* enter define mode
      iret = nf_create('atmos_r30_out.nc', NF_CLOBBER, ncid72)
      call check_err(iret)
ccc   print *, ' just after creating atmos_r30_out.nc '

* define dimensions
      iret = nf_def_dim(ncid72, 'longitude', lon_len, lon_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid72, 'latitude',  lat_len, lat_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid72, 'levels',  klevlen, klevdim)
      call check_err(iret)
      iret = nf_def_dim(ncid72, 'levelsp',  klvplen, klvpdim)
      call check_err(iret)
      iret = nf_def_dim(ncid72, 'nwaves',  nwavlen, nwavdim)
      call check_err(iret)
      iret = nf_def_dim(ncid72, 'mwaves',  mwavlen, mwavdim)
      call check_err(iret)

* define variables
      lon_dims(1) = lon_dim
      iret = nf_def_var(ncid72, 'longitude', NF_REAL,lon_rank,lon_dims,
     &lon_id)
      call check_err(iret)
      lat_dims(1) = lat_dim
      iret = nf_def_var(ncid72, 'latitude', NF_REAL, lat_rank,lat_dims,
     &lat_id)
      call check_err(iret)
      klevdims(1) = klevdim
      iret = nf_def_var(ncid72, 'levels', NF_REAL, klevrank, klevdims,
     &klevid)
      call check_err(iret)
      klvpdims(1) = klvpdim
      iret = nf_def_var(ncid72, 'levelsp', NF_REAL, klvprank, klvpdims,
     &klvpid)
      call check_err(iret)
      nwavdims(1) = nwavdim
      iret = nf_def_var(ncid72, 'nwaves', NF_REAL, nwavrank, nwavdims,
     &nwavid)
      call check_err(iret)
      mwavdims(1) = mwavdim
      iret = nf_def_var(ncid72, 'mwaves', NF_REAL, mwavrank, mwavdims,
     &mwavid)
      call check_err(iret)

      ir2d_dims(2) = lat_dim
      ir2d_dims(1) = lon_dim

      is3d_dims(3) = klevdim
      is2d_dims(2) = nwavdim
      is3d_dims(2) = nwavdim
c *2 is for complex variables
      is2d_dims(1) = mwavdim
      is3d_dims(1) = mwavdim

      iret = nf_def_var(ncid72, 'vor'   , NF_DOUBLE,is3d_rank,is3d_dims,
     &idvor)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'div'   , NF_DOUBLE,is3d_rank,is3d_dims,
     &iddiv)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'temp'  , NF_DOUBLE,is3d_rank,is3d_dims,
     &idtemp)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'psln'  , NF_DOUBLE,is2d_rank,is2d_dims,
     &idpsln)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'qmix'  , NF_DOUBLE,is3d_rank,is3d_dims,
     &idqmix)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'vorm'  , NF_DOUBLE,is3d_rank,is3d_dims,
     &idvorm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'divm'  , NF_DOUBLE,is3d_rank,is3d_dims,
     &iddivm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'tempm' , NF_DOUBLE,is3d_rank,is3d_dims,
     &idtmpm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'pslnm' , NF_DOUBLE,is2d_rank,is2d_dims,
     &idpsm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'qmixm' , NF_DOUBLE,is3d_rank,is3d_dims,
     &idqmxm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'zs'    , NF_DOUBLE,is2d_rank,is2d_dims,
     &idzs)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'radten', NF_DOUBLE,is3d_rank,is3d_dims,
     &idrdtn)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'ts'    , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idts)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'slwdp' , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idswdp)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'precip', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idprcp)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'evapor', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idevpr)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'snowfl', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsnfl)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'soilm' , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsolm)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'runoff', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idrunf)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'snwdpt', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsndp)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'snwmlt', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsnml)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'calb'  , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idcalb)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'sublim', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsubl)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'ycemlt', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idycmt)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'kflags', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idkflg)
      call check_err(iret)
      
      klevdims(1) = klevdim
      klvpdims(1) = klvpdim

      iret = nf_def_var(ncid72, 'tmean' , NF_DOUBLE, klevrank, klevdims,
     &idtmn)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'tmnshr', NF_DOUBLE, klvprank, klvpdims,
     &idtmsr)
      call check_err(iret)

      iret = nf_def_var(ncid72, 'qao'   , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idqao)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'taux'  , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idtaux)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'tauy'  , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idtauy)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'per'   , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idper)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'sno'   , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsno)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'to'    , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsst)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'yce'   , NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idyce)
      call check_err(iret)
      iret = nf_def_var(ncid72, 'snwrnf', NF_DOUBLE,ir2d_rank,ir2d_dims,
     &idsnrf)
      call check_err(iret)


* assign attributes
      iret = nf_put_att_text(ncid72, lon_id, 'long_name', 9,
     &'longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, lon_id, 'units', 12,
     &'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, lat_id, 'long_name', 8,
     &'latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, lat_id, 'units', 13,
     &'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, klevid, 'long_name', 8,
     &'klevels ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, klevid, 'units', 13,
     &'levels       ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, klvpid, 'long_name', 8,
     &'klevelsp')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, klvpid, 'units', 13,
     &'levels       ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, nwavid, 'long_name', 8,
     &'nwave   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, nwavid, 'units', 13,
     &'waves        ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, mwavid, 'long_name', 8,
     &'mwave   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, mwavid, 'units', 13,
     &'waves        ')
      call check_err(iret)

      iret = nf_put_att_text(ncid72, idvor , 'long_name', 27,
     &'vorticity spectral         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idvor , 'units', 9, 'cm**2/sec')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, iddiv , 'long_name', 27,
     &'divergence spectral        ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, iddiv , 'units', 9, 'cm**2/sec')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtemp, 'long_name', 27,
     &'temperature spectral       ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtemp, 'units', 9, 'C        ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idpsln, 'long_name', 27,
     &'sfc pressure spectral      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idpsln, 'units',11, 'dynes/cm**2')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqmix, 'long_name', 27,
     &'mixing ratio spectral      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqmix, 'units',11, 'gm/gm      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idvorm, 'long_name', 27,
     &'vorticity spectral tau-1   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idvorm, 'units', 9, 'cm**2/sec')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, iddivm, 'long_name', 27,
     &'divergence spectral tau-1  ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, iddivm, 'units', 9, 'cm**2/sec')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmpm, 'long_name', 27,
     &'temperature spectral tau-1 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmpm, 'units', 9, 'C        ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idpsm , 'long_name', 27,
     &'sfc pressure spectral tau-1')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idpsm , 'units',11, 'dynes/cm**2')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqmxm, 'long_name', 27,
     &'mixing ratio spectral tau-1')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqmxm, 'units',11, 'gm/gm      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idzs  , 'long_name', 27,
     &'zstar        spectral      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idzs  , 'units',11, 'cm         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idrdtn, 'long_name', 27,
     &'radition ten spectral      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idrdtn, 'units',11, 'C/sec      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idts  , 'long_name', 27,
     &'skin temperature           ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idts  , 'units',11, 'K          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idswdp, 'long_name', 27,
     &'SLWDP (see code for defin) ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idswdp, 'units',11, 'cal/cm**2  ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idprcp, 'long_name', 27,
     &'Precipitation              ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idprcp, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idevpr, 'long_name', 27,
     &'evaporation                ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idevpr, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnfl, 'long_name', 27,
     &'snowfall                   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnfl, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idrunf, 'long_name', 27,
     &'runoff                     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idrunf, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsndp, 'long_name', 27,
     &'snow depth                 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsndp, 'units',11, 'cm         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsolm, 'long_name', 27,
     &'soil moisture              ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsolm, 'units',11, 'cm         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnml, 'long_name', 27,
     &'snowmelt                   ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnml, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idcalb, 'long_name', 27,
     &'surface albedo             ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idcalb, 'units',11, 'percent    ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsubl, 'long_name', 27,
     &'Sublimation                ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsubl, 'units',11, 'cm/day     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idycmt, 'long_name', 27,
     &'Ycemlt (see code for defin)')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idycmt, 'units',11, 'cm/sec     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idkflg, 'long_name', 27,
     &'Packed cloud flags         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idkflg, 'units',11, 'none       ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmn , 'long_name', 27,
     &'tmean (see code for defin) ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmn , 'units',11, 'K          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmsr, 'long_name', 27,
     &'tmnshr(see code for defin) ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtmsr, 'units',11, 'K          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqao , 'long_name', 27,
     &'total heat flux into ocean ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idqao , 'units',11, 'cal/cm**2/s')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtaux, 'long_name', 27,
     &'East-west wind stress      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtaux, 'units',11, 'dyne/cm**2 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtauy, 'long_name', 27,
     &'north-sth wind stress      ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idtauy, 'units',11, 'dyne/cm**2 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idper , 'long_name', 27,
     &'liq. water flux into ocean ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idper , 'units',11, 'cm/sec     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsno , 'long_name', 27,
     &'snowfall on ocean          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsno , 'units',11, 'cm/sec     ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsst , 'long_name', 27,
     &'SST on atm grid            ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsst , 'units',11, 'K          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idyce , 'long_name', 27,
     &'sea ice thickness          ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idyce , 'units',11, 'cm         ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnrf, 'long_name', 27,
     &'frz. water flux into ocean ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, idsnrf, 'units',11, 'cm/sec     ')
      call check_err(iret)

ccc   realval(1) = -999
ccc   iret = nf_put_att_real(ncid72, ivbl_id, '_FillValue', NF_REAL, 1,
ccc  &realval)
ccc   call check_err(iret)
      iret = nf_put_att_text(ncid72, NF_GLOBAL, 'source', 40,
     &'Atmosphere restart file                 ')
      call check_err(iret)
      iret = nf_put_att_text(ncid72, NF_GLOBAL, 'title', 16,
     &'atmos_r30.nc    ')
      call check_err(iret)

c put out ctlblk as attributes
      iret = nf_put_att_text (ncid72, NF_GLOBAL,'expnam1', 8, expnr(1))
      call check_err(iret)
      iret = nf_put_att_text (ncid72, NF_GLOBAL,'expnam2', 8, expnr(2))
      call check_err(iret)
c     no integer*8 - convert int*8 to real and store
      eps = 0.01
      r4word = ktau
      r4word = r4word + eps
      iret = nf_put_att_real (ncid72, NF_GLOBAL,'rktau'   , NF_REAL, 1, 
     &                                           r4word)
      call check_err(iret)
      r8word = ifcst
      r8word = r8word + eps
      iret = nf_put_att_double(ncid72, NF_GLOBAL,'rifcst'  ,NF_DOUBLE,1, 
     &                                           r8word)
      call check_err(iret)
      r8word = ifcstm
      r8word = r8word + eps
      iret = nf_put_att_double(ncid72, NF_GLOBAL,'rifcstm' ,NF_DOUBLE,1, 
     &                                           r8word)
      call check_err(iret)
      i4byte = idt
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'idt'     , NF_INT , 1, 
     &                                           i4byte)
      call check_err(iret)
      i4byte = ls
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ls'      , NF_INT , 1, 
     &                                          i4byte)
      call check_err(iret)
      i4byte = jrow
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'jrow'    , NF_INT , 1, 
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'day'    , NF_DOUBLE,
     &                                          1, day)
      call check_err(iret)
      r4word = ktaust
      r4word = r4word + eps
      iret = nf_put_att_real (ncid72, NF_GLOBAL,'rktaust' , NF_REAL, 
     &                                        1, r4word)
      call check_err(iret)
      i4byte = karray(1)
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'karray1',  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = karray(2)
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'karray2',  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = karray(3)
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'karray3',  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = karray(4)
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'karray4',  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = karray(5)
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'karray5',  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intrs
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intrs'  ,  NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intpr
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intpr'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intzm
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intzm'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intrd 
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intrd'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intdk
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intdk'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = iday
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'iday'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_text (ncid72, NF_GLOBAL,'month'  , 8, month)
      call check_err(iret)
      i4byte = iyear
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'iyear'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intbal
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intbal' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = intsw
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'intsw'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = nin
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'nin'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = lts
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'lts'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = igmt
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'igmt'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = ltc
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ltc'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = kanal
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'kanal'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sprec'  , NF_DOUBLE,
     &                                          1, sprec)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sevap'  , NF_DOUBLE,
     &                                          1, sevap)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sqbar'  , NF_DOUBLE,
     &                                          1, sqbar)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'stbar'  , NF_DOUBLE,
     &                                          1, stbar)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sgvft'  , NF_DOUBLE,
     &                                          1, sgvft)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'scaht'  , NF_DOUBLE,
     &                                          1, scaht)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sgrad'  , NF_DOUBLE,
     &                                          1, sgrad)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'fct'    , NF_DOUBLE,
     &                                          1, fct)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'fcq'    , NF_DOUBLE,
     &                                          1, fcq)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'diss'   , NF_DOUBLE,
     &                                          1, diss)
      call check_err(iret)
      i4byte = 0
      if (balnce) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'balnce' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (prnt) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'prnt'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (first) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'first'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (rsrt) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'rsrt'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (seasnl) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'seasnl' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (lrad) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'lrad'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (implct) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'implct' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (resetm) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'resetm' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (resetp) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'resetp' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (ldisk) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ldisk'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (kdisk) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'kdisk'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (map) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'map'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (sgmt) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'sgmt'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (jumps) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'jumps'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (resetd) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'resetd' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (resets) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'resets' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (resetn) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'resetn' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (bndry) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'bndry'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (hybrid) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'hybrid' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (diurnl) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'diurnl' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      r4word = ktauo
      r4word = r4word + eps
      iret = nf_put_att_real (ncid72, NF_GLOBAL,'ktauo'  , NF_REAL,
     &                                        1, r4word)
      call check_err(iret)
      r8word = ifcsto
      r8word = r8word + eps
      iret = nf_put_att_double(ncid72, NF_GLOBAL,'ifcsto' , NF_DOUBLE,
     &                                        1, r8word)
      call check_err(iret)
      i4byte = idto
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'idto'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = kpero
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'kpero'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = mono
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'mono'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      r4word = oyear
      iret = nf_put_att_real (ncid72, NF_GLOBAL,'oyear'  , NF_REAL,
     &                                        1, r4word)
      call check_err(iret)
      i4byte = ioprt
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ioprt'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (lmon) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'lmon'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'apo'    , NF_double,
     &                                          1, apo)
      call check_err(iret)
      i4byte = monend
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'monend' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = mon
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'mon'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = kper
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'kper'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'ap'     , NF_double,
     &                                          1, ap)
      call check_err(iret)
      i4byte = kpermo
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'kpermo' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = ito
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ito'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = its
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'its'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = iys
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'iys'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (lrs) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'lrs'    , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = iseg
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'iseg'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      iret = nf_put_att_text (ncid72, NF_GLOBAL,'leg'    , 8, leg)
      call check_err(iret)
      i4byte = ih
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'ih'     , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (swamp) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'swamp'  , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (mlayer) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'mlayer' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = 0
      if (qflx) i4byte = 1
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'qflx'   , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)
      i4byte = nmonth
      iret = nf_put_att_int  (ncid72, NF_GLOBAL,'nmonth' , NF_INT , 1,
     &                                           i4byte)
      call check_err(iret)


      endif
cnetcdf


c write output some info from control block
      write(6,400) ktau, iday, month, iyear, iunit
c
  400 format(' At ktau = ',i10,' The day month year is ',i3,a8,i10,
     * ' The control block is written on unit ',i3,'.')
c
      else
c
c  write control block to formatted unit
      write(iunit,410) expnr,           ktau,  ifcst, ifcstm,
     *                   idt,     ls,   jrow,    day, ktaust,
     *                karray,
     *                 intrs,  intpr,  intzm,  intrd,  intdk,
     *                  iday,  month,  iyear, intbal,  intsw,
     *                   nin,    lts,   igmt,    ltc,  kanal,
     *                 sprec,  sevap,  sqbar,  stbar,  sgvft,
     *                 scaht,  sgrad,    fct,    fcq,   diss,
     *                balnce,   prnt,  first,   rsrt, seasnl,
     *                  lrad, implct, resetm, resetp,  ldisk,
     *                 kdisk,    map,   sgmt,  jumps, resetd,
     *                resets, resetn,  bndry, hybrid, diurnl,
     *                 ktauo, ifcsto,   idto,  kpero,   mono,
     *                 oyear,  ioprt,   lmon,    apo, monend,
     *                   mon,   kper,     ap, kpermo,    ito,
     *                   its,    iys,    lrs,   iseg,    leg,
     *                    ih,  swamp, mlayer,   qflx, nmonth,
     *            narray

      endif
c
c
  410 format(' The control block: ',
     *       /,1x,2a10,6i10,e10.4,i10,/,1x,10i10,/,1x,i10,a10,8i10,
     *       /,1x,10e10.4,
     *       /,1x,10l10,/,1x,10l10,/,1x,5i10,e10.4,i10,l10,e10.4,i10,
     *       /,1x,2i10,e10.4,i10,3i10,l10,i10,a10,/,1x,i10,3l10,6i10,
     *       /,1x,10i10)
c
      return
      end
