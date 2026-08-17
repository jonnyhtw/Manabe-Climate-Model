      subroutine wrtatm(iunit)

cnetcdf
      use netcdf
cnetcdf
c
c This subroutine writes the atmosphere restart data from
c the designated unit passed as an argument.
c
#include "parm.h"
#include "parmjx.h"
#include "parmkp.h"
#include "parmout.h"
#include "parmdsk.h"
c
      complex vor  , div  , t  , psln
      complex vorm , divm , tm , psm
      complex vordt, divdt, tdt, psdt
      complex qmix , qmixm, qmxdt
      complex ucos , vcos , psdx, psdy, zs, ek
#ifdef rad
      complex radten
#endif
c
#include "anal1.h"
c
#include "anal2.h"
c
#include "cloudp.h"
c
#include "grdfld.h"
c
#include "iohist.h"
c
#include "nxstps.h"
c
#include "qa.h"
c
#include "radten.h"
c
#include "sfharm.h"
c
#include "surfo.h"
c
#include "water.h"
c
#include "wimplt.h"
c
cnetcdf
     
#include "netcdf.inc"

#include "ncstuff_atm.h"

#include "lats.h"

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
      real*4  klvp(klvplen)
* attribute vectors
      real*4  realval(1)

      real*8 rkflg(ix,il)

      if (iunit .ne. 72) then

cnetcdf

      write (iunit) vor
      write (iunit) div
      write (iunit) t
      write (iunit) psln
      write (iunit)  qmix
      write (iunit) vorm
      write (iunit) divm
      write (iunit) tm
      write (iunit) psm
      write (iunit)  qmixm
      write (iunit) zs
#ifdef rad
      write (iunit) radten
#endif
      write (iunit) ts
#ifdef rad
      write (iunit) slwdp
#endif
#ifdef moist
      write (iunit) precip
      write (iunit) evapor
      write (iunit) snowfl
#endif
#ifdef hydrol
      write (iunit) soilm
      write (iunit) runoff
      write (iunit) snwdpt
      write (iunit) snwmlt
#endif
#ifdef rad
      write (iunit) calb
#endif
#ifdef moist
      write (iunit) sublim
#endif
#ifdef hydrol
      write (iunit) ycemlt
#endif
#ifdef prdcld
      write (iunit) kflags
#ifdef covcld
      write (iunit) kfgcov,crx
#endif
#endif
#ifdef impcit
      write (iunit) tmean,tmnshr
#endif
#if defined cpldleg3 || defined cpldleg1
      write (iunit) qao
      write (iunit) taux
      write (iunit) tauy
      write (iunit) per
      write (iunit) sno
      write (iunit) to
      write (iunit) yce
      write (iunit) snwrnf
# ifdef windbc
      write (iunit) winds
      write (iunit) wind3
# endif
#endif
      write (iunit) samta,samsw,samlw,samto,samyce,samqao,samsno
#if defined cpldleg3 || defined cpldleg1
     *           ,sampme
#endif
     *           , toan,ycean, qoan, swan, rlwan,  taan, snoan
#if defined cpldleg3 || defined cpldleg1
     *           ,pmean
#endif
cnetcdf
      else
c  put out attributes first
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samta'  , NF_DOUBLE,
     *                                          1, samta)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samsw'  , NF_DOUBLE,
     *                                          1, samsw)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samlw'  , NF_DOUBLE,
     *                                          1, samlw)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samto'  , NF_DOUBLE,
     *                                          1, samto)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samyce' , NF_DOUBLE,
     *                                          1, samyce)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samqao' , NF_DOUBLE,
     *                                          1, samqao)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'samsno' , NF_DOUBLE,
     *                                          1, samsno)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'sampme' , NF_DOUBLE,
     *                                          1, sampme)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'toan'   , NF_DOUBLE,
     *                                          1, toan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'ycean'  , NF_DOUBLE,
     *                                          1, ycean)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'qoan'   , NF_DOUBLE,
     *                                          1, qoan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'swan'   , NF_DOUBLE,
     *                                          1, swan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'rlwan'  , NF_DOUBLE,
     *                                          1, rlwan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'taan'   , NF_DOUBLE,
     *                                          1, taan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'snoan'  , NF_DOUBLE,
     *                                          1, snoan)
      call check_err(iret)
      iret = nf_put_att_double (ncid72, NF_GLOBAL,'pmean'  , NF_DOUBLE,
     *                                          1, pmean)
      call check_err(iret)

* leave define mode
      iret = nf_enddef(ncid72)
      call check_err(iret)

c compute alongs
      aintvl = 360.0/float(ix)
      do 24 i=1,ix
       lon(i) = (i-1)*aintvl
ccc    if (i.le.10) print *, i,lon(i)
   24 continue
c
      do 23 j=1,il
       lat(j) = ylat(j)
ccc    if (j.le.10) print *, j,lat(j)
   23 continue

      do k=1,klevlen
       klev(k) = k
      end do

      do k=1,klvplen
       klvp(k) = k
      end do

      do m=1,mwavlen
       mwav(m) = m
      end do

      do n=1,nwavlen
       nwav(n) = n
      end do

      iret = nf_put_var_real(ncid72, lon_id, lon)
      call check_err(iret)
      iret = nf_put_var_real(ncid72, lat_id, lat)
      call check_err(iret)
      iret = nf_put_var_real(ncid72, klevid, klev)
      call check_err(iret)
      iret = nf_put_var_real(ncid72, klvpid, klvp)
      call check_err(iret)
      iret = nf_put_var_real(ncid72, nwavid, nwav)
      call check_err(iret)
      iret = nf_put_var_real(ncid72, mwavid, mwav)
      call check_err(iret)

c   put restart data
      iret = nf_put_var_double(ncid72, idvor , vor)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, iddiv , div)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtemp, t)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idpsln, psln)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idqmix, qmix)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idvorm, vorm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, iddivm, divm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtmpm, tm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idpsm , psm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idqmxm, qmixm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idzs  , zs)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idrdtn, radten)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idts  , ts)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idswdp, slwdp)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idprcp, precip)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idevpr, evapor)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsnfl, snowfl)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsolm, soilm)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idrunf, runoff)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsndp, snwdpt)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsnml, snwmlt)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idcalb, calb)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsubl, sublim)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idycmt, ycemlt)
      call check_err(iret)
      rkflg = kflags
      rkflg = rkflg + 0.01
      iret = nf_put_var_double(ncid72, idkflg, rkflg)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtmn , tmean)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtmsr, tmnshr)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idqao , qao)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtaux, taux)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idtauy, tauy)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idper , per)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsno , sno)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsst , to)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idyce , yce)
      call check_err(iret)
      iret = nf_put_var_double(ncid72, idsnrf, snwrnf)
      call check_err(iret)

c  close netcdf file
      iret = nf_close(ncid72)
      call check_err(iret)

      endif
cnetcdf
c
      return
      end
