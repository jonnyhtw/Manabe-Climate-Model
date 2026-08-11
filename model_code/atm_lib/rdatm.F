      subroutine rdatm(iunit)

cnetcdf
      use netcdf
cnetcdf
c

c
c This subroutine reads the atmosphere restart data from
c the designated unit passed as an argument.
c
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
      integer  ivarid

      real*8 rkflg(ix,il)

      if (iunit .ne. 71) then

cnetcdf

      read (iunit) vor
      read (iunit) div
      read (iunit) t
      read (iunit) psln
      read (iunit)  qmix
      read (iunit) vorm
      read (iunit) divm
      read (iunit) tm
      read (iunit) psm
      read (iunit)  qmixm
      read (iunit) zs
#ifdef rad
      read (iunit) radten
#endif
      read (iunit) ts
#ifdef rad
      read (iunit) slwdp
#endif
#ifdef moist
      read (iunit) precip
      read (iunit) evapor
      read (iunit) snowfl
#endif
#ifdef hydrol
      read (iunit) soilm
      read (iunit) runoff
      read (iunit) snwdpt
      read (iunit) snwmlt
#endif
#ifdef rad
      read (iunit) calb
#endif
#ifdef moist
      read (iunit) sublim
#endif
#ifdef hydrol
      read (iunit) ycemlt
#endif
#ifdef prdcld
      read (iunit) kflags
#ifdef covcld
      read (iunit) kfgcov
#endif
#endif
#ifdef impcit
      read (iunit) tmean,tmnshr
#ifdef tmean_warm
      do k=1,kx
       tmean (k) = 300.0
       tmnshr(k) = 0.0
      end do
      tmnshr(kp) = 0.0
#endif
#endif
#if defined cpldleg3 || defined cpldleg1
      read (iunit) qao
      read (iunit) taux
      read (iunit) tauy
      read (iunit) per
      read (iunit) sno
      read (iunit) to
      read (iunit) yce
      read (iunit) snwrnf
# ifdef windbc
      read (iunit) winds
      read (iunit) wind3
# endif
#endif
      read (iunit) samta,samsw,samlw,samto,samyce,samqao,samsno
#if defined cpldleg3 || defined cpldleg1
     *           ,sampme
#endif
     *           , toan,ycean, qoan, swan, rlwan,  taan, snoan
#if defined cpldleg3 || defined cpldleg1
     *           ,pmean
#endif
      else
cnetcdf

c get attributes first
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samta'  , samta)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samsw'  , samsw)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samlw'  , samlw)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samto'  , samto)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samyce' , samyce)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samqao' , samqao)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'samsno' , samsno)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sampme' , sampme)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'toan'   , toan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'ycean'  , ycean)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'qoan'   , qoan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'swan'   , swan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'rlwan'  , rlwan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'taan'   , taan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'snoan'  , snoan)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'pmean'  , pmean)
      call check_err(iret)

c now get data
      iret = nf_inq_varid(ncid71, "vor", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, vor)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "div", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, div)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "temp", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, t)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "psln", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, psln)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "qmix", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, qmix)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "vorm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, vorm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "divm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, divm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "tempm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, tm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "pslnm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, psm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "qmixm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, qmixm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "zs", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, zs)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "radten", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, radten)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "ts", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, ts)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "slwdp", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, slwdp)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "precip", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, precip)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "evapor", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, evapor)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "snowfl", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, snowfl)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "soilm", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, soilm)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "runoff", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, runoff)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "snwdpt", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, snwdpt)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "snwmlt", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, snwmlt)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "calb", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, calb)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "sublim", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, sublim)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "ycemlt", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, ycemlt)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "kflags", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, rkflg)
      call check_err(iret)
      kflags = rkflg 
      iret = nf_inq_varid(ncid71, "tmean", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, tmean)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "tmnshr", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, tmnshr)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "qao", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, qao)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "taux", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, taux)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "tauy", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, tauy)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "per", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, per)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "sno", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, sno)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "to", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, to)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "yce", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, yce)
      call check_err(iret)
      iret = nf_inq_varid(ncid71, "snwrnf", ivarid)
      call check_err(iret)
      iret = nf_get_var_double(ncid71, ivarid, snwrnf)
      call check_err(iret)

c  close netcdf file
      iret = nf_close(ncid71)
      call check_err(iret)

      endif
cnetcdf
c
      return
      end
