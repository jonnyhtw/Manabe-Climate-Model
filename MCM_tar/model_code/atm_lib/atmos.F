      subroutine atmos
c
c   subroutines called: symdump,timing,fftint(fftkkg),cpoly(parmtr),
c                       indata,parmtr,march,uvphrm,getime,passa,energy,
c                       linear,tderiv(linear),passb,robtq(linear),zmint,
c                       spint,rstart(indata),lgndre1(parmtr)
c                       (also swapan,zeroac(swapan),gmtan(swapan) for
c                        diurnal models)
c
c   atmos is the main program in the supersource spectral model.
c   see documentsupersrc (rjs) for descriptions of: the parameters in
c   the model, variables in the common blocks, listing of online
c   integrals, flow chart of the model, and descriptions of the
c   subroutines.
c
c
#include "parm.h"
#include "parmjx.h"
#include "parmkp.h"
#include "parms.h"
#include "parmncc.h"
#include "parmout.h"
#include "parmdsk.h"
#ifdef fixcld
      parameter ( lx=nc )
#endif
#ifdef prdcld
      parameter ( lx=kx+2 )
#endif

      parameter ( k2=kx+2 ) 
      parameter ( lpmx=isc*mtrun+jtrun+1 )
      parameter ( ktaua=0 )
c
      complex vor   , div  , t    , psln, ucos , vcos
      complex vorm  , divm , tm   , psm , ucosm, vcosm
      complex vordt , divdt, tdt  , psdt
      complex qmix  , qmixm, qmxdt
      complex ps    , psdx , psdy ,  zs, ek
      complex  aift , bift , ekift, zaift, zbift
      complex  tift , tuift, tvift, psift, spift
      complex  qeift, quift, qvift
#ifdef rad
      complex  radten
#endif
#ifdef dirnal
      complex  qdot , somg  , qdift, ogift
#endif
      complex  rzs
cSGI
      integer*8 isecday, nstep8
      integer*8 nday
cSGI
cnetcdf
      integer iret, ncid, var_id,start(2), count(2)

#include "netcdf.inc"

cnetcdf
c
c
#ifdef dirnal
      dimension expno(2)
#endif
      dimension anala(awrds),analb(bwrds)
c
#include "absrb.h"
c
#include "albedo.h"
c
#include "albsea.h"
c
#include "anal2.h"
c
#include "asblck.h"
c
#if defined cpldleg3 || defined cpldleg1
#include "ascoma.h"
#endif
      common / atmtmp /  rsrt2, cnvly, gmrad, strad, tbar, qbar
      common / atmtmp / nsteps, ktaub, ktauc
#ifdef dirnal
      common / atmtmp /   irts, ngmt
#endif
cSGI
      common / atmtmp2/   isecday
cSGI
c
#include "auxatm.h"
c
#include "cabd.h"
c
#include "cexpnam.h"
c
#include "cld.h"
c
c
#include "cld1.h"
c
#include "cloudp.h"
c
#include "co2lw.h"
c
#include "ctlblk.h"
c
#include "dispat.h"
c
#include "dxordy.h"
c
#include "establ.h"
c
#include "gauss.h"
c
#include "grdfld.h"
c
#include "grid.h"
c
#include "gridw.h"
c
#include "guvm.h"
c
#include "h2olw.h"
c
#include "ift.h"
c
#include "iohist.h"
c
cSGI
cinclude "matrix.h"
cSGI
c
#include "intgrl.h"
c
#include "latdat.h"
c
#include "maimc.h"
c
#include "mrstr.h"
c
#include "nxstps.h"
c
#include "onlint.h"
c
#include "param.h"
c
#include "planck.h"
c
#include "polyls.h"
c
#include "pstar.h"
c
#include "qa.h"
c
#include "radten.h"
c
#include "rdintl.h"
c
#include "sfharm.h"
c
#include "sftuvm.h"
c
#include "sgwgt1.h"
c
#include "sqdot.h"
c
#include "sstcof.h"
c
#include "surfc.h"
c
#include "surfo.h"
c
#include "sw.h"
c
#include "swnuin.h"
c
#include "swnuot.h"
c
#include "swtabl0.h"
c
#include "timeck.h"
c
#include "trans.h"
c
#include "tsmoth.h"
c
#ifdef ice_dynamics
c
c variables in ice model used in atmospheric model (cfi)
c
#include "varice.h"
#endif
c
      dimension pk(ix,il)
c
#include "water.h"
c
#include "wimplt.h"
c
#include "zrise.h"
c
#include "zsrvar.h"
c
      equivalence ( anala(1),   vfu(1,1))
#ifdef moist
      equivalence ( analb(1), convr(1,1))
#endif
#ifdef dry
      equivalence ( analb(1), convt(1,1))
#endif
c
      namelist / indat1 / irt, ndays
      namelist / indat2 / alpha, dt, grow, rk0
      namelist / indat3 / ktaust, intpr, intrs, intzm,
     *                    intrd , intdk, intsr, ioprt
c
c  definitions for the namelist variables:
c     irt     = initial timestep number
c     ndays   = number of days to be integrated.
c     alpha   = robert smoothing constant used in tderiv and robtq.
c               alpha is stored in common /tsmoth/
c     dt      = timestep length in seconds.  its value is mult by 5.0
c               in parmtr for implicit models.
c               mod(86400,5*ifix(dt)) must equal zero. the following 64
c               numbers are the only numbers which satisfy that
c               condition :
c                 17280, 8640, 5760, 4320, 3456, 2880, 2160, 1920, 1728,
c                  1440, 1152, 1080,  960,  864,  720,  640,  576,  540,
c                   480,  432,  384,  360,  320,  288,  270,  240,  216,
c                   192,  180,  160,  144,  135,  128,  120,  108,   96,
c                    90,   80,   72,   64,   60,   54,   48,   45,   40,
c                    36,   32,   30,   27,   24,   20,   18,   16,   15,
c                    12,   10,    9,    8,    6,    5,    4,    3,    2,
c                     1.
c               the values that have been used successfully are in the
c               range of 180 to 540. 540 was used in a t30 version of
c               the model. if the selected value of "dt" does not work,
c               the variable "alpha" may need to be changed. if the
c               value of "dt" is changed, it is highly recommended that
c               the value of ktau be changed so that the day number of
c               the model can be calculated by dividing ktau by the
c               number of time steps in a model day.
c     grow      the grow mountains option.  if grow is true mountains
c               are grown linearly until day 20 of the experiment.
c               (see day in the ctlblk.)  after that grow has no effect.
c               if grow is false, the input mountains are used as is.
c     rk0       rk0 is the tunable constant for the gravity wave drag
c               parmeterization.  if rk0 equals 0.0, most of the
c               gravity wave drag code is skipped.  if rk0 is non-
c               zero, the surface elevation variances are read in
c               from unit 85.
c
c     ktaust  = timestep to stop run.
c     intpr   = mod with ktau;interval to calculate and print integrals
c               from enregy. (in timesteps)
c     intrs   = mod with ktau; interval to call rstart, write restart
c               tape in indata. (in timesteps)
c     intzm   = mod with ktau; interval to calculate and print zonal
c               means in zmint. (in timesteps)
c     intrd   = mod with ifcst; interval to do radiation in impt.
c               (in seconds)
c     intdk   = mod with ifsct; interval to do full long-wave radiation
c               update. (in days, later mult by 86400)
c     intsr   = mod with ktau; interval to do short-wave rad update.
c               (in timesteps) 't' (diurnal) models only.
c     ioprt   = mod with ktau; interval to print stuff from tsc.
c               (in timesteps) air-sea models only.
c
cSGI 
c  default namelist variables
c     data irt     /     0 /
c     data ndays   /    30 /
cifdef r15
c     data dt      / 360.0 /
cendif
cifdef r21
c     data dt      / 240.0 /
cendif
cifdef r30
c     data dt      / 180.0 /
cendif
cifdef t21
c     data dt      / 360.0 /
cendif
cifdef t42
c     data dt      / 180.0 /
cendif
cifdef r15
c     data alpha   / 0.020 /
cendif
cifdef r21
c     data alpha   / 0.030 /
cendif
cifdef r30
c     data alpha   / 0.040 /
cendif
cifdef t21
c     data alpha   / 0.020 /
cendif
cifdef t42
c     data alpha   / 0.040 /
cendif
c     data grow    / .true.   /
c     data rk0     / 0.0      /
c     data intsr   /        9 /
cSGI
c
#ifdef dirnal
      data expno /8hsundials,8h        /
      data lg / 51 /
      data lh / 52 /
#endif
c
c   <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <> <>
c
c
      if(.not.la) go to 1
cSGI
c  default namelist variables
      irt     =     0
      ndays   =    30
#ifdef r15
      dt      = 360.0
#endif
#ifdef r21
      dt      = 240.0
#endif
#ifdef r30
      dt      = 180.0
#endif
#ifdef t21
      dt      = 360.0
#endif
#ifdef t42
      dt      = 180.0
#endif
#ifdef r15
      alpha   = 0.020
#endif
#ifdef r21
      alpha   = 0.030
#endif
#ifdef r30
      alpha   = 0.040
#endif
#ifdef t21
      alpha   = 0.020
#endif
#ifdef t42
      alpha   = 0.040
#endif
      grow    = .true.
      rk0     = 0.0
      intsr   =        9
      ioprt   = 99999999

      isecday = 86400
cSGI

c
c  if not first time in atmos this run, go to top of timestep loop.
c
c  call initialization programs
c
c     call symdump
c
      call timing
c
#ifdef internal_fft
      call fftrig
#endif
#ifdef cray_fft
      call fftint(ix)
#endif
c
#ifdef triang
      call initrc
#endif
c
      call cpoly
c
      call lgndre1
c
      call prescl
c
#ifdef dirnal
      ltmm =ltm
      ngmt =0
      adski=adsk
      bdski=bdsk
      lgi  =lg
      lhi  =lh
      adsko=adsk+2
      bdsko=bdsk+2
      lgo  =lg+2
      lho  =lh+2
c
      rewind 79
      write (79) ngmt,ngmt,ltmm,ngmt
      rewind 79
#endif
#if defined cpldleg3 || defined cpldleg1
c
c   read array for surface water (or ice) runoff into the ocean
c
ccc   read (irunof) ijrnf
cread in netcdf data
      nf_read = 0
      iret = nf_open('ijrnf_r30.nc',nf_read, ncid)
      call check_err(iret)
      start(1) = 1
      start(2) = 1
      count(1) = ix
      count(2) = il
      iret = nf_inq_varid(ncid, "ijrnf", var_id)
      call check_err(iret)
      iret = nf_get_vara_int (ncid, var_id, start, count, ijrnf)
      call check_err(iret)
      iret = nf_close(ncid)
      call check_err(iret)
cnetcdf

#endif
c
c   read annual mean bare-surface albedo array
c
ccc   read (ialb) alb
c read in netcdf data
      iret = nf_open('alb_kukla_r30.nc',nf_read, ncid)
      call check_err(iret)
      start(1) = 1
      start(2) = 1
      count(1) = ix
      count(2) = il
      iret = nf_inq_varid(ncid, "alb", var_id)
      call check_err(iret)
      iret = nf_get_vara_real (ncid, var_id, start, count, alb)
      call check_err(iret)
      iret = nf_close(ncid)
      call check_err(iret)
cnetcdf

c
c
c  read starting timestep (irt)  and number of days
c    to be run (ndays).
      read (5,indat1)
      write(6,indat1)
c  read robert constant (alpha) and timestep length (dt)
      read (5,indat2)
      write(6,indat2)
c
c  read in variance of obs heights
c
      if (rk0 .ne. 0.0) then
ccc     read (izsvar)  zsvar
c read in netcdf data
        nf_read = 0
        iret = nf_open('zstarvar_r30.nc',nf_read, ncid)
        call check_err(iret)
        start(1) = 1
        start(2) = 1
        count(1) = ix
        count(2) = il
        iret = nf_inq_varid(ncid, "zsvar", var_id)
        call check_err(iret)
        iret = nf_get_vara_real (ncid, var_id, start, count, zsvar)
        call check_err(iret)
 
      endif
c
c  read in data from restart tape
      call indata
c
#ifdef ice_dynamics
c
c read in ice restart (cfi)
c
      call indati
#endif
c
c  store the experiment name for the control block
      expnr(1) = expn1
      expnr(2) = expn2
c
#ifdef dirnal
      ltc =1
      irts=irt
      call zeroac
c
      rewind 79
      read (79) ngmt
      rewind 79
#endif
#ifdef impcit
      implct=.true.
#endif
c
      call parmtr
c
#ifdef fixcld
      call cldht
#endif
c
#ifdef consrv_ps
c
c  In conserving ln(ps) ifdef, call intcnp to
c  save off initial values of psln, psm, dtx2
c
      call intcnp
c
#endif
c
#if !defined cpldleg3 && !defined cpldleg1
c
c q-flux and fixed-sst models in original model (cfi)
c
#ifdef old_qflux_fixedsst
      if (.not.swamp) then
#ifdef mixed_layer
          if (qflx) call qcoef
#else
          call coef
#endif
      endif
#endif
c
#endif
c
      tprec = sprec
      tevap = sevap
      gvft  = sgvft
      caht  = scaht
      gmrad = sgrad
c
#ifdef dirnal
      nin=0
      nou=2
#endif
c
      first  = .true.
      idt    = dt
      nsteps = 86400/idt
cSGI
      nstep8 = nsteps
cSGI
      rsrt2  = 1./sqrt(2.)
      resetd = .true.
#ifdef dirnal
      igmt   = nsteps/ltm
#endif
c
#ifdef changdt
c
c   The following can be used to reset the timestep if desired
c
ccc   if (ktau .eq. ???????) then
ccc     ifcold = ifcst
ccc     ktold  = ktau
ccc     iyrold = iyear
c   compute new ifcst and ktau:
ccc     ifcst = ???????
ccc     ktau  = ifcst/idt
ccc     iyear = ???????
c   note that resetm needs irt
ccc     irt   = ktau
ccc     irtsav= ktau
# ifdef dirnal
ccc     irts=irt
# endif
ccc     write(6,9006) ifcold, ifcst, ktold, ktau, iyrold, iyear
ccc   endif
#endif
c
c   read ending timestep, print, restart, zonalmean, radiation intervals
c
c  below is a cluge to determine if intdk were input
c  from namelist
      intdks = intdk
c
      read (5,indat3)
      write(6,indat3)
c
c  reset ktaust iff ndays nonzero
      if (ndays .ne. 0) ktaust = irt + ndays*nsteps
c
c  reset intdk iff intdk input from namelist
      if (intdks .ne. intdk) intdk = intdk * intrd
#ifdef dirnal
      intsw = 86400 * intsr
#endif
      ktaub = ktaua + intrs
#ifndef no_reset_tmean
      ktauc = ktaua + 100 * nsteps
#else
      ktauc = 0
#endif
      write (6,140) ktaust,intpr,intrs,intzm,intrd,intdk,intsr
      write (6,141) ioprt
 141  format('   ioprt =',i10,/)
c
c
c
c   ************************ begin timestep loop *******************
c
c
   1  continue
c
ccccc day=float(ifcst) * (1.0/86400.0)
c   use nday to compute day
      nday   = ifcst / isecday
      rem    = ifcst - nday * isecday
      day    = float(nday) + rem / 86400.0
#ifdef seanal
      isecsr = mod(ifcst,isecday)
      ihour  = isecsr/3600
      itemp  = isecsr-3600*ihour
      imin   = itemp/60
      isec   = itemp-60*imin
#endif
      prnt   = mod(ktau,intpr).eq.0
      tau    = ktau
c
c  set logical control variables
      resetm = .false.
      rsrt   = mod(ktau+1   ,intrs ).eq.0
      resetp = mod(ifcst    ,intbal).eq.0
      resetn = mod(ifcst    ,intbal).eq.0
      balnce = mod(ifcst+idt,intbal).eq.0
#ifdef rad
      lrad   = mod(ifcst+idt,intrd ).eq.0
      ldisk  = mod(ifcst+idt,intdk ).eq.0
      if (ktau   .eq.      0) lrad   = .true.
#endif
#ifdef impcit
      if (ktau   .eq.    irt) resetm = .true.
#endif
#ifdef rad
      if (ktau   .eq.  ktaua) lrad   = .true.
      if (resetd .and.  lrad) ldisk  = .true.
      if (lrad )              resetd = .false.
#endif
#ifdef dirnal
      resets =           (ktau-irts).lt.nsteps
      kdisk  = mod(ifcst+idt, intsw).le.86400
      sgmt   = mod(   ktau+1,  igmt).eq.0
      ltc    = mod(   ktau+1,nsteps)/igmt
      if (resets.and.lrad) kdisk = .true.
      if (ltc   .eq.    0) ltc   = ltm
#endif
#ifdef cpldleg3
      if (  prnt) write(6,9107)
#endif
      if (  prnt) call timing
#ifdef seanal
      if   (lrad) call march
#endif
      if (balnce) prnt=.true.
c
c   set surface boundary data from ocean
      if(kper.eq.0) then
c   zero boundary data
#if defined cpldleg3 || defined cpldleg1
        write(6,9181) ktau
#endif
        do 181 j=1,il
         do 181 i=1,ix
          qao(i,j) = 0.0
          sno(i,j) = 0.0
#ifdef ice_dynamics
c
c restart net heat fluxes over ice and open water (cfi)
c
	  qfice(i,j) = 0.0
          qfow(i,j) = 0.0
c
c restart wind stress vectors (cfi)
c
          tauxic(i,j) = 0.
          tauyic(i,j) = 0.
          windvx(i,j) = 0.
          windvy(i,j) = 0.
c
#endif
#if defined cpldleg3 || defined cpldleg1
           per(i,j) = 0.0
           taux(i,j) = 0.0
           tauy(i,j) = 0.0
# ifdef windbc
          winds(i,j) = 0.0
          wind3(i,j) = 0.0
# endif
#endif
  181   continue
      endif
#ifdef cpldleg1
      if(la .and. kpermo .ne. 0) then
        write(6,9085) ktau, kpermo
        stop
      endif
      if(kpermo .eq. 0) then
c   zero arrays for storing monthly boundary data
        do 185 j=1,il
          do 185 i=1,ix
            qaot(i,j) = 0.0
            pert(i,j) = 0.0
            snot(i,j) = 0.0
           tauxt(i,j) = 0.0
           tauyt(i,j) = 0.0
# ifdef windbc
          windst(i,j) = 0.0
          wind3t(i,j) = 0.0
# endif
  185   continue
        write(6,9090) ktau
      endif
#endif
c
      if (mod(ifcst,isecday).eq.0 .or. la .or. kper.eq.0) then
c
#if defined cpldleg3 || defined cpldleg1
c
c   Load surface boundary data for atmos grid  
c
# ifdef annlmn
        do 4020 j=1,il
         do 4020 i=1,ix
          to (i,j) = bca(i,j,1,1)
          yce(i,j) = bca(i,j,2,1)
4020    continue
# endif
c
# ifdef seanal
# ifdef cpldleg1
c  compute months needed for model with seasonal forcing.
c  this assumes boundary conditions are monthly averages that are
c  treated as valid for the midpoint of the month
c
c   *****    ***  definitions  ***           ******
c   *  wmahd = weight for the month ahead         *
c   *  wmbhd = weight for the month behind        *
c   *  mbhd  = month behind                       *
c   *  mahd  = month ahead                        *
c   *****                                    ******
c
        secmon = (365.0 * 86400.0) / 12.0
        mnthfl = (ifcst + 0.5*secmon) / secmon
        wmahd  = (ifcst + 0.5*secmon - mnthfl * secmon) / secmon
        wmbhd  = 1.0 - wmahd
        mbhd   = mod(mnthfl+11,12) + 1
        mahd   = mod(mnthfl,12) + 1
# endif
c
c  load to and yce with new boundary data
c
        do 4021 j=1,il
         do 4021 i=1,ix
# ifdef cpldleg3
           to (i,j) = bca(i,j,1,1)
           yce(i,j) = bca(i,j,2,1)
# endif
# ifdef cpldleg1
           to (i,j) = bca(i,j,1,mahd)*wmahd + bca(i,j,1,mbhd)*wmbhd
           yce(i,j) = bca(i,j,2,mahd)*wmahd + bca(i,j,2,mbhd)*wmbhd
# endif
 4021   continue
# endif
#endif
c
c  zero global integral vbls
        tto = 0.0
        tyc = 0.0
        twt = 0.0
c
c   load ts with to where there is no ice or land
c   also compute global mean ocn temp and sea ice integrals.
        istpls = 0
        do 405 j1=1,il
         do 182 i=1,ix
          if (ts(i,j1) .gt. 0.0) then
            zin(i) = 1.0
          else
            zin(i) = 0.0
          endif
  182    continue
c
c   convert water equivalent ice to ice in passa and passb
         do 183 i=1,ix
          if (yce(i,j1) .gt. 0.0) then
            pk(i,j1) = 1.0
          else
            pk(i,j1) = 0.0
          endif
  183    continue
c
        ttoj = 0.0
        tycj = 0.0
        seaj = 0.0
        ls   = min0(j1,2*iy+1-j1)
        wt   = 0.5*weight(ls)*cosqr(ls)*float(3-nhem)
c
#ifdef cpldleg3
c   check land sea mask between atm and ocn models
        if (.not. la) then
          do 4030 i=1,ix
           if (zin(i) .eq. smask(i,j1)) then
             istpls = 1
             write(6,407)
             write(6,408) i,j1,zin(i),smask(i,j1)
           endif
 4030     continue
        endif
#endif
c
        do 403 i=1,ix
#ifdef cpldleg3
c   Convert ocean surface temperature to degrees Kelvin
          to(i,j1) = to(i,j1) + tfra
#endif
         ts(i,j1) = abs(ts(i,j1))
         pk(i,j1) = max(pk(i,j1),zin(i))
         ts(i,j1) = to(i,j1)*(1.0-pk(i,j1))
     *            + ts(i,j1)*     pk(i,j1)
c   now pk=1 implies land or ice  and pk=0 sea with no ice
  403   continue
c
c   set ts negative over sea
        do 404 i=1,ix
         ts(i,j1) = ts(i,j1)*(zin(i)-1.0)+ts(i,j1)*zin(i)
  404   continue
c
c  zonally sum ocn temp and sea ice
        do 406 i=1,ix
         ttoj = ttoj +  to(i,j1)*(1.0-zin(i))
         tycj = tycj + yce(i,j1)*(1.0-zin(i))
         seaj = seaj +           (1.0-zin(i))
  406   continue
c
c  compute global mean to and yce
        tto = tto + ttoj*rix*wt
        tyc = tyc + tycj*rix*wt
        twt = twt + seaj*rix*wt
  405   continue
c
#ifdef cpldleg3
c
        if (istpls .eq. 1 ) call aquit (407)
c
        tto = tto/twt - tfra
#endif
#ifndef cpldleg3
        tto = tto/twt
#endif
        tyc = tyc/twt
c
        write(6,171) ktau, tto, tyc
#if defined seanal && defined cpldleg1
        write(6,173) mahd, mbhd, wmahd, wmbhd
#endif
      endif
c
c  store ocn temp and sea ice when lrad is true.
c  note: tto and tyc are computed on the first timestep of the day,
c        lrad is usually true on the last timestep of the day.
      if (lrad) then
        samto  = samto  + tto
        samyce = samyce + tyc
      endif
c
      la=.false.
c
c   reset tmean and tmnshr during the first ktauc timesteps
#ifdef impcit
      if (ktau .le. ktauc) then
        if (mod(ktau,nstep8).eq.0) then
          resetm=.true.
          do 75 k=1,kx
           tmean(k)=rsrt2*real(t(1,1,k))+tfra
 75       continue
          do 73 k=2,kx
           tmnshr(k)=tmean(k)-tmean(k-1)
 73       continue
          tmnshr(1 )=0.0
          tmnshr(kp)=0.0
          write(6,711) ktau,tmean,tmnshr
        endif
      endif
#endif
c
c   change topography gradually at the start of a new experiment
c   also see grmult,impvar and passa to raise topography
      if (grow) then
        rr   = dt/1728000.0
        stts = 0.0
        zrat = amin1(rr*(tau-stts),1.0)
        do 97 n=1,jx
         do 97 m=1,mx
          rzs(m,n,1)=zrat*zs(m,n,1)
 97     continue
        if (prnt .and. rr*(tau-stts).le.1.0)
     *      write(6,9997) day, tau, zrat
      endif
c
c  call uvphrm which in spectral space converts vor and div into
c  ucos and vcos.
c
      call uvphrm
c
      call getime(cutime,tuvp)
c
c   initialize variables used for timing model components
      thxp = 0.0
      tgrm = 0.0
      thex = 0.0
      tpas = 0.0
      tder = 0.0
      t2ff = 0.0
c
#ifndef dirnal
c   update ktau temporarily for control block written on adsk & bdsk
        if (rsrt) then
          ktau=ktau+1
c
          call wrtctl(adsk, 0)
          call wrtctl(bdsk, 0)
c
          ktau=ktau-1
        endif
#endif
c
c
c  begin first transformation to grid point data for most of physics
c
      call passa
c
      if(prnt) then
        write(6,107) expnr, tau  , day
     *             , ihour, imin , isec , iday , month, iyear
        call energy
      endif
c
      call getime(cutime,eltime)
c
c  call linear which compute the horizontal diffusion in spectral space.
c
      call linear
c
c  call tderiv which performs a leap frog time step
c
      call tderiv
c
      call getime(cutime,eltime)
      tder=tder+eltime
c
c  second transformation to grid domain to do moist processes follows
c
      call passb
c
c  call rotq which time filters the data using a robert filter
c
      call robtq
c
#ifdef rad
c   compute global mean radiation tendency
      if (lrad) then
        radscl=86400.*cpp*og
        call spint ( radten, ps, gmrad, radscl )
        write(6,92) gmrad
      endif
#endif
c
#ifdef dirnal
      if(sgmt) then
       call gmtan
      endif
#endif
c
      call getime(cutime,eltime)
      t2ff=t2ff+eltime
c
      if(prnt) then
        write(6,108) tprec,train,tsnow,tevap,thtfx
      endif
c
      if (balnce) then
        write(6,337) gvft
        write(6,338) caht
#ifdef moist
        write(6,340) gmadj
        write(6,336) gvfr
#endif
        write(6,333) tuvp,thxp,tgrm,thex,tder,t2ff,tpas
      endif
c
c   print zonal mean quantities at each latitude
      if(mod(ktau,intzm).eq.0) call zmint
c
      sprec = tprec
      sevap = tevap
      sgvft = gvft
      scaht = caht
c
c   compute normalization terms for temperature and mixing ratio
      if (balnce) then
        psbar=rsrt2*real(ps(1,1,1))
        rpsbr=1./psbar
        call spint ( t  , ps, tbar, rpsbr )
        call spint (qmix, ps, qbar, og    )
        if (ktau .le. ktaub) then
          stbar=tbar
#ifdef moist
          sqbar=qbar
#endif
        else
          cnvly=grav/(cpp*psbar)
          strad=0.0
          strad=float(intbal)*sgrad/86400.
          fcq  =1.0
#ifdef moist
          if (qbar.ne.0.) fcq=(sqbar+sevap-sprec)/qbar
#endif
          sfct=sqrt(2.)*(cnvly*(sgvft+scaht+strad)+stbar-tbar)
          do 49 k=1,kx
           tm(1,1,k)=cmplx(real(tm(1,1,k))+sfct,0.0)
           t (1,1,k)=cmplx(real(t (1,1,k))+sfct,0.0)
 49       continue
c
          do 98 k=1,kx
            do 98 n=1,jx
              do 98 m=1,mx
                qmix (m,n,k) = fcq * qmix (m,n,k)
                qmixm(m,n,k) = fcq * qmixm(m,n,k)
   98     continue
c
          fct=sfct/sqrt(2.)
          write(6,344) fct,fcq
          call spint ( t  , ps, stbar, rpsbr )
          call spint (qmix, ps, sqbar, og    )
        endif
      endif
#ifdef rad
      sgrad=gmrad
#endif
c  load samta for anals code when lrad is .true.
      if (lrad) then
        rpsbr = 1.0/(rsrt2*real(ps(1,1,1)))
        call spint (t,ps,tabar,rpsbr)
        samta = samta+tabar+tfra
      endif
c
c  increment time step counter and forecast time in seconds
      ktau   = ktau+1
      kper   = kper+1
#ifdef cpldleg1
      kpermo = kpermo+1
#endif
      ifcstm = ifcst
      ifcst  = ifcst+idt
      nday   = ifcst/86400
      rem    = ifcst-nday*86400
      day    = float(nday)+rem/86400.0
      rsrt   = mod(ktau,intrs).eq.0
#ifdef cpldleg1
# ifdef seanal
c  compute month of this timestep for boundary data computation
      lmon   = .true.
      call march
      lmon   = .false.
# endif
#endif
#ifdef dirnal
      sgmt   = mod(ktau+1,igmt).eq.0
      if (mod(ifcst,86400).eq.0) rewind kdsk
#endif
c
      if (mod(ifcst,iseg) .eq. 0) then
#if !defined cpldleg3 && !defined cpldleg1
        call ascomm
#endif
#if defined cpldleg3 || defined cpldleg1
        call bcocn
#endif
      endif
c
      if ((day+float(idt)/86400.0) .ge. ap ) then
        call anals
        ap = ap + 365.0/12.0
      endif
c
c   compute time averaged hydrology variables over the restart interval
      if(rsrt) then
        rpfac=float(idt*intrs)/86400.
         pfac=1.0/rpfac
        write(6,360) rpfac
        do 72 jl=1,il
         do 72 i=1,ix
#ifdef moist
          precip(i,jl) = pfac*precip(i,jl)
          snowfl(i,jl) = pfac*snowfl(i,jl)
          evapor(i,jl) = pfac*evapor(i,jl)
          runoff(i,jl) = pfac*runoff(i,jl)
          snwmlt(i,jl) = pfac*snwmlt(i,jl)
          sublim(i,jl) = pfac*sublim(i,jl)
# if defined cpldleg3 || defined cpldleg1
          snwrnf(i,jl) = pfac*snwrnf(i,jl)
# endif
#endif
          ycemlt(i,jl) = pfac*ycemlt(i,jl)
   72   continue
c
c  call rstart (entry point in indata) to write out data to history
c   file.
c
        call rstart
c
      endif
c
#ifdef dirnal
      if(mod(ifcst,86400).eq.0) then
        ngmt=ngmt+1
        rewind 79
        write (79) ngmt,ktau,ltmm,igmt
        rewind 79
        write(6,159) ngmt
        nin  =iabs(nin-2)
        nou  =iabs(nin-2)
        adski=adsk+nin
        adsko=adsk+nou
        bdski=bdsk+nin
        bdsko=bdsk+nou
        lgi  =lg+nin
        lgo  =lg+nou
        lhi  =lh+nin
        lho  =lh+nou
        ltc  =0
        rewind adski
        rewind adsko
        rewind bdski
        rewind bdsko
        rewind lgi
        rewind lgo
        rewind lhi
        rewind lho
      else
        if(sgmt) then
          narray(1)   =ngmt+1
          ltc         =ltc+1
          call wrtctl(adsko, 0)
          call wrtctl(bdsko, 0)
          call wrtctl(lgo, 0)
          call wrtctl(lho, 0)
          read (adski)
          read (bdski)
          read (lgi  )
          read (lhi  )
        endif
      endif
#endif
c
c   return to main program at end of atmos day
      if(mod(ifcst,iseg) .eq. 0) go to 777
      if (ktau.gt.ktaust) go to 77
        first=.false.
        go to 1
c
  777 continue
      return
c
c
c
c   end of timestep loop
c
  77  continue
      iunit  = 6
      ifrmat = 1
      call wrtctl(iunit, ifrmat)
c
#ifdef dirnal
      if (nin.ne.0) call swapan
#endif
c
      call timing
c
      stop 77
c
#ifdef rad
   92 format(' truncated global mean radiation tendency=',e14.7,
     *   ' langleys/day'/)
#endif
  105 format('   the input starting timestep equals',i10)
  107 format (/' expnr = ',2a8,'  tau = ',f9.0,'  day = ',f9.3,
     *       3x,i2,1h:,i2,1h:,i2,2x,i2,1x,a4,1x,i4)
  108 format(' precip=',e14.7,3x,'rain  =',e14.7,3x,'snow  =',e14.7,
     *  3x,'evapor=',e14.7,3x,'heatfx=',e14.7/)
  140 format('  options for this run -',/,
     *       '   ktaust=',i10,/,'   intpr =',i10,/,'   intrs =',i10,/,
     *       '   intzm =',i10,/,'   intrd =',i10,/,'   intdk =',i10,/,
     *       '   intsr =',i10)
#ifdef dirnal
  159 format(' diurnal gmt averaging count updated to ',i3,' days'/)
#endif
  171 format(/,' ocean boundary data loaded at ktau = ',i10,
     *         '  to = ',e15.8,'  yce = ',e15.8)
#if defined seanal && defined cpldleg1
  173 format(' month ahead and behind =',2I3,
     *       ' weight ahead and behind =',2F7.3)
#endif
  333 format(' times for model components=',7e14.7/)
#ifdef moist
  336 format(' global mean vfr=',e14.7/)
#endif
  337 format(' global mean vft=',e14.7,' langleys'/)
  338 format(' global mean convective adjustment heating=',e14.7,
     *' langleys'/)
  340 format(' global mean negative mixing ratio=',e14.7/)
  344 format(' temperature normalization term=',e14.7,
     *       ' mixing ratio normalization factor=',e14.7/)
  360 format(' hydrology for restart averaged over ',f5.2,' days')
#ifdef cpldleg3
  407 format(//,' land-sea masks do not match',/,
     *          ' the atm and ocn grids do not match.')
  408 format(/,' the bad points are i=',i5,' jrow=',I5,/,
     *         ' zin and smask equal',2e16.8)
#endif
#ifdef impcit
  711 format(' tmean reset at ktau=',i10,7(/,1x,10f10.3))
#endif
#ifdef changdt
9006  format(//,' *** changing ifcst from ',i12,' to ',i12,/,
     $          ' *** changing  ktau from ',i12,' to ',i12,/,
     $          ' *** changing  year from ',i12,' to ',i12,/)
#endif
#ifdef cpldleg1
9085  format(//,' kpermo must be zero at the beginning of the run;'
     *          '  ktau, kpermo=',2i10,//)
9090  format(' zero arrays for storing monthly boundary data for ocean'
     *       ' at ktau=',i10)
#endif
#if defined cpldleg3 || defined cpldleg1
9107  format(//)
9181  format(' zero arrays for storing boundary data for ocean at',
     *       ' ktau=',i10)
#endif
9997  format(' ***changing topography gradually at day=',f8.3,
     *       '  tau=',f8.0,'  zratio=',e13.5)
c
      end
