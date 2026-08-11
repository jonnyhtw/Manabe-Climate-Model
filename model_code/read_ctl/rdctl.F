      subroutine rdctl(iunit)
c
c  This subroutine reads the control block
c
#include "ctlblk.h"

cnetcdf
#include "netcdf.inc"

#include "ncstuff_atm.h"

      real*4 r4word
      real*8 r8word
      integer*4 i4byte
      integer*4 iret

cnetcdf
c
c
c
#ifndef old_input_data
cnetcdf
      if (iunit .ne. 71) then
cnetcdf
      read (iunit)     expnr,           ktau,  ifcst, ifcstm,
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

c open netcdf unit
      iret = nf_open('atmos_r30_in.nc',0, ncid71)
      call check_err(iret)
c get ctlblk data from attributes
      iret = nf_get_att_text (ncid71, NF_GLOBAL,'expnam1',  expnr(1))
      call check_err(iret)
      iret = nf_get_att_text (ncid71, NF_GLOBAL,'expnam2',  expnr(2))
      call check_err(iret)
      iret = nf_get_att_real (ncid71, NF_GLOBAL,'rktau'   , r4word)
      call check_err(iret)
      ktau = r4word
      iret = nf_get_att_double(ncid71, NF_GLOBAL,'rifcst'  , r8word)
      call check_err(iret)
      ifcst = r8word
      iret = nf_get_att_double(ncid71, NF_GLOBAL,'rifcstm' , r8word)
      call check_err(iret)
      ifcstm = r8word
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'idt'     , i4byte)
      call check_err(iret)
      idt = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ls'      , i4byte)
      call check_err(iret)
      ls = i4byte 
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'jrow'    , i4byte)
      call check_err(iret)
      jrow = i4byte
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'day'   , day)
      call check_err(iret)
      iret = nf_get_att_real (ncid71, NF_GLOBAL,'rktaust' , r4word)
      call check_err(iret)
      ktaust = r4word
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'karray1',  i4byte)
      call check_err(iret)
      karray(1) = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'karray2',  i4byte)
      call check_err(iret)
      karray(2) = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'karray3',  i4byte)
      call check_err(iret)
      karray(3) = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'karray4',  i4byte)
      call check_err(iret)
      karray(4) = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'karray5',  i4byte)
      call check_err(iret)
      karray(5) = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intrs'  ,  i4byte)
      call check_err(iret)
      intrs = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intpr'  ,  i4byte)
      call check_err(iret)
      intpr = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intzm'  ,  i4byte)
      call check_err(iret)
      intzm = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intrd'  ,  i4byte)
      call check_err(iret)
      intrd = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intdk'  ,  i4byte)
      call check_err(iret)
      intdk = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'iday'   ,  i4byte)
      call check_err(iret)
      iday = i4byte
      iret = nf_get_att_text (ncid71, NF_GLOBAL,'month'  ,  month)
      call check_err(iret)
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'iyear'  ,  i4byte)
      call check_err(iret)
      iyear = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intbal' ,  i4byte)
      call check_err(iret)
      intbal = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'intsw'  ,  i4byte)
      call check_err(iret)
      intsw = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'nin'    ,  i4byte)
      call check_err(iret)
      nin = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'lts'    ,  i4byte)
      call check_err(iret)
      lts = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'igmt'   ,  i4byte)
      call check_err(iret)
      igmt = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ltc'    ,  i4byte)
      call check_err(iret)
      ltc = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'kanal'  ,  i4byte)
      call check_err(iret)
      kanal = i4byte
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sprec',  sprec)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sevap',  sevap)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sqbar',  sqbar)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'stbar',  stbar)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sgvft',  sgvft)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'scaht',  scaht)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'sgrad',  sgrad)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'fct'  ,  fct)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'fcq'  ,  fcq)
      call check_err(iret)
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'diss' ,  diss)
      call check_err(iret)
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'balnce' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        balnce = .false.
      else
        balnce = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'prnt'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        prnt = .false.
      else
        prnt = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'first'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        first = .false.
      else
        first = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'rsrt'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        rsrt = .false.
      else
        rsrt = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'seasnl' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        seasnl = .false.
      else
        seasnl = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'lrad'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        lrad = .false.
      else
        lrad = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'implct' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        implct = .false.
      else
        implct = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'resetm' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        resetm = .false.
      else
        resetm = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'resetp' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        resetp = .false.
      else
        resetp = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ldisk'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        ldisk = .false.
      else
        ldisk = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'kdisk'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        kdisk = .false.
      else
        kdisk = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'map'    ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        map = .false.
      else
        map = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'sgmt'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        sgmt = .false.
      else
        sgmt = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'jumps'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        jumps = .false.
      else
        jumps = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'resetd' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        resetd = .false.
      else
        resetd = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'resets' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        resets = .false.
      else
        resets = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'resetn' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        resetn = .false.
      else
        resetn = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'bndry'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        bndry = .false.
      else
        bndry = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'hybrid' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        hybrid = .false.
      else
        hybrid = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'diurnl' ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        diurnl = .false.
      else
        diurnl = .true.
      endif
      iret = nf_get_att_real (ncid71, NF_GLOBAL,'ktauo'  ,  r4word)
      call check_err(iret)
      ktauo = r4word
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'ifcsto' ,  r8word)
      call check_err(iret)
      ifcsto = r8word
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'idto'   ,  i4byte)
      call check_err(iret)
      idto = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'kpero'  ,  i4byte)
      call check_err(iret)
      kpero = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'mono'   ,  i4byte)
      call check_err(iret)
      mono = i4byte
      iret = nf_get_att_real (ncid71, NF_GLOBAL,'oyear'  ,  r4word)
      call check_err(iret)
      oyear = r4word
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ioprt'  ,  i4byte)
      call check_err(iret)
      ioprt = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'lmon'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        lmon = .false.
      else
        lmon = .true.
      endif
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'apo'  ,  apo)
      call check_err(iret)
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'monend' ,  i4byte)
      call check_err(iret)
      monend = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'mon'    ,  i4byte)
      call check_err(iret)
      mon = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'kper'   ,  i4byte)
      call check_err(iret)
      kper = i4byte
      iret = nf_get_att_double (ncid71, NF_GLOBAL,'ap'   ,  ap)
      call check_err(iret)
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'kpermo' ,  i4byte)
      call check_err(iret)
      kpermo = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ito'    ,  i4byte)
      call check_err(iret)
      ito = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'its'    ,  i4byte)
      call check_err(iret)
      its = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'iys'    ,  i4byte)
      call check_err(iret)
      iys = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'lrs'    ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        lrs = .false.
      else
        lrs = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'iseg'   ,  i4byte)
      call check_err(iret)
      iseg = i4byte
      iret = nf_get_att_text (ncid71, NF_GLOBAL,'leg'    ,  leg)
      call check_err(iret)
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'ih'     ,  i4byte)
      call check_err(iret)
      ih = i4byte
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'swamp'  ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        swamp = .false.
      else
        swamp = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'mlayer' , i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        mlayer = .false.
      else
        mlayer = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'qflx'   ,  i4byte)
      call check_err(iret)
      if (i4byte .eq. 0) then
        qflx = .false.
      else
        qflx = .true.
      endif
      iret = nf_get_att_int  (ncid71, NF_GLOBAL,'nmonth' ,  i4byte)
      call check_err(iret)
      nmonth = i4byte
      endif

cnetcdf
cSGI temp
      call wrtctl(6,1)
cSGI
#endif
#ifdef old_input_data
      read (iunit)     expnr(1),        ktau,  ifcst, ifcstm,
     *                   idt,     ls,   jrow,    day, ktaust,
     *                karray,
     *                 intrs,  intpr,  intzm,  intrd,  intdk,
     *                  iday,  month(1:4),  iyear, intbal,  intsw,
     *                   nin,    lts,   igmt,    ltc,  kanal,
     *                 sprec,  sevap,  sqbar,  stbar,  sgvft,
     *                 scaht,  sgrad,    fct,    fcq,   diss,
     *                balnce,   prnt,  first,   rsrt, seasnl,
     *                  lrad, implct, resetm, resetp,  ldisk,
     *                 kdisk,    map,   sgmt,  jumps, resetd,
     *                resets, resetn,  bndry, hybrid, diurnl,
     *                 ktauo, ifcsto,   idto,  kpero,   mono,
     *                 oyear,  ioprt, ktaubd,   apoa,   apob,
     *                   mon,   kper,    apa,    apb,    ito,
     *                   its,    iys,    lrs,   iseg,    leg,
     *                    ih,  swamp, mlayer,   qflx, nmonth 
#endif
#ifdef change_name
c
c quick fix of names for the simulation (cfi)
c
      expnr(1) = 'r15mix   '
      expnr(2) = 'qflux    '
#endif
c
c  write to output some info from control block
      write(6,400) ktau, iday, month, iyear, iunit
c
  400 format(' At ktau = ',i10,' The day month year is ',i3,a8,i10,
     * ' The control block is read on unit ',i3,'.')
c
      return
      end
