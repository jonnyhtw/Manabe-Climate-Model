      character*8 expnr, month, leg
cSGI
      integer*8  ktau,  ifcst, ifcstm
      integer*8   idt,     ls,   jrow,  ktaust
      integer*8  karray
      integer*8 intrs,  intpr,  intzm,  intrd,  intdk
      integer*8  iday,  iyear, intbal,  intsw
      integer*8   nin,    lts,   igmt,    ltc,  kanal
      integer*8 ktauo, ifcsto,   idto,  kpero,   mono
      integer*8 ioprt, monend
      integer*8   mon,   kper, kpermo,    ito
      integer*8   its,    iys,   iseg
      integer*8    ih, nmonth
      integer*8  narray
c
      logical*8 balnce,   prnt, first,   rsrt, seasnl, lrad, implct
      logical*8 resetm, resetp, ldisk,  kdisk,    map, sgmt
      logical*8 resets, resetn, bndry, hybrid, diurnl
      logical*8  jumps, resetd,   lrs,  swamp, mlayer, qflx
      logical*8 lmon
c
      integer*8  ktauo2, ifcst2,  idto2, kpero2,  mono2
      integer*8  ito2,   its2
cSGI
c
      common / ctlblk /   expnr(2),           ktau,  ifcst, ifcstm
      common / ctlblk /        idt,     ls,   jrow,    day, ktaust
      common / ctlblk /  karray(5)
      common / ctlblk /      intrs,  intpr,  intzm,  intrd,  intdk
      common / ctlblk /       iday,  month,  iyear, intbal,  intsw
      common / ctlblk /        nin,    lts,   igmt,    ltc,  kanal
      common / ctlblk /      sprec,  sevap,  sqbar,  stbar,  sgvft
      common / ctlblk /      scaht,  sgrad,    fct,    fcq,   diss
      common / ctlblk /     balnce,   prnt,  first,   rsrt, seasnl
      common / ctlblk /       lrad, implct, resetm, resetp,  ldisk
      common / ctlblk /      kdisk,    map,   sgmt,  jumps, resetd
      common / ctlblk /     resets, resetn,  bndry, hybrid, diurnl
      common / ctlblk /      ktauo, ifcsto,   idto,  kpero,   mono
      common / ctlblk /      oyear,  ioprt,   lmon,    apo, monend
      common / ctlblk /        mon,   kper,     ap, kpermo,    ito
      common / ctlblk /        its,    iys,    lrs,   iseg,    leg
      common / ctlblk /         ih,  swamp, mlayer,   qflx, nmonth
      common / ctlblk / narray(15)
c
c ocean control block (another copy of some of the above)
c This is a temporary "fix" to the problem that ocean variables
c are updated after the control block is written in rstart(indata)
c The real fix is to have a seperate control block for the ocean.
c
      common / ocnblk /  ktauo2, ifcst2,  idto2, kpero2,  mono2
      common / ocnblk /  oyear2,   apo2,   ito2,   its2
