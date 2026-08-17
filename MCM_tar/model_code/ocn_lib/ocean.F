      subroutine ocean
c
c=======================================================================
c
c                        GFDL Modular Ocean Model
c
c
c             A primitive equation ocean model implemented by
c
c                R.C. Pacanowski   K. Dixon    A. Rosati
c
c                                at the
c
c              Geophysical Fluid Dynamics Laboratory /NOAA
c                         Princeton, NJ. 08542
c
c                         based on the work of
c
c         Kirk Bryan: A numerical method for the study of the  of the
c         circulation world ocean: 1969, J. Computat. Phys 4 347-376
c
c                                 and
c
c         the invaluable work of Mike Cox & Burt Semtner on earlier
c                       fortran implementations.
c
c
c        The GFDL Modular Ocean Model (acronym MOM) is a three
c        dimensional primitive equation ocean model intended  to be
c        a flexible tool useful for ocean and coupled air-sea modeling
c        applications over a wide range of space & time scales.
c        It is also intended to run efficiently on scalar, vector &
c        parallel architectures. The programming approach is
c        modular and additions to this model are encouraged to follow
c        this approach. Additional modules will be added with time and
c        new versions of this model will be released when ready. Please
c        refer to the READ_ME file for more details until the user's
c        guide is available.
c
c        Requirements:
c
c        Standard fortran 77 is used (except for use of namelist)
c        The preprocessor "cpp" (available on systems using "c" or UNIX)
c
c        Please refer questions or problems to:
c
c        R.C.P. e-mail: rcp@gfdl.princeton.edu    phone (609) 452-6527
c        K.D.   e-mail: kd@gfdl.princeton.edu     phone (609) 452-6574
c        A.R.   e-mail: ar@gfdl.princeton.edu     phone (609) 452-6584
c
c        Acknowlegments:
c
c        We wish to thank ...
c
c        Chris Kerr of Cray Research for his help in getting MOM to
c        autotask & multitask correctly.
c
c        Disclaimer:
c
c        MOM is an ocean modeling research tool at GFDL.
c        Others may freely use it but we assume no responsibility
c        for problems or incorrect use of MOM. To address this, we
c        have provided extensive diagnostic capability including
c        Global energetics & Regional term balances for the tracer and
c        momentum equations. It is left to the user to satisfy
c        (him/her)self that a particular configuration is working
c        correctly.
c
c
c=======================================================================
c
      external blkdta
c
#include "param.h"
#include "accel.h"
#include "cdiag.h"
#include "chmix.h"
#ifdef inland
#include "cinlnd.h"
#endif
#ifdef isopycmix
#include "cisop.h"
#endif
#include "coord.h"
#if defined fourfil || defined firfil
# include "cpolar.h"
#endif
#include "cprnts.h"
#include "cregin.h"
#ifdef ppvmix
# include "cppmix.h"
#endif
#include "crelax.h"
#include "ctask.h"
#include "ctmngr.h"
#include "cvbc.h"
#include "cvmix.h"
#include "emode.h"
#include "grdvar.h"
#ifdef selfdoc
#include "histjr.h"
#endif
#include "index.h"
#include "iounit.h"
#include "levind.h"
#include "scalar.h"
#ifdef multitasking
# include "cshrbf.h"
# include "timelv.h"
#else
# include "slabs.h"
#endif
#ifdef seaice
#include "seayce.h"
#endif
#if defined cpldleg2 && defined saltnorm
# include "conblk.h"
# include "sltdat.h"
      equivalence (sltsum,konblk(100))
#endif
#include "switch.h"
#include "temp.h"
#include "versno.h"
#ifdef isopycmix
c
      dimension ahhd(km)
c
#endif
#ifdef rigidlid
c
      dimension kmz(imt,jmt)
c
#endif
#ifndef readrmsk
c
      dimension ysr(nhreg), yer(nhreg), ylimr(nhreg+1)
c
c    "ysr", "yer" & "ylimr" are used to set up default horizontal
c                           regional masks (zonal bands)
c
#endif
c
      namelist /eddy/   am, ah, fkpm, fkph, bvdc, bvvc, cdbot
#if defined constvmix && defined blvmix
     $,                 afkph, dfkph, sfkph, zfkph
#endif
#ifdef tcvmix
     $,                 aq
#endif
#ifdef isopycmix
     $,                 ahisop, slmxr
#endif
#ifdef ppvmix
     $,                 fricmx, wndmix
#endif
     $,                 aidif
      namelist /tsteps/ dtts, dtuv, dtsf
      namelist /parms/  acor, mxscan, sor, crit
      namelist /contrl/ init, days, dgnstc, nmix, eb, restrt
     $,                 tsi, ntasks, ncon, travg, snaps
#ifdef kdtermb
     $,                 avterm
#endif
#ifdef worldocn
     $,                 histt
#endif
#if defined cpldleg3 || defined cpldleg2
     $,                 nsteps, pintro, prdiag
#endif
#ifdef islands
# ifdef oldcongrad
congrad
      namelist /ibox/ aisi, aiei, ajsi, ajei
# else
      namelist /iland/ alonis, alatis
# endif
#endif
#if defined cfctracers || defined ccl4tracer
      namelist /nmcfcs/ tryr0
#endif
#if defined cpldleg3 || defined cpldleg2
c
c     Varible "lo" is true the first time ocean is called during
c     any execution of a coupled model run or the first time called
c     during a "leg 2" ocean only model run as part of a coupled
c     model experiment.
c     When "lo" is true, all ocean initialization procedures found
c     before reaching the timestepping loop are performed (whether
c     as the very first timestep of a integration or as the first
c     timestep when picking up a run from a restart tape).
c     When "lo" is false, the model skips down to the timestepping
c     loop and integrates the ocean model "nsteps" timesteps before
c     ending.
c     In a scheme in which the coupled model executable is run in
c     one month increments, and the ocean-atmosphere coupling is done
c     once per day, (the ocean model integrates one model day per call)
c     the following holds.  Variable "lo" is true only on the first day
c     of a month, and "nsteps" would equal the number of ts timesteps
c     in a day. For all other days in the the month, "lo" is false.
c     In ocean-only models, "lo" is always true.  This is because the
c     ocean model is called only once per run, with each run being of
c     one month duration. Variable "nsteps" then is set to the a very
c     large number (9.e9) and the variable "days" is used to end the
c     integration after 1 month.
c     Variable "lo" is set in the program airsea which calls ocean.
c
#include "asblck.h"
#endif
c
c=======================================================================
c     begin introductory section
c=======================================================================
c
#if defined cpldleg3 || defined cpldleg2
      if( .not. lo) go to 3333
#endif
      open (stdin,file='ocean.in')
c
c     write MOM version info
c
      write (stdout,'(/25x,a80/)') model
c
c-----------------------------------------------------------------------
c     initialize various quantities
c-----------------------------------------------------------------------
c
      pi     = c4*atan(c1)
      radian = c360/(c2*pi)
      omega  = pi/43082.0
c
c-----------------------------------------------------------------------
c     read in run parameters
c-----------------------------------------------------------------------
c
      read  (stdin, contrl)
      write (stdout,contrl)
      read  (stdin, eddy)
      write (stdout,eddy)
      read  (stdin, tsteps)
      write (stdout,tsteps)
      read  (stdin, parms)
      write (stdout,parms)
#ifdef islands
# ifdef oldcongrad
congrad
      read (stdin, ibox)
      write(stdout,ibox)
# else
      read (stdin, iland)
      write(stdout,iland)
# endif
#endif
#if defined cfctracers || defined ccl4tracer
      read (stdin, nmcfcs)
      write(stdout,nmcfcs)
#endif
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c     print out tracer diffusion coefficients in isopycnal case
c-----------------------------------------------------------------------
c
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9101) am
      write(stdout,9102) ahisop, ah, slmxr
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
c
c-----------------------------------------------------------------------
c     set up model grids in x, y, & z
c-----------------------------------------------------------------------
c
      call grids
#if defined blhmix && defined worldocn
c
c     set lateral diffusivity for tracers to vary with k-level
c     as was done for old cyber runs. here, set ah to unity and
c     load diffusion coefficients into "hdcsck"
c
# ifdef isopycmix
      ahisop = c1
# else
      ah = c1
# endif
#ifdef cyberAH
      ahb = 500.
      ahm = 526.1
      epn = -.002
      do 80 k=1,km
        epn1  = epn * zt(k) / 100.0
        epn2  = exp(epn1)
        epn3  = ahm * epn2
        hdcsck(k) = (ahb + epn3)*10000.0*(5.0 - (k-1)/(4.0*(km-1)/km))
# ifdef isopycmix
        ahhd  (k) = ah + ahisop*hdcsck(k)
# endif
#else
c
c   Brian & Lewis(1979) horizontal diffusivity
ccccc ahs = 1.8e+3
ccccc ahb = 0.5e+3
      ahs = 5.0e+3
      ahb = 1.0e+3
ctd 1/7/97
      ahs = 1.5e+3
      ahb = 0.65e+3
ctd 1/7/97
      do 80 k=1,km
        hdcsck(k) = ahb + (ahs - ahb)*exp(-zt(k)/50000.0)
        hdcsck(k) = hdcsck(k) * 1.0e+4
# ifdef isopycmix
        ahhd  (k) = ah + ahisop*hdcsck(k)
# endif
#endif
c
# ifdef noAH
        hdcsck(k) = c0
# endif
# ifdef sameAH
        hdcsck(k) = ah
# endif
80    continue
c
#endif
#if defined consthmix && !defined isopycmix
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9105) ah
      write(stdout,9106) hdcsck
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
#ifdef isopycmix
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9105) ah, ahisop
      write(stdout,9106) hdcsck
      write(stdout,9115) ahhd
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
#if defined constvmix && defined blvmix
c
c     compute depth dependent vertical diffusion coefficients for
c     tracers using the relationship of Bryan and Lewis
c     (results are independent of base value "fkph")
c
      fkph = c1
      do 90 k=1,km
        vdcsck(k) = (afkph + (dfkph/pi)*(atan(sfkph*(zw(k) - zfkph))))
# ifdef sameFKPH
        vdcsck(k) = c1
# endif
90    continue
#endif
c
#ifdef constvmix
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9104) fkph
      write(stdout,9106) vdcsck
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
c
#if defined fourfil || defined firfil
c
c     compute sin and cos values for vector correction before filtering
c
      fx =  1.0e-10
      fxa = dxt(1)/radius
c
      do 110 i=2,imum1
        fxb = fxa*float(i-2)
        spsin(i) = sin(fxb)
        spcos(i) = cos(fxb)
        if (abs(spsin(i)) .lt. fx) spsin(i) = c0
        if (abs(spcos(i)) .lt. fx) spcos(i) = c0
110   continue
c
      spsin(1) = c0
      spcos(1) = c0
      spsin(imu) = c0
      spcos(imu) = c0
c
c     set up model indices for filtering high latitudes
c
      jfrst  = indp (rjfrst, yt, jmt)
      jft0   = indp (rjft0, yt, jmt)
      jft1   = indp (rjft1, yt, jmt)
      jft2   = indp (rjft2, yt, jmt)
      jfu0   = indp (rjfu0, yu, jmt)
      jfu1   = indp (rjfu1, yu, jmt)
      jfu2   = indp (rjfu2, yu, jmt)
      jskpt  = jft2-jft1
      jskpu  = jfu2-jfu1
      njtbft = (jft1-jfrst+1)+(jmtm1-jft2+1)
      njtbfu = (jfu1-jfrst+1)+(jmtm1-jfu2+1)
      if (njtbft .gt. jmtfil .or. njtbfu .gt. jmtfil) then
        write (stdout,9599) njtbft, njtbfu
        write (stderr,9599) njtbft, njtbfu
        stop '>ocean 1'
      endif
#  ifdef firfil
c
c     set "numflt" and "numflu" to filter more at higher latitudes
c     using "jft0" and "jfu0" as reference latitude rows.
c
      numfmx = imt * p25
      refcos = cst(jft0)
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9501) refcos, yt(jft0)
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
      do 170 j=jfrst,jmt
        if ((j .le. jft1 .or. j .ge. jft2) .and. j .ge. jfrst) then
          jj = j - jfrst + 1
          if (j .ge. jft2) jj = jj - jskpt + 1
          numflt(jj) = max(1,int(refcos/cst(j)))
          if (numflt(jj) .gt. numfmx) numflt(jj) = numfmx
# if defined cpldleg3 || defined cpldleg2
          if (pintro) then
# endif
          write(stdout,9502) j, jj, numflt(jj), yt(j)
# if defined cpldleg3 || defined cpldleg2
          endif
# endif
        endif
170   continue
      refcos = csu(jfu0)
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write(stdout,9503) refcos, yu(jfu0)
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
      do 180 j=jfrst,jmt
        if ((j .le. jfu1 .or. j .ge. jfu2) .and. j .ge. jfrst) then
          jj = j - jfrst + 1
          if (j .ge. jfu2) jj = jj - jskpu + 1
          numflu(jj) = max(1,int(refcos/csu(j)))
          if (numflu(jj) .gt. numfmx) numflu(jj) = numfmx
# if defined cpldleg3 || defined cpldleg2
          if (pintro) then
# endif
          write(stdout,9502) j, jj, numflu(jj), yu(j)
# if defined cpldleg3 || defined cpldleg2
          endif
# endif
        endif
180   continue
c
#  endif
#endif
c-----------------------------------------------------------------------
c     list code options chosen by user for this model run
c-----------------------------------------------------------------------
c
      call docmnt (stdout, .false., .false., .false., .false.)
c
c-----------------------------------------------------------------------
c     open the disk datasets
c-----------------------------------------------------------------------
c
      nbuf = 2
#ifdef crayio
c  nblk=512 because minimum block length is 2 (=512/512+1) in call wopen.
c  nblk=2  gives: UT009 - WARNING - BUFFERS REQUESTED .LT. MINIMUM,
c  DEFAULT USED (default is 2).
      call ostart (kontrl, nvar11, 512, 1)
#else
      call ostart (kontrl, nvar11, nvar11, 1)
#endif
      call ostart (kflds, nkflds*nwds, nwds, 1)
      call ostart (labs(1), jmt*nslab, nslab, nbuf)
      call ostart (labs(2), jmt*nslab, nslab, nbuf)
#ifdef twoslabdisks
      if (ntlev .eq. 3) then
        call ostart (labs(3), jmt*nslab, nslab, nbuf)
      endif
#else
      call ostart (labs(3), jmt*nslab, nslab, nbuf)
#endif
c
c
c=======================================================================
c     end introductory section
c=======================================================================
c
c
c
c-----------------------------------------------------------------------
c     is this a start from initial conditions or a restart?
c     in either case, set the number of vertical levels on the "t" grid
c-----------------------------------------------------------------------
c
cnetcdf
ccc   open (iorest, access='SEQUENTIAL', form='UNFORMATTED')
cnetcdf
c
      if (init) then
        call ocn1st
      else
        call rdrest
#if defined diskless && defined t90_1st
c       when reading in a C90 restart on the T90, one must read in the
c       integer kmt field separately and then put it out using a call
c       to oputi.  the read is accomplished by the ordkmt entry in odam.
c
        call ordkmt (iorest, kmt)
c
# ifdef f90type
        call oputi (kflds, nwds, 5*nwds+1, kmt)
# else
        write(stdout,8166)
        write(stderr,8166)
        stop 'F90TYP'
8166    format(//' ##### ERROR ##### f90type option must be'/
     $           ' ##### defined when t90_1st is defined'//)
# endif
#endif
c
# ifdef f90type
        call ogeti (kflds, nwds, 5*nwds+1, kmt)
# else
        call oget (kflds, nwds, 5*nwds+1, kmt)
# endif
      endif
c
c   verify that kmt=0 at j=jmt
      kmtj82 = 0
      do 8282 i=1,imt
        kmtj82 = kmtj82 + kmt(i,jmt)
8282  continue
      if (kmtj82 .ne. 0) then
        write(stdout,9982) (kmt(iw,jmt),iw=1,imt)
9982    format(/,' kmt is not 0 at j=jmt; kmt =',/,(1x,40i3))
        stop
      endif
c
c-----------------------------------------------------------------------
c     compute number of vertical levels on the "u" grid
c-----------------------------------------------------------------------
c
      do 800 j=1,jmt
        kmu(imt,j) = 0
800   continue
c
      do 900 i=1,imt
        kmu(i,jmt) = 0
900   continue
c
      do 1000 j=1,jmtm1
        do 990 i=1,imtm1
         kmu(i,j) = min (kmt(i,j), kmt(i+1,j), kmt(i,j+1), kmt(i+1,j+1))
990     continue
1000  continue
#ifdef cyclic
c
c     set cyclic conditions
c
      do 1100 j=1,jmt
        kmu(imt,j) = kmu(2,j)
1100  continue
#endif
#if defined worldocn && defined bumpkmu
# if defined fourfil || defined firfil
C
C SEARCH FOR PLACES IN THE OCEAN WHICH HAVE ONLY TWO VELOCITY
C    POINTS IN THE EAST-WEST DIRECTION AND ARE ALSO FILTERED BY THE
C    MODEL.  WHEN FOUND, MODIFY THE kmu VALUES AT THOSE POINTS IN
C    THE FOLLOWING MANNER TO AVOID COMPUTATIONAL INSTABILITIES.
C
#  if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#  endif
      write(stdout,9722) JFU1, JFU2
9722  format(' JFU1, JFU2 = ',2i8)
      write(stdout,9731)
9731  FORMAT(/,1X,'the following kmu values were modified to help'
     *       ,' ensure stable computations at velocity points',//)
#  if defined cpldleg3 || defined cpldleg2
      endif
#  endif
      DO 737 J=1,JMTM1
        IF ((J .LE. JFU1) .OR. (J .GE. JFU2)) THEN
          DO 738 I = 2, IMTM1
            KX  = kmu(I,J)
            KXX = KX
            IM1 = I - 1
            IM2 = I - 2
            IP1 = I + 1
            IP2 = I + 2
            IF (I .EQ. 2) IM2 = IMTM2
            IF (I .EQ. IMTM1) IP2 = 3
            IF (KX .NE. 0) THEN
              LCT = 0
              DO 739 K = 1, KM
                IF (kmu(IP1,J).GE.KX .AND. kmu(IM1,J).GE.KX) GO TO 741
                IF (kmt(IP1,J).GE.KX .AND. kmu(IP2,J).GE.KX) GO TO 741
                IF (kmt(IM1,J).GE.KX .AND. kmu(IM2,J).GE.KX) GO TO 741
                KX  = KX - 1
                LCT = 1
  739         CONTINUE
  741         CONTINUE
#  if defined cpldleg3 || defined cpldleg2
              if (pintro) then
#  endif
              IF (LCT .EQ. 1) write(stdout,9742) I,J,KXX,KX
 9742         FORMAT(1X,'I=',I4,2X,'J=',I4,2X,'OLD kmu=',I4,2X,
     *                 'NEW kmu=',I4)
#  if defined cpldleg3 || defined cpldleg2
              endif
#  endif
              kmu(I,J) = KX
            ENDIF
  738     CONTINUE
       ENDIF
  737 CONTINUE
c
#  ifdef cyclic
      do 6543 j=1,jmt
        kmu( 1 ,j) = kmu(IMTM1,j)
        kmu(IMT,j) = kmu(  2  ,j)
 6543 continue
#  endif
# endif
#endif
#ifdef symmetry
c
c     set symmetry conditions
c
      do 1200 i=1,imt
        kmu(i,jmt) = kmu(i,jmtm2)
1200  continue
#endif
#ifdef inland
c
c-----------------------------------------------------------------------
c     initialize inland sea communications
c-----------------------------------------------------------------------
c
      call inlnd1
#endif
c
c---------------------------------------------------------------------
c     set the horizontal regional masks and names to be used when
c     computing averages on the "t" grid in subroutine region.
c     also set the vertical regional masks and names for use (along
c     with the horizontal ones) in term balance calculations for
c     tracers & momentum
c---------------------------------------------------------------------
c
c     open unit for tracer averages
c
      if (iotavg .ne. stdout) open (iotavg, access='SEQUENTIAL'
     $,                             form='FORMATTED')
#ifdef testnorth
      write(stdout,123) iotrm2
123   format(' WILL PUT TERM BALANCES TO UNIT iotrm2 =',i3)
c
c    read in masks ("normsk") that divide the ocean into Atlantic/Artic
c    and Indo-Pacific sectors for meridional transport calculations
c    from file "northmasks" (unit ionrth)
c
      write(stdout,9149) ionrth
9149  format(//' READY TO READ IN northmasks AS UNIT',i3,/)
      open (ionrth,file='northmasks')
      call nor1st (ionrth)
      close (ionrth)
#endif
#ifdef readrmsk
c
c    read in horizontal & vertical regional masks ("mskhr" & "mskvr")
c    and names ("hregnm" & "vregnm") on file "regionmasks" (unit iormsk)
c
      open (iormsk,file='regionmasks')
      call reg1st (iormsk, .true., .false., .false., .false., .true.)
      close (iormsk)
#else
c
c     set up the horizontal regions:
c     specify "mskhr" and "hregnm" values rather than reading them in
c     from a file.  (arbitrarily defaulted below for 5 zonal bands)
c
      xlimr  = xt(2)
      xstart = xt(indp (xlimr, xt, imt))
      xlimr  = xt(imtm1)
      xend   = xt(indp (xlimr, xt, imt))
c
      ylimr(1) = -90.0
      ylimr(2) = -60.0
      ylimr(3) = -23.0
      ylimr(4) =  23.0
      ylimr(5) =  60.0
      ylimr(6) =  90.0
c
      do 185 nr=1,nhreg
        jsr = indp (ylimr(nr), yt, jmt)
        if (ylimr(nr) .ge. yt(jsr)) jsr = jsr + 1
        jer = indp (ylimr(nr+1), yt, jmt)
        if (ylimr(nr+1) .lt. yt(jer)) jer = jer - 1
        ysr(nr) = yt(jsr)
        yer(nr) = yt(jer)
        write (hregnm(nr),9301) xstart, xend, ysr(nr), yer(nr)
185   continue
c
      do 200 j=1,jmt
        jmask = 0
        do 190 nr=1,nhreg
          if (yt(j).ge.ysr(nr) .and. yt(j).le.yer(nr)) jmask = nr
190     continue
        do 195 i=1,imt
          mskhr(i,j) = 0
          if (kmt(i,j) .gt. 0 .and. xt(i) .ge. xstart
     $       .and. xt(i) .le. xend)  mskhr(i,j) = jmask
195     continue
#ifdef cyclic
        mskhr(1,j)   = 0
        mskhr(imt,j) = 0
#endif
200   continue
c
c     set up the vertical regions:
c     specify "mskvr" and "vregnm" values. also the starting & ending
c     levels for the vertical regions. (arbitrarily defaulted as shown)
c
      do 270 k=1,km
        mskvr(k) = 0
270   continue
      llvreg(1,1) = 1
      llvreg(1,2) = 4
      llvreg(2,1) = 5
      llvreg(2,2) = 9
      do 290 n=1,nvreg
        ks = llvreg(n,1)
        ke = llvreg(n,2)
        itop = nint((zw(ks)-dzt(ks))*0.01)
        ibot = nint(zw(ke)*0.01)
        write (vregnm(n),'(a6,i4,a5,i4,a1)') "  dpt: ",itop, "m to "
     $,                                      ibot, "m"
        do 280 k=ks,ke
          mskvr(k) = n
280     continue
290   continue
#endif
c
c     print out regional mask info
c
      call reg1st (iotavg, .true., .false., .false., .true., .false.)
c
c---------------------------------------------------------------------
c     compute the surface area and volume of the ocean regions
c     (values used in subroutine region are done in terms of meters,
c     rather than centimeters)
c---------------------------------------------------------------------
c
      areag = c0
      volgt  = c0
c
      do 300  k=1,km
        volgk(k) = c0
300   continue
c
      do 380 n=0,numreg
        areat(n) = c0
        areau(n) = c0
        volt(n)  = c0
        volu(n)  = c0
380   continue
      do 400 mask=1,nhreg
        areab(mask) = c0
        volbt(mask) = c0
        do 390  k=1,km
          volbk(mask,k) = c0
390     continue
400   continue
c
      do 500 j=2,jmtm1
        do 490 i=2,imtm1
          mask = mskhr(i,j)
          kz = kmt(i,j)
          if (kz .gt. 0) then
            boxat = cst(j) * dxt(i) * dyt(j)
            if (kmu(i,j) .ne. 0) then
              boxau = csu(j) * dxu(i) * dyu(j)
            else
              boxau = c0
            endif
            areat(0) = areat(0) +  boxat
            areau(0) = areau(0) +  boxau
            if (mask .gt. 0) areab(mask)  = areab(mask) + boxat * 1.e-4
            do 480 k=1,kz
              if (mask .gt. 0) then
                volbk(mask,k) = volbk(mask,k) + boxat * dzt(k) * 1.e-6
              endif
              dvolt   = boxat * dzt(k)
              if (kmu(i,j) .ge. k) then
                dvolu   = boxau * dzt(k)
              else
                dvolu   = c0
              endif
              volt(0) = volt(0) +  dvolt
              volu(0) = volu(0) +  dvolu
              n = nhreg*(mskvr(k)-1) + mask
              if (n .gt. 0 .and. mask .gt. 0) then
                volt(n) = volt(n) +  dvolt
                volu(n) = volu(n) +  dvolu
                do 470 nv=1,nvreg
                  ks = llvreg(nv,1)
                  if (k .eq. ks) then
                    areat(n) = areat(n) +  boxat
                    areau(n) = areau(n) +  boxau
                  endif
470             continue
              endif
480         continue
          endif
490     continue
500   continue
c
c     prevent possible divides by zero in subroutine diag2 by replacing
c     regions with zero areas or volumes with very small values
c
      do 550 n=0,numreg
        if (areat(n) .eq. c0) areat(n) = epsln
        if (areau(n) .eq. c0) areau(n) = epsln
        if (volt(n)  .eq. c0) volt(n)  = epsln
        if (volu(n)  .eq. c0) volu(n)  = epsln
550   continue
c
      do 600 mask=1,nhreg
        do 590 k=1,km
          volbt(mask) = volbt(mask) + volbk(mask,k)
          volgk(k) = volgk(k) + volbk(mask,k)
590     continue
        areag = areag + areab(mask)
        volgt = volgt + volbt(mask)
600   continue
c
      call reg1st (iotavg, .true., .true., .true., .false., .false.)
c
      area   = c0
      volume = c0
      ocnp   = 0
c
c     this comment directive turns off autotasking on the YMP
c     for the following loop
c
cfpp$ noconcur l
      do 700 j=2,jmtm1
        do 690 i=2,imtm1
          if (kmt(i,j) .gt. 0) then
            area   = area + cst(j)*dxt(i)*dyt(j)
            volume = volume + cst(j)*dxt(i)*dyt(j)*zw(kmt(i,j))
            ocnp   = ocnp + float(kmt(i,j))
          endif
690     continue
700   continue
c
      sum = c0
      do 790 n=1,numreg
        sum = sum + volt(n)
790   continue
      sum    = 100.0*sum/volume
      pctocn = 100.0*ocnp/float((imt-2)*(jmt-2)*km)
      diffa  = 100.0 * (c1 - (area - 10000.0*areag)/area)
c
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
      write (stdout,9341) area, volume
      write (stdout,9342) diffa, numreg, sum, pctocn
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
c
      if (init) then
c
c---------------------------------------------------------------------
c       print map of "kmt" levels. (to print "kmu", set ikmu=1)
c---------------------------------------------------------------------
c
        ikmu = 0
        if (ikmu .ne. 0) then
          write (stdout,9401)
        else
          write (stdout,9402)
        endif
        do 1300 ibk=1,imt,40
          isp  = ibk
          iept = ibk + 40 - 1
          if(iept.gt.imt) iept=imt
          write (stdout,'(/, 4x, 40i3)') (ii, ii=isp,iept)
          do 1290 jrev=1,jmt
            j=jmt-jrev+1
            if (ikmu .ne. 0) then
              write (stdout,'(1x,i3, 40i3)')
     $               j, (kmu(i,j),i=isp,iept)
            else
              write (stdout,'(1x,i3, 40i3)')
     $               j, (kmt(i,j),i=isp,iept)
            endif
1290      continue
1300    continue
      endif
#ifdef islands
# ifdef oldcongrad
congrad
c
c     convert island latitudes & longitudes to nearest model indices
c
      do 1009 l=1,nisle
        isis(l) = indp (aisi(l), xt, imt)
        ieis(l) = indp (aiei(l), xt, imt)
        jsis(l) = indp (ajsi(l), yt, jmt)
        jeis(l) = indp (ajei(l), yt, jmt)
        if (isis(l) .eq. 1)   isis(l) = 2
        if (ieis(l) .eq. imt) ieis(l) = imtm1
1009  continue
# else
c
c-----------------------------------------------------------------------
c     find the island perimeter paths
c-----------------------------------------------------------------------
c
      call iperim (kmz)
# endif
#endif
c
c---------------------------------------------------------------------
c     compute depths and reciprocal depths
c---------------------------------------------------------------------
c
      do 1400 j=1,jmt
        do 1390 i=1,imt
          hr(i,j) = c0
          h(i,j)  = c0
          if (kmu(i,j) .ne. 0) then
            hr(i,j) = c1/zw(kmu(i,j))
            h (i,j) = zw(kmu(i,j))
          endif
1390    continue
1400  continue
#ifdef rigidlid
c
c-----------------------------------------------------------------------
c     compute an array to indicate "interior" stream function grid boxes
c-----------------------------------------------------------------------
c
      do 1500 j=1,jmt
        kmz(1,j)   = 0
        kmz(imt,j) = 0
1500  continue
c
      do 1600 i=1,imt
        kmz(i,1)   = 0
        kmz(i,jmt) = 0
1600  continue
c
      do 1700 j=2,jmtm1
        do 1690 i=2,imu
          kmz(i,j) = min(kmu(i-1,j-1), kmu(i,j-1), kmu(i-1,j), kmu(i,j))
     $
1690    continue
1700  continue
#  ifdef cyclic
c
c     set cyclic conditions
c
      do 1800 j=1,jmt
        kmz(1,j) = kmz(imtm1,j)
1800  continue
#  endif
# ifdef neptune
c---------------------------------------------------------------------
c     compute the streamfunction "pnept" which represents topographic
c     stress in large scale ocean models ("the Neptune Effect", after
c     Greg Holloway, personal communication).
c     The transport streamfunction  tendency is "corrected" toward
c     "pnept" via appended relaxation ("slnept" is the scale length)
c---------------------------------------------------------------------
c
      slnept = 10.0e5
      do 1850 j=1,jmt
        f = c2*omega*sine(j)
        do 1840 i=1,imt
          pnept(i,j) = -f*slnept*zw(kmz(i,j))*zw(kmz(i,j))
          unept(i,j) = c0
          vnept(i,j) = c0
1840    continue
1850  continue
c
c     calculate depth independent velocity components that arise
c     from the Neptune effect
c
      do 1870 j=1,jmtm1
        do 1860 i=2,imtm1
            unept(i,j) = -(pnept(i,j+1)-pnept(i,j))*dyur(j)*
     $                     min(hr(i-1,j),hr(i,j))
            vnept(i,j) =  (pnept(i+1,j+1)-pnept(i,j+1))*dxur(i)*
     $                     min(hr(i,j+1),hr(i,j))*cstr(j+1)
1860    continue
# ifdef cyclic
        unept(1,j) =   unept(imtm1,j)
        unept(imt,j) = unept(2,j)
        vnept(1,j) =   vnept(imtm1,j)
        vnept(imt,j) = vnept(2,j)
# endif
1870  continue
c
# endif
c
c---------------------------------------------------------------------
c     compute longitudinal start & end indices for tracer, velocity &
c     stream function calculations
c---------------------------------------------------------------------
c
      do 1900 l=1,lseg
        do 1890 j=1,jmt
          isz(j,l) = 0
          iez(j,l) = 0
          ist(j,l) = 0
          iet(j,l) = 0
          isu(j,l) = 0
          ieu(j,l) = 0
1890    continue
1900  continue
      lsegp = lseg + 1
#ifdef skipland
c
c     indices for tracer calculations over ocean points
c
      do 2000 j=2,jmtm1
        l=1
        do 1990 i=2,imtm1
          if (kmt(i-1,j) .eq. 0 .and. kmt(i,j) .ne. 0) ist(j,l) = i
          if (i .eq. 2 .and. kmt(1,j) .ne. 0)          ist(j,l) = i
          if (kmt(i,j) .ne. 0 .and. kmt(i+1,j) .eq. 0) iet(j,l) = i
          if (i .eq. imtm1 .and. kmt(i+1,j) .ne. 0)    iet(j,l) = i
          if (kmt(i,j) .ne. 0 .and. kmt(i+1,j) .eq. 0) l=l+1
          if (l .gt. lsegp) then
            write (stderr,9499) ' tracer         ',j, lsegp
            stop ' ocean 2'
          endif
1990    continue
2000  continue
c
c     indices for velocity calculations over ocean points
c
      do 2100 j=2,jmtm1
        l=1
        do 2090 i=2,imum1
          if (kmu(i-1,j) .eq. 0 .and. kmu(i,j) .ne. 0) isu(j,l) = i
          if (i .eq. 2 .and. kmu(1,j) .ne. 0)          isu(j,l) = i
          if (kmu(i,j) .ne. 0 .and. kmu(i+1,j) .eq. 0) ieu(j,l) = i
          if (i .eq. imum1 .and. kmu(i+1,j) .ne. 0)    ieu(j,l) = i
          if (kmu(i,j) .ne. 0 .and. kmu(i+1,j) .eq. 0) l=l+1
          if (l .gt. lsegp) then
            write (stderr,9499) ' velocity       ',j, lsegp
            stop ' ocean 3'
          endif
2090    continue
2100  continue
#endif
c
c     indices for stream function calculation
c
      do 2200 j=3,jscan
        l=1
        do 2190 i=2,imum1
          if (kmz(i-1,j) .eq. 0. .and. kmz(i,j) .ne. 0.) isz(j,l)=i
          if (i .eq. 2 .and. kmz(imt,j) .ne. 0.) isz(j,l)=2
          if (kmz(i,j) .ne. 0. .and. kmz(i+1,j) .eq. 0.) iez(j,l)=i
          if (i .eq. imum1 .and. kmz(i+1,j) .ne. 0.) iez(j,l)=i
          if (kmz(i,j) .ne. 0. .and. kmz(i+1,j) .eq. 0.) l=l+1
          if (l .gt. lsegp) then
            write (stderr,9499) ' stream function',j, lsegp
            stop ' ocean 4'
          endif
2190    continue
2200  continue
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
#ifdef skipland
c
      write (stdout,9451) ' tracer '
c
      do 2300 jrev=1,jmt
        j = jmt - jrev + 1
        write (stdout,9461) j, (ist(j,l),iet(j,l),l=1,lseg)
2300  continue
c
      write (stdout,9451) ' velocity '
c
      do 2400 jrev=1,jmt
        j = jmt - jrev + 1
        write (stdout,9461) j, (isu(j,l),ieu(j,l),l=1,lseg)
2400  continue
      endif
#endif
c
      write (stdout,9451) ' stream function '
c
      do 2500 jrev=1,jmt
        j = jmt - jrev + 1
        write (stdout,9461) j, (isz(j,l),iez(j,l),l=1,lseg)
2500  continue
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#endif
#ifdef fourfil
c
c---------------------------------------------------------------------
c     find and print start & end indices for filtering
c---------------------------------------------------------------------
c
# if defined cpldleg3 || defined cpldleg2
c
      if (pintro) then
# endif
      write (stdout,9551)
      if (lsegf.gt.11) write (stdout,9552)
      write (stdout,9553)
# if defined cpldleg3 || defined cpldleg2
      endif
c
# endif
      call findex (kmt, jmtfil, km, jft1, jft2, imt, istf, ietf)
# if defined cpldleg3 || defined cpldleg2
c
      if (pintro) then
# endif
      write (stdout,9554)
# if defined cpldleg3 || defined cpldleg2
      endif
c
# endif
      call findex (kmu, jmtfil, km, jfu1, jfu2, imu, isuf, ieuf)
# if defined cpldleg3 || defined cpldleg2
c
      if (pintro) then
# endif
      write (stdout,9555)
# if defined cpldleg3 || defined cpldleg2
      endif
c
# endif
      call findex (kmz, jmtfil, 1, jfu1, jfu2, imt, iszf, iezf)
#endif
c
c---------------------------------------------------------------------
c     print out ts timestep multipliers
c---------------------------------------------------------------------
c
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write (stdout,9601) (dtxcel(k),k=1,km)
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
#ifdef isopycmix
c
c-----------------------------------------------------------------------
c     store the square root of the tracer timestep acceleration values
c     into variable "dtxsqr" for use in isopycnal mixing
c-----------------------------------------------------------------------
c
      do 2550 k=1,km
       dtxsqr(k) = sqrt(dtxcel(k))
2550  continue
#endif
c
c-----------------------------------------------------------------------
c     read in control block & set relative time in years for tmngr
c-----------------------------------------------------------------------
c
# ifdef f90type
#  ifdef gtime
      call ogeti (kontrl, nvar11, 1, itt)
#  else
      if (nvar11 .NE. 2) then
        write (stdout,9919) nvar11
        write (stderr,9919) nvar11
        stop 'NVAR11'
      else
        call ogeti (kontrl, 1, 1, itt)
        call oget  (kontrl, 1, 2, years)
      endif
9919  format (/' ### ERROR ### problem with F90TYPE option in ocean.F'/
     $         ' ############# code expected nvar11 = 2 but got ',i8,/)
#  endif
# else
      call oget (kontrl, nvar11, 1, itt)
# endif
#ifdef gtime
c     relyr is set in subroutine gtmngr
#else
      relyr = years
#endif
c
c-----------------------------------------------------------------------
c     compute permuting disc indicators and read in 2 levels of
c     stream function.
c-----------------------------------------------------------------------
c
      ndisk  = mod(itt+0,3) + 1
      ndiska = mod(itt+1,3) + 1
c
      call oget (kflds, nwds, (ndisk-1)*nwds+1, p(1,1,2))
      call oget (kflds, nwds, (ndiska-1)*nwds+1, p(1,1,1))
c
c-----------------------------------------------------------------------
c     initialize various things
c-----------------------------------------------------------------------
c
# ifdef sfctopog
c
c      to minimize errors in sea level relation the first time thru
c      read in a previous 2-D field of model sea level anomalies
c
cSGI
cc      open (iosealev, file='sealev.dat', access='SEQUENTIAL',
cc   $                  form='UNFORMATTED')
cSGI
cnetcdf
ccc     read(iosealev) sealev
cread in netcdf data
      iret = nf_open('sealevel_r30_in.nc',nf_read, ncid)
      call check_err(iret)
      iret = nf_inq_varid(ncid, "sealev", var_id)
      call check_err(iret)
      iret = nf_get_var_double (ncid, var_id, sealev)
      call check_err(iret)
cnetcdf
c
# endif
      do 2600 j=1,jmt
        do 2590 i=1,imt
          ztd(i,j) = c0
          zu(i,j)  = c0
          zv(i,j)  = c0
2590    continue
2600  continue
c
c     initialize the coriolis factor
c
      do 2660 j=1,jmt
        fcor(j) = c2*omega*sine(j)
2660  continue
c
c     initialize the diffusion factors
c
#ifdef biharmonic
      asubh = sqrt(ah)
      asubm = sqrt(am)
#endif
#ifdef consthmix
      asubh = ah
      asubm = am
#endif
#ifdef nlhmix
      asubh = c1
      asubm = c1
#endif
      do 2670 j=2,jmtm1
        bbt(j) = asubh*c8*cstr(j)*cstr(j)
        cct(j) = asubh*csu(j  )*dyur(j  )*dytr(j)*cstr(j)
        ddt(j) = asubh*csu(j-1)*dyur(j-1)*dytr(j)*cstr(j)
c
        bbu(j) = asubm*c8*csur(j)*csur(j)
        ccu(j) = asubm*cst(j+1)*dytr(j+1)*dyur(j)*csur(j)
        ddu(j) = asubm*cst(j  )*dytr(j  )*dyur(j)*csur(j)
        ggu(j) = asubm*(c1-tng(j)*tng(j))/(radius*radius)
        hhu(j) = asubm*c2*sine(j)/(radius*csu(j)*csu(j))
2670  continue
c
c-----------------------------------------------------------------------
c     set up the division of slabs into chunks by means of starting
c     and ending indices.  there may be one or more slabs per chunk.
c     each chunk will be a task to be completed by a single processor.
c     when "multitasking", multiple processors are deployed against the
c     available tasks. tasks may be done in any order but within a task
c     slabs are handled sequentially (in order).
c-----------------------------------------------------------------------
c
      nchunk = min(ntasks,jmt-2)
      inc   = max((nchunk+jmt-2)/nchunk, 1)
      jxtra = inc*nchunk - (jmt-1)
      jtemp = 1
# if defined cpldleg3 || defined cpldleg2
      if (pintro) then
# endif
      write (stdout,'(/,1x,"Tasking configuration: ")')
      write (stdout,'(/,1x,a6,4x,a10,4x,a10,4x,a10)') "task #"
     $,        "chunk size","start row","end row"
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
      do 2680 ntask=1,nchunk
        jsrowb(ntask) = jtemp
        jerowb(ntask) = min(jmt-1,jsrowb(ntask)+inc)
        if (ntask .le. jxtra) jerowb(ntask) = jerowb(ntask)-1
        jtemp         = jerowb(ntask)
        ilen          = jerowb(ntask)-jsrowb(ntask)+1
# if defined cpldleg3 || defined cpldleg2
        if (pintro) then
# endif
         write (stdout,'(1x,i5,5x,i10,4x,i10,4x,i10)')
     $              ntask, ilen, jsrowb(ntask), jerowb(ntask)
# if defined cpldleg3 || defined cpldleg2
      endif
# endif
2680  continue
c
c-----------------------------------------------------------------------
c     do all consistency checks now
c-----------------------------------------------------------------------
c
      call checks
c
#ifdef timing
c    measure cp & wall clock time neglecting overhead of "ocean"
c
      cps   = second()
      walls = timef()
c
#endif
c-----------------------------------------------------------------------
c     open units for term balances & saving data
c-----------------------------------------------------------------------
c
      open (iosnap, access='SEQUENTIAL', form='UNFORMATTED')
      if (iotrmb .ne. stdout) open (iotrmb, access='SEQUENTIAL'
     $,                             form='FORMATTED')
c
#ifdef worldocn
# ifdef seaice
      iceflg = .TRUE.
# else
      iceflg = .FALSE.
# endif
#endif
c
# ifdef saltnorm
c   Salt normalization term, sltsum, is saved from the previous run
c    in konblk(100)
c   totslt is computed in ocn1st, which is called only
c    on the initial run. Compute it here for subsequent runs.
      sinit = 0.0346786
      totslt = sinit - 0.035
# endif
c
#ifdef selfdoc
c     set logical flag for writing of cheat_sheet file to TRUE so
c     that it may be written just once per model execution
c
      lcheat1st = .TRUE.
c
#endif
#if defined cpldleg3 || defined cpldleg2
3333  continue
c
#endif
c-----------------------------------------------------------------------
c     start the time step loop
c-----------------------------------------------------------------------
c
#if defined cpldleg3 || defined cpldleg2
      do 3400 loop=1,nsteps
# ifdef cpldleg3
        loop1 = .false.
        if (loop .eq. 1) loop1 = .true.
# endif
        endseg = .false.
        if (loop .eq. nsteps) endseg = .true.
#else
      do 3400 loop=1,9999999
#endif
c
c-----------------------------------------------------------------------
c       update timestep counter, set corresponding model time, and set
c       all time dependent logical switches to determine program flow
c-----------------------------------------------------------------------
#ifdef gtime
c
c      use an updated time manager developed for use in G group
c      long term climate integrations
c
        call gtmngr (dtts)
c
ccc     write(stdout,1462) stamp, eoday, loday
1462    format(' >> ',a32,' eoday=',L4,' loday=',L4)
c
#else
c
c      use the MOM 1 time manager
c
        call tmngr (dtts)
c
#endif
c-----------------------------------------------------------------------
c       should a copy of the data be saved for off line analysis?
c-----------------------------------------------------------------------
c
        if (snapts) then
          write (iosnap) itt, dtts, relyr, stamp
          write (iosnap) xt, yt, zt
        endif
c
c-----------------------------------------------------------------------
c       adjust various quantities for normal/mixing timesteps
c-----------------------------------------------------------------------
c
        mxpas2 = .false.
        eots   = .true.
        if (mixts) then
          if (eb) eots = .false.
          c2dtts = dtts
          c2dtuv = dtuv
          c2dtsf = dtsf
          do 2700 j=1,jmt
            do 2690 i=1,imt
              p(i,j,2) = p(i,j,1)
2690        continue
2700      continue
        else
          c2dtts = c2*dtts
          c2dtuv = c2*dtuv
          c2dtsf = c2*dtsf
        endif
c
c-----------------------------------------------------------------------
c       update  pointers (indices) to tau-1, tau, & tau+1 discs
c-----------------------------------------------------------------------
c
#if defined diskless
#  ifdef multitasking
        ndiskb = mod(itt+2,ntlev) + 1
        ndisk  = mod(itt  ,ntlev) + 1
        ndiska = mod(itt+1,ntlev) + 1
#  else
        ndiskb = mod(itt+1,ntlev) + 1
        ndisk  = mod(itt  ,ntlev) + 1
        ndiska = ndiskb
#  endif
#else
# ifdef twoslabdisks
        ndiskb = mod(itt+1,ntlev) + 1
        ndisk  = mod(itt  ,ntlev) + 1
        ndiska = ndiskb
# else
        ndiskb = mod(itt+2,ntlev) + 1
        ndisk  = mod(itt  ,ntlev) + 1
        ndiska = mod(itt+1,ntlev) + 1
# endif
#endif
c
2800    continue
c
c-----------------------------------------------------------------------
c       initialize various quantities for on line integral analysis
c-----------------------------------------------------------------------
c
        if (prntsi) then
          do 2900 j=1,jmt
            do 2890 k=0,km
              ektot(k,j) = c0
              do 2880 m=1,nt
                dtabs(k,m,j) = c0
2880          continue
2890        continue
2900      continue
        endif
c
c-----------------------------------------------------------------------
c       the memory slab window is "ntau" time levels deep.
c       set time pointers (indices) for the memory slab window.
c       (note: these pointers (indices) are distinct from the disc
c        pointers (indices) .)
c       nc = current tau, np = tau + 1, nm = tau - 1
c       note: np is currently not used since tau+1 data is stored in
c       the buffer slab "ta, ua, va"
c-----------------------------------------------------------------------
c
        nm = mod(itt+1,ntau) + 1
        nc = mod(itt  ,ntau) + 1
        np = mod(itt+2,ntau) + 1
c
#ifdef seaice
c-----------------------------------------------------------------------
c     read fields of sea ice on disk into the proper timelevel positions
c-----------------------------------------------------------------------
c
        if (mixts) then
          call oget (kflds, nwds, (5+ndisk )*nwds+1, yc(1,1,nm))
        else
          call oget (kflds, nwds, (5+ndiskb)*nwds+1, yc(1,1,nm))
        endif
        call oget (kflds, nwds, (5+ndisk )*nwds+1, yc(1,1,nc))
c
#endif
c
c-----------------------------------------------------------------------
c       cycle chunks of slabs from disk into the memory slab window.
c       perform calculations then cycle chunks back to disk. on tracer
c       averaging or diagnostic time steps, allow only one chunk of
c       "jmt-2" slabs ... otherwise arbitrary chunk sizes are allowed.
c-----------------------------------------------------------------------
c
#ifdef worldocn
        if (diagts .or. tavgts .or. snapts .or. histts) then
#else
        if (diagts .or. tavgts .or. snapts) then
#endif
          ntask = 0
          call step (ntask)
        else
#ifdef multitasking
c
c         the parallel region starts with this loop.
c         everything called from "step.F" is within the parallel region
c
cmic$ doall
#endif
          do 3000 ntask=1,nchunk
            call step (ntask)
3000      continue
        endif
#ifdef seaice
c
c       put tau+1 sea ice thickness onto disk
c
        call oput (kflds, nwds, (5+ndiska)*nwds+1, yca)
c
#endif
c
c-----------------------------------------------------------------------
c       print on line integral analysis
c-----------------------------------------------------------------------
c
        if (eots .and. prntsi) then
          do 3100 j=2,jmtm1
            do 3070 k=km,1,-1
              ektot(0,1) = ektot(0,1) + ektot(k,j)
3070        continue
            do 3090 m=1,nt
              do 3080 k=km,1,-1
                dtabs(0,m,1) = dtabs(0,m,1) + dtabs(k,m,j)
3080          continue
3090        continue
3100      continue
          ektot(0,1) = ektot(0,1)/volume
          do 3200 m=1,nt
            dtabs(0,m,1)=dtabs(0,m,1)/volume
3200      continue
          write (stdout,9602) itt, stamp, ektot(0,1), dtabs(0,1,1)
     $,                      dtabs(0,2,1), mscan
#ifdef testrelx
     $,                      resmax, ires, jres
            write(stdout,3223) itt, wresmx(1), wresmx(mscan)
3223        format(2x,'relax info written in ocean.F at itt=',i10,
     $             5x,'wresmx=',2e12.3)
            write(iorelx) itt, ires, jres, resmax, crtp, mxscan, mscan,
     $                    (wresmx(iwre),iwre=1,mxscan)
            write(iorelx) itt, stamp, ektot(0,1)
#endif
        endif
c
c-----------------------------------------------------------------------
c       construct the vorticity
c-----------------------------------------------------------------------
c
        call vort
c
c-----------------------------------------------------------------------
c       solve for the new stream function
c-----------------------------------------------------------------------
c
#ifdef rigidlid
# ifdef congrad5pt
c       use 5 point laplacian
c
        call congr5
# endif
# ifdef congrad9pt
c       use 9 point laplacian allowing for more accurate reconstrucion
c       of the rigid lid surface pressure
c
        call congr9
# endif
# ifdef oldrelax
c       use sequential over-relaxation to solve the laplacian equation
c       the same as in the codes of Cox (1984) and Semtner (1974)
c
        call relax
# endif
# ifdef hypergrid
c       use sequential over-relaxation to solve the laplacian along
c       diagonals to allow for complete vectorization & parallelization
c
        call hyper
# endif
# ifdef neptune
c       apply a correction to the streamfunction in accordance with the
c       neptune effect (stnept = timescale for restoring streamfunction
c       toward "pnept")
c
        stnept = 60.0*86400.0
c
        do 3220 j=1,jmt
          do 3210 i=1,imt
            p(i,j,1) = p(i,j,1) + (pnept(i,j)-p(i,j,1))/stnept*dtsf
3210      continue
3220    continue
c
# endif
c
c-----------------------------------------------------------------------
c       should a copy of the data be saved for off line analysis?
c-----------------------------------------------------------------------
c
        if (eots .and. snapts) then
          do 3250 j=1,jmt
            do 3240 i=1,imt
              res(i,j) = p(i,j,1)
3240        continue
3250      continue
          write (iosnap) res
        endif
#endif
c
#ifdef worldocn
c
c-----------------------------------------------------------------------
c       write history data for analysis when "histts" is true
c-----------------------------------------------------------------------
c
        if (histts .and. eots) then
# ifdef sfctopog
c
c         call routine to compute surfpx and surfpy for history file
c         (as done in diag2 routine on diagnostic timesteps)
cccc      call calc_psgrad
c
          call calc_sea_level
# endif
          ktau = itt
          call histw (ktau)
        endif
#endif
c-----------------------------------------------------------------------
c       calculate & print the remaining diagnostics outside the parallel
c       region
c-----------------------------------------------------------------------
c
c
        if (eots .and. diagts) call diag2
#ifdef kdtermb
c
c       compute time averaged regional term balances and write them out
c
        if (eots .and. termav) call diag3
#endif
c
c-----------------------------------------------------------------------
c       if this is the end of the 1st pass of an euler backward timestep
c       set the input disc units so that the proper levels are fetched
c       on the next pass. the output for the 2nd pass will be placed on
c       the "ndiska" unit. return to the top to do the 2nd pass.
c-----------------------------------------------------------------------
c
        if (mixts .and. eb) then
          eots   = .true.
          mixts  = .false.
          mxpas2 = .true.
          ndiskx = ndiskb
          ndiskb = ndisk
          ndisk  = ndiska
          ndiska = ndiskx
          go to 2800
        endif
#if !defined diskless || defined multitasking
# if !defined twoslabdisks
c
c-----------------------------------------------------------------------
c       if this is the end of the 2nd pass of an euler backward timestep
c       transfer data which was temporarily written to "ndiskx" to its
c       final destination (the original "ndiska")
c-----------------------------------------------------------------------
c
        if (mxpas2) then
          ndiska = ndisk
          ndisk  = ndiskb
c
          do 3300 j=2,jmtm1
            call oget (labs(ndiskx), nslab, (j-1)*nslab+1, ta)
            call oput (labs(ndiska), nslab, (j-1)*nslab+1, ta)
3300      continue
# ifdef seaice
c
          call oget (kflds, nwds, (5+ndiskx)*nwds+1, yca)
          call oput (kflds, nwds, (5+ndiska)*nwds+1, yca)
# endif
        endif
# endif
#endif
c
c-----------------------------------------------------------------------
c       after everything is complete for this time step, update the
c       stream functions & control block on disk. also, set "first"
c       to .false.
c-----------------------------------------------------------------------
c
        years = relyr
        call oput (kflds, nwds, mod(itt+1,3)*nwds+1, p(1,1,1))
# ifdef f90type
#  ifdef gtime
        call oputi (kontrl, nvar11, 1, itt)
#  else
        call oputi (kontrl, 1, 1, itt)
        call oput  (kontrl, 1, 2, years)
#  endif
# else
        call oput (kontrl, nvar11, 1, itt)
# endif
        if (loop .eq. 1) first = .false.
c
        if (last) go to 3401
3400  continue
3401  continue
#ifdef timing
c
c     use cray specific timing calls
c
      cpe   = second()
      walle = timef()
      write (stdout,'(/,1x, a16, f9.3, a12, f9.3, a19)')
     $      'timing summary: ', cpe - cps,  ' cp seconds,',
     $                (walle - walls)/1000.0, ' wall clock seconds'
      write (stdout,'(16x, a1, g14.7, a42,/)')
     $    '(',(cpe-cps)/(loop*imt*jmt*km),
     $' cp seconds per grid point per time step)'
#endif
c
#if defined cpldleg3 || defined cpldleg2
      if ( .not. last) return
#endif
c
c-----------------------------------------------------------------------
c     documentation info at end of run
c-----------------------------------------------------------------------
c
      call docmnt (stdout, .true., .true., .true., .true.)
c
c-----------------------------------------------------------------------
c     save restart at end of run if restrt is true
c-----------------------------------------------------------------------
c
      if (restrt) call wrrest
c
c-----------------------------------------------------------------------
c     close all units
c-----------------------------------------------------------------------
c
      call oclose (kontrl)
      call oclose (kflds)
      call oclose (labs(1))
      call oclose (labs(2))
#ifdef twoslabdisks
      if (ntlev .eq. 3) then
       call oclose (labs(3))
      endif
#else
      call oclose (labs(3))
#endif
      close (unit=stdin,status='keep')
      if (stderr .ne. stdout) close (unit=stderr,status='keep')
      if (iotavg .ne. stdout) close (unit=iotavg,status='keep')
      if (iorest .ne. stdout) close (unit=iorest,status='keep')
      if (iotrmb .ne. stdout) close (unit=iotrmb,status='keep')
#ifdef worldocn
      close (unit=intrnl,status='keep')
      close (unit=iohist,status='keep')
#endif
      write (stdout,'(1x,a10,1x,a32)') 'MOMdone at',stamp
      close (unit=stdout,status='keep')
c
#if defined cpldleg3 || defined cpldleg2
      return
#else
      stop
#endif
c
#if defined isopycmix || defined worldocn
9101  format(/' am = ',e12.6,' cm**2/sec for horizontal mixing of ',
     $'momentum')
9102  format(/' ahisop = ',e12.6,' along isopyncal tracer mixing '
     $,  '(cm**2/sec) '/' ah = ',e12.6,' cm**2/sec for background '
     $,  'horizontal tracer diffusion'/' slmxr = ',e12.6,' to limit '
     $,  'slope used in computing mixing tensors')
#endif
#ifdef constvmix
9104  format(/' fkph = base for vertical mixing of tracers (cm**2/sec)'
     $ ,2x,e12.6)
#endif
#if defined consthmix && !defined isopycmix
9105  format(/' ah = ',e12.6,' cm**2/sec for horizontal mixing of ',
     $'tracers')
#endif
#ifdef isopycmix
9105  format(/' ah = ',e12.6,' cm**2/sec for extra background ',
     $'horizontal mixing of tracers',/
     $' ahisop = ',e12.6,' cm**2/sec base for isopycnal diffusion ',
     $'coefficient for tracers')
#endif
9106  format(' k-level scale factors '/,(8(2x,e12.6)))
9115  format(' k-level scale factors for extra background horizontal',
     $' mixing of tracers',/,8(2x,e12.6))
9301  format ('lon: ',f5.1,' => ',f5.1,'  lat: ',f5.1,' => ',f5.1,2i7)
9341  format (//,'  Regional & Global ocean statistics:'
     $,/,'  the total ocean surface area (t grid) =',1pe15.8,'cm**2'
     $,/,'  the total ocean volume (t grid)       =',1pe15.8,'cm**3')
9342  format ('  the horizontal regional masks cover',f8.3
     $, '% of the total ocean surface area.'/
     $, '  there are ', i6, ' regions over which tracer & '
     $, 'momentum balances will be computed',/,'  accounting for '
     $, f6.2,'% of the total ocean volume.'/
     $, 1x,f6.2,'% of the grid points lie within the ocean.'/)
9401  format(/t50,'number of levels on "u,v" grid')
9402  format(/t50,'number of levels on "t,s" grid')
9451  format (/' ==== start and end indices for',a17,'====')
#if defined cpldleg3 || defined cpldleg2
9461  format (' j=',i3,5x,10(2i4,2x))
#else
9461  format (' j=',i3,5x,5(2i5,10x))
#endif
9499  format (/' error => lseg too small for',a15,' indices'
     $        /'          j =',i5,'  lseg + 1 =',i8)
#ifdef firfil
9501  format(/' firfil reference cosine for tracers =',e12.6,' (',
     $       f7.3,' deg)'/'     jrow      jj numflt(jj)  latitude')
9502  format(1x,3i8,6x,f7.3)
9503  format(/' firfil reference cosine for velocities =',e12.6,' (',
     $       f7.3,' deg)'/'     jrow      jj numflu(jj)  latitude')
#endif
#if defined fourfil  || defined firfil
9551  format (/' ==== start and end indices for fourier filtering ====')
9552  format (' only 11 sets of indices fit across the page.',
     $       '  others will not be printed.'/)
9553  format (/,' == filtering indices for t,s ==')
9554  format (/,' == filtering indices for u,v ==')
9555  format (/,' == filtering indices for stream function ==')
9599  format (/,' error => jmtfil must be >= max(njtbft,njtbfu)',
     $        /,'          njtbft=',i8,' njtbfu=',i8)
#endif
9601  format(/,' "dtxcel(km)" tracer timestep multipliers:',/,10(f9.3))
9602  format (1x,'ts=',i9, 1x, a32, ', ke=', 1pe13.6,' dtemp=',1pe13.6
#ifdef testrelx
     $,        ' dsalt=',1pe13.6,' scans=',i3
     $,        ' res=',e8.2, 2i3)
#else
     $,        ' dsalt=',1pe13.6,' scans=',i3)
#endif
c
      end
