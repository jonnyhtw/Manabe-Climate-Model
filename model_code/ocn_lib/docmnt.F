            subroutine docmnt (nunit, wrdoc, wrc0, wrc1, wrc2)
      implicit real*8 (a-h,o-z)
c
c=======================================================================
c
c     This code arranges data to be written out as a form of
c     documentation summarizing the model characteristics that
c     uniquely define this run.
c     In order to reap the full benefits of this feature, the user will
c     need to review this routine when configuring a model run or adding
c     new options to the code.
c     We envision this feature as being beneficial to those setting up
c     new model runs, comparing different model runs, and for analysis
c     purposes (much of the information used in setting up mapping
c     procedures and other analysis products are contained in the
c     data written by this routine).
c
c     Five namelist writes are available.  Which ones are written
c     is determined from the 4 logical variables passed into the
c     subroutine (wrdoc, wrc0, wrc1, & wrc2).  The unit to which
c     they are written is set by "nunit".
c
c     If "wrdoc" is true, then the "ndocs" & "docums" namelists
c     are written...
c       The 1st, "ndocs" simply contains one variable (ndoc); the
c     dimension of "docum(ndoc)" found in the 2nd namelist write.
c       The 2nd, "docums" consists of a series of 60 character-long
c     lines of information the user should use to list the variables
c     in the subsequent namelist writes, and to describe the model
c     run in general.
c
c     If "wrc0" is true, then the "cntrl0" namelist is written...
c       Namelist "cntrl0" contains the 2 variables needed to restart
c     the model (itt, years).
c
c       The "cntrl1" and "cntrl2" namelists contain information
c     describing the model's configuration and specifications.
c
c     If "wrc1" is true, then the "cntrl1" namelist is written...
c       "cntrl1" contains variables (non-dimensioned) that will be
c     common to most, if not all, incarnations of the model code.
c
c     If "wrc2" is true, then the "cntrl2" namelist is written...
c       "cntrl2" is intended to contain the larger (dimensioned)
c     data that will help identify the model specifics to users
c     and analysis programs alike.
c
c=======================================================================
c
#include "param.h"
#include "accel.h"
#if defined cpldleg3 || defined cpldleg2
#include "cexpnam.h"
#endif
#include "chmix.h"
#ifdef inland
#include "cinlnd.h"
#endif
#ifdef isopycmix
#include "cisop.h"
#endif
#include "coord.h"
#include "cregin.h"
#include "ctmngr.h"
#include "cvmix.h"
#include "docnam.h"
#include "index.h"
#include "scalar.h"
#ifdef restorst
#include "slabs.h"
#endif
#include "switch.h"
#include "versno.h"
c
c   "ndocum" is the number of 60 character lines of documentation in
c      docum(ndocum)" describing the model written by namelist.
c   "ifdmax" is the maximum number of 12 character ifdef options (in
c      (dimension of "cifdef") that may be listed as being set in the
c       control block (both ifdmax & the 12-character length can be
c       user specified)
c   "nifdef" is the actual number of ifdef options that are set.
c
      parameter (ndocum = 22)
      parameter (ifdmax = 50)
c
      character * 60 docum(ndocum)
#if !defined cpldleg3 && !defined cpldleg2
      character * 60 expnam
#endif
      character * 12 cifdef(ifdmax), blnk12
c
c   the logical variables "wrdoc, wrc0, wrc1, &  wrc2" let the user
c   set which namelists are to be written to unit "nunit".
c
      logical wrdoc, wrc0, wrc1, wrc2
c
      namelist /ndocs / ndoc
c
      namelist /docums/ docum
c
      namelist /cntrl0/
     $     itt, years
c
      namelist /cntrl1/
     $      expnam, tyear, tmonth, tday, thour, tmin, tsec, relyr
     $     ,time0, julian, yrnorm, yrleap, nleap
     $     ,dtts, dtuv, dtsf, init, nmix, eb, acor, am, ah, fkpm,fkph
#if defined constvmix && defined blvmix
     $     ,fkph, afkph, dfkph, sfkph, zfkph   
#endif
#ifdef isopycmix
     $     ,ahisop, slmxr
#endif
     $     ,id, jd, kd, ntrace, lsegs, nisl, area, volume
     $     ,nhorr, nverr, ntreg, ifdmx, nifdef
#ifdef fourfil
     $     ,lsegsf, jfrst, jft0, jft1, jft2, jfu0, jfu1, jfu2
#endif
#ifdef firfil
     $     ,jfrst, jft0, jft1, jft2, jfu0, jfu1, jfu2
#endif
#ifdef restorst
     $     ,rests
#endif
#ifdef inland
     $     ,ninld
#endif
c
      namelist /cntrl2/
     $     xt, yt, xu, yu, zt
     $     ,trname, cifdef, hregnm, mskhr, vregnm, mskvr
#ifdef islands
     $     ,alonis, alatis
#endif
     $     ,dtxcel
#ifdef firfil
     $     ,numflt, numflu
#endif
#ifdef constvmix
     $     ,vdcsci, vdcscj, vdcsck, vvcsci, vvcscj, vvcsck
#endif
#if defined constmix || defined isopycmix
     $     ,hdcsck, hvcsck
#endif
#ifdef inland
     $     ,iinlnd, jinlnd, kinlnd, vinlnd, binlnd
#endif
c
      save nifdef, cifdef
c
      data (docum(n),n=1,ndocum)
     1   /'namelist  cntrl0 = itt, years                               '
     2   ,'namelist  cntrl1 =                                          '
     3   ,' expnam, tyear, tmonth, tday, thour, tmin, tsec, relyr      '
     4   ,',time0, julian, yrnorm, yrleap, nleap                       '
     5   ,',dtts, dtuv, dtsf, init, nmix, eb, acor, am, ah, fkpm, fkph '
#if defined constvmix && defined blvmix
     6   ,',afkph, dfkph, sfkph, zfkph                                 '
#else
     6   ,' ==> constvmix && blvmix not set <===                       '
#endif
#ifdef isopycmix
     7   ,',ahisop, slmxr                                              '
#else
     7   ,'==> isopycmix not set <===                                  '
#endif
     8   ,',id, jd, kd, ntrace, lsegs, nisl, area, volume              '
     9   ,',nhorr, nverr, ntreg, ifdmx, nifdef                         '
#ifdef fourfil
     $   ,',lsegsf, jfrst, jft0, jft1, jft2, jfu0, jfu1, jfu2          '
#endif
#ifdef firfil
     $   ,',jft0, jft1, jft2, jfu0, jfu1, jfu2                         '
#endif
#ifndef firfil
#  ifndef fourfil
     $   ,' ==> no filtering set <===                                  '
#  endif
#endif
#ifdef restorst
     1   ,',rests                                                      '
#else
     1   ,'==> restorst not set <===                                   '
#endif
     2   ,'namelist  cntrl2 =                                          '
     3   ,'xt(id), yt(jd), xu(id), yu(jd), zt(kd)                      '
     4   ,'trname, cifdef, hregnm, mskhr, vregnm, mskvr                '
#ifdef islands
     5   ,',alonis(nisl), alatis(nisl)                                 '
#else
     5   ,' ==> islands ifdef not set <==                              '
#endif
     6   ,',dtxcel(km)                                                 '
#ifdef firfil
     7   ,',numflt(jmtfil), numflu(jmtfil)                             '
#else
     7   ,' ==> firfil not set  <===                                   '
#endif
#ifdef constvmix
     8   ,',vdcsci(id), vdcscj(jd), vdcsck(kd)                         '
     9   ,',vvcsci(id), vvcscj(jd), vvcsck(kd)                         '
#else
     8   ,' ==> constvmix not set  <===                                '
     9   ,' ==> constvmix not set  <===                                '
#endif
#if defined constvmix || defined isopycmix
     $   ,',hdcsck(kd), hvcsck(kd)                                     '
#else
     $   ,'  ==> neither consthmix nor isopycmix set  <===             '
#endif
#ifdef inland
     1   ,',iinlnd(ninld,2), jinlnd(ninld), kinlnd(ninld,2)            '
     2   ,',vinlnd(ninld), binlnd(ninld,2)                             '
#else
     1   ,' ==> inland not set  <===                                   '
     2   ,' ==> inland not set  <===                                   '
#endif
     $    /
      data blnk12 /'            '/
c
c     The experiment name is defaulted to the MOM version info here
c
#if defined cpldleg3 || defined cpldleg2
c   For leg3 or leg2 of coupled models, the experiment name is loaded
c    into expn1, expn2 in  program airsea
      expnam( 1: 8) = expn1
      expnam( 9:16) = expn2
      expnam(17:60) = ' '
#else
      expnam = model
#endif
c
c     User specified tracer names are placed into "trname" here.
c
      do 100 m=1,nt
        trname(m) = '**unknown***'
100   continue
      trname(1) = 'potentl temp'
      trname(2) = 'salinity    '
#ifdef cybertra
      trname(3) = 'tracera     '
#endif
#ifdef seaice
      ycname = 'sea ice     '
#endif
#ifdef worldocn
c
# ifdef agetracers
c     set idealized age tracer name used in GFDL world ocean model
c
        trname(ntage1) = 'age tracer  '
# endif
# ifdef cfctracers
c     set chlorofluorocarbon tracer names used in GFDL world ocean model
c
        trname(ntcfc1) = 'cfc-11      '
        trname(ntcfc2) = 'cfc-12      '
# endif
# ifdef ccl4tracer
c     set carbon tetrachloride tracer name
c
        trname(ntccl4) = 'carbon-tet  '
# endif
#endif
c
c-----------------------------------------------------------------------
c     list options that are in effect on first timestep
c-----------------------------------------------------------------------
c
      if (first) then
c
#if defined cpldleg3 || defined cpldleg2
      if (pintro) then
#endif
        write (stdout,9001)
        write (stdout,9002)
        write (stdout,9003)
        write (stdout,9004)
        write (stdout,9005)
        write (stdout,9006)
#ifdef worldocn
        write (stdout,9007)
        write (stdout,9008)
#endif
c
c     set ifdef list for control block info.
c     (maximum count set in parameter "ifdmax")
c
        blnk12 = '            '
        do 200 n=1,ifdmax
          cifdef(n) = blnk12
200     continue
c
        nifdef = 0
        cifdef(1) = '0 ifdefs set'
c
        write (stdout,'(/,a45)')
     $   ' The following "ifdef" options are in effect:'
c
c     >> optimization options <<
#ifdef multitasking
        nifdef = nifdef + 1
        cifdef(nifdef) = 'multitasking'
#endif
#ifdef timing
        nifdef = nifdef + 1
        cifdef(nifdef) = 'timing      '
#endif
c
c     >> grid options <<
#ifdef cyclic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyclic      '
#endif
#ifdef skipland
        nifdef = nifdef + 1
        cifdef(nifdef) = 'skipland    '
#endif
#ifdef symmetry
        nifdef = nifdef + 1
        cifdef(nifdef) = 'symmetry    '
#endif
c
c     >> filtering options <<
#ifdef firfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'firfil      '
        write (stdout,*) ' note ==> firfil not fully checked out yet'
#endif
#ifdef fourfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fourfil     '
#endif
c
c     >> external mode velocity options <<
#ifdef freesfc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'freesfc     '
        write (stdout,*) ' stop ==> free surface not implemented yet'
        if (nifdef .lt. 10000) stop '2 docum'
#endif
#ifdef rigidlid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'rigidlid    '
#  ifdef sfctopog
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sfctopog    '
#  endif
#  ifdef congrad5pt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'congrad5pt  '
#  endif
#  ifdef congrad9pt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'congrad9pt  '
#  endif
#  ifdef hypergrid
        nifdef = nifdef + 1
        cifdef(nifdef) = 'hypergrid   '
#  endif
#  ifdef islands
        nifdef = nifdef + 1
        cifdef(nifdef) = 'islands     '
#  endif
#  ifdef neptune
        nifdef = nifdef + 1
        cifdef(nifdef) = 'neptune     '
#  endif
#  ifdef oldrelax
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldrelax    '
#  endif
#endif
c
c     >> vertical mixing scheme options <<
#ifdef constvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'constvmix   '
#endif
#ifdef implicitvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'implicitvmix'
#endif
#ifdef ppvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ppvmix      '
        write (stdout,*) ' warning ==> ppvmix not fully checked out yet'
#endif
#ifdef tcvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'tcvmix      '
        write (stdout,*) ' warning ==> tcvmix not fully checked out yet'
#endif
#ifdef leq
        nifdef = nifdef + 1
        cifdef(nifdef) = 'leq         '
        write (stdout,*) ' warning ==> leq not fully checked out yet'
#endif
#ifdef lalg
        nifdef = nifdef + 1
        cifdef(nifdef) = 'lalg        '
        write (stdout,*) ' warning ==> lalg not fully checked out yet'
#endif
c
c     >> horizontal mixing scheme options <<
#ifdef biharmonic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'biharmonic  '
        write (stdout,*) ' warning ==> biharmonic not fully checked'
#endif
#ifdef consthmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'consthmix   '
#endif
#ifdef isopycmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'isopycmix   '
#endif
#ifdef nlhmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'nlhmix  '
        write (stdout,*) ' warning ==> nlhmix not fully checked'
#endif
c
c     >> i/o options <<
#ifdef crayio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'crayio      '
#endif
#ifdef diskless
        nifdef = nifdef + 1
        cifdef(nifdef) = 'diskless    '
#endif
#ifdef fio
        nifdef = nifdef + 1
        cifdef(nifdef) = 'fio         '
#endif
#ifdef f90type
        nifdef = nifdef + 1
        cifdef(nifdef) = 'f90type     '
#endif
#ifdef twoslabdisks
        nifdef = nifdef + 1
        cifdef(nifdef) = 'twoslabdisks'
#endif
#ifdef t90_1st
        nifdef = nifdef + 1
        cifdef(nifdef) = 't90_1st     '
#endif
c
c     >> miscellaneous options <<
#ifdef kdtermb
        nifdef = nifdef + 1
        cifdef(nifdef) = 'kdtermb     '
#endif
#ifdef kdtravg
        nifdef = nifdef + 1
        cifdef(nifdef) = 'kdtravg     '
#endif
#ifdef keepterms
        nifdef = nifdef + 1
        cifdef(nifdef) = 'keepterms   '
#endif
#ifdef nohilats
        nifdef = nifdef + 1
        cifdef(nifdef) = 'nohilats    '
#endif
#ifdef readrmsk
        nifdef = nifdef + 1
        cifdef(nifdef) = 'readrmsk  '
#endif
#ifdef restorst
        nifdef = nifdef + 1
        cifdef(nifdef) = 'restorst    '
#endif
#ifdef testcfl
        nifdef = nifdef + 1
        cifdef(nifdef) = 'testcfl     '
#endif
# ifdef cflj
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cflj        '
# endif
#ifdef testnorth
        nifdef = nifdef + 1
        cifdef(nifdef) = 'testnorth   '
#endif
c
c     >>  options developed for GFDL world ocean runs  <<
#ifdef worldocn
        nifdef = nifdef + 1
        cifdef(nifdef) = 'worldocn    '
#endif
#ifdef agetracers
        nifdef = nifdef + 1
        cifdef(nifdef) = 'agetracers  '  
#endif
#ifdef amazon
        nifdef = nifdef + 1
        cifdef(nifdef) = 'amazon      '  
#endif
#ifdef bctest
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bctest      '
#endif
#ifdef blhmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'blhmix      '  
#endif
#ifdef blvmix
        nifdef = nifdef + 1
        cifdef(nifdef) = 'blvmix      '  
#endif
#ifdef bumpkmu
        nifdef = nifdef + 1
        cifdef(nifdef) = 'bumpkmu     '  
#endif
#ifdef ccl4tracer
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ccl4tracer   '  
#endif
#ifdef cfctracers
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cfctracers  '  
#endif
#ifdef cjbcfcs
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cjbcfcs     '  
#endif
#ifdef gptcfcs
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gptcfcs     '  
#endif
#ifdef cpldleg2
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cpldleg2     '
#endif
#ifdef cpldleg3
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cpldleg3     '
#endif
#ifdef cyberAH 
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyberAH     '  
#endif
#ifdef cyberfil
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cyberfil    '  
#endif
#ifdef cybertra
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cybertra    '  
#endif
#ifdef cybervbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cybervbc    '
#endif
#ifdef gtime
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gtime       '
# ifdef gotogtime    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'gotogtime   '  
# endif
# ifdef cpldinit    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'cpldinit    '  
# endif
#endif
#ifdef inland
        nifdef = nifdef + 1
        cifdef(nifdef) = 'inland      '  
#endif
#ifdef january
        nifdef = nifdef + 1
        cifdef(nifdef) = 'january     '  
#endif
#ifdef kdtest
        nifdef = nifdef + 1
        cifdef(nifdef) = 'kdtest      '  
#endif
#ifdef km12
        nifdef = nifdef + 1
        cifdef(nifdef) = 'km12        '  
#endif
#ifdef km18
        nifdef = nifdef + 1
        cifdef(nifdef) = 'km18        '  
#endif
#ifdef km25
        nifdef = nifdef + 1
        cifdef(nifdef) = 'km25        '  
#endif
#ifdef noAH    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'noAH        '  
#endif
#ifdef noadvection  
        nifdef = nifdef + 1
        cifdef(nifdef) = 'noadvection '  
#endif
#ifdef noclinic
        nifdef = nifdef + 1
        cifdef(nifdef) = 'noclinic    '  
#endif
#ifdef noramp  
        nifdef = nifdef + 1
        cifdef(nifdef) = 'noramp      '  
#endif
#ifdef notracer
        nifdef = nifdef + 1
        cifdef(nifdef) = 'notracer    '  
#endif
#ifdef ocnr15
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ocnr15      '  
#endif
#ifdef ocnr30
        nifdef = nifdef + 1
        cifdef(nifdef) = 'ocnr30      '  
#endif
#ifdef oldcongrad
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldcongrad  '  
#endif
#ifdef old_dtxcel
        nifdef = nifdef + 1
        cifdef(nifdef) = 'old_dtxcel  '  
#endif
#ifdef old_dzt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'old_dzt     '  
#endif
#ifdef oldstate
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldstate    '
#endif
#ifdef oldtime
        nifdef = nifdef + 1
        cifdef(nifdef) = 'oldtime     '  
#endif
#ifdef onetoone
        nifdef = nifdef + 1
        cifdef(nifdef) = 'onetoone    '  
#endif
#ifdef r15vbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'r15vbc      '
#endif
#ifdef saltnorm
        nifdef = nifdef + 1
        cifdef(nifdef) = 'saltnorm    '
#endif
#ifdef sameAH    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sameAH      '  
#endif
#ifdef sameFKPH    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'sameFKPH    '  
#endif
#ifdef seaice
        nifdef = nifdef + 1
        cifdef(nifdef) = 'seaice      '
#endif
#ifdef seasonal
        nifdef = nifdef + 1
        cifdef(nifdef) = 'seasonal    '
#endif
#ifdef selfdoc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'selfdoc     '
#endif
#ifdef testrelx    
        nifdef = nifdef + 1
        cifdef(nifdef) = 'testrelx    '  
#endif
#ifdef windbc
        nifdef = nifdef + 1
        cifdef(nifdef) = 'windbc      '
#endif
#ifdef writekmt
        nifdef = nifdef + 1
        cifdef(nifdef) = 'writekmt    '  
#endif
#ifdef yearzero
        nifdef = nifdef + 1
        cifdef(nifdef) = 'yearzero    '  
#endif
#ifdef zonalavg
        nifdef = nifdef + 1
        cifdef(nifdef) = 'zonalavg    '  
#endif
        nn = (nifdef+4)/5
        do 300 i=1,nn
          ns = 1+(5*(i-1))
          ne = ns + 4
          if (ne .gt. nifdef) ne = nifdef
          write (stdout,9101) (cifdef(n),n=ns,ne)
300     continue
        write (stdout,'(a21,/)') ' end of options list.'
#if defined cpldleg3 || defined cpldleg2
      endif
#endif
      endif
c
c-----------------------------------------------------------------------
c     Values set in parameter statements will be stored in local
c     variable names, to allow them to be written out as namelist.
c-----------------------------------------------------------------------
c
      if (wrdoc) then
        ndoc = ndocum
        write (nunit,ndocs)
        write (nunit,docums)
      endif
c
      if (wrc0) then
        write (nunit,cntrl0)
      endif
c
      if (wrc1) then
        id = imt
        jd = jmt
        kd = km
        ntrace = nt
        lsegs = lseg
        nisl = nisle
        ifdmx = ifdmax
        nhorr = nhreg
        nverr = nvreg
        ntreg = numreg
#ifdef fourfil
        lsegsf = lsegf
#endif
#ifdef inland
        ninld = ninlnd
#endif
        write (nunit,cntrl1)
      endif
c
      if (wrc2) then
        write (nunit,cntrl2)
      endif
c
      return
c
9001  format(/' The following "ifdef" options are available:'
     $  /' =========================================================='
     $ //'     >> optimization options <<'
     $ //' multitasking=allows parallelism at the latitude row level'
     $  /'               (for cray only ...not fully optimized yet!)'
     $  /'               (answers check but no ssread/sswrite yet)'
     $  /' skipland   = do calculations only over ocean segments'
     $  /'              (otherwise ... do from 2 ..imt-1)'
     $  /' timing     = gives cp & wall clock for calls to step'
     $  /'               (for cray only)'
     $ //'     >> grid options <<'
     $ //' cyclic     = east-west cyclic conditions are set'
     $  /'              (otherwise assume solid walls in east & west)'
     $  /' symmetry   = symmetry conditions set at northern boundary'
     $  /'              (otherwise assume a solid wall)'
     $ //'     >> filtering options <<'
     $ //' firfil     = use finite impulse response filter at high lats'
     $  /' fourfil    = perform Fourier filtering at high latitudes'
     $   )
9002  format(/'     >> external mode velocity options <<'
     $ //' freesfc    = free surface model (no relaxation needed)'
     $  /'               (not implemented yet!)'
     $  /' rigidlid   = use relaxation method for barotropic mode'
     $  /'  sfctopog  = change selfdoc history tape to include info'
     $  /'              for off-line calc of local sea level variations'
     $  /'              due ocean circulation (need selfdoc option)'
     $  /'  congrad5pt= preconditioned conjugate gradient scheme'
     $  /'                (recommended over oldrelax & hypergrid)'
     $  /'  congrad9pt= preconditioned conjugate gradient scheme'
     $  /'                (the most accurate poisson solver.)'
     $  /'                (uses a 9pt laplacian. Allows for a more)'
     $  /'                (accurate surface pressure reconstruction)'
     $  /'  hypergrid = successive over-relaxation on diagonals'
     $  /'  islands   = include islands in relaxation topography'
     $  /'  neptune   = represent topographic stress via Neptune Effect'
     $  /'  oldrelax  = original sequential relaxation'
     $   )
9003  format(/'     >> vertical mixing scheme options <<'
     $ //' constvmix  = constant vertical mixing'
     $  /' gargettmix = Gargett vertical mixing'
     $  /'               (not implemented yet!)'
     $  /' implicitvmix= implicit treatment of vertical mixing'
     $  /'               (otherwise assume explicit treatment)'
     $  /' ppvmix     = Pacanowski-Philander vertical mixing'
     $  /'               (not fully checked out yet. user beware !)'
     $  /' tcvmix     = turbulence closure vertical mixing'
     $  /'               (not fully checked out yet. user beware !)'
     $  /'   lalg     = macro length scale from algebraic relation'
     $  /'               (not fully checked out yet. user beware !)'
     $  /'   leq      = macro length scale from differential eqn'
     $  /'               (not fully checked out yet. user beware !)'
     $  /' windvmix   = modified Krauss-Turner wind-forced mixing'
     $  /'               (not implemented yet!)'
     $   )
9004  format(/'     >> horizontal mixing scheme options <<'
     $ //' biharmonic = bi-harmonic horizontal mixing del**4'
     $  /'               (not fully checked out yet. user beware !)'
     $  /' consthmix  = constant horizontal mixing del**2'
     $  /' isopycmix  = orient tracer diffusion along isopycnals'
     $  /' nlhmix     = non-linear horizontal mixing'
     $  /'               (not fully checked out yet. user beware !)'
     $   )
9005  format(/'     >> i/o options <<'
     $ //' crayio     = getwa & putwa i/o for now'
     $  /' diskless   = simulate disk storage using memory array'
     $  /' fio        = fortran direct access i/o'
     $  /' f90type    = use ogeti and oputi i/o for integers'
     $  /' twoslabdisks = save disk space used by slabs by having'
     $  /'                2 rather than 3 disk units (tau+1 and tau-1'
     $  /'              timelevels share one disk unit) '
     $  /' t90_1st    = read restart data on the Cray T90 that was'
     $  /'              written on the Cray C90'
     $   )
9006  format(/'     >> miscellaneous options <<'
     $ //' kdtermb    = write regional tracer term balance info for '
     $  /'              use in graphical analysis program (kd)'
     $  /' kdtravg    = write regional tracer average info for   '
     $  /'              use in graphical analysis program (kd)'
     $  /' keepterms  = use arrays (instead of statement functions) to'
     $  /'              retain terms in momentum and tracer equations'
     $  /' nohilats   = no advective metric terms in the momentum'
     $  /'              equations.(ie: tropical model domains)'
     $  /' readrmsk   = read in region info from file "regionmasks"'
     $  /' restorst   = restore surface tracers to prescribed values'
     $  /'              via newtonain damping'
     $  /' testcfl    = monitors velocities & checks for cfl violations'
     $  /'    cflj    = print u, v, w for all k at row j for cfl'
     $  /'              violations'
     $  /' testnorth  = compute basin northward transports for   '
     $  /'              globe, Atlantic&Arctic, Indo-Pacific     '
     $  /' =========================================================='
     $   )
9007  format(/'     >>  options developed for GFDL world ocean runs  <<'
     $ //' worldocn   = make available options developed for GFDL'
     $  /'              world ocean runs'
     $  /'  agetracers= carry an idealized age tracer'
     $  /'  amazon    = mix amazon river outflow in 3 dimenstions'
     $  /'  bctest    = on first timestep, compute zonal avg of sfc bc'
     $  /'  blhmix    = make horizontal diffusion coeff vary with'
     $  /'              depth using the relationship of Bryan & Lewis'
     $  /'              (does not affect viscosity and must be used'
     $  /'              with constvmix or isopycmix defined)'
     $  /'  blvmix    = make vertical diffusion coefficient vary with'
     $  /'              depth using the relationship of Bryan & Lewis'
     $  /'              (does not affect vertical viscosity and must'
     $  /'              be used with constvmix defined)'
     $  /'  bumpkmu   = modify kmu values in filter areas'
     $  /'  ccl4tracer= carry a carbon tetrachloride tracer'
     $  /'  cfctracers= carry cfc-11 and cfc-12 as tracers'
     $  /'   cjbcfcs  = use atmospheric time histories of CJ Beegle'
     $  /'   gptcfcs  = use atmospheric time histories of GP Thiele'
     $  /'  cpldleg2  = code changes for leg2 of coupled model runs'
     $  /'  cpldleg3  = code changes for leg3 of coupled model runs'
     $  /'  cyberAH   = set AHisop as in cyber code'
     $  /'  cyberfil  = set filtering lats to match previous Cyber runs'
     $  /'  cybertra  = specify tracera as in cyber code'
     $  /'  cybervbc  = use surface boundary conditions from Cyber runs'
     $  /'  gtime     = use KDs time manager developed for G-group'
     $  /'              (see 5 key assumptions listed in gtmngr.F)'
     $  /'   gotogtime= read in old non-gtime restart file but run'   
     $  /'              using gtime and write out gtime type restart'
     $  /'   cpldinit = resets time indicators used with gtime option'
     $  /'              to initial conditions when starting a coupled'
     $  /'              model from a pre-existing restart file'
     $  /'  inland    = enable communication with inland seas'
     $  /'              NOTE: see tracer.F for variations'
     $  /'  january   = used for perpetual january test case'
     $   )
9008  format(
     $   '  kdtest    = used for various tests of ocn only code'
     $  /'  km12      = use 12-level depths from previous cyber runs'
     $  /'  km18      = use 18-level depths'
     $  /'   old_dzt  = dzt for km18 levels, dzt(1)=32m'
     $  /'  km25      = use 25-level depths'
     $  /'  noAH      = set coefficients to zero for isopycmix'
     $  /'  noadvection= ignore advection of tracers'
     $  /'  noclinic  = set surface momentum fluxes to zero'
     $  /'  noramp    = ignore ramp for surf boundary conditions'
     $  /'  notracer  = set surface tracer fluxes to zero'
     $  /'  ocnr15    = match ocean grid with rhomboidal 15 spectral'
     $  /'              grid for atmospheric climate model 2:1'
     $  /'  ocnr30    = match ocean grid with rhomboidal 30 spectral'
     $  /'              grid for atmospheric climate model 2:1'
     $  /'  oldcongrad= congrad from dpr version 1 (no seed)'
     $  /'  old_dtxcel= timestep multipliers for max dtxcel=10'
     $  /'  oldstate  = replicate eq of state of Cyber runs'
     $  /'  oldtime   = replicate calendar cycle of Cyber runs'
     $  /'  onetoone  = set dtxcel variable to 1.0 all all k-levels'
     $  /'              (default for worldocn runs accelerate to 10x)'
     $  /'  r15vbc    = use bc dixon made for YMP runs'
     $  /'  saltnorm  = enables salt normalization; does not'
     $  /'              work when dtxcel varies with depth'
     $  /'  sameAH    = set AHisop = bkgrd AH at all k-levels'
     $  /'  sameFKPH  = set constant vertical tracer mixing'
     $  /'  seaice    = include computation of sea ice'
     $  /'  seasonal  = include seasonal cycle by cycling thru monthly'
     $  /'              surface boundary conditions (default, no cycle)'
     $  /'  selfdoc   = self documenting history data for analysis '
     $  /'  testrelx  = store and print info for relaxation scheme'
     $  /'  windbc    = include wind and wind**3 surface boundary'
     $  /'              conditions from atmos model for wind mixing'
     $  /'  writekmt  = write the kmt array on stdout in topog.F'
     $  /'  yearzero  = set year0 = 0 (not 1) to agree with prev runs'
     $  /'  zonalavg  = compute and printout tracers and u, v & w as'
     $  /'              zonal averages rather than by row'
     $  /' =========================================================='
     $   )
9101  format(5(3x,a12))
c
      end
