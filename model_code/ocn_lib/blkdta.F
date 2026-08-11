      block data blkdta
      implicit real*8 (a-h,o-z)
#include "param.h"
#include "accel.h"
#if defined worldocn && defined ocnr30
# ifdef amazon
# include "camazon.h"
# endif
#endif
cSGI
#if defined worldocn
# ifdef testnorth
# include "cdiag.h"
# endif
#endif
#include "cstate.h"
cSGI
#ifdef testcfl
# include "ccfl.h"
#endif
#include "chmix.h"
#if defined cpldleg3 || defined cpldleg2
# ifdef inland
# include "cinlnd.h"
# endif
#endif
#ifdef isopycmix
#include "cisop.h"
#endif
#include "coord.h"
#include "cprnts.h"
#ifdef ppvmix
#include "cppmix.h"
#endif
#ifdef tcvmix
#include "ctcmix.h"
#endif
#include "ctask.h"
#include "ctmngr.h"
#include "cvmix.h"
#include "grdvar.h"
#include "iounit.h"
#if defined cpldleg3 || defined cpldleg2
# include "histjr.h"
# include "index.h"
#else
#if defined firfil || defined fourfil
#include "pfil.h"
#endif
#endif
#include "scalar.h"
#ifdef seaice
#include "seayce.h"
#endif
#ifdef restorst
#include "slabs.h"
#endif
#include "switch.h"
#include "versno.h"
c
c     specify current MOM version info
c
      data model/
     $' GFDL Modular Ocean Model (MOM) version 1.0 ... Dec 1, 1990  '/
c
c     initialize time at initial conditions. see "ctmngr.h"
c
#ifdef worldocn
# ifdef yearzero
c   set year0 to conform to previous runs
      data year0  /0/
# else
      data year0  /1/
# endif
#else
      data year0  /0/
#endif
#ifdef gtime
      data nvar11 /4/
#else
      data nvar11 /2/
c
      data month0 /1/
      data day0   /1/
      data hour0  /0/
      data minut0 /0/
      data sec0   /0/
#endif
c
#if defined cpldleg3 || defined cpldleg2
      data months / 'jan.','feb.','mar.','apr.','may ','june',
     $              'july','aug.','sep.','oct.','nov.','dec.'/
#endif
c
#ifndef gtime
c     specify julian calendar
c
      data julian /.true./
#endif
c
c     initialize all quantities that define the grid. see "coord.h"
c
#ifdef worldocn
c     for models to match grids of climate groups spectral models,
c     variables stlat & stlon are set in grids.F
c
#else
      data stlon, stlat / -4.0, -90.0 /
#endif
c
c     specifications for "nxres" distinct longitudinal domains (deg)
c
#ifdef worldocn
c     for r15 & r30 cases, variables xmax, xmin, xwid & idir
c     are set at the beginning of grids.F
c     for r15 & r30 cases, variables ymax, ymin, ywid & jdir
c     are not used, due to gaussian latitudes
c
#else
      data xmax /4.0/
      data xmin /4.0/
      data xwid /368.0/
      data idir /1/
c
c     specifications for "nyres" distinct latitudinal domains (deg)
c
      data ymax /3.0/
      data ymin /3.0/
      data ywid /180.0/
      data jdir /1/
#endif
c
#include "thick.h"
c
c     initialize some switches. see "switch.h"
c
      data init /.false./, first/.true./, tsi /1.0/, dgnstc/365.25/
      data travg /365.25/, snaps/365.25/
#if defined cpldleg3 || defined cpldleg2
      data pintro /.true./
#endif
#ifdef worldocn
      data hist0 /.true./, histt/1.0e+6/

cSGI
# ifdef testnorth
      data (nornam(ns),ns=0,2)
     $     /'  Global Total  ', ' Atlantic/Arctic', '  Indo/Pacific  '/

# endif
#endif
c
c     initialize some vertical mixing coefficients. see "cvmix.h"
c
      data  bvdc, bvvc / 0.00134, 0.0134/
      data fkph, fkpm / 20.0, 1.0/
c
#ifdef constvmix
      data vdcsci / imt*1.0 /
      data vdcscj / jmt*1.0 /
      data vdcsck / km *1.0 /
      data vvcsci / imt*1.0 /
      data vvcscj / jmt*1.0 /
      data vvcsck / km *1.0 /
# ifdef blvmix
ccccc Ahv range of 0.15 to 1.3, crossover at 2500m
ccccc data afkph, dfkph, sfkph, zfkph /0.725, 1.22, 4.5e-5, 2500.0e2/
c     Brian & Lewis values,
c      with Ahv range of 0.3 to 1.3, crossover at 2500m.
ctd 1/7/97      data afkph, dfkph, sfkph, zfkph /0.8, 1.05, 4.5e-5, 2500.0e2/
       data afkph, dfkph, sfkph, zfkph /0.725, 1.22, 4.5e-5, 2500.0e2/
ctd 1/7/97
# endif
#endif
c
c     initialize some horizontal mixing coefficients. see "chmix.h"
c
      data am, ah /1.e9, 2.e7/
#if defined consthmix || defined isopycmix
      data hdcsck / km*1.0 /
      data hvcsck / km*1.0 /
#endif
c
c     set bottom drag coefficient. see "scalar.h"
c
      data cdbot /0.0/
#ifdef testcfl
c
c     factor by which the CFL criteria must be exceeded in order to
c     stop MOM
c
      data cflcrt /1.5/
#endif
#if defined implicitvmix || defined isopycmix
c
c     default to fully implicit vertical mixing. zero gives explicit
c
      data aidif /1.0/
#else
      data aidif /0.0/
#endif
#ifdef multitasking
c
c     initialize the number of tasks for cycling slabs. see "ctask.h"
c     "ntasks" should be set to the number of processors. default is
c     for an 8 processor cray ymp.
c
      data ntasks /8/
#else
c
c     uni-tasking case. one task is the default. "ntasks" can be set
c     greater than 1 but this leads to redundant calculations at the
c     interface between tasks.
c
      data ntasks /1/
#endif
c
c     set type of time mixing. see "switch.h"
c
      data eb /.true./
c
c     set number of time steps between time mixing. see "switch.h"
c
      data nmix /17/
c
c     set number of convective passes in "tracer.F". see "scalar.h"
c
#if defined cpldleg3 || defined cpldleg2
# ifdef km12
      data ncon /4/
# endif
# ifdef km18
      data ncon /7/
# endif
#else
      data ncon /1/
#endif
c
c     default starting & ending coordinates for matrix prints during
c     diagnostic time steps (lats & lons in degrees, depths in cm)
c     see "cprnts.h"
c
#ifndef zonalavg
      data prlat /-60.0, 0.0, 27.0, 55.0, 6*100.0/
#endif
#ifdef zonalavg
      data prlat /nlatpr*100.0/
#endif
      data prslon /180.0, 180.0, 180.0, 180.0, 6*0.0/
      data prelon /250.0, 250.0, 250.0, 250.0, 6*0.0/
      data prsdpt /nlatpr*0.0/
      data predpt /nlatpr*6000.e2/
      data slatxy, elatxy, slonxy, elonxy /-90.0, 90.0, 3.0, 357.0/
c
c     default the i/o units. see "iounit.h"
c
c     (note labs units must be sequential for the "diskless" option)
c
      data kontrl /11/, kflds /12/, labs /13, 14, 15/
      data iorest /21/, restrt /.false./
#ifdef kdtermb
c     specify fortran i/o units used in time averaging of term balances
      data iointb /69/
      data iotrm2 /79/
#endif
# ifdef testnorth
c     unit for northward transport sector masks
      data ionrth /82/
# endif
#ifdef kdtravg
      data iormsk /55/, iotavg /stdout/, iotavg2 /87/,
     $     iotrmb/stdout/, iosnap/22/
#else
      data iormsk /55/, iotavg /stdout/, iotrmb/stdout/, iosnap/22/
#endif
#ifdef worldocn
      data iokmt /49/
# if defined cpldleg3 || defined cpldleg2
      data iobc  /35/
      data iobcw /36/
# else
      data iobc /78/
# endif
      data iohist /26/
      data intrnl /45/
      data iochet /91/
      data iodumy /90/
#endif
#ifdef sfctopog
      data iosealev /107/
#endif
#if defined cpldleg3 || defined cpldleg2
      data iobca /74/
#endif
#if defined cpldleg2 || defined cpld_restore
      data ioobs /34/
#endif
#ifdef testrelx
      data iorelx /89/
#endif
c
c     initialize physical constants. see "scalar.h"
c
      data grav /980.6/, radius /6370.e5/
c
#if defined firfil || defined fourfil
c
c     set latitudes used in filtering of tracer and velocity fields
c     see "pfil.h"
c
#if defined cpldleg3 || defined cpldleg2
# ifdef ocnr15
# ifdef cyberfil
c     set filter ranges to match previous cyber runs
c
      data rjfrst /-78.0/
      data rjft0, rjft1, rjft2 /-46.7, -51.1, 51.1/
      data rjfu0, rjfu1, rjfu2 /-62.2, -66.7, 66.7/
# else
      data rjfrst /-78.0/
      data rjft0, rjft1, rjft2 /-64.5, -68.9, 68.9/
      data rjfu0, rjfu1, rjfu2 /-66.7, -71.1, 71.1/
# endif
# endif
# ifdef ocnr30
      data rjfrst /-78.0/
      data rjft0, rjft1, rjft2 /-66.0, -68.2, 68.2/
      data rjfu0, rjfu1, rjfu2 /-67.1, -69.3, 69.3/
# endif
#else
      data rjfrst /-81.0/
      data rjft0, rjft1, rjft2 /-51.0, -57.0, 57.0/
      data rjfu0, rjfu1, rjfu2 /-54.0, -60.0, 60.0/
#endif
c
#endif
#ifdef firfil
c     "numflt" specifies the number of filter applications to use
c              for the tracers
c     "numflu" specifies the number of filter applications to use
c              for velocities and vorticity
c
      data numflt /jmtfil*1/
      data numflu /jmtfil*1/
c
#endif
c
c     "dtxcel" specifies the depth dependent dtts timestep multipliers
c      see "accel.h"
c
#if defined worldocn && !defined onetoone
# ifdef km12
c      worldocn 12-level runs defaulted to acceleration at level 6 and
c      below reaching maximum of 10.0 at level 12
c
       data dtxcel /5*1.0, 1.2, 2.3, 4.0, 7.0, 8.5, 9.4, 10.0/
# endif
# ifdef km18
#  ifdef old_dtxcel
c      worldocn 18-level runs defaulted to acceleration at level 9 and
c      below reaching maximum of 10.0 at level 18
c
       data dtxcel /8*1.0, 1.2, 1.5, 2.2, 3.4, 5.0,
     $                     7.1, 8.3, 9.1, 9.6, 10.0/
c
#  else
c      worldocn 18-level runs defaulted to acceleration at level 15 and
c      below reaching maximum of 5.0 at level 18
c
ctd 8/19/96 1/7/97       data dtxcel /14*1.0, 2.4, 4.0, 4.5, 5.0/
       data dtxcel / 10*1.000,
     $    1.146, 1.369, 1.697,  2.185, 2.921, 4.042, 2*5.000 /
ctd 8/19/96 1/7/97
c
#  endif
# endif
#else
c      no acceleration of ts timestep with depth
c
       data dtxcel /km*1.0/
#endif
c
#ifdef restorst
c
c      time scale for restoring surface tracers to prescribed values
c      see "setvbc.F"
c
       data rests /50.0/
#endif
#if defined constvmix && defined implicitvmix
      data vvclim, vdclim / 1.0e6, 1.0e6/
#endif
#ifdef ppvmix
c
c     initialize constants for pacanowski & philander vertical mixing
c     see "cppmix.h"
c
      data  wndmix / 10.0/
      data fricmx  / 50.0/
# ifdef implicitvmix
c
c     simulate convective adjustment with large mixing coefficient
c     limits
c
      data vvclim, vdclim  / 1.0e6, 1.0e6/
# else
c
c     in regions of gravitational instability set mixing limits to the
c     maximum consistent with the "cfl" criterion. convective adjustment
c     will also act on the instability.
c
      data vvclim, vdclim / 50.0, 50.0/
# endif
#endif
#ifdef isopycmix
c
c     initialize constants for isopycnal oriented mixing of tracers
c     see "cisop.h"
c
      data slmxr  / 100.0/
      data ahisop / 2.0e7/
#endif
#ifdef tcvmix
c
c     data initialize constants for turbulence closure vertical mixing
c     see "ctcmix.h"
c
      data a1, b1, a2, b2, cc1 /0.92, 16.6, 0.74, 10.1, 0.08/
      data e1, e2, e3 /1.8, 1.33, 1.0/
      data aq /2.0e7/
# ifdef implicitvmix
      data vvclim, vdclim, vdqlim /1.e6, 1.e6, 1.e6/
# else
      data vvclim, vdclim, vdqlim /1.e1, 1.e1, 1.e1/
# endif
      data alpha /.1/
      data vk /0.40/
      data sq /0.2/
      data small /1.e-6/
      data big /1.e15/
      data deps /1.e-10/
      data c100 /100.0/
      data p48 /0.48/, p032 /0.032/, c15/15.0/, p1/0.1/, p4/0.4/
      data p154 /.154/
      data gh /imtkm*0.154/
      data gmp1 /imtkm*0.154/
      data sm /imtkm*0.39/
      data sh /imtkm*0.49/
      data rhoh2o /1.025/
#endif
#ifdef inland
c
c     set locations for inland mixing scheme
c     values are extremely model dependent (see "cinlnd.h")
c
# ifdef worldocn
#  if defined ocnr15
c     set mediterranean inland mixing for ocean model matched with
c     the r15 spectral model
c
      data iinlnd /2, 95/
      data jinlnd /30/
c     data vinlnd /1.50e12/
      data vinlnd /1.75e12/
cc    data vinlnd /0.2522555e+16/
#  endif
#  if defined ocnr30
#    ifdef amazon
c
c     specify grid range over which amazon river mixing will extend
      data kamzn  /   2 /
      data iamzns / 168 /
      data iamzne / 171 /
      data jamzns /  41 /
      data jamzne /  44 /
#    endif
c     set mediterranean inland mixing
c        and Hudson Bay inland miximg for ocean model matched with
c        the r30 spectral model
c
cc    data iinlnd /2, 191, 153, 159/
      data iinlnd /2, 153, 191, 159/
      data jinlnd /58, 69/
      data vinlnd /1.75e12, 0.15e12/
#  endif
#  if defined km12
c     set mediterranean inland mixing for 12-level model
      data kinlnd / 1, 5/
#  endif
#  if defined km18
c     set mediterranean inland mixing
c        and Hudson Bay inland miximg for 18-level model
cc    data kinlnd / 1, 8, 1, 7/
      data kinlnd / 1, 1, 8, 7/
#  endif
# else
      data iinlnd /ninlnd*1,ninlnd*1/
      data jinlnd /ninlnd*1/
      data kinlnd /ninlnd*1,ninlnd*1/
      data vinlnd /ninlnd*0.0/
# endif
#endif
#ifdef worldocn
# ifdef gtime
c     set the number of days per month
c
      data dayspermon /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c                       J   F   M   A   M   J   J   A   S   O   N   D
c
# else
c     set the number of days per month in "blkdta" rather than "tmngr"
c     (similar to an improvement to be made in the next version of MOM)
c
      data dpm /31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0,
     $          31.0, 30.0, 31.0/
c
# endif
#endif
#ifdef seaice
c
c     ycejam = sea ice jam criterion [cm], above this value advection
c              of sea ice may be set to zero
c     ahyc   = horizontal sea ice diffusion coefficient [cm**2/sec]
c     freeze = temperature [deg C] at which water freezes into sea ice
c
      data ycejam /400.0/
      data ahyc   /1.0e7/
      data altht  /80.0 /
      data freeze /-2.0 /
#endif
#ifdef islands
# ifdef ocnr30
c
      data alonis / 20.0, 172.0, 135.0,  47.0, 280.0, 344.0, 18.0/
      data alatis /-75.0, -43.0, -25.0, -18.0,  21.0,  65.5, 78.0/
# endif
#endif
c
cSGI moved from state.F
#include "dncoef.h"
cSGI
      end
