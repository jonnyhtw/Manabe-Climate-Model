c====================== include file "switch.h" ========================
c
c     all time dependent decisions are made in time manager tmngr.F
c     and communicated elsewhere in the model via logical switches .
c
c-----------------------------------------------------------------------
c     Inputs: (set through namelist)
c-----------------------------------------------------------------------
c
c  -------------------------
c  Floating Point Inputs:
c  -------------------------
c     dgnstc  = number of days between diagnostic calculations:
c               global energetics & regional momentum & tracer term
c               balances
c     travg   = number of days between regional tracer averages
c               (under horizontal regions)
#ifdef worldocn
c     histt   = number of days between history writes
c
c     tnt     =  number of days between calculation of tracer northward
c                transport.
c
#endif
c
c     tsi     = number of days between printing of time step info
c     days    = number of days to integrate
c
c     snaps   = number of days between taking instantaneous snapshots
c               (saving data) for off line analysis.
c               symbolically ... the data saved is:
c
c                 write (iosnap) itt, dtts, relyr, stamp
c                 write (iosnap) xt, yt, zt
c                 do ... jrow=2,jmtm1
c                   write (iosnap) smf, stf, slab(jrow)
c                 continue
c                 write (iosnap) stream function
#ifdef kdtermb
c     avterm  = number of days between time averaging term balance
c               diagnostics and writing them out (instantaneous values
c               that go into the average are "dgnstc" days apart)
#endif
c
c
c     acor    = (>0, 0) = (implicit, explicit) treatment of coriolis
c               term
c     eb      = (false,true) implies (forward, euler backward) mixing
c               time step
c
#if defined cpldleg3 || defined cpldleg2
c     nsteps  = integrate the ocean model "nsteps" timesteps before ending.
c               In a scheme in which the coupled model executable is run in
c               one month increments, and the ocean-atmosphere coupling is done
c               once per day, (the ocean model integrates 1 model day per call),
c               "nsteps" would equal the number of ts timesteps in 1 day.
c               In ocean-only models, the ocean model is called only once per
c               run, with each run being of 1 month duration. Variable "nsteps"
c               then is set to a very large number (9.e9) and the variable
c               "days" is used to end the integration after 1 month.
c     pintro  = (false,true) = (do not, do) execute the write statements in the
c               introductory section of subroutine ocean.
c     prdiag  = number of days between printing of diagnostic info;
c               must be a multiple of dgnstc.
#endif
c  -------------------------
c  Integer Inputs:
c  -------------------------
c     nmix    = number of time steps between time step mixing
c               to damp leap frog time splitting
c
c  -------------------------
c  Logical Inputs:
c  -------------------------
c     init    = (false,true) to indicate if this run is a
c               (restart, start from initial conditions)
c     restrt  = (false,true) = (do not, do) write a restart at the end
c               of this run
c
c-----------------------------------------------------------------------
c     outputs: (set in tmngr.F)
c-----------------------------------------------------------------------
c
c  -------------------------
c  Logical Outputs:
c  -------------------------
c     the following logical switches are set within tmngr every
c     time step based on the above requirements.
c
c     diagts  = (false,true) = (do not, do) diagnostic calculation
c                                          based on dgnstc
c     tavgts  = (false,true) = (do not, do) tracer averages
c                                          based on travg
c     snapts  = (false,true) = (do not, do) save a snapshot of the data
c                                          based on snaps
c     mixts   = (false,true) = (do not, do) mixing time step
c                                          based on nmix
c     mxpas2  = (false,true) = (do not, do) 2nd pass of a mixing
c                                          time step (if eb = true.)
c                                          (set in ocean.F)
c     eots    = end of a time step. always true except for first
c               pass of an euler backward time step
c
c     prntsi  = (false,true) = (do not ,do) print time step info
c                                          based on ntsi
c     last    = (false,true) = (do not, do) stop the integration when
c                                          this time step completes
c                                          based on days
c     first   = (false,true) =  when it is (not the first, the first)
c                               time step of a run
c     eom     = true on the 1st timestep of a new month,otherwise false
#ifdef gtime
c     eoday    = true when within one-half timestep of the end of a day,
c                otherwise false (end of the month)
c     eomonth  = true when within a half timestep of the end of a month
c     eoyear   = true when within a half timestep of the end of a year
c     eodecade = true when within a half timestep of the end of a decade
c     loday    = true when it is the last timestep of a day, otherwise
c                false  (last of the day)
c     lomonth  = true when it is the last timestep of a month
c     loyear   = true when it is the last timestep of a year
c     lodecade = true when it is the last timestep of a decade
c
#endif
#ifdef kdtermb
c     termav  = (false,true) = (do not, do) time average term balance
c                              diagnostics and write them out
c                              based on avterm
# if defined cpldleg3 || defined cpldleg2
c     prterm  = (false,true) = (do not, do) print the Regional Term Balances for
c                              Momentum terms and Tracer terms for each region
c                              based on (mod(month,6).eq.1).
# endif
#endif
#ifdef worldocn
c     histts  = (false,true) = (do not, do ) write history data for
c                                          analysis: based on histt
c     hist0   = (false,true) =  when it is (not the first, the first)
c                               write of history data of a run
c     tntts   = (false,true) = (do not, do) calculation of tracer
c               northward transport on this time step. based on tnt.
c
#endif
#if defined cpldleg3 || defined cpldleg2
c     endmon  = true on the last time step of the month
c               otherwise ... false
c     endseg  = true on the last time step of the segment (=1 day)
c               otherwise ... false
c     prdgts  = (false,true) = (do not, do) print diagnostic info
c                               on a diagts time step. based on prdiag
#endif
#ifdef cpldleg3
c     loop1   = true on the  first timestep of the day
c               otherwise ... false
#endif
c
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c
      logical diagts, mixts, mxpas2, prntsi
      logical eb, last, init, first
      logical restrt, snapts, tavgts, eom, eots
#ifdef gtime
      logical eoday, eomonth, eoyear, eodecade
      logical loday, lomonth, loyear, lodecade
#endif
#ifdef worldocn
      logical histts, hist0, tntts
#endif
#if defined cpldleg3 || defined cpldleg2
      logical endmon, endseg
      logical pintro, prdgts
#endif
#ifdef cpldleg3
      logical loop1
#endif
#ifdef kdtermb
      logical termav
# if defined cpldleg3 || defined cpldleg2
      logical prterm
# endif
#endif
c
c-----------------------------------------------------------------------
c   put logical variables in switchl common block
c-----------------------------------------------------------------------
c
      common /switchl/ diagts, mixts, mxpas2, prntsi
      common /switchl/ eb, last, init, first
      common /switchl/ restrt, snapts, tavgts, eom, eots
#ifdef gtime
      common /switchl/ eoday, eomonth, eoyear, eodecade
      common /switchl/ loday, lomonth, loyear, lodecade
#endif
#ifdef worldocn
      common /switchl/ histts, hist0, tntts
#endif
#if defined cpldleg3 || defined cpldleg2
      common /switchl/ endmon, endseg
      common /switchl/ pintro, prdgts
#endif
#ifdef cpldleg3
      common /switchl/ loop1
#endif
#ifdef kdtermb
      common /switchl/ termav
# if defined cpldleg3 || defined cpldleg2
      common /switchl/ prterm
# endif
#endif
c
c-----------------------------------------------------------------------
c   put integer variables in switchi common block
c-----------------------------------------------------------------------
c
      common /switchi/ nmix
#if defined cpldleg3 || defined cpldleg2
      common /switchi/ nsteps
#endif
c
c-----------------------------------------------------------------------
c   put floating point variables in switchr common block
c-----------------------------------------------------------------------
c
      common /switchr/ dgnstc, tsi, days, travg, snaps, acor
#ifdef kdtermb
      common /switchr/ avterm
#endif
#ifdef worldocn
      common /switchr/ histt, tnt
#endif
#if defined cpldleg3 || defined cpldleg2
      common /switchr/ prdiag
#endif
c
