c====================== include file "iounit.h" ========================
c
c     i/o units & related variables
c
c     ndiskb  = disk pointer (index) used for tau-1 variables
c               = (1,2,3) based on itt
c     ndisk   = disk pointer (index) used for tau   variables
c               = (1,2,3) based on itt
c     ndiska  = disk pointer (index) used for tau+1 variables
c               = (1,2,3) based on itt
c     kontrl  = disk unit used for control block
c     kflds   = disk unit used for two dimensional fields
c     labs    = disk units for slabs (one for each time level)
c     iorest  = unit for history file (restarts)
c     iotavg  = unit for tracer averages
#ifdef kdtermb
c     iointb = internal unit for instantaneous term balances to be
c              written before time averaging
c     iotrm2 = external unit for writing time averaged term balances
c              (output to be used in kd's graphical analysis program)
#endif
#ifdef testnorth
c     ionrth = unit for reading in masks ("normsk") to divide
c              the ocean into Atlantic-Arctic and Indo-Pacific sectors
c              for meridional transport calculations
#endif
#ifdef kdtravg
c     iotavg2 = unit for unformatted tracer averages used in kd's
c               graphical analysis program
#endif
c     iormsk  = unit for regional masks (horizontal & vertical)
c     iotrmb  = unit for tracer & momentum term balances
c     iosnap  = unit for writing instantaneous snapshots of variables
#ifdef worldocn
c
c     i/o units used in GFDL world ocean models
c
c     iokmt   = unit for reading model depths
c     iobc    = unit for reading temperature, salinity, and wind stress
c               boundary conditions
c     iohist  = history file unit number used by keith et al, IEEE data
c               in 'smmaps' input form.
c     intrnl  = internal unit for history data arrangement
c     iochet  = unit for cheat sheet list output
c     iodumy  = dummy internal unit for cheat sheet list arrangement
#endif
#ifdef sfctopog
c     iosealev = unit for 2-D sealev restart field
#endif
#if defined cpldleg3 || defined cpldleg2
c     iobca   = unit for reading temperature and salinity surface flux      
c               correction terms, terr and serr    
c     iobcw   = unit for writing temperature and sea ice
c               boundary conditions for the atmosphere model
#endif
#if defined cpldleg2 || defined cpld_restore
c     ioobs   = unit for reading observed temperature, salinity, and
c               sea ice for computing terr and serr flux adjustments
#endif
#ifdef testrelx
c     iorelx  = unit for storing info about scans for the relaxation
c                method
#endif
c
      common /iounit/ ndiskb, ndisk, ndiska, kontrl, kflds, labs(3)
     $,               iorest, iotavg, iormsk, iotrmb, iosnap
#ifdef kdtermb
     $,               iointb, iotrm2
#endif
#ifdef testnorth
     $,               ionrth
#endif
#ifdef kdtravg
     $,               iotavg2
#endif
#ifdef sfctopog
     $,               iosealev
#endif
#ifdef worldocn
     $,               iohist, intrnl, iochet, iodumy
     $,               iokmt, iobc
# if defined cpldleg3 || defined cpldleg2
     $,               iobca, iobcw  
# endif
# if defined cpldleg2 || defined cpld_restore
     $,               ioobs       
# endif
#endif
#ifdef testrelx
     $,               iorelx
#endif
c
c     nvar11 = number of variables being written to unit 11 kontrol
      common /iouniti/ nvar11
c
