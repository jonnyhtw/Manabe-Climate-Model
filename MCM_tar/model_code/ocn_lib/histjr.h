c====================== include file "histjr.h" ========================
#ifdef selfdoc
c     sets parameters used to determine the number of history variables
c     to be written to unit "iohist" and contains common blocks which
c     contain info specific to the selfdoc option.
c
c     NOTE: modifications of this routine should be made only after
c           reviewing the way in which it works, to ensure that
c           understanding that the modified code will work for the
c           various pre-exiting options.
c
c     NOTE: this code assumes that the ifdef seaice option is selected
c           only during leg2 or leg3 of an air sea project
c
c-----------------------------------------------------------------------
c
      parameter ( iblk0=(2+km*(nt+3)+2+nt+2+km+1) )
c
c
# ifdef cfctracers
      parameter ( ibcfcs = 14 )
# else
      parameter ( ibcfcs = 0 )
# endif
# ifdef ccl4tracer
      parameter ( ibccl4 = 7 )
# else
      parameter ( ibccl4 = 0 )
# endif
# if defined cfctracers || defined  ccl4tracer
#  if defined cpldleg2 || defined cpldleg3
      parameter ( iblgxc = 1 )
#  else
      parameter ( iblgxc = 2 )
#  endif
# else
      parameter ( iblgxc = 0 )
# endif
c
      parameter ( iblktr = ibcfcs + ibccl4 + iblgxc )
c
c
# if defined cpldleg2 || defined cpldleg3
# ifdef seaice
      parameter ( iblegs = 6 )
# else
      parameter ( iblegs = 2 )
# endif
# else
      parameter ( iblegs = 0 )
# endif
# if defined cpldleg2 || defined cpld_restore
      parameter ( ibleg2 = 8
# ifdef seasonal
     $                   + 2
# endif
     $           )
# else
      parameter ( ibleg2 = 0 )
# endif
# if defined cpldleg3 && !defined cpld_restore
      parameter ( ibleg3 = 0 )
# else
      parameter ( ibleg3 = 0 )
# endif
c
      parameter ( iblkas = iblegs + ibleg2 + ibleg3 )
c
c
      parameter ( iblkxt = 0
#ifdef sfctopog
     $                   + 7
#endif
     $           )
c
c
      parameter ( iblk = iblk0 + iblktr + iblkas + iblkxt )
c
c-----------------------------------------------------------------------
c
c     explanation of "iblk" value and the pointwise data arrangement:
c     NOTE: modifications of this routine should be made only after
c           reviewing the way in which it works, to ensure that
c           understanding that the modified code will work for the
c           various pre-exiting options.
c
c     "iblk" is the total number of history variables associated with
c            each i,j point.
c
c     "iblk0" is the minimum number of variables to be written out in
c            a plain vanilla run.  additional variables are added
c            depending upon the ifdef options selected, and are
c            accounted for by iblktr, iblkas and iblkxt, below.
c
c     "iblktr" is the number of history  variables corresponding to
c             gaseous ocean tracers being carried in the model (set in
c             response to the "cfctracers" and "ccl4tracer" options)
c     "ibccl4" & "ibcfcs" are the number of additional boundary
c             condition fields written out for carbon-tet and the
c             chlorofluorocarbons, respectively
c     "iblgxc" is the number of history variables added to handle the
c             variables in surface gas exchange rate calculations
c
c     "iblkas" is the number of history variables to be written out that
c             correspond to GFDL air-sea coupled model project runs.
c             (either the "leg2" potion of the run or the full "leg3"
c             portion)
c     "iblegs" is the number of variables added to the history file that
c             are common to both the leg2 or leg3 portions of air-sea
c             runs, such as some sea ice info.
c     "ibleg2" & "ibleg3" are the number of additional variables written
c             out specifically for either the leg2 or leg3 portions of
c             the air sea runs, respectively.
c
c     "iblkxt" is the number of extra history variables to be written
c             out that do not fall under the iblk0, iblktr, or iblkas
c             headings.  it can also be used for newly created variables
c
c-----------------------------------------------------------------------
c
c     the iblk0 variables:
c     * the first two locations are always given to the number of ocean
c         levels at that gridpoint on the ts & uv grids (.le. "km")
c     * the next km*(nt+3) positions are filled with the tracers,
c         u & v (total velocity) & w (from tracer) at each of the "km"
c         levels (including land points).
c     * the next 2 positions are filled with the external mode u & v
c     * the next nt positions contain the surface tracer fluxes for
c         each of the "nt" tracers
c     * the next nt positions contain the surface tracer fluxes for
c         each of the "nt" tracers
c     * the next 2 positions contain the surface momentum fluxes
c         otherwise refered to as the x and y components of the wind
c         stress field
c     * the next km positions contain the ocean convection flags
c     * the very last variable position of the base portion (the iblk0
c         position) of the data vector is given to the barotropic
c         transport streamfunction
c
c-----------------------------------------------------------------------
c
c     "varnam" & "varunt" contain the names and units of measure of the
c          variables written to the history file
c     "varhgr" & "varvgr" contain the 2 character flags identifying
c          the horizontal and vertical grids on which the variables lie
c          horizontal - "ts" for tracer grid ("xt" & "yt")
c                       "uv" for velocity grid ("xv" & "yv")
c          vertical   - "zt" for midpoints of model layers ("zt")
c                       "zw" for layer interfaces ("zw")
c     "klevar" is the k-level appropriate for plotting the data
c          (use zw(0) for air-sea interface)
c     "nvchk" is used to check that the variable number being set for
c          the cheat sheet matches the data position in the pointwise
c          arrangment of history data within the history data slab
c
      character * 32 stampb, stampe
      character * 16 varnam, varunt
      character * 4 cmonth
      character * 2 varhgr, varvgr
c
c
      common /histrw/ knthst
      common /histrc/ stampb, stampe, cmonth
      common /histrc/ varnam(iblk), varunt(iblk)
      common /histrc/ varhgr(iblk), varvgr(iblk)
      common /histri/ klevar(iblk)
      common /histrw/ nvchk
c
      logical lcheat1st
      common /histrl/ lcheat1st
#endif
c
c    jhists = starting j-row for data to be written to history file
c             for off-line analysis
c    jhiste = ending j-row for data to be written to history file
c             for off-line analysis
c
      common /histrw/  jhists, jhiste
c
c-----------------------------------------------------------------------
