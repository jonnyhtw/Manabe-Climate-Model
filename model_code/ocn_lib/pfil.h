c====================== include file "pfil.h" ==========================
c
c     define latitudinal domain for filtering
c
c     filter t to yield equiv dx at lat rjft0 from lat rjfrst to rjft1
c       and lat rjft2 to yt(jmtm1)
c     filter u to yield equiv dx at lat rjfu0 from lat rjfrst to rjfu1
c       and lat rjfu2 to yt(jmtm2)
c
c     lsegf  = max number of longitudinal strips for filtering
c     jmtfil = max number of latitudes to be filtered
c     numflt = specifies the number of filter applications to use
c              for the tracers when using the finite impulse
c              response filter
c     numflu = specifies the number of filter applications to use
c              for velocities and vorticity when using the finite 
c              impulse response filter
c
#ifdef worldocn
c     parameters for GFDL world ocean models
c
# ifdef ocnr15
c     parameters matching those of previous r15 cyber runs
c
      parameter (lsegf=8,  jmtfil=20)
# endif
# ifdef ocnr30
c
      parameter (lsegf=8, jmtfil=22)
# endif
#else
c     parameters for MOM sample case
c
      parameter (lsegf=5, jmtfil=20)
#endif
c
      common /cfilr/ rjfrst, rjft0, rjft1, rjft2, rjfu0, rjfu1, rjfu2
      common /cfil/   jfrst,  jft0,  jft1,  jft2,  jfu0,  jfu1,  jfu2
     $,               jskpt, jskpu, njtbft, njtbfu
#ifdef firfil
     $,               numflt(jmtfil), numflu(jmtfil)
#endif
c
