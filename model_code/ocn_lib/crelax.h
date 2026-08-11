c====================== include file "crelax.h" ========================
c
c     variables used in relaxation. see "relax.F" & "hyper.F"
c
c     inputs:  ( set in namelist)

c     mxscan  = max number of allowable scans
c     crit    = criterian used to stop relaxation
c     sor     = sucessive over-relaxation constant
c
c     output:  (from relax or hyper)
c
c     mscan   = actual number of scans taken
c
#ifdef testrelx
      parameter (mxwres = 500)
#endif
      common /crelax/  mxscan, mscan, crit, sor
#ifdef testrelx
     $,                resmax, crtp, ires, jres, wresmx(mxwres)
#endif
c
