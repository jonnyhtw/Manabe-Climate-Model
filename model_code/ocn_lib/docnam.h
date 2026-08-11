c====================== include file "docnam.h" ========================
c
c    info from docum.F that can be used elsewhere.
#ifdef kdtravg
c    the experiment name description "expnam" is loaded in docum
#endif
#if defined cpldleg3 || defined cpldleg2
c    For leg3 or leg2 of coupled models, read the experiment name 
c     in program airsea
#endif
c    user specified tracer names are place into "trname" in docum
c
#if defined cpldleg3 || defined cpldleg2
      character*60 expnam
#endif
      character*12 trname
#ifdef seaice
      character*12 ycname
#endif
c
#if defined cpldleg3 || defined cpldleg2
      common /docnam/ expnam
#endif
      common /docnam/ trname(nt)
#ifdef seaice
     $,               ycname
#endif
c
