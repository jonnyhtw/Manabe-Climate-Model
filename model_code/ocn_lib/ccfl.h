c====================== include file "ccfl.h" ==========================
c
c
c     cflcrt  = factor by which the cfl criteria must be exceeded in
c               order to stop MOM  (see blkdta.F)
c     cflum   = zonal velocity which comes closest to its cfl criteria
c     cflup   = percent of cfl criteria reached by "cflum"
c     icflu   = "i" coordinate of "cflum"
c     jcflu   = "j" coordinate of "cflum"
c     kcflu   = "k" coordinate of "cflum"
c
c     cflvm   = meridional velocity which comes closest to its cfl
c               criteria
c     cflvp   = percent of cfl criteria reached by "cflvm"
c     icflv   = "i" coordinate of "cflvm"
c     jcflv   = "j" coordinate of "cflvm"
c     kcflv   = "k" coordinate of "cflvm"
c
c     cflwm   = vertical velocity which comes closest to its cfl
c               criteria
c     cflwp   = percent of cfl criteria reached by "cflwm"
c     icflw   = "i" coordinate of "cflwm"
c     jcflw   = "j" coordinate of "cflwm"
c     kcflw   = "k" coordinate of "cflwm"
c
      common /ccfl/ cflup, cflum, icflu, jcflu, kcflu
     $,             cflvp, cflvm, icflv, jcflv, kcflv
     $,             cflwp, cflwm, icflw, jcflw, kcflw
     $,             cflcrt
c
