c====================== include file "ccflj.h" ==========================
#ifdef cflj
c
c
c     cfluj(j,k)   = maximun value of cfl criteria for the u-component at 
c                    row j, level k
c     cflvj(j,k)   = maximun value of cfl criteria for the v-component at 
c                    row j, level k
c     cflwj(j,k)   = maximun value of cfl criteria for the w-component at 
c                    row j, level k
c     imaxu(j,k)   = "i" coordinate of "cfluj" at level k
c     imaxv(j,k)   = "i" coordinate of "cflvj" at level k
c     imaxw(j,k)   = "i" coordinate of "cflwj" at level k
c
      common / cflcom / cfluj(jmt,km), cflvj(jmt,km), cflwj(jmt,km)
      common / cflcom / imaxu(jmt,km), imaxv(jmt,km), imaxw(jmt,km)
c
#endif
