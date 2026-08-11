      subroutine getime(acc,elaps)
c
#ifdef wall_clock_time
      integer time
#endif
      data alst/0.0/,init/0/,carry/0.0/
c
#ifdef cray_clk
      acc=second()
#endif
#ifdef wall_clock_time
c  the line below does not work on a sparc (sun 4)
      acc = time()
c     acc = 0
#endif
      acc   = acc-carry
      elaps = acc-alst
      alst  = acc
      if(init .ne. 0)  go to 100
      elaps = 0.0
      carry = acc
      acc   = 0.0
      alst  = 0.0
      init  = 1
  100 continue
      return
      end
