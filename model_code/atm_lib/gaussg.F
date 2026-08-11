      subroutine gaussg(nzero,f,wt)
c
c     this routine calculates the cosine of the colatitudes 'f' and
c     the weights 'wt' , for the gaussian quadrature with 'nzero' points
c     between pole and equator
c
      parameter ( nzermx = 100 )
c
      dimension f(nzero), wt(nzero)
c
      if (nzero .gt. nzermx) then
        write (6,60)
        stop 555
      endif
c
c  xlim is an critical number for the interation below
c   It is machine precision dependent. If the value is too small
c   the interation will NOT converge.
c   10-12 => 64 bit machine, 10-6 => 32 bit machine
#ifndef a_32_bit_machine
      xlim  = 1.0e-12
#endif
#ifdef a_32_bit_machine
      xlim  = 1.0e-6
#endif
      ir    = nzero + nzero
      pi    = 4.0 * atan(1.0)
      fi    = ir
      fi1   = fi + 1.0
      piov2 = 0.5 * pi
      fn    = piov2 / nzero
c
      do 10 i=1,nzero
        wt(i) = i - 0.5
   10 continue
c
      do 20 i=1,nzero
        f(i) = sin(wt(i) * fn + piov2)
   20 continue
c
      dn  = fi  / sqrt(4.0 * fi  * fi  - 1.0)
      dn1 = fi1 / sqrt(4.0 * fi1 * fi1 - 1.0)
      a   = dn1 * fi
      b   = dn  * fi1
      irp = ir + 1
      irm = ir - 1
c
      do 2 i=1,nzero
c
    5   continue
c
        call ordleg(g , f(i), ir )
        call ordleg(gm, f(i), irm)
        call ordleg(gp, f(i), irp)
        gt    = (f(i) * f(i) - 1.0) / (a * gp - b * gm)
        ftemp = f(i) - g * gt
        gtemp = f(i) - ftemp
        f(i)  = ftemp
        if (abs(gtemp) .gt. xlim) go to 5
    2 continue
c
      do 6 i=1,nzero
        a = 2.0 * (1.0 - f(i) * f(i))
        call ordleg(b,f(i),irm)
        b     = b * b * fi * fi
        wt(i) = a * (fi - 0.5) / b
    6 continue
c
   60 format(// ' execution halted in gaussg.increase size of array',
     *          ' maximum')
c
      return
      end
