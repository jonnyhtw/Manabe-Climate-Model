      subroutine ordleg(sx,coa,ir)
c
      irpp = ir + 1
      irppm = irpp - 1
      delta =   acos(coa)
      sqr2= sqrt(2.)
c
      theta=delta
      c1=sqr2
      do 20 n=1,irppm
      fn=n
      fn2=fn+fn
      fn2sq=fn2*fn2
      c1=c1* sqrt(1.000-1.000/fn2sq)
   20 continue
c
      n=irppm
      ang=fn*theta
      s1=0.
      c4=1.0
      a=-1.0
      b=0.
      n1=n+1
      do 27 kk=1,n1,2
      k=kk-1
      if (k.eq.n) c4=0.5*c4
      s1=s1+c4* cos(ang)
      a=a+2.0
      b=b+1.0
      fk=k
      ang=theta*(fn-fk-2.0)
      c4=(a*(fn-b+1.000)/(b*(fn2-a)))*c4
   27 continue
      sx=s1*c1
c
      return
      end
