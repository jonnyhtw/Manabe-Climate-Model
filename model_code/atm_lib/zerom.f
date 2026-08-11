      subroutine zerom(vor,kmx)
c
c  zerom zeroes any component of a complex data field (vor) of length
c   kmx which is smaller in magnitude than 10**(-20).
c
      dimension vor(1)
c
      kmxt=2*kmx
      do 1 k=1,kmxt
      if(abs(vor(k)).lt.1.e-20) vor(k)=0.0
  1   continue
      return
      end
