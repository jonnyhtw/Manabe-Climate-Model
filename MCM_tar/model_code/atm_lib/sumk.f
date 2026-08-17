      subroutine sumk(ain, aout, imax1, imax2, imax3,imax4)
c
      dimension ain (imax1,imax2,imax3)
      dimension aout(imax1,imax2,imax3,imax4)
c
cfpp$ noconcur R
      if (imax3 .eq. 1) then
        do 10 ls=1,imax4
         do 10 n=1,imax2
          do 10 m=1,imax1
           ain(m,n,1) = ain(m,n,1) + aout(m,n,1,ls)
   10   continue
      else
        do 20 ls=1,imax4
         do 20 k=1,imax3
          do 20 n=1,imax2
           do 20 m=1,imax1
            ain(m,n,k) = ain(m,n,k) + aout(m,n,k,ls)
   20   continue
      endif
c
      return
      end
