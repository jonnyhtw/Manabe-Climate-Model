      subroutine mapit (iunit, idim, jdim, sclfctr, addfctr, arrayin)
      implicit real*8 (a-h,o-z)
c
      dimension arrayin(idim,jdim)
      dimension prarray(192,80)
c
      write(iunit,50) sclfctr, addfctr
      nperline = 24
      is = 1
c
      do j=1,jdim
        do i=1,idim
          prarray(i,j) = arrayin(i,j)*sclfctr + addfctr
        end do
      end do
c
10    continue
c
      ie = is + nperline
      if (ie .gt. idim) ie=idim
      write(iunit,90) is,ie
c
      do j=1,jdim
        jj=jdim+1-j
        write(iunit,60) jj, (prarray(ii,jj),ii=is,ie)
      end do
c
      if (ie .NE. idim) then
       is = is + nperline
       goto 10
      endif
c
c
      return
50    format('========== multiply scalefactor =',e12.6,
     $       6x,' add factor =',e12.6)
cc60    format(i4,1x,30(1x,f4.1))
60    format(i4,1x,30(f5.1))
90    format(6x,i4,116('='),i4)
      end
