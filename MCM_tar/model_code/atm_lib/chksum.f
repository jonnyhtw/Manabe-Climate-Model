      subroutine chksum(array,n)
c
c subroutine chksum is called from energy and computes a checksum
c  of the global integrals.  the subroutine closely follows a
c  version developed by elaine jacobsen (cdc) and is described
c  in cuab user memo #3 (april 22, 1983).
c
c  modified to just do integer sum 21 nov 89....rjs
c
c     dimension array(n)
cSGI
cc    integer array(n)
      integer*8 array(n), mm
cSGI
      mm = 0
cSGI  m  = 0
      do 20 i=1,n
c       mm = mm + q8sextb(array(i),32,32) + q8sextb(array(i),32,0)
c       m = xor(m,array(i))
        mm = mm + array(i)
   20 continue
      l = m + mm
      write(6,100) l
c
  100 format(' checksum=',2x,z20/)
c
      return
      end
