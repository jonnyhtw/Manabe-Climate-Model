      subroutine aspace (iunit)
c
  100 format (1x)
c
    1 read (iunit,100,end=99)
      go to 1
c
   99 backspace iunit
      return
      end
