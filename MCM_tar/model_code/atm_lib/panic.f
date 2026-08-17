      subroutine panic( m )
c
      character  m*(*)
c
c historical gfdl routines
c
      write(6,100) m
100   format(' ',a50)
c
      call  aquit( 20)
c
      end
