      subroutine hexpft(svble,gvble,kmx)
c
c    this subroutine is called by passa and passb.  hexpft calls scalar
c    and fftg.
c  definitions:
c               svble = input spectral variable
c               gvble = output grid pt variable
c               kmx   = number of horizonal fields in vor to be
c                       transformed.
c
#include "parm.h"
#include "parmjx.h"
c
      complex svble
      complex ftvar
      complex ftvark
      complex ftvar4
      complex ftvar1
      complex ftdum
c
      common / fftstr / ftvar ( mx, kx , nhem )
c   the array ftdum makes fftstr larger and avoids loader messages
      common / fftstr / ftdum ( mx, kx , nhem ,6)
c
      dimension gvble ( ix, kmx, nhem ), svble(mx,jx,kmx)
      dimension ftvark( mx, kx , nhem )
      dimension ftvar4( mx, nhem,  4  )
      dimension ftvar1( mx,  1 , nhem )
      dimension ftvara( mx, nhem ,1 )
      dimension ftvarb( mx, nhem ,1 )
      dimension ftvarc( mx, nhem ,1 )
      dimension ftvard( mx, nhem ,1 )
c
      equivalence ( ftvar (1,1,1), ftvark (1,1,1) )
      equivalence ( ftvar (1,1,1), ftvar4 (1,1,1) )
      equivalence ( ftvar (1,1,1), ftvar1 (1,1,1) )
      equivalence ( ftvar4(1,1,1), ftvara (1,1,1) )
      equivalence ( ftvar4(1,1,2), ftvarb (1,1,1) )
      equivalence ( ftvar4(1,1,3), ftvarc (1,1,1) )
      equivalence ( ftvar4(1,1,4), ftvard (1,1,1) )
c
c .....................................................
c
      if (kmx.eq.kx) then
        call scalar ( ftvark, svble, kmx )
        call fftg   ( gvble, ftvark, kmx*nhem )
      else if (kmx.eq.4) then
c
c  kmx equals four contiguous 1-level vbls; psdx,psdy,psln,zs
        call scalar ( ftvara, svble(1,1,1),  1  )
        call scalar ( ftvarb, svble(1,1,2),  1  )
        call scalar ( ftvarc, svble(1,1,3),  1  )
        call scalar ( ftvard, svble(1,1,4),  1  )
        call fftg   ( gvble, ftvar4,   4*nhem )
      else if (kmx.eq.1) then
        call scalar ( ftvar1, svble, kmx )
        call fftg   ( gvble, ftvar1, kmx*nhem )
      else
        write(6,100)
        call aquit (100)
        stop 100
      endif
c
  100 format('  kmx.ne.kx or 4 or 1 in hexpft.')
c
      return
      end
