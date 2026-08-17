c
      integer awrds, bwrds
      integer adsko, bdsko, adski, bdski
      integer lgo, lgi, lho, lhi
      common / dsk / adsko, adski, bdsko, bdski
      common / lglh / lgo, lgi, lho, lhi
      parameter ( awrds=ix * ivf*kx )
      parameter ( bwrds=ix * ( icv*kx + irm) )
      parameter ( kwrds=ix*(kp+2) )
c
