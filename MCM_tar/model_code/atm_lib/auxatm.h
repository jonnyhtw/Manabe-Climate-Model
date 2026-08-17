      common / auxatm / vfu   (ix,kx), vfv   (ix,kx), vft   (ix,kx)
#ifdef moist
      common / auxatm / vfr   (ix,kx), convr (ix,kx)
#endif
      common / auxatm / convt (ix,kx)
#ifdef rad
      common / auxatm / zsnit (ix)   , ztufxt(ix)   , zsng  (ix)
      common / auxatm / zufxnt(ix)   , zshfx (ix)   , zaltnt(ix)
# ifdef cloud_forcing
      common / auxatm / zsnitc (ix)  , ztufxtc(ix)  , zsngc  (ix)
      common / auxatm / zufxntc(ix)  , zplalb (ix)  , zplalbc(ix)
# endif
#endif
#ifdef moist
      common / auxatm / fillmr(ix,kx)
#endif
      common / auxatm / dadadj(ix,kx)
#ifdef moist
      common / auxatm / relhum(ix,kx)
#endif
