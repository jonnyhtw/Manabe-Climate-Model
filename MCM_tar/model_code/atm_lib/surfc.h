      common / surfc / zin (ix), cd   (ix,lmx), gmzkx(ix)
#ifdef moist
      common / surfc / ess (ix), dm   (ix), rms  (ix)
      common / surfc / dw  (ix), fluxw(ix)
#endif
#ifdef rad
      common / surfc / hlt (ix), pkice(ix)
#endif
#ifdef moist
      common / surfc / ralb(ix), salb (ix)
#endif
#ifdef moist
      common / surfc / ekd (ix), fkd  (ix)
#endif
#ifdef moist
      common / surfc / ekdr(ix), fkdr (ix)
#endif
